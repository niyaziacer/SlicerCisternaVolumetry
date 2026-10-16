"""Shared cisterna segmentation core.

This module is imported by both the CLI tool (`cisterna_analysis.py`) and the
3D Slicer extension (`CisternaMagnaVolumetry`) - both interfaces use the SAME
code path so results stay consistent. A copy of this file is vendored inside
the Slicer extension (see
`CisternaMagnaVolumetry/CisternaMagnaVolumetryLib/cisterna_core.py`) -
if you change the logic, update BOTH copies.

Method summary:
- CSF mask source: native_tissues (label=1) minus ventricular labels -> "external CSF"
- Anatomical landmarks: native_structures (4th ventricle=11, brainstem=35)
- Cistern extraction: mm-based ROI boxes + morphological filtering + connected
  component selection
- VolBrain's atlas has NO cistern-level labels (see repo README) - this is why a
  landmark-based heuristic is used instead of an atlas shortcut.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
from nibabel.orientations import (
    apply_orientation,
    axcodes2ornt,
    io_orientation,
    ornt_transform,
)
from scipy import ndimage as ndi


VENTRICULAR_LABELS = (4, 11, 49, 50, 51, 52)

# The ROI boxes in make_cistern_masks() rely on axis index (0=P-A, 1=S-I, 2=L-R),
# calibrated against this axis order. Different VolBrain outputs (different
# scanner/protocol, e.g. 7T vs 1.5T) can come in a different axis orientation
# (e.g. L,P,S instead of P,S,R) - so we normalize to this common orientation
# before analysis, and reorient outputs back to the original input orientation
# when saving.
CANONICAL_AXCODES = ("P", "S", "R")

CISTERN_LABEL_MAP = {"cisterna_magna": 1, "cisterna_ambiens": 2, "cisterna_interpeduncularis": 3}
CISTERN_COLORS = {
    "cisterna_magna": (0.85, 0.10, 0.10),
    "cisterna_ambiens": (0.95, 0.85, 0.10),
    "cisterna_interpeduncularis": (0.10, 0.75, 0.85),
}


@dataclass
class Landmark:
    name: str
    bbox_min: np.ndarray
    bbox_max: np.ndarray
    centroid: np.ndarray


def bbox_and_centroid(mask: np.ndarray, name: str) -> Landmark:
    idx = np.argwhere(mask)
    if idx.size == 0:
        raise ValueError(f"Landmark not found: {name}")
    return Landmark(
        name=name,
        bbox_min=idx.min(axis=0),
        bbox_max=idx.max(axis=0),
        centroid=idx.mean(axis=0),
    )


def clamp_range(start: int, end: int, max_size: int) -> Tuple[int, int]:
    s = max(0, start)
    e = min(max_size - 1, end)
    if s > e:
        s, e = e, s
    return s, e


def box_mask(shape: Tuple[int, int, int], x: Tuple[int, int], y: Tuple[int, int], z: Tuple[int, int]) -> np.ndarray:
    m = np.zeros(shape, dtype=bool)
    m[x[0] : x[1] + 1, y[0] : y[1] + 1, z[0] : z[1] + 1] = True
    return m


def label_components(mask: np.ndarray) -> Tuple[np.ndarray, int]:
    structure = np.ones((3, 3, 3), dtype=bool)
    return ndi.label(mask, structure=structure)


def select_components_by_seeds(
    candidate_mask: np.ndarray,
    seeds: List[np.ndarray],
    max_components_per_seed: int = 1,
) -> np.ndarray:
    labeled, num = label_components(candidate_mask)
    if num == 0:
        return np.zeros_like(candidate_mask, dtype=bool)

    objects = ndi.find_objects(labeled)
    comp_centroids = {}
    comp_sizes = {}
    for i, slc in enumerate(objects, start=1):
        if slc is None:
            continue
        comp_vox = np.argwhere(labeled == i)
        comp_centroids[i] = comp_vox.mean(axis=0)
        comp_sizes[i] = int(comp_vox.shape[0])

    selected_labels = set()
    for seed in seeds:
        distances = []
        for c in comp_centroids:
            d = np.linalg.norm(comp_centroids[c] - seed)
            distances.append((d, -comp_sizes[c], c))
        distances.sort()
        for _, _, comp_id in distances[:max_components_per_seed]:
            selected_labels.add(comp_id)

    if not selected_labels:
        largest = max(comp_sizes.items(), key=lambda kv: kv[1])[0]
        selected_labels.add(largest)

    out = np.isin(labeled, list(selected_labels))
    return out


def clean_mask(mask: np.ndarray, open_iter: int = 1, close_iter: int = 1) -> np.ndarray:
    structure = np.ones((3, 3, 3), dtype=bool)
    m = ndi.binary_opening(mask, structure=structure, iterations=open_iter)
    m = ndi.binary_closing(m, structure=structure, iterations=close_iter)
    return m


def make_cistern_masks(
    external_csf: np.ndarray,
    struct_data: np.ndarray,
    voxel_sizes: Tuple[float, float, float] = (1.0, 1.0, 1.0),
) -> Dict[str, np.ndarray]:
    """external_csf and struct_data must be in CANONICAL_AXCODES axis order
    (see get_reorientation / Reorientation)."""
    shape = external_csf.shape
    # voxel_sizes must be in CANONICAL_AXCODES order (P-A, S-I, L-R), mm/voxel.
    vx0, vx1, vx2 = voxel_sizes

    # ROI offsets below are defined in mm (calibrated against a 1mm isotropic
    # reference dataset) and converted to voxel counts using voxel_sizes. This
    # keeps the boxes at the same PHYSICAL size across resolutions (e.g. 7T
    # 0.7mm isotropic) - using a fixed voxel count would shrink the box by
    # ~30% on 0.7mm data.
    def mm0(mm: float) -> int:
        return int(round(mm / vx0))

    def mm1(mm: float) -> int:
        return int(round(mm / vx1))

    def mm2(mm: float) -> int:
        return int(round(mm / vx2))

    fourth = bbox_and_centroid(struct_data == 11, "4th_ventricle")
    brainstem = bbox_and_centroid(struct_data == 35, "brainstem")
    midline = int(round(fourth.centroid[2]))

    # Cisterna magna ROI (posterior-inferior neighborhood of the 4th ventricle,
    # toward the foramen magnum).
    # NOTE: lower bound = fourth.bbox_min[0] - cisterna magna spans the same P-A
    # range as the 4th ventricle but lies inferior to it (confirmed against a
    # manual mask: i=134 vs ventricle bbox_min[0]=133, 1-voxel difference, on
    # the 1mm reference dataset). Does not cross the brainstem's posterior
    # boundary (bbox_max[0]), so it doesn't bleed into the prepontine/
    # premedullary space (brainstem bbox: 114-142, this value is typically >=133).
    x_magna = clamp_range(int(fourth.bbox_min[0]), int(fourth.bbox_max[0]) + mm0(45), shape[0])
    y_magna = clamp_range(int(brainstem.bbox_min[1]) - mm1(20), int(fourth.centroid[1]) + mm1(10), shape[1])
    z_magna = clamp_range(midline - mm2(30), midline + mm2(30), shape[2])
    roi_magna = box_mask(shape, x_magna, y_magna, z_magna)
    seed_magna = np.array([x_magna[0] + mm0(8), int(fourth.centroid[1]) - mm1(10), midline])

    # Cisterna interpeduncularis ROI (anterior-inferior midline of the midbrain)
    x_inter = clamp_range(int(brainstem.bbox_min[0]) - mm0(35), int(brainstem.bbox_min[0]) + mm0(8), shape[0])
    y_inter = clamp_range(int(brainstem.bbox_max[1]) - mm1(30), int(brainstem.bbox_max[1]) + mm1(18), shape[1])
    z_inter = clamp_range(midline - mm2(14), midline + mm2(14), shape[2])
    roi_inter = box_mask(shape, x_inter, y_inter, z_inter)
    # NOTE: seed P-A anchored at brainstem.bbox_min[0] itself (not -10mm in front of it).
    # Anatomical reasoning: the interpeduncular fossa sits essentially AT the anterior
    # brainstem surface, between the cerebral peduncles - confirmed against a manual mask
    # on the 1mm reference dataset (manual centroid P-A=115.1 vs brainstem.bbox_min[0]=114,
    # ~1-voxel difference). The old -10mm offset placed the seed close to a spurious small
    # CSF blob instead (10.1 vs 11.1 voxel distance - a fragile near-tie that picked the
    # wrong connected component, missing the true structure entirely, Dice=0).
    seed_inter = np.array([int(brainstem.bbox_min[0]), int(brainstem.bbox_max[1]) - mm1(8), midline])

    # Cisterna ambiens ROI (bilateral, wraps the lateral midbrain surface).
    # Per anatomical reference (Wikipedia "Ambient cistern"; laterally bounded by
    # the mesial temporal lobe - uncus/parahippocampal gyrus/dentate gyrus/fornix
    # fimbria - and medially by the lateral midbrain surface), the ambient cistern
    # is CONTINUOUS anteriorly with the interpeduncular cistern and posteriorly
    # with the quadrigeminal cistern.
    # NOTE: x_amb lower bound changed from `centroid[0]-5mm` to `bbox_min[0]`
    # (same anterior anchor as the interpeduncularis seed) so the box reaches the
    # anterior junction with the interpeduncular cistern instead of starting
    # ~25mm posterior to it. Upper bound reduced from `bbox_max[0]+28mm` to
    # `+15mm`: the quadrigeminal cistern junction is a short posterior transition,
    # not a deep extension - the old +28mm value let the candidate region reach
    # into cerebellar exterior/white matter and occipital lingual gyrus territory
    # (found via neighborhood-label analysis on the 7T case).
    #
    # TRIED AND REVERTED: widening y_amb's superior bound (+20mm -> +45mm) and
    # removing the z midline gap, to deliberately fold the quadrigeminal cistern
    # into cisterna_ambiens (Wikipedia/Gray's: "some sources define the
    # quadrigeminal cistern as part of the ambient cistern" - no universally
    # agreed dividing line). The wider box DID contain the anatomically correct
    # merged bilateral+midline structure as its largest component (case 1:
    # 11114 vox spanning both hemispheres + midline; case 2: 22248 vox, same
    # pattern) - but `select_components_by_seeds()` picks the component
    # NEAREST to the seed regardless of size, and the wider box also creates
    # many small spurious CSF specks. On case 2, two tiny specks (88 + 27
    # voxels) beat the true 22248-voxel structure by only 5-9 voxels of
    # centroid distance, collapsing the result to 115 voxels; on case 1, the
    # correct component itself absorbed extra cortical sulcal CSF (10257 vox,
    # touching posterior cingulate gyrus/cerebellum exterior). Unstable in
    # opposite directions on the two cases - reverted rather than shipped
    # without a robust (e.g. largest-component-within-radius, not
    # nearest-regardless-of-size) selection rule and validation data to test it
    # against. See docs/validation for the full writeup.
    x_amb = clamp_range(int(brainstem.bbox_min[0]), int(brainstem.bbox_max[0]) + mm0(15), shape[0])
    y_amb = clamp_range(int(brainstem.bbox_max[1]) - mm1(28), int(brainstem.bbox_max[1]) + mm1(20), shape[1])
    z_amb_l = clamp_range(midline - mm2(38), midline - mm2(10), shape[2])
    z_amb_r = clamp_range(midline + mm2(10), midline + mm2(38), shape[2])
    roi_amb = box_mask(shape, x_amb, y_amb, z_amb_l) | box_mask(shape, x_amb, y_amb, z_amb_r)
    seed_amb_l = np.array([int(brainstem.centroid[0]), int(brainstem.bbox_max[1]) - mm1(8), midline - mm2(18)])
    seed_amb_r = np.array([int(brainstem.centroid[0]), int(brainstem.bbox_max[1]) - mm1(8), midline + mm2(18)])

    candidates = {
        "cisterna_magna": clean_mask(external_csf & roi_magna),
        "cisterna_interpeduncularis": clean_mask(external_csf & roi_inter),
        "cisterna_ambiens": clean_mask(external_csf & roi_amb),
    }

    masks = {
        "cisterna_magna": select_components_by_seeds(candidates["cisterna_magna"], [seed_magna]),
        "cisterna_interpeduncularis": select_components_by_seeds(candidates["cisterna_interpeduncularis"], [seed_inter]),
        "cisterna_ambiens": select_components_by_seeds(candidates["cisterna_ambiens"], [seed_amb_l, seed_amb_r]),
    }

    # Resolve overlaps: priority interpeduncularis > ambiens > magna (reduces
    # incorrect double-labeling near the midline).
    inter = masks["cisterna_interpeduncularis"]
    amb = masks["cisterna_ambiens"] & (~inter)
    mag = masks["cisterna_magna"] & (~inter) & (~amb)

    masks["cisterna_interpeduncularis"] = inter
    masks["cisterna_ambiens"] = amb
    masks["cisterna_magna"] = mag

    return masks


def volume_from_mask(mask: np.ndarray, voxel_volume_mm3: float) -> Dict[str, float]:
    vox = int(mask.sum())
    mm3 = vox * voxel_volume_mm3
    cm3 = mm3 / 1000.0
    return {"voxels": vox, "volume_mm3": mm3, "volume_cm3": cm3}


class Reorientation:
    """Holds the forward/backward transform between the input's original axis
    orientation and CANONICAL_AXCODES, plus voxel sizes reordered into
    canonical axis order."""

    def __init__(self, affine: np.ndarray, target_axcodes: Tuple[str, str, str] = CANONICAL_AXCODES):
        orig_ornt = io_orientation(affine)
        target_ornt = axcodes2ornt(target_axcodes)
        self.to_canonical = ornt_transform(orig_ornt, target_ornt)
        self.from_canonical = ornt_transform(target_ornt, orig_ornt)
        self._source_axis_for_canonical = [int(self.to_canonical[i, 0]) for i in range(3)]

    def canonical_voxel_sizes(self, voxel_sizes: Tuple[float, float, float]) -> Tuple[float, float, float]:
        return tuple(float(voxel_sizes[i]) for i in self._source_axis_for_canonical)

    def to_canon(self, data: np.ndarray) -> np.ndarray:
        return apply_orientation(data, self.to_canonical)

    def from_canon(self, data: np.ndarray) -> np.ndarray:
        return apply_orientation(data, self.from_canonical)


def run_cistern_analysis(
    t1_data: np.ndarray,
    tissues_data: np.ndarray,
    structures_data: np.ndarray,
    affine: np.ndarray,
    voxel_sizes: Tuple[float, float, float],
) -> Tuple[Dict[str, np.ndarray], np.ndarray, Dict[str, Dict[str, float]], Dict[str, float]]:
    """Full analysis: orientation normalization + cistern masks + volumes.

    Inputs (t1_data, tissues_data, structures_data) must be in their ORIGINAL
    axis order (nibabel or Slicer array order doesn't matter, as long as
    `affine` is defined consistently with that order). Output masks are
    returned in that same original order.

    Returns: (cistern_masks, external_csf, volumes, total_csf)
    """
    reorient = Reorientation(affine)
    tissues_canon = reorient.to_canon(tissues_data.astype(np.int16))
    structures_canon = reorient.to_canon(structures_data.astype(np.int16))
    canonical_voxel_sizes = reorient.canonical_voxel_sizes(voxel_sizes)

    csf_tissue_canon = tissues_canon == 1
    ventricular_canon = np.isin(structures_canon, VENTRICULAR_LABELS)
    external_csf_canon = csf_tissue_canon & (~ventricular_canon)

    cistern_masks_canon = make_cistern_masks(external_csf_canon, structures_canon, canonical_voxel_sizes)

    cistern_masks = {
        name: reorient.from_canon(mask.astype(np.int16)).astype(bool) for name, mask in cistern_masks_canon.items()
    }
    external_csf = reorient.from_canon(external_csf_canon.astype(np.int16)).astype(bool)

    voxel_volume_mm3 = float(np.prod(voxel_sizes))
    volumes = {name: volume_from_mask(mask, voxel_volume_mm3) for name, mask in cistern_masks.items()}
    total_csf = volume_from_mask(external_csf, voxel_volume_mm3)

    return cistern_masks, external_csf, volumes, total_csf
