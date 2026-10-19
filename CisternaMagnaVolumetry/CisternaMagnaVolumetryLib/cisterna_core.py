"""Shared cisterna segmentation core.

This module is imported by both the CLI tool (`cisterna_analysis.py`) and the
3D Slicer extension (`CisternaMagnaVolumetry`) - both interfaces use the SAME
code path so results stay consistent. A copy of this file is vendored inside
the Slicer extension (see
`CisternaMagnaVolumetry/CisternaMagnaVolumetryLib/cisterna_core.py`) -
if you change the logic, update BOTH copies.

Method summary (v2, October 2026):
- CSF mask source: native_tissues (label=1) minus ventricular labels -> "external CSF"
- Anatomical landmarks from native_structures (vol2Brain / Neuromorphometrics labels):
  brainstem=35, 4th ventricle=11, ventral DC=61/62, cerebellum=38-41, vermis=71-73
- The midbrain level is found from the brainstem itself: the pontomesencephalic
  junction is the highest slice where the anterior brainstem edge jumps forward
  (the pons). Cisterns are then defined relative to the midbrain surface:
  interpeduncular = midline CSF in front of the midbrain floor (not below the
  junction, so the prepontine cistern is excluded; not more than 8 mm in front of
  the upper pons, so the suprasellar cistern is excluded); ambiens = lateral CSF
  within 8 mm of the midbrain surface at midbrain level; magna = CSF below the
  cerebellum and behind the medulla, near the midline.
- The previous box-based method is kept as make_cistern_masks_v1() for comparison.
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

CISTERN_LABEL_MAP = {"cisterna_magna": 1, "cisterna_ambiens": 2, "cisterna_interpeduncularis": 3, "cisterna_cruralis": 4, "cisterna_quadrigeminalis": 5}
CISTERN_COLORS = {
    "cisterna_magna": (0.85, 0.10, 0.10),
    "cisterna_ambiens": (0.95, 0.85, 0.10),
    "cisterna_interpeduncularis": (0.10, 0.75, 0.85),
    "cisterna_cruralis": (0.30, 0.80, 0.30),
    "cisterna_quadrigeminalis": (0.75, 0.35, 0.85),
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


def make_cistern_masks_v1(
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



# --- v2 (landmark/surface-based) ------------------------------------------------
IP_DOWN, IP_HALFW = 0, 8          # interpeduncular: mm below the junction, half-width (mm)
AMB_DOWN = 6                      # ambiens: mm below the pontomesencephalic junction (infratentorial part)
QUAD_UP, QUAD_DEPTH, QUAD_HALFW = 12, 15, 12   # quadrigeminal cistern beside the tectum (mm)
QUAD_DEPTH_UP, QUAD_HALFW_UP = 30, 7             # above the brainstem (pineal region, vein of Galen)
CRURAL_FRAC = 0.4                 # anterior fraction of the midbrain A-P extent taken as the cerebral peduncle
MESIAL_TEMPORAL = (31, 32, 47, 48, 116, 117, 170, 171)   # amygdala, hippocampus, entorhinal area, parahippocampal gyrus
BS, V4, VDC, CBL, VERMIS = 35, 11, (61, 62), (38, 39, 40, 41), (71, 72, 73)

def _largest_touching(cand, touch):
    lab, n = ndi.label(cand, np.ones((3,3,3), bool))
    if n == 0: return np.zeros_like(cand)
    ids = np.unique(lab[touch & cand]); ids = ids[ids > 0]
    if len(ids) == 0: return np.zeros_like(cand)
    sizes = ndi.sum(cand, lab, ids)
    keep = ids[sizes >= max(20, 0.05 * sizes.max())]
    return np.isin(lab, keep)

def make_cistern_masks(ext, st, vs=(1.0, 1.0, 1.0), info=None):
    """Landmark/surface-based cistern masks (v2). ext and st must be in CANONICAL_AXCODES order (P, S, R);
    vs are voxel sizes in that order (mm)."""
    v0, v1, v2 = vs; sh = ext.shape
    bs = st == BS; mb_like = bs | np.isin(st, VDC)
    if not bs.any():
        raise ValueError("Landmark not found: brainstem (label 35)")
    if not (st == V4).any():
        raise ValueError("Landmark not found: 4th_ventricle (label 11)")
    # midline: brainstem centre (L-R) in its upper half - more central than the 4th ventricle in tilted heads
    q = np.argwhere(bs); q = q[q[:, 1] >= np.median(q[:, 1])]
    mid = int(round(q[:, 2].mean()))
    # --- pontomesencephalic junction: anterior edge of the brainstem jumps forward at the pons
    Ss = np.where(bs.any(axis=(0, 2)))[0]; top = Ss.max()
    ant = {S: np.argwhere(bs[:, S, :])[:, 0].min() for S in Ss}
    a_mid = np.median([ant[S] for S in Ss if S > top - round(10 / v1)])
    pons = [S for S in Ss if ant[S] <= a_mid - 6 / v0]
    pmj = max(pons) if pons else top - round(18 / v1)
    s_hi = top + round(3 / v1)
    s_mid = pmj + 1                                   # lowest midbrain slice
    # the ambient cistern also extends below the tentorium (lower part, superior cerebellar artery) beside the
    # upper pons: the lateral slab starts AMB_DOWN mm below the junction
    s_lo = pmj + 1 - round(AMB_DOWN / v1)
    slab = np.zeros(sh, bool); slab[:, s_lo:s_hi + 1, :] = True
    # in-plane distance (mm) to the midbrain (brainstem ∪ ventral DC) within the slab
    mb = mb_like & slab
    d_mb = np.full(sh, 1e9)
    for S in range(s_lo, s_hi + 1):
        if mb[:, S, :].any():
            d_mb[:, S, :] = ndi.distance_transform_edt(~mb[:, S, :], sampling=(v0, v2))
    P = np.arange(sh[0])[:, None, None] * np.ones(sh); X = (np.arange(sh[2])[None, None, :] - mid) * v2 * np.ones(sh)
    # anterior / posterior edges of the midbrain per slice (mm in P index units)
    a_edge = np.full(sh[1], np.nan); p_edge = np.full(sh[1], np.nan); b_ant = np.full(sh[1], np.nan)
    for S in range(s_lo, s_hi + 1):
        q = np.argwhere(mb[:, S, :])
        if len(q): a_edge[S], p_edge[S] = q[:, 0].min(), q[:, 0].max()
        r = np.argwhere(bs[:, S, :])
        if len(r): b_ant[S] = r[:, 0].min()
    A = a_edge[None, :, None] * np.ones(sh); Pe = p_edge[None, :, None] * np.ones(sh); Ba = b_ant[None, :, None] * np.ones(sh)
    # --- interpeduncular: midline CSF in front of the brainstem floor, between the peduncles, up to 5 mm in front of them
    # S range: from 3 mm below the pontomesencephalic junction (pontomesencephalic sulcus) to the top of the brainstem;
    # P range: from 8 mm in front of the pons' upper anterior surface back to the fossa floor (midbrain anterior surface)
    S_ = np.arange(sh[1])[None, :, None] * np.ones(sh)
    a_lim = ant[pmj] - 8 / v0
    cand_i = ext & (S_ >= pmj - round(IP_DOWN / v1)) & (S_ <= s_hi) & (np.abs(X) <= IP_HALFW) & (P >= a_lim) & (P <= a_mid + 1)
    inter = _largest_touching(cand_i, cand_i & (np.abs(X) <= 3))
    # --- lateral cisterns (crural + ambiens): CSF between the midbrain and the mesial temporal lobe
    # (uncus, parahippocampal gyrus, hippocampus/fimbria). A voxel belongs here if it is within 8 mm of the
    # midbrain and not clearly on the temporal side of the gap (d_mb <= d_temporal + 3 mm).
    tl = np.isin(st, MESIAL_TEMPORAL) & slab
    d_t = np.full(sh, 1e9)
    for S in range(s_lo, s_hi + 1):
        if tl[:, S, :].any():
            d_t[:, S, :] = ndi.distance_transform_edt(~tl[:, S, :], sampling=(v0, v2))
    lat = ext & slab & ~inter & (np.abs(X) >= 5) & (d_mb <= 8) & (d_mb <= d_t + 3)
    # A-P division (Rhoton): crural cistern lies beside the cerebral peduncle (P2A segment of the PCA), the ambient
    # cistern beside the rest of the midbrain up to the lateral end of the quadrigeminal plate (P2P); behind that
    # is the quadrigeminal cistern, which is not included. The peduncle is approximated as the anterior
    # CRURAL_FRAC of the midbrain's A-P extent in each slice.
    # In the upper slices the ventral DC (hypothalamus) reaches far forward; do not let it pull the anterior edge
    # of the "midbrain" more than 13 mm in front of the brainstem's own anterior surface.
    A2 = np.fmax(A, Ba - 13 / v0)
    split = A2 + CRURAL_FRAC * (Pe - A2)
    cand_c = lat & (P >= A2 - 5 / v0) & (P < split) & (S_ >= s_mid)
    split = np.where(S_ >= s_mid, split, -np.inf)      # below the midbrain everything lateral is ambiens
    cand_a = lat & (P >= split) & (P <= Pe)
    crur = _largest_touching(cand_c, cand_c & (d_mb <= 2))
    amb = _largest_touching(cand_a, cand_a & (d_mb <= 2))
    # --- quadrigeminal cistern (cisterna venae magnae cerebri): midline CSF behind the tectum, from the lowest
    # midbrain slice up to QUAD_UP mm above the top of the brainstem (pineal region, below the splenium), at most
    # QUAD_DEPTH mm behind the posterior surface of the midbrain and within QUAD_HALFW mm of the midline.
    # Laterally it is continuous with the ambient cistern; voxels already in the ambient cistern stay there.
    post = np.full(sh[1], np.nan)
    for S in range(sh[1]):
        r = np.argwhere(bs[:, S, :])
        if len(r): post[S] = r[:, 0].max()
    last = post[top]
    post = np.where(np.isnan(post), last, post)
    Qp = post[None, :, None] * np.ones(sh)
    S_ = np.arange(sh[1])[None, :, None] * np.ones(sh)
    # beside the tectum the cistern is wide and shallow (it continues into the superior cerebellar cistern
    # further back, which is not included); above the brainstem (pineal region, vein of Galen) it is narrow and
    # runs back under the splenium
    upper = S_ > top
    halfw = np.where(upper, QUAD_HALFW_UP, QUAD_HALFW)
    front = np.where(upper, Qp + 4 / v0, Qp - 2 / v0)
    depth = np.where(upper, QUAD_DEPTH_UP, QUAD_DEPTH)
    cand_q = ext & ~amb & ~crur & (S_ >= s_mid) & (S_ <= top + round(QUAD_UP / v1)) & (np.abs(X) <= halfw) \
        & (P >= front) & (P <= Qp + depth / v0)
    seed_q = cand_q & (np.abs(X) <= 3) & (P <= Qp + 6 / v0) & (S_ <= top)
    quad = _largest_touching(cand_q, seed_q)
    # --- magna: CSF below the cerebellum and behind the medulla
    cbl = np.isin(st, CBL + VERMIS)
    any_c = cbl.any(axis=1); low_c = np.where(any_c, cbl.argmax(axis=1), sh[1])        # lowest cerebellar S per (P,R) column
    below_c = np.arange(sh[1])[None, :, None] <= (low_c[:, None, :] + 1)
    v4 = np.argwhere(st == V4); v4_low = v4[:, 1].min()
    bs_post = np.full(sh[1], np.nan)
    for S in range(sh[1]):
        r = np.argwhere(bs[:, S, :])
        if len(r): bs_post[S] = r[:, 0].max()
    bs_post = np.where(np.isnan(bs_post), np.nanmax(bs_post), bs_post)
    Bp = bs_post[None, :, None] * np.ones(sh)
    S_ = np.arange(sh[1])[None, :, None] * np.ones(sh)
    cand_m = ext & below_c & (P >= Bp - 2 / v0) & (np.abs(X) <= 15) & (S_ <= v4_low + 8 / v1)
    seed_m = cand_m & (np.abs(X) <= 4)
    mag = _largest_touching(cand_m, seed_m)
    if info is not None: info.update(pmj=int(pmj), top=int(top), mid=mid, v4_low=int(v4_low))
    return {"cisterna_magna": mag & ~inter & ~amb & ~crur & ~quad, "cisterna_ambiens": amb & ~inter,
            "cisterna_cruralis": crur & ~inter & ~amb, "cisterna_interpeduncularis": inter,
            "cisterna_quadrigeminalis": quad & ~amb & ~inter}


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
