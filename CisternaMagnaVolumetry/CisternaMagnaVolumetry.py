"""3D Slicer scripted module: CisternaMagnaVolumetry.

Computes cisterna magna, cisterna ambiens, and cisterna interpeduncularis
volumes from VolBrain segmentation output (native_tissues + native_structures),
adds the result as a colored/labeled segmentation to the scene, and displays a
volume table.

Usage: first load the T1, native_tissues*.nii.gz, and native_structures*.nii.gz
files into Slicer (drag and drop / Add Data). Then, in this module, pick the
corresponding volumes from the three dropdowns and click "Calculate".

The analysis logic lives in `CisternaMagnaVolumetryLib/cisterna_core.py` (vendored copy;
the source of truth is `src/cisterna_core.py` in the CLI repo - keep both
copies in sync). See the CLI repo README for method, limitations, and
validation results.
"""

import logging
import os
import sys

import numpy as np
import vtk

import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
    ScriptedLoadableModuleWidget,
)
from slicer.util import VTKObservationMixin

try:
    import qt
    import ctk
except ImportError:
    # Not running inside Slicer (e.g. importing this file for a syntax check).
    qt = None
    ctk = None


DISPLAY_NAMES = {
    "cisterna_magna": "Cisterna magna",
    "cisterna_ambiens": "Cisterna ambiens",
    "cisterna_interpeduncularis": "Cisterna interpeduncularis",
    "cisterna_cruralis": "Cisterna cruralis",
    "cisterna_quadrigeminalis": "Cisterna quadrigeminalis",
}
CISTERN_ORDER = ["cisterna_magna", "cisterna_interpeduncularis", "cisterna_cruralis", "cisterna_ambiens", "cisterna_quadrigeminalis"]


class CisternaMagnaVolumetry(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Cisterna Magna Volumetry"
        self.parent.categories = ["Segmentation"]
        self.parent.dependencies = []
        self.parent.contributors = ["Niyazi (VolBrain/7T cisterna volumetry project)"]
        self.parent.helpText = (
            "Computes the volumes of the cisterna magna, interpeduncularis, cruralis, ambiens and "
            "quadrigeminalis from either VolBrain/vol2Brain output (native_tissues + "
            "native_structures) or a SynthSeg 2.0 segmentation with cortical parcellation "
            "(SlicerSynthSeg, FreeSurfer mri_synthseg --parc, or turkbrain.net). The method is "
            "landmark-based (VolBrain's atlas has no cistern-level labels): the cisterns are "
            "defined relative to the midbrain and cerebellum surfaces (method v2, Oct 2026). "
            "See the README for the validation."
        )
        self.parent.acknowledgementText = ""


class CisternaMagnaVolumetryWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        # --- Input section ---
        parametersCollapsibleButton = ctk.ctkCollapsibleButton()
        parametersCollapsibleButton.text = "Input"
        self.layout.addWidget(parametersCollapsibleButton)
        parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

        self.sourceCombo = qt.QComboBox()
        self.sourceCombo.addItem("VolBrain / vol2Brain (native_tissues + native_structures)", "volbrain")
        self.sourceCombo.addItem("SynthSeg 2.0 + cortical parcellation", "synthseg")
        parametersFormLayout.addRow("Segmentation source:", self.sourceCombo)

        self.t1Selector = self._makeVolumeSelector(
            "The T1 volume that was segmented (reference geometry; with SynthSeg also used to recover "
            "thin CSF around the midbrain)."
        )
        parametersFormLayout.addRow("T1 volume:", self.t1Selector)

        self.tissuesSelector = self._makeVolumeSelector("VolBrain native_tissues_job*.nii.gz")
        self.tissuesLabel = qt.QLabel("native_tissues:")
        parametersFormLayout.addRow(self.tissuesLabel, self.tissuesSelector)

        self.structuresSelector = self._makeVolumeSelector("VolBrain native_structures_job*.nii.gz")
        self.structuresLabel = qt.QLabel("native_structures:")
        parametersFormLayout.addRow(self.structuresLabel, self.structuresSelector)

        self.synthsegSelector = self._makeVolumeSelector(
            "SynthSeg 2.0 result run WITH cortical parcellation: a label volume (*_synthseg.nii.gz) or the "
            "segmentation created by the SlicerSynthSeg module.",
            nodeTypes=["vtkMRMLLabelMapVolumeNode", "vtkMRMLScalarVolumeNode", "vtkMRMLSegmentationNode"],
        )
        self.synthsegLabel = qt.QLabel("SynthSeg result:")
        parametersFormLayout.addRow(self.synthsegLabel, self.synthsegSelector)

        self.sourceHint = qt.QLabel(
            "Make the SynthSeg segmentation with SlicerSynthSeg (model: SynthSeg 2.0 + cortical "
            "parcellation), with FreeSurfer (<code>mri_synthseg --i T1.nii.gz --o seg.nii.gz --parc</code>) "
            "or on turkbrain.net. Without parcellation the parahippocampal/entorhinal limits are missing."
        )
        self.sourceHint.setWordWrap(True)
        self.sourceHint.setStyleSheet("color: gray; font-size: 11px;")
        parametersFormLayout.addRow(self.sourceHint)

        self.sourceCombo.connect("currentIndexChanged(int)", self.onSourceChanged)
        self.applyButton = qt.QPushButton("Calculate")
        self.applyButton.toolTip = "Calculate cistern volumes and add them as a segmentation to the scene."
        parametersFormLayout.addRow(self.applyButton)

        self.statusLabel = qt.QLabel("")
        self.statusLabel.setWordWrap(True)
        parametersFormLayout.addRow(self.statusLabel)

        # --- Results section ---
        resultsCollapsibleButton = ctk.ctkCollapsibleButton()
        resultsCollapsibleButton.text = "Results"
        self.layout.addWidget(resultsCollapsibleButton)
        resultsLayout = qt.QVBoxLayout(resultsCollapsibleButton)

        self.resultsTable = qt.QTableWidget()
        self.resultsTable.setColumnCount(4)
        self.resultsTable.setHorizontalHeaderLabels(["Cistern", "Voxels", "Volume (mm³)", "Volume (cm³)"])
        self.resultsTable.horizontalHeader().setStretchLastSection(True)
        self.resultsTable.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)
        self.resultsTable.setRowCount(0)
        resultsLayout.addWidget(self.resultsTable)

        noteLabel = qt.QLabel(
            "Note: landmark-based method (v2). Neither VolBrain nor SynthSeg has cistern labels; the cisterns "
            "are cut out of the CSF mask relative to the midbrain and cerebellum. Validated against manual "
            "masks on 2 cases only. The cisterna magna ends at the lower edge of the segmentation "
            "(about the foramen magnum); VolBrain may crop it higher. Do not mix VolBrain and SynthSeg "
            "values in one study. Check every result visually."
        )
        noteLabel.setWordWrap(True)
        noteLabel.setStyleSheet("color: gray; font-size: 11px;")
        resultsLayout.addWidget(noteLabel)

        self.layout.addStretch(1)

        self.applyButton.connect("clicked(bool)", self.onApplyButton)
        saved = qt.QSettings().value("CisternaMagnaVolumetry/source", "volbrain")
        idx = self.sourceCombo.findData(saved)
        self.sourceCombo.setCurrentIndex(idx if idx >= 0 else 0)
        self.onSourceChanged()

    def source(self):
        return self.sourceCombo.itemData(self.sourceCombo.currentIndex) or "volbrain"

    def onSourceChanged(self, index=None):
        syn = self.source() == "synthseg"
        for w in (self.tissuesLabel, self.tissuesSelector, self.structuresLabel, self.structuresSelector):
            w.setVisible(not syn)
        for w in (self.synthsegLabel, self.synthsegSelector, self.sourceHint):
            w.setVisible(syn)
        qt.QSettings().setValue("CisternaMagnaVolumetry/source", self.source())

    def _makeVolumeSelector(self, tooltip, nodeTypes=None):
        selector = slicer.qMRMLNodeComboBox()
        selector.nodeTypes = nodeTypes or ["vtkMRMLScalarVolumeNode", "vtkMRMLLabelMapVolumeNode"]
        selector.selectNodeUponCreation = False
        selector.addEnabled = False
        selector.removeEnabled = False
        selector.noneEnabled = False
        selector.showChildNodeTypes = False
        selector.setMRMLScene(slicer.mrmlScene)
        selector.setToolTip(tooltip)
        return selector

    def onApplyButton(self):
        t1Node = self.t1Selector.currentNode()
        syn = self.source() == "synthseg"
        if syn:
            synNode = self.synthsegSelector.currentNode()
            if not (t1Node and synNode):
                slicer.util.errorDisplay("Please select the T1 volume and the SynthSeg result.")
                return
        else:
            tissuesNode = self.tissuesSelector.currentNode()
            structuresNode = self.structuresSelector.currentNode()
            if not (t1Node and tissuesNode and structuresNode):
                slicer.util.errorDisplay("Please select T1, native_tissues, and native_structures volumes.")
                return

        self.statusLabel.text = "Calculating..."
        self.applyButton.enabled = False
        slicer.app.processEvents()

        try:
            if self.logic is None:
                self.logic = CisternaMagnaVolumetryLogic()
            if syn:
                volumes = self.logic.runSynthSeg(t1Node, synNode)
            else:
                volumes = self.logic.run(t1Node, tissuesNode, structuresNode)
        except Exception as e:
            logging.exception("CisternaMagnaVolumetry computation failed")
            slicer.util.errorDisplay(f"Error: {e}")
            self.statusLabel.text = "An error occurred - see the Python console for details."
        else:
            self._populateResultsTable(volumes)
            self.statusLabel.text = "Done. Segmentation 'Cisterna volumes' added to the scene."
        finally:
            self.applyButton.enabled = True

    def _populateResultsTable(self, volumes):
        extra = volumes.get("_extra", {})
        extraRows = [("Lateral ventricles", "lateral_ventricles_mm3"), ("3rd ventricle", "third_ventricle_mm3"),
                     ("4th ventricle", "fourth_ventricle_mm3"), ("Extra-ventricular CSF (total)", "extraventricular_csf_mm3"),
                     ("Intracranial volume", "intracranial_mm3")]
        extraRows = [r for r in extraRows if r[1] in extra]
        self.resultsTable.setRowCount(len(CISTERN_ORDER) + len(extraRows))
        for row, name in enumerate(CISTERN_ORDER):
            v = volumes[name]
            self.resultsTable.setItem(row, 0, qt.QTableWidgetItem(DISPLAY_NAMES[name]))
            self.resultsTable.setItem(row, 1, qt.QTableWidgetItem(str(v["voxels"])))
            self.resultsTable.setItem(row, 2, qt.QTableWidgetItem(f"{v['volume_mm3']:.1f}"))
            self.resultsTable.setItem(row, 3, qt.QTableWidgetItem(f"{v['volume_cm3']:.3f}"))
        for i, (label, key) in enumerate(extraRows):
            row = len(CISTERN_ORDER) + i
            self.resultsTable.setItem(row, 0, qt.QTableWidgetItem(label))
            self.resultsTable.setItem(row, 1, qt.QTableWidgetItem(""))
            self.resultsTable.setItem(row, 2, qt.QTableWidgetItem(f"{extra[key]:.1f}"))
            self.resultsTable.setItem(row, 3, qt.QTableWidgetItem(f"{extra[key] / 1000:.3f}"))


class CisternaMagnaVolumetryLogic(ScriptedLoadableModuleLogic):
    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self._ensureDependencies()

    def _ensureDependencies(self):
        try:
            import nibabel  # noqa: F401
            import scipy  # noqa: F401
        except ImportError:
            proceed = slicer.util.confirmOkCancelDisplay(
                "This module requires the 'nibabel' and 'scipy' packages. Install them now?\n"
                "(You may need to restart Slicer after installation.)"
            )
            if proceed:
                slicer.util.pip_install("nibabel scipy")
            else:
                raise RuntimeError("nibabel/scipy are not installed, cannot compute.")

        moduleDir = os.path.dirname(__file__)
        if moduleDir not in sys.path:
            sys.path.insert(0, moduleDir)

    def run(self, t1Node, tissuesNode, structuresNode):
        from CisternaMagnaVolumetryLib.cisterna_core import run_cistern_analysis, CISTERN_COLORS

        t1_array = slicer.util.arrayFromVolume(t1Node)
        tissues_array = slicer.util.arrayFromVolume(tissuesNode)
        structures_array = slicer.util.arrayFromVolume(structuresNode)

        if tissues_array.shape != t1_array.shape or structures_array.shape != t1_array.shape:
            raise ValueError(
                "T1, native_tissues, and native_structures volumes must have the same shape "
                f"(T1={t1_array.shape}, tissues={tissues_array.shape}, structures={structures_array.shape}). "
                "Make sure you selected files from the same VolBrain job."
            )

        affine = self._arrayAffineFromVolumeNode(t1Node)
        spacing = t1Node.GetSpacing()  # in (i, j, k) order (VTK/Slicer convention)
        voxel_sizes = (spacing[2], spacing[1], spacing[0])  # convert to array (k,j,i) order

        cistern_masks, external_csf, volumes, total_csf = run_cistern_analysis(
            t1_array.astype("float64"),
            tissues_array.astype("int16"),
            structures_array.astype("int16"),
            affine,
            voxel_sizes,
        )

        self._createSegmentation(t1Node, cistern_masks, CISTERN_COLORS)

        volumes = dict(volumes)
        volumes["total_external_csf"] = total_csf
        return volumes

    # SynthSeg segment names -> FreeSurfer label values (names as written by SynthSeg 2.0, by the SynthSeg 1.0
    # table of SlicerSynthSeg, and by FreeSurferColorLUT; compared after normalisation)
    _SYNTHSEG_NAME_TO_LABEL = {
        "brain stem": 16, "4th ventricle": 15, "3rd ventricle": 14, "csf": 24,
        "left ventral dc": 28, "right ventral dc": 60, "left ventraldc": 28, "right ventraldc": 60,
        "left cerebellum white matter": 7, "left cerebellum cortex": 8,
        "right cerebellum white matter": 46, "right cerebellum cortex": 47,
        "left hippocampus": 17, "right hippocampus": 53, "left amygdala": 18, "right amygdala": 54,
        "left lateral ventricle": 4, "right lateral ventricle": 43,
        "left inferior lateral ventricle": 5, "right inferior lateral ventricle": 44,
        "left inf lat vent": 5, "right inf lat vent": 44,
        "ctx lh entorhinal": 1006, "ctx rh entorhinal": 2006,
        "ctx lh parahippocampal": 1016, "ctx rh parahippocampal": 2016,
    }

    @staticmethod
    def _norm(name):
        return " ".join(name.lower().replace("-", " ").replace("_", " ").split())

    def _synthsegArray(self, synNode, refNode):
        """SynthSeg label array on the grid of refNode. A segmentation node (from SlicerSynthSeg) is converted
        back to FreeSurfer label values via the segment names; other segments get a neutral non-zero value."""
        if synNode.IsA("vtkMRMLSegmentationNode"):
            seg = synNode.GetSegmentation()
            arr = None
            found = set()
            for i in range(seg.GetNumberOfSegments()):
                sid = seg.GetNthSegmentID(i)
                name = self._norm(seg.GetSegment(sid).GetName())
                value = self._SYNTHSEG_NAME_TO_LABEL.get(name, 99)
                m = slicer.util.arrayFromSegmentBinaryLabelmap(synNode, sid, refNode)
                if arr is None:
                    arr = np.zeros(m.shape, dtype=np.int32)
                arr[m > 0] = value
                if value != 99:
                    found.add(value)
            if 24 not in found or 16 not in found:
                raise ValueError("The segmentation has no 'csf' or 'brain-stem' segment. Run SynthSeg 2.0 "
                                 "(not 1.0) with cortical parcellation.")
            return arr
        return slicer.util.arrayFromVolume(synNode).astype(np.int32)

    def runSynthSeg(self, t1Node, synNode):
        from CisternaMagnaVolumetryLib.cisterna_core import run_cistern_analysis_synthseg, CISTERN_COLORS

        if synNode.IsA("vtkMRMLSegmentationNode"):
            refNode = t1Node
        else:
            refNode = synNode          # SynthSeg writes a 1 mm grid; the T1 is resampled onto it if needed
        seg = self._synthsegArray(synNode, refNode)
        if 24 not in np.unique(seg):
            raise ValueError("No CSF label (24) in the SynthSeg result: run SynthSeg 2.0 with cortical "
                             "parcellation (SynthSeg 1.0 has no extra-ventricular CSF label).")
        if not np.isin(seg, (1006, 1016, 2006, 2016)).any():
            logging.warning("SynthSeg result has no cortical parcellation; only hippocampus and amygdala are "
                            "used as the lateral limit of the crural/ambient cisterns.")
        t1 = slicer.util.arrayFromVolume(t1Node)
        resampled = None
        if t1.shape != seg.shape:
            resampled = slicer.modules.volumes.logic().ResampleVolumeToReferenceVolume(t1Node, refNode)
            t1 = slicer.util.arrayFromVolume(resampled)
        affine = self._arrayAffineFromVolumeNode(refNode)
        spacing = refNode.GetSpacing()
        voxel_sizes = (spacing[2], spacing[1], spacing[0])
        masks, external_csf, volumes, extra = run_cistern_analysis_synthseg(seg, affine, voxel_sizes, t1)
        if resampled is not None:
            slicer.mrmlScene.RemoveNode(resampled)
        self._createSegmentation(refNode, masks, CISTERN_COLORS)
        slicer.util.setSliceViewerLayers(background=t1Node)
        volumes = dict(volumes)
        volumes["_extra"] = extra
        return volumes

    @staticmethod
    def _arrayAffineFromVolumeNode(volumeNode):
        """slicer.util.arrayFromVolume() returns a numpy array indexed in
        (k,j,i) order (the row-major numpy equivalent of VTK/ITK's
        fastest-varying-I convention). This function returns the 4x4
        (nibabel-style) affine that maps that array axis order to RAS."""
        ijkToRas = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(ijkToRas)
        m = np.array([[ijkToRas.GetElement(r, c) for c in range(4)] for r in range(4)])
        # array axes (a0=k, a1=j, a2=i) -> ijk (i,j,k): P encodes this permutation
        P = np.array(
            [
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        return m @ P

    def _createSegmentation(self, referenceVolumeNode, cistern_masks, colors):
        segmentationNode = slicer.mrmlScene.GetFirstNodeByName("Cisterna volumes")
        if segmentationNode is None:
            segmentationNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "Cisterna volumes")
        else:
            segmentationNode.GetSegmentation().RemoveAllSegments()
        segmentationNode.CreateDefaultDisplayNodes()
        segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(referenceVolumeNode)

        for name in CISTERN_ORDER:
            mask = cistern_masks[name]
            segmentId = segmentationNode.GetSegmentation().AddEmptySegment(
                name, DISPLAY_NAMES[name], colors[name]
            )
            slicer.util.updateSegmentBinaryLabelmapFromArray(
                mask.astype("uint8"), segmentationNode, segmentId, referenceVolumeNode
            )

        segmentationNode.GetDisplayNode().SetVisibility(True)
        slicer.util.setSliceViewerLayers(background=referenceVolumeNode)


class CisternaMagnaVolumetryTest(ScriptedLoadableModuleTest):
    """Automated test skeleton. Since a real test requires actual VolBrain
    sample data (patient data, not shared in this repo), this only verifies
    that the module's dependencies import cleanly. Real validation lives in
    the CLI repo's docs/validation/ (two cases)."""

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_ModuleImports()

    def test_ModuleImports(self):
        self.delayDisplay("Testing cisterna_core import")
        moduleDir = os.path.dirname(__file__)
        if moduleDir not in sys.path:
            sys.path.insert(0, moduleDir)
        from CisternaMagnaVolumetryLib import cisterna_core  # noqa: F401

        self.assertTrue(hasattr(cisterna_core, "run_cistern_analysis"))
        self.delayDisplay("Test passed")
