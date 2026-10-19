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
}
CISTERN_ORDER = ["cisterna_magna", "cisterna_ambiens", "cisterna_interpeduncularis"]


class CisternaMagnaVolumetry(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Cisterna Magna Volumetry"
        self.parent.categories = ["Segmentation"]
        self.parent.dependencies = []
        self.parent.contributors = ["Niyazi (VolBrain/7T cisterna volumetry project)"]
        self.parent.helpText = (
            "Computes cisterna magna, cisterna ambiens, and cisterna interpeduncularis "
            "volumes from VolBrain native_tissues/native_structures output. The method is "
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
        parametersCollapsibleButton.text = "Input (VolBrain volumes)"
        self.layout.addWidget(parametersCollapsibleButton)
        parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

        self.t1Selector = self._makeVolumeSelector(
            "T1 volume that was the input to VolBrain (used as reference geometry)."
        )
        parametersFormLayout.addRow("T1 volume:", self.t1Selector)

        self.tissuesSelector = self._makeVolumeSelector("VolBrain native_tissues_job*.nii.gz")
        parametersFormLayout.addRow("native_tissues:", self.tissuesSelector)

        self.structuresSelector = self._makeVolumeSelector("VolBrain native_structures_job*.nii.gz")
        parametersFormLayout.addRow("native_structures:", self.structuresSelector)

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
            "Note: landmark-based method (v2) - VolBrain's atlas has no cistern-level "
            "labels. Compared with manual masks on 2 cases: interpeduncular Dice 0.51, "
            "cisterna magna Dice 0.60. Main limitation: VolBrain's brain mask crops the "
            "image at the level of the foramen magnum, so part of the cisterna magna is "
            "never in the CSF mask. Check every result visually."
        )
        noteLabel.setWordWrap(True)
        noteLabel.setStyleSheet("color: gray; font-size: 11px;")
        resultsLayout.addWidget(noteLabel)

        self.layout.addStretch(1)

        self.applyButton.connect("clicked(bool)", self.onApplyButton)

    def _makeVolumeSelector(self, tooltip):
        selector = slicer.qMRMLNodeComboBox()
        selector.nodeTypes = ["vtkMRMLScalarVolumeNode", "vtkMRMLLabelMapVolumeNode"]
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
        self.resultsTable.setRowCount(len(CISTERN_ORDER))
        for row, name in enumerate(CISTERN_ORDER):
            v = volumes[name]
            self.resultsTable.setItem(row, 0, qt.QTableWidgetItem(DISPLAY_NAMES[name]))
            self.resultsTable.setItem(row, 1, qt.QTableWidgetItem(str(v["voxels"])))
            self.resultsTable.setItem(row, 2, qt.QTableWidgetItem(f"{v['volume_mm3']:.1f}"))
            self.resultsTable.setItem(row, 3, qt.QTableWidgetItem(f"{v['volume_cm3']:.3f}"))


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
