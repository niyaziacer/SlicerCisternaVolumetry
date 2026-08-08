cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(CisternaMagnaVolumetry)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/niyaziacer/SlicerCisternaVolumetry")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Niyazi")
set(EXTENSION_DESCRIPTION "Computes cisterna magna, cisterna ambiens, and cisterna interpeduncularis volumes from VolBrain native_tissues/native_structures output and adds them as a segmentation to the scene.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/niyaziacer/SlicerCisternaVolumetry/main/CisternaMagnaVolumetry.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/niyaziacer/SlicerCisternaVolumetry/main/figure.%20cisternas%20WhatsApp%20Image%202026-08-06%20at%2012.02.02.jpeg")
set(EXTENSION_DEPENDS "NA")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(CisternaMagnaVolumetry)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
