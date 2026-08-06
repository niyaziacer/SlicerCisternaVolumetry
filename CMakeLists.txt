cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(CisternaMagnaVolumetry)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/niyaziacer/CisternaMagnaVolumetry")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Niyazi")
set(EXTENSION_DESCRIPTION "Computes cisterna magna, cisterna ambiens, and cisterna interpeduncularis volumes from VolBrain native_tissues/native_structures output and adds them as a segmentation to the scene.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/niyaziacer/CisternaMagnaVolumetry/main/CisternaMagnaVolumetry.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/niyaziacer/CisternaMagnaVolumetry/main/screenshot.jpg")
set(EXTENSION_DEPENDS "NA")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(CisternaMagnaVolumetry)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
