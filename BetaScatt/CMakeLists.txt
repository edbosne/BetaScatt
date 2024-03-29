# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(BetaScatt)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup ROOT include directories and compile definitions
# Setup include directory for this project
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(BetaScatt BetaScatt.cc ${sources} ${headers})
target_link_libraries(BetaScatt ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
message(STATUS "Copying files ${destDir}")
set(USER_SCRIPTS
  Ranlux.conf
#  file
  )

foreach(_script ${USER_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

set(destDir ${PROJECT_BINARY_DIR}/Test)
set(srcDir ${PROJECT_SOURCE_DIR}/Test)
message(STATUS "Copying Test directory ${destDir}")
make_directory(${destDir})

file(GLOB testFiles RELATIVE ${srcDir} ${srcDir}/*)
foreach(testFiles ${testFiles})
    set(srcTemplatePath ${srcDir}/${testFiles})
    if(NOT IS_DIRECTORY ${srcTemplatePath})
        configure_file(
                ${srcTemplatePath}
                ${destDir}/${testFiles}
                COPYONLY)
    endif(NOT IS_DIRECTORY ${srcTemplatePath})
endforeach(testFiles)

set(destDir ${PROJECT_BINARY_DIR}/isotope_spectrums)
set(srcDir ${PROJECT_SOURCE_DIR}/isotope_spectrums)
message(STATUS "Copying Test directory ${destDir}")
make_directory(${destDir})

file(GLOB testFiles RELATIVE ${srcDir} ${srcDir}/*)
foreach(testFiles ${testFiles})
    set(srcTemplatePath ${srcDir}/${testFiles})
    if(NOT IS_DIRECTORY ${srcTemplatePath})
        configure_file(
                ${srcTemplatePath}
                ${destDir}/${testFiles}
                COPYONLY)
    endif(NOT IS_DIRECTORY ${srcTemplatePath})
endforeach(testFiles)


#----------------------------------------------------------------------------
# Add definitions for debug and compatibility
#
#add_definitions(-DG4ANALYSIS_USE -DGEOMETRY_DEBUG -DGDEBUG)
add_definitions(-DG4ANALYSIS_USE -DGDEBUG)





