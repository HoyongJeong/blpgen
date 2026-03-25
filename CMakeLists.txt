################################################################################
###
###   CMakeLists.txt for blpgen.
###
### - Authors: Hoyong   Jeong   (hoyong5419@korea.ac.kr)
###            Antonino Cannavo (  acannavo@bnl.gov    )
###
################################################################################



###-----------------------------------------------------------------------------
### Setup the project
###-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(BLPGEN)



###-----------------------------------------------------------------------------
### Setup include directory for this project
###-----------------------------------------------------------------------------
include_directories(${PROJECT_SOURCE_DIR}/include)



###-----------------------------------------------------------------------------
### Find ROOT package
###-----------------------------------------------------------------------------
find_package(ROOT REQUIRED EG)



###-----------------------------------------------------------------------------
### Setup ROOT include directories and compile definitions.
###-----------------------------------------------------------------------------
include(${ROOT_USE_FILE})



###-----------------------------------------------------------------------------
### Setup PLUTO using user's environment variable "PLUTOSYS"
### You can edit directory path depending on your system.
###-----------------------------------------------------------------------------
link_directories($ENV{PLUTOSYS}/lib)
include_directories($ENV{PLUTOSYS}/include)



###-----------------------------------------------------------------------------
### Locate sources and headers for this project
### NB: headers are included so they will show up in IDEs
###-----------------------------------------------------------------------------
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)



###-----------------------------------------------------------------------------
### Add the executable, and link it to the libraries
###-----------------------------------------------------------------------------
add_executable(blpgen main.cc ${sources} ${headers})
target_link_libraries(blpgen ${ROOT_LIBRARIES} Pluto)



###-----------------------------------------------------------------------------
### For internal Geant4 use - but has no effect if you build this blpgen
### standalone
###-----------------------------------------------------------------------------
add_custom_target(BLPGEN DEPENDS blpgen)
