
cmake_minimum_required(VERSION 2.8)

project(pp_WW_fully_leptonic CXX)

# Find MoMEMta
find_package(MoMEMta CONFIG REQUIRED)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable -std=c++11 -O3 -fcx-fortran-rules -fcx-limited-range")


file(GLOB_RECURSE SOURCES "*.cc")


include_directories("include")
include_directories("SubProcesses/P1_Sigma_sm_uux_epvemumvmx")


add_library(me_pp_WW_fully_leptonic SHARED ${SOURCES})
target_link_libraries(me_pp_WW_fully_leptonic momemta::momemta)

