# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(dwptst)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Build executable------------------------------------------------------------
add_executable(dwptst dwptst.C)
target_link_libraries(dwptst Garfield::Garfield)

add_executable(dwp_np_2e13 dwp_np_2e13.C)
target_link_libraries(dwp_np_2e13 Garfield::Garfield)

add_executable(dwp_np_6e13 dwp_np_6e13.C)
target_link_libraries(dwp_np_6e13 Garfield::Garfield)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(dwptst OpenMP::OpenMP_CXX)
endif()

# ---Copy all files locally to the build directory-------------------------------

