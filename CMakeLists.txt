cmake_minimum_required(VERSION 3.16)
project(EdgesCGALWeightedDelaunayND LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

add_executable(EdgesCGALWeightedDelaunayND src/EdgesCGALWeightedDelaunayND.cpp)

# Link CGAL; for Eigen, some distros don't export Eigen3::Eigen target.
if(TARGET Eigen3::Eigen)
  target_link_libraries(EdgesCGALWeightedDelaunayND PRIVATE CGAL::CGAL Eigen3::Eigen)
else()
  target_link_libraries(EdgesCGALWeightedDelaunayND PRIVATE CGAL::CGAL)
  if(EIGEN3_INCLUDE_DIR)
    target_include_directories(EdgesCGALWeightedDelaunayND PRIVATE ${EIGEN3_INCLUDE_DIR})
  endif()
endif()

target_compile_options(EdgesCGALWeightedDelaunayND PRIVATE -O3 -DNDEBUG)