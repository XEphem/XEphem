cmake_minimum_required(VERSION 3.10)

project(libip C)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(SRC_FILES 
  explodegif.c
  fits.c
  fsmatch.c
  gaussfit.c
  lstsqr.c
  median.c
  sqr.c
  stars.c
  stats.c
  walk.c
  wcs.c)

add_library(ip STATIC ${SRC_FILES})
target_link_libraries(ip astro)
target_include_directories(ip PUBLIC
                           "${CMAKE_CURRENT_BINARY_DIR}"
                           "${CMAKE_CURRENT_SOURCE_DIR}")
