cmake_minimum_required(VERSION 3.10)

project(libz C)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(SRC_FILES 
  adler32.c
  compress.c
  crc32.c
  deflate.c
  inffast.c
  inflate.c
  inftrees.c
  testzlib.c
  trees.c
  uncompr.c
  zutil.c)

add_library(z STATIC ${SRC_FILES})
target_include_directories(z PUBLIC
                           "${CMAKE_CURRENT_BINARY_DIR}"
                           "${CMAKE_CURRENT_SOURCE_DIR}")
