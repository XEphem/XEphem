cmake_minimum_required(VERSION 3.10)

project(xephem VERSION 4.1.0 LANGUAGES C)
set(DISPLAY_NAME "XEphem")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Required packages
find_package(Motif REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(X11 REQUIRED)

# Optional packages
find_package(JPEG)
find_package(PNG)
find_package(ZLIB)

# Define installation directories
if((NOT WIN32) AND (NOT NATIVE_OSX_APP))
  include(GNUInstallDirs)
endif()
if(NATIVE_OSX_APP)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.app/Contents")
  set(CMAKE_INSTALL_FULL_LOCALEDIR "${CMAKE_INSTALL_PREFIX}/Resources/locale")
  set(CMAKE_INSTALL_BINDIR "MacOS")
  set(XEPHEM_DATADIR "${CMAKE_INSTALL_PREFIX}/Resources")
elseif(WIN32)
  set(CMAKE_INSTALL_BINDIR ".")
  set(CMAKE_INSTALL_DATAROOTDIR ".")
  set(CMAKE_INSTALL_DATADIR ".")
  set(CMAKE_INSTALL_FULL_LOCALEDIR "locale")
  set(CMAKE_INSTALL_LOCALEDIR "${CMAKE_INSTALL_FULL_LOCALEDIR}")
  set(XEPHEM_DATADIR "${CMAKE_INSTALL_DATADIR}")
else()
  set(XEPHEM_DATADIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
endif()

# Build libraries
add_subdirectory(libastro)
add_subdirectory(libip)
if(NOT JPEG_FOUND)
  add_subdirectory(libjpegd)
endif()
add_subdirectory(liblilxml)
if(NOT PNG_FOUND)
  add_subdirectory(libpng)
endif()
if(NOT ZLIB_FOUND)
  add_subdirectory(libz)
endif()

# Build main program
add_subdirectory(GUI/xephem)

# Install
install(TARGETS xephem
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES ${CMAKE_BINARY_DIR}/xephem.1
        DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
install(DIRECTORY GUI/xephem/auxil
        DESTINATION ${XEPHEM_DATADIR})
install(DIRECTORY GUI/xephem/catalogs
        DESTINATION ${XEPHEM_DATADIR})
install(DIRECTORY GUI/xephem/fifos
        DESTINATION ${XEPHEM_DATADIR})
install(DIRECTORY GUI/xephem/fits
        DESTINATION ${XEPHEM_DATADIR})
install(DIRECTORY GUI/xephem/gallery
        DESTINATION ${XEPHEM_DATADIR})
install(DIRECTORY GUI/xephem/help
        DESTINATION ${XEPHEM_DATADIR})
install(DIRECTORY GUI/xephem/lo
        DESTINATION ${XEPHEM_DATADIR})
