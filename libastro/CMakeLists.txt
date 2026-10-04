cmake_minimum_required(VERSION 3.10)

project(libastro C)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(SRC_FILES 
  aa_hadec.c
  aberration.c
  actan.c
  airmass.c
  anomaly.c
  ap_as.c
  atlas.c
  auxil.c
  bdl.c
  chap95.c
  chap95_data.c
  circum.c
  comet.c
  constel.c
  dbfmt.c
  deep.c
  deltat.c
  earthsat.c
  eq_ecl.c
  eq_gal.c
  formats.c
  helio.c
  jupmoon.c
  libration.c
  magdecl.c
  marsmoon.c
  misc.c
  mjd.c
  moon.c
  mooncolong.c
  moonnf.c
  nutation.c
  obliq.c
  parallactic.c
  parallax.c
  plans.c
  plmoon.c
  plshadow.c
  precess.c
  reduce.c
  refract.c
  rings.c
  riset.c
  riset_cir.c
  satmoon.c
  sdp4.c
  sgp4.c
  sphcart.c
  sun.c
  thetag.c
  twobody.c
  umoon.c
  utc_gst.c
  vsop87.c
  vsop87_data.c)

add_library(astro STATIC ${SRC_FILES})
target_include_directories(astro PUBLIC
                           "${CMAKE_CURRENT_BINARY_DIR}"
                           "${CMAKE_CURRENT_SOURCE_DIR}")
