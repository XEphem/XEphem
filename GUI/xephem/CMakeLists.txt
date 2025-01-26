cmake_minimum_required(VERSION 3.10)

set(SRC_FILES 
  aavso.c
  annotmenu.c
  broadcast.c
  calmenu.c
  closemenu.c
  compiler.c
  coordsmenu.c
  datamenu.c
  db.c
  dbmenu.c
  earthmap.c
  earthmenu.c
  fallbacks.c
  favmenu.c
  formats.c
  fsmenu.c
  gallerymenu.c
  glance.c
  gsc.c
  gscnet.c
  helpmenu.c
  homeio.c
  hznmenu.c
  imregmenu.c
  indimenu.c
  jpeg2pm.c
  jupmenu.c
  listmenu.c
  mainmenu.c
  marsmenu.c
  marsmmenu.c
  moonmenu.c
  moviemenu.c
  msgmenu.c
  netmenu.c
  objmenu.c
  obslog.c
  patchlevel.c
  plot_aux.c
  plotmenu.c
  preferences.c
  progress.c
  ps.c
  query.c
  rotated.c
  satmenu.c
  saveres.c
  scope.c
  sites.c
  skybinary.c
  skyeyep.c
  skyfifos.c
  skyfiltmenu.c
  skyfits.c
  skyhist.c
  skyip.c
  skylist.c
  skytoolbar.c
  skyviewmenu.c
  solsysmenu.c
  splash.c
  srchmenu.c
  sunmenu.c
  time.c
  tips.c
  trailmenu.c
  ucac.c
  uranusmenu.c
  usno.c
  versionmenu.c
  webdbmenu.c
  xe2.c
  xe3.c
  xephem.c
  xmisc.c)

set(TARGET_LIB astro ip lilxml m png ssl Xm Xt Xext Xmu X11 z)
if(NOT JPEG_FOUND)
  set(TARGET_LIB ${TARGET_LIB} jpegd)
else()
  set(TARGET_LIB ${TARGET_LIB} jpeg)
endif()

add_executable(xephem ${SRC_FILES})
target_link_libraries(xephem ${TARGET_LIB})

add_custom_command(TARGET xephem
                   POST_BUILD
                   COMMAND nroff -man ${CMAKE_CURRENT_SOURCE_DIR}/xephem.man > xephem.1
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                   COMMENT "Building manpage"
                   VERBATIM)
