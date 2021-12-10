# $Header$

# AUTOMAKE_OTPIONS =
# ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

bin_PROGRAMS = xephem

xephem_SOURCES =  \
	aavso.c         \
	annotmenu.c     \
	broadcast.c     \
	calmenu.c       \
	closemenu.c     \
	compiler.c      \
	coordsmenu.c    \
	datamenu.c      \
	db.h db.c            \
	dbmenu.c        \
	dm.h \
	earthmap.c      \
	earthmenu.c     \
	fallbacks.c     \
	favmenu.c       \
	formats.c       \
	fsmenu.c        \
	gallerymenu.c   \
	glance.c        \
	gsc.c           \
	gscnet.c        \
	helpmenu.c      \
	homeio.c        \
	hznmenu.c       \
	indiapi.h indimenu.c      \
	imregmenu.c     \
	jpeg2pm.c       \
	jupmenu.c       \
	listmenu.c      \
	mainmenu.c      \
	map.h \
	marsmenu.c      \
	marsmmenu.c     \
	moonmenu.c      \
	moviemenu.c     \
	msgmenu.c       \
	net.h       \
	netmenu.c       \
	objmenu.c       \
	obslog.c        \
	patchlevel.h patchlevel.c    \
	plot.h plot_aux.c      \
	plotmenu.c      \
	preferences.c   \
	progress.c      \
	ps.h ps.c            \
	query.c         \
	rotated.h rotated.c       \
	satmenu.c       \
	saveres.c       \
	scope.c         \
	sites.h sites.c         \
	skybinary.c     \
	skyeyep.h skyeyep.c       \
	skyfifos.c      \
	skyfiltmenu.c   \
	skyfits.c       \
	skyhist.h skyhist.c       \
	skyip.h skyip.c         \
	skylist.h skylist.c       \
	skytoolbar.h skytoolbar.c    \
	skyviewmenu.c   \
	solsysmenu.c    \
	splash.c        \
	srchmenu.c      \
	sunmenu.c       \
	time.c          \
	tips.c          \
	trails.h trailmenu.c     \
	tools/lx200xed/liblx200.h \
	uranusmenu.c    \
	ucac.c          \
	usno.c          \
	versionmenu.c   \
	webdbmenu.c     \
	xe2.c           \
	xe3.c           \
	xephem.h xephem.c        \
	xmisc.c


xephem_CFLAGS = -D_GNU_SOURCE -I. -I../../libastro -I../../libip -I../../liblilxml -I/usr/include/libpng16
xephem_CPPFLAGS = -D_GNU_SOURCE -I. -I../../libastro -I../../libip -I../../liblilxml -I/usr/include/libpng16
xephem_LDFLAGS =  -L../../libastro -L../../libip -L../../liblilxml
xephem_LDADD = -lastro -lip -llilxml

# AM_CFLAGS = -I. -Iinclude

# AM_LIBTOOLFLAGS=--silent

dist_man_MANS = xephem.1

# Skipping the tools/xephemdb tools/xedb tools/simpleIndi which will lead to install problems with debuild
nobase_dist_pkgdata_DATA = xephem.man $(wildcard auxil/*) $(wildcard gallery/*) \
	$(wildcard fits/*) $(wildcard help/*.html) $(wildcard help/png/*) $(wildcard lo/*) \
	$(wildcard tools/*.pl) $(wildcard tools/*.awk) $(wildcard tools/*.c) $(wildcard tools/lx200xed/*) \
	$(wildcard catalogs/*) $(wildcard fifos/[MR]*)

BUILT_SOURCES = xephem.1

# DIST_SUBDIRS = auxil fifos gallery fits help lo tools catalogs
SUBDIRS = fifos

xephem.1: xephem.man
	nroff -man $? > $@
