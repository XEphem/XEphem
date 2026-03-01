# Sample Makefile for xephem v4.1.0

# Sample compile and link flags for a few systems. Find your system or similar
# and uncomment as required. The default is for linux or Mac OS X systems.
# The only real challenge usually is to set the MOTIFI and MOTIFL macros.

# The basic idea to build XEphem is to first go compile the libraries in
# ../../lib* then come back here and compile and link all the *.c files into
# one executable, xephem.

# These -I and -L flags point to the supporting XEphem libraries
LIBINC = -I../../libastro -I../../libip -I../../liblilxml -I../../libjpegd -I../../libpng -I../../libz
LIBLNK = -L../../libastro -L../../libip -L../../liblilxml -L../../libjpegd -L../../libpng -L../../libz
LIBLIB = -lastro -lip -llilxml -ljpegd -lpng -lz

# MOTIFI is the directory containing the Xm directory of include files.
# MOTIFL is the directory containing the libXm.a library.
# MOTIF can be used to set one common directory if appropriate.
# On 64 bit systems, try MOTIFI = /usr/include and MOTIFL = /usr/lib64
ifdef MOTIF
MOTIFI = $(MOTIF)
MOTIFL = $(MOTIF)
endif
ifndef MOTIFI
MOTIFI = /usr/include
endif
ifndef MOTIFL
MOTIFL = /usr/lib
endif

# For Linux and Apple OS X
# /opt/homebrew is for modern Apple Silicon; the other paths, for older macOS.

UNAME_S = $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  PLATI = -I/opt/homebrew/include \
          -I/usr/local/opt/openssl@3/include \
          -I/opt/X11/include
  PLATL = -L/opt/homebrew/lib \
          -L/usr/local/opt/openssl@3/lib \
          -L/opt/X11/lib
endif

CC = gcc
CLDFLAGS = -g
CFLAGS = $(LIBINC) $(CLDFLAGS) -O2 -Wall -I$(MOTIFI) $(PLATI)
LDFLAGS = $(LIBLNK) $(CLDFLAGS) -L$(MOTIFL) $(PLATL)
XLIBS = -lXm -lXt -lXext -lXmu -lX11
LIBS = $(XLIBS) $(LIBLIB) -lm -lssl

# Variables installation for Linux
prefix := /usr/local
bindir = $(prefix)/bin
datadir = $(prefix)/share
mandir = $(datadir)/man
export datadir

# static linking on Apple using X11 libs from ports
# CC = gcc
# OLI = /opt/local/include
# OLL = /opt/local/lib
# CLDFLAGS = -g
# CFLAGS = $(LIBINC) $(CLDFLAGS) -O2 -Wall -I$(MOTIFI) -I$(OLI)
# LDFLAGS = $(LIBLNK) $(CLDFLAGS) -L$(MOTIFL) -L$(OLL)
# XLIBS = -lXm $(OLL)/libXp.a $(OLL)/libXt.a $(OLL)/libSM.a $(OLL)/libICE.a $(OLL)/libXext.a $(OLL)/libXmu.a $(OLL)/libX11.a $(OLL)/libxcb.a $(OLL)/libXau.a $(OLL)/libXdmcp.a
# LIBS = $(XLIBS) $(LIBLIB) -lm -lssl

# for ppc Apple OS X to make universal (i386 and ppc binary), requires
# Xcode 2.2, /Developers/SDK/MacOSX10.4u.sdk and all libs, including libXm,
# built similarly
# CLDFLAGS = -arch i386 -arch ppc
# CFLAGS = $(LIBINC) $(CLDFLAGS) -O2 -Wall -isysroot /Developer/SDKs/MacOSX10.4u.sdk -I$(MOTIFI) -I/usr/X11R6/include
# LDFLAGS = $(LIBLNK) $(CLDFLAGS) -Wl,-syslibroot,/Developer/SDKs/MacOSX10.4u.sdk -L$(MOTIFL) -L/usr/X11R6/lib

# MKS Enterprise toolkit
# MOTIFI = -I/usr/include
# MOTIFL = -L/usr/lib
# CC = cc
# CLDFLAGS =
# CFLAGS = $(LIBINC) $(CLDFLAGS) -O -I$(MOTIFI)
# LDFLAGS = $(LIBLNK) $(CLDFLAGS) -L$(MOTIFL) -W/subsystem:windows -s
# XLIBS = -lXm -lXt -lX11
# LIBS = $(XLIBS) $(LIBLIB) -lm

# SVR4 derivatives:
# MOTIFI = -I/usr/include
# MOTIFL = -L/usr/lib
# CC = cc
# CLDFLAGS = 
# CFLAGS = $(LIBINC) $(CLDFLAGS) -O -I$(MOTIFI)
# LDFLAGS = $(LIBLNK) $(CLDFLAGS) -L$(MOTIFL) 
# XLIBS = -lXm -lXt -lX11
# LIBS = $(XLIBS) $(LIBLIB) -lsocket -lnsl -lc -lm /usr/ucblib/libucb.a
# Note: if you get regex undefined, add -lgen to the end of LIBS.

# Solaris:
# Motif stuff might also be in /usr/dt/share/{include,lib}.
# be sure /opt/SUNWspro/bin and /usr/ccs/bin are in your PATH ahead of /usr/ucb.
# MOTIFI = -I/usr/dt/include
# MOTIFL = -L/usr/dt/lib -R/usr/dt/lib
# CC = cc
# CLDFLAGS =
# CFLAGS = $(LIBINC) $(CLDFLAGS) -O -I$(MOTIFI) -I/usr/openwin/include
# LDFLAGS = $(LIBLNK) $(CLDFLAGS) -L$(MOTIFL) -L/usr/openwin/lib
# XLIBS = -lXm -lXt -lX11
# LIBS = $(XLIBS) $(LIBLIB) -lm -lsocket -lnsl

# HP-UX
# MOTIFI = -I/usr/include
# MOTIFL = -L/usr/lib
# CC = cc
# CLDFLAGS =
# CFLAGS = $(LIBINC) $(CLDFLAGS) -Aa -fast -I$(MOTIFI)
# LDFLAGS = $(LIBLNK) $(CLDFLAGS) -L$(MOTIFL)

# AIX
# MOTIFI = -I/usr/include
# MOTIFL = -L/usr/lib
# CC = xlc
# CLDFLAGS =
# CFLAGS = $(LIBINC) $(CLDFLAGS) -O2 -qlanglvl=ansi -qarch=com -qmaxmem=16384 -I$(MOTIFI)
# LDFLAGS = $(LIBLNK) $(CLDFLAGS) -L$(MOTIFL)

.PRECIOUS:	xephem

INCS =	db.h dm.h indiapi.h map.h net.h patchlevel.h plot.h ps.h \
	rotated.h sites.h skyeyep.h skyhist.h skyip.h skylist.h skytoolbar.h \
	trails.h xephem.h

OBJS =			\
	aavso.o		\
	annotmenu.o	\
	broadcast.o	\
	calmenu.o	\
	closemenu.o	\
	compiler.o	\
	coordsmenu.o	\
	datamenu.o	\
	db.o		\
	dbmenu.o	\
	earthmap.o	\
	earthmenu.o	\
	fallbacks.o	\
	favmenu.o	\
	formats.o	\
	fsmenu.o	\
	gallerymenu.o	\
	glance.o	\
	gsc.o		\
	gscnet.o	\
	helpmenu.o	\
	homeio.o	\
	hznmenu.o	\
	indimenu.o	\
	imregmenu.o	\
	jpeg2pm.o	\
	jupmenu.o	\
	listmenu.o	\
	mainmenu.o	\
	marsmenu.o	\
	marsmmenu.o	\
	moonmenu.o	\
	moviemenu.o	\
	msgmenu.o	\
	netmenu.o	\
	objmenu.o	\
	obslog.o	\
	patchlevel.o	\
	plot_aux.o	\
	plotmenu.o	\
	preferences.o	\
	progress.o	\
	ps.o		\
	query.o		\
	rotated.o	\
	satmenu.o	\
	saveres.o	\
	scope.o		\
	sites.o		\
	skybinary.o	\
	skyeyep.o	\
	skyfifos.o	\
	skyfiltmenu.o	\
	skyfits.o	\
	skyhist.o	\
	skyip.o		\
	skylist.o	\
	skytoolbar.o	\
	skyviewmenu.o	\
	solsysmenu.o	\
	splash.o	\
	srchmenu.o	\
	sunmenu.o	\
	time.o		\
	tips.o		\
	trailmenu.o	\
	uranusmenu.o	\
	ucac.o		\
	usno.o		\
	versionmenu.o	\
	webdbmenu.o	\
	xe2.o		\
	xe3.o		\
	xephem.o	\
	xmisc.o

all: xephem xephem.1

xephem: libs $(INCS) $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

xephem.1: xephem.man
	cp $? $@

libs:
	$(MAKE) -C ../../libastro
	$(MAKE) -C ../../libip
	$(MAKE) -C ../../libjpegd
	$(MAKE) -C ../../liblilxml
	$(MAKE) -C ../../libpng
	$(MAKE) -C ../../libz

clean:
	rm -fr *.o ../../lib*/*.[ao]

clobber: clean
	rm -fr xephem xephem.1

# Installation
install: install-binary install-man install-data install-desktopfiles

install-binary:
	install -Dm 755 xephem -t $(DESTDIR)$(bindir)

install-man:
	install -Dm 644 xephem.1 -t $(DESTDIR)$(mandir)/man1

install-data:
	install -dm 755 $(DESTDIR)$(datadir)/xephem
	@for d in auxil catalogs fifos fits gallery help lo; do \
	  echo "cp -R $$d $(DESTDIR)$(datadir)/xephem"; \
	  cp -R $$d $(DESTDIR)$(datadir)/xephem; \
	done

install-desktopfiles:
	$(MAKE) -C ../../desktopfiles DESTDIR=$(DESTDIR) install
