# $Header$

# AUTOMAKE_OTPIONS =
# ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

SUBDIRS = libastro libip liblilxml GUI/xephem

# make perl files executable
all-local:
	find . \( -name "*.pl" -o -name "*.awk" \) -exec chmod +x {} \;

if DEBIAN
install-data-local:
	$(MKDIR_P) -m 755 $(DESTDIR)/etc/X11/app-defaults
	echo "XEphem.ShareDir: ${datarootdir}/$(PACKAGE_NAME)" > $(DESTDIR)/etc/X11/app-defaults/XEphem
else
install-data-local:
	$(MKDIR_P) -m 755 $(DESTDIR)/${datarootdir}/X11/app-defaults
	echo "XEphem.ShareDir: ${datarootdir}/$(PACKAGE_NAME)" > $(DESTDIR)/${datarootdir}/X11/app-defaults/XEphem
endif
