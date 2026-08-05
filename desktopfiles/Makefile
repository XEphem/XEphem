# Makefile for manage desktop files
# Requires utility 'convert' of the package ImageMagick.

CONVERT = convert

SRCPNGICON = ../GUI/xephem/XEphem.png
PNGICONS = \
	32x32/apps/xephem.png \
	48x48/apps/xephem.png \
	64x64/apps/xephem.png \
	128x128/apps/xephem.png
DESKTOPFILE = xephem.desktop
METAINFOFILE = xephem.metainfo.xml

iconsdir = $(DESTDIR)$(datadir)/icons/hicolor
applicationsdir = $(DESTDIR)$(datadir)/applications
metainfodir = $(DESTDIR)$(datadir)/metainfo


all:

install: install-icons install-desktopfile install-metainfo

install-icons: $(PNGICONS)

$(PNGICONS): $(SRCPNGICON)
	@d=$(iconsdir)/$(shell dirname $@); \
	s=$(shell echo $@ | sed -e 's/[a-z]*\([0-9]*\).*/\1/'); \
	echo "install -m 755 -d $$d"; \
	install -m 755 -d $$d; \
	echo "convert -geometry $${s}x$${s} -depth 8 -background none $< $(iconsdir)/$@"; \
	$(CONVERT) -geometry $${s}x$${s} -depth 8 -background none $< $(iconsdir)/$@

install-desktopfile:
	install -Dm 644 $(DESKTOPFILE) -t $(applicationsdir)

install-metainfo:
	install -Dm 644 $(METAINFOFILE) -t $(metainfodir)
