# Make the fifos used by xephem, assuming the default configuration.
#
# on non-POSIX systems, you make fifos using mknod as follows:
#    mknod <file> p

# Debian/Ubuntu has a policy not to allow creation of special files in makefiles,
# so we simply don't support them on Debian/Ubuntu (condiditional on DEBIAN)

# Another problem: some RPM builder (Fedora?) call make instal twice, once in BUILDROOT and
#  a second time in BUILD, so the second time the mkfifo will fail because the
# fifos are already there. Workaround is a conditional mkfifo.
if !DEBIAN
install-data-local:
	$(MKDIR_P) -m 755 $(DESTDIR)/${datarootdir}/xephem/fifos
	cd $(DESTDIR)/${datarootdir}/xephem/fifos &&  if [ ! -e xephem_loc_fifo ] ; then mkfifo -m 0644 xephem_loc_fifo xephem_db_fifo xephem_in_fifo xephem_out_fifo ; fi

uninstall-local:
	rm -rf $(DESTDIR)/${datarootdir}/xephem/fifos
endif
