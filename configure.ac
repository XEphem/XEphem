#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([xephem], [4.1.0], [mathar@mpia.de])
AC_CONFIG_SRCDIR([GUI])
# AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
# AC_CONFIG_HEADERS([config.h])

LT_INIT([disable-shared])

AC_PROG_CC
AC_PROG_CPP


# warning: order in the following libraries list may matter!
AC_CHECK_LIB([xml2], [xmlMalloc])
AC_CHECK_LIB([X11], [XFillArc])
AC_CHECK_LIB([Xmu], [XmuAreaCopy])
AC_CHECK_LIB([Xext], [XextFindDisplay])
AC_CHECK_LIB([Xt], [XtDisplay])
AC_CHECK_LIB([Xm], [XmCreateText])
AC_CHECK_LIB([z], [adler32])
AC_CHECK_LIB([ssl], [SSL_read])
AC_CHECK_LIB([jpeg], [jpeg_huff_decode])
AC_CHECK_LIB([png16], [png_read_png])
AC_CHECK_LIB([png], [png_read_png])
AC_CHECK_LIB([m], [cos])
# the following apparently needed on older HP-UX, perhaps -lxnet too?
AC_CHECK_LIB([crypt], [crypt])
AC_CHECK_LIB([crypto], [RSA_sign])
# the following apparently not needed on must systems but perhaps on MACs?
AC_CHECK_LIB([SM], [SmcOpenConnection])
AC_CHECK_LIB([ICE], [IceOpenConnection])
AC_CHECK_LIB([xcb], [xcb_free_pixmap])
AC_CHECK_LIB([Xau], [XauReadAuth])
AC_CHECK_LIB([Xdmcp], [XdmcpFill])

# AC_CHECK_HEADERS([png.h X11/X.h Xm/Xm.h])

# we need a way to prevent fifos being installed by automake on Debian/Ubuntu
# and to create the app-defaults in a different place of the file system on Debian/Ubuntu
DEBIAN=m4_esyscmd_s(uname -v |fgrep Ubuntu |wc -l)
AM_CONDITIONAL([DEBIAN], [test x$DEBIAN = x1])

# AC_SUBST([DEBIAN])

AC_CONFIG_FILES([ Makefile libastro/Makefile libip/Makefile liblilxml/Makefile GUI/xephem/Makefile GUI/xephem/fifos/Makefile ])
AC_OUTPUT
