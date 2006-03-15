# $Id: Makefile.am,v 1.1 2006/03/15 14:18:20 mtw Exp $

#pkgconfigdir = $(libdir)/pkgconfig

#lib_LIBRARIES = liblattice.a
bin_PROGRAMS = treekin

treekin_SOURCES = main.c       \
		  globals.c    globals.h   \
                  err.c        err.h       \
                  mxmesch.c    matrix.h    \
                  schur.c      \
                  calc.c       calc.h      \
                  mxccm.c      mxccm.h     \
                  barparser.c  barparser.h \
	          exp_matrix.c exp_matrix.h

#pkginclude_HEADERS = lattice.h      \
#                     ll_coor.h      \
#                     ll_hash_util.h \
#                     ll_queue.h     \
#                     ll_SAW_bit.h   \
#                     ll_SAW_utils.h

#pkgconfig_DATA = liblattice.pc

INCLUDES = @GLIB_CFLAGS@

LDADD = @GLIB_LIBS@

EXTRA_DIST = VERSION.txt
# End of file
