ggo_verbose = $(ggo_verbose_@AM_V@)
ggo_verbose_ = $(ggo_verbose_@AM_DEFAULT_V@)
ggo_verbose_0 = @echo "  GGO      $@";

SUBDIRS = man

SUFFIXES = _cmdline.c _cmdline.h .ggo

bin_PROGRAMS = treekin

treekin_LDADD = $(mpack_LIBS)

treekin_SOURCES = \
        mxccm.cpp mxccm.h \
        globals.cpp globals.h \
        barparser.cpp barparser.h bardata.h\
        calcpp.cpp  calcpp.h \
        calc.cpp calc.h \
        exp_matrix.cpp  exp_matrix.h \
        treekin_cmdline.c treekin_cmdline.h \
        treekinCastableTypes.h main.cpp


GENGETOPT_CMDL = \
        treekin_cmdline.c \
        treekin_cmdline.h

GENGETOPT_FILES = treekin.ggo

EXTRA_DIST = \
        ${GENGETOPT_FILES} \
        ${GENGETOPT_CMDL} \
        VERSION.txt


if TREEKIN_AM_SWITCH_HAS_GENGETOPT

## Below are the rules how to create the test suite executables
# how we make the commandline interface parser
.ggo_cmdline.c:
	$(ggo_verbose)$(GENGETOPT) -i $<

_cmdline.c_cmdline.h :

endif

globals.cpp:  treekin_cmdline.c

# remove all generated man pages
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)
