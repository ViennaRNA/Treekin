dist_man_MANS = treekin.1

bin_PROGRAMS = treekin

treekin_SOURCES = main.c\
                  globals.c globals.h \
                  calc.c  calc.h \
                  mxccm.c mxccm.h \
                  barparser.c barparser.h \
                  exp_matrix.c  exp_matrix.h \
                  treekin_cmdline.c treekin_cmdline.h \
                  calcpp.cpp  calcpp.h

EXTRA_DIST = treekin.ggo treekin_cmdline.c treekin_cmdline.h

globals.c:  treekin_cmdline.h

if MAKE_CMDLINE_PARSER
treekin_cmdline.h treekin_cmdline.c:  treekin.ggo
	@GENGETOPT@ -i treekin.ggo
endif



# End of file
