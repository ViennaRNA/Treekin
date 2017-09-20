dist_man_MANS = treekin.1

bin_PROGRAMS = treekin

treekin_LDADD=$(mpack_LIBS)

treekin_SOURCES = mxccm.cpp mxccm.h \
                  globals.cpp globals.h \
                  barparser.cpp barparser.h bardata.h\
                  calcpp.cpp  calcpp.h \
                  calc.cpp calc.h \
                  exp_matrix.cpp  exp_matrix.h \
                  treekin_cmdline.c treekin_cmdline.h \
                  treekinCastableTypes.h main.cpp
                  
                  

EXTRA_DIST = treekin.ggo treekin_cmdline.c treekin_cmdline.h

globals.cpp:  treekin_cmdline.h

if MAKE_CMDLINE_PARSER
treekin_cmdline.h treekin_cmdline.c:  treekin.ggo
	@GENGETOPT@ -i treekin.ggo
endif



# End of file
