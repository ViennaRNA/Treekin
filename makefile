# $Id: makefile,v 1.4 2003/08/09 15:58:19 mtw Exp $
CC      = gcc
SCRDIR  = .
INCDIR  = -I$(SCRDIR) -I/scratch/mtw/meschach-1.2b
CDEBUG  = -g3
COPTIM  = -O -march=i686
CFLAGS  = -Wall $(INCDIR) $(CDEBUG) -DDMALLOC
LIBS    = -lm -ldmalloc /scratch/mtw/meschach-1.2b/meschach.a 
EXE     = treekin
OBJS    = \
	 main.o\
	 calc.o\
	 mxccm.o\
	 globals.o\
	 barparser.o

%.o : $(SCRDIR)%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(EXE):	$(OBJS)
	$(CC) -o $(EXE) $(CFLAGS) $(OBJS) $(LIBS)

clean:
	rm -f  *.o core

veryclean: clean
	rm -f $(EXE)

new:   clean %.o 

