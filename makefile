# $Id: makefile,v 1.2 2003/07/16 12:27:01 mtw Exp $
CC      = gcc
SCRDIR  = .
INCDIR  = -I$(SCRDIR) -I/scratch/mtw/meschach
CDEBUG  = -g3
COPTIM  = -O3 -march=i686
CFLAGS  = -Wall $(INCDIR) $(CDEBUG)
LIBS    = -lm /scratch/mtw/meschach/meschach.a
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

