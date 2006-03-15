# $Id: makefile,v 1.16 2006/03/15 11:08:15 mtw Exp $
CC      = gcc
SCRDIR  = .
INCDIR  = -I$(SCRDIR)
CDEBUG  = -g3 -O
COPTIM  = -O3 -ffast-math
CARCH   = -march=pentium4
CFLAGS  = -Wall -pedantic $(INCDIR) $(COPTIM) $(CARCH)
LIBS    = -lm 
EXE     = treekin
OBJS    = \
	 main.o\
	 err.o mxmesch.o\
	 schur.o\
	 calc.o\
	 mxccm.o\
	 globals.o\
	 barparser.o\
	 exp_matrix.o	

%.o : $(SCRDIR)%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(EXE):	$(OBJS)
	$(CC) -o $(EXE) $(CFLAGS) $(OBJS) $(LIBS)

clean:
	rm -f  *.o core

veryclean: clean
	rm -f $(EXE)

new:   clean %.o 

