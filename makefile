# $Id: makefile,v 1.15 2003/10/23 10:47:28 mtw Exp $
CC      = gcc
SCRDIR  = .
INCDIR  = -I$(SCRDIR)
CDEBUG  = -g3 -O
COPTIM  = -O3 -march=i686
CFLAGS  = -Wall -pedantic $(INCDIR) $(COPTIM)
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

