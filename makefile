# $Id: makefile,v 1.12 2003/09/25 13:51:01 mtw Exp $
CC      = gcc
SCRDIR  = .
INCDIR  = -I$(SCRDIR) -I/scratch/mtw/meschach-1.2b
CDEBUG  = -g3 -O
COPTIM  = -O3 -march=i686
CFLAGS  = -Wall $(INCDIR) $(CDEBUG)
LIBS    = -lm /scratch/mtw/meschach-1.2b/meschach.a 
EXE     = treekin
OBJS    = \
	 main.o\
	 barparser.o\
	 calc.o\
	 mxccm.o\
	 globals.o\
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

