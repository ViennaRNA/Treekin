# $Id: makefile,v 1.8 2003/09/11 11:45:30 mtw Exp $
CC      = gcc
SCRDIR  = .
INCDIR  = -I$(SCRDIR) -I/scratch/mtw/meschach-1.2b
CDEBUG  = -g3
COPTIM  = -O3 -march=i686
CFLAGS  = -Wall $(INCDIR) $(CDEBUG) 
LIBS    = -lm /scratch/mtw/meschach-1.2b/meschach.a 
EXE     = treekin
OBJS    = \
	 main.o\
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

