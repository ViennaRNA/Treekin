/* main.c */
/* Last changed Time-stamp: <2003-07-16 19:14:07 mtw> */
/* static char rcsid[] = "$Id: main.c,v 1.3 2003/07/16 17:14:56 mtw Exp $"; */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include "calc.h"
#include "mxccm.h"  
#include "barparser.h"
#include "globals.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

int in_nr, n;

/* ================ main begins here ================ */
int main (int argc, char **argv) {

  TypeBarData *Data;
  InData *InD;
  double *U, *S, *p0, *p8, *R;
  int  dim;
  
  /* arrays  */ 
  double *Energies;  /* array 4 the energies of all entries from subopt */
                     /* needed for calculation of p8 in full process */
  int *lmin_nr_so;   /* array containing info which lmin is which # in subopt */
  int *assoc_gradbas;/* array 4 associated gradient basins of each # of subopt */ 
  
  parse_commandline(argc, argv);
 
  if(opt.method == 'F'){             /* full process */
    if(opt.absrb != 0) {
      fprintf(stderr, "absorbing states not yet implemented in full process\n");
      exit(999);
    }
    dim = ParseInfile(opt.INFILE, &InD, &Energies, &lmin_nr_so, &assoc_gradbas);
    MxInit (dim);
    U = MxMethodeFULL(InD);          
    p8 = MxEqDistrFULL (Energies);
    S = MxSymmetr (U, p8);
    p0 = MxStartVec ();
    MxIterate_effective_lmins (p0, p8, S, lmin_nr_so, assoc_gradbas);
    MxMemoryCleanUp();
    
    if (opt.pini != NULL) free(opt.pini);
    free(p8);
    free(S);
    free(U);
    free(InD);
    free(Energies);
    free(lmin_nr_so);
    free(assoc_gradbas);
  }
  else if(opt.method == 'F'){
    if(opt.absrb != 0) {
      fprintf(stderr, "absorbing states not yet implemented for input matrix case\n");
      exit(999);
    }
    dim = ParseBarfile (opt.INFILE, &Data);
    ParseRatesFile(&R, dim);
    MxInit (dim);

  }
  else {                         /* tree process */
    dim = ParseBarfile (opt.INFILE, &Data);
    MxInit (dim);
    U = MxBar2Matrix (Data);
    p0 = MxStartVec ();
    p8 = MxEqDistr (Data);
    if(opt.absrb > 0) MxEVnonsymMx(U, &S);
    else S = MxSymmetr (U, p8);
    MxIterate (p0, S);
    MxMemoryCleanUp();
    
    if (opt.pini != NULL) free(opt.pini);
    free(p8);
    free(S);
    free(U);
    free(Data);
    /* saddles freen !!! */
  }
  return 0;
}

/* End of file */



