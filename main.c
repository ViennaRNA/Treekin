/* main.c */
/* Last changed Time-stamp: <2003-09-22 18:52:44 mtw> */
/* static char rcsid[] = "$Id: main.c,v 1.8 2003/09/23 16:29:55 mtw Exp $"; */

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
    dim = ParseInfile(opt.INFILE, &InD, &Energies, &lmin_nr_so, &assoc_gradbas);
    MxInit (dim);
    U = MxMethodeFULL(InD);
    p0 = MxStartVec ();
    p8 = MxEqDistrFULL (Energies);
    if(opt.absrb) MxEVnonsymMx(U, &S);
    else S = MxSymmetr (U, p8);
    MxIterate_FULL (p0, p8, S, assoc_gradbas, lmin_nr_so[0]);
    MxMemoryCleanUp();
    
    if (opt.pini != NULL) free(opt.pini);
    free(U);free(S);free(p8);
    free(InD);
    free(Energies);
    free(lmin_nr_so);
    free(assoc_gradbas);
  }
  else {                         /* tree/rates  process */
    dim = ParseBarfile (opt.INFILE, &Data);
    if(opt.method == 'I')  ParseRatesFile(&R, dim);
    MxInit (dim);
    U = MxBar2Matrix (Data, R);
    p0 = MxStartVec ();
    p8 = MxEqDistr (Data);
    if(opt.matexp){
      MxExponent(p0,p8,U);
    }
    else{
      if(opt.absrb) MxEVnonsymMx(U, &S);
      else S = MxSymmetr (U, p8);
      MxIterate (p0, p8, S);
      free(S);free(p8);
    }
    MxMemoryCleanUp();
    
    if (opt.pini != NULL) free(opt.pini);
    free(U);
    free(Data);
    /* saddles freen !!! */
  }
  return 0;
}

/* End of file */



