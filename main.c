/* main.c */
/* Last changed Time-stamp: <2003-09-22 18:52:44 mtw> */
/* static char rcsid[] = "$Id: main.c,v 1.11 2003/09/26 14:41:19 mtw Exp $"; */

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

/* ============= main begins here ============== */
int main (int argc, char **argv) {

  TypeBarData *Data;
  InData *InD;
  double *U, *S, *p0, *p8, *R;
  int  dim,lmins;
  
  parse_commandline(argc, argv);
 
  if(opt.method == 'F'){         /* full process */
    dim = ParseInfile(opt.INFILE, &InD, &lmins);
    MxInit (dim);
    U = MxMethodeFULL(InD);
    p0 = MxStartVec ();
    p8 = MxEqDistrFULL (E);
    if(opt.absrb) MxEVnonsymMx(U, &S);
    else S = MxSymmetr (U, p8);
    MxIterate_FULL (p0, p8, S, lmins);
    MxMemoryCleanUp();
    
    if (opt.pini != NULL) free(opt.pini);
    free(U);free(S);free(p8);
    free(InD);
  }
  else {                   /* tree/rates  process */
    dim = ParseBarfile (opt.INFILE, &Data);
    if(opt.method == 'I')  ParseRatesFile(&R, dim);
    MxInit (dim);
    U = MxBar2Matrix (Data, R);
    p0 = MxStartVec ();
    p8 = MxEqDistr (Data);
    if(opt.matexp) MxExponent(p0,p8,U);
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
  }
  return 0;
}

/* End of file */



