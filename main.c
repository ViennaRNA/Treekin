/* main.c */
/* Last changed Time-stamp: <2003-12-02 16:01:42 mtw> */
/* static char rcsid[] = "$Id: main.c,v 1.15 2005/06/21 10:08:31 mtw Exp $"; */

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

  BarData *Data;
  InData *InD;
  double *U, *S, *p0, *p8, *R = NULL;
  int  dim,fpt;

  fpt = 0; /* first passage time */
  parse_commandline(argc, argv);
 
  if(opt.method == 'F')      /* full process */
    dim = ParseInfile(opt.INFILE, &InD);
  else
    dim = ParseBarfile (opt.INFILE, &Data);
  if(opt.method == 'I')  ParseRatesFile(&R, dim);
  MxInit (dim);
  if(opt.method == 'F')
    U = MxMethodeFULL(InD);
  else
    U = MxBar2Matrix (Data, R);
  p0 = MxStartVec ();
  if(opt.method == 'F') 
    p8 = MxEqDistrFULL (E);
  else
    p8 = MxEqDistr (Data);
  if(fpt) MxFPT(U, p8);
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
  if(opt.method == 'F') free(InD);
  else free(Data);
  return 0;
}

/* End of file */



