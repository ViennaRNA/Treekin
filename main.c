/*=================================================================*/
/*=   main.c                                                      =*/
/*=   main file for treekin                                       =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-11-10 13:04:09 mtw>          =*/
/*=   $Id: main.c,v 1.21 2006/11/10 13:51:39 mtw Exp $            =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include "calc.h"
#include "mxccm.h"  
#include "barparser.h"
#include "globals.h"

int
main (int argc, char **argv)
{
  BarData *Data=NULL;
  double *U=NULL, *S=NULL, *p0=NULL, *p8=NULL, *R=NULL;
  int  dim;

  parse_commandline(argc, argv);
 
  if(opt.method == 'F')      /* full process */
    dim = ParseInfile(opt.INFILE, &R);
  else
    dim = ParseBarfile (opt.INFILE, &Data);
  if(opt.method == 'I'){
    if (opt.binrates)
      dim =  MxBinRead(&R);
    else
      ParseRatesFile(&R, dim);
  }
  MxInit (dim);

  p0 = MxStartVec ();
  if(opt.method == 'F') 
    p8 = MxEqDistrFULL (E);
  else
    p8 = MxEqDistr (Data);
  
  if(opt.method == 'F')
    U = MxBar2Matrix (NULL, R);
  else
    U = MxBar2Matrix (Data, R);
  
  if(opt.fpt){
    MxFPT(U, p8);
    MxFirstPassageTime(U, p8);
  }

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
  if(Data != NULL) free(Data);
  return 0;
}

/* End of file */



