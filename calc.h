/* calc.h */
/* Last changed Time-stamp: <2003-09-10 14:40:58 mtw> */
/*  static char rcsid[] = "$Id: calc.h,v 1.6 2003/09/25 13:51:01 mtw Exp $"; */

#ifndef _CALC_H_
#define _CALC_H_
#include "barparser.h"

void    MxInit (int d);
double *MxBar2Matrix (TypeBarData *Data, double *);
double *MxEqDistr (TypeBarData *Data);
double *MxEqDistrFULL (SubInfo *E);
double *MxMethodeFULL(InData *InD);
double *MxMethodeINPUT (TypeBarData *Data, double *);
double *MxSymmetr (double *U, double *PI);
double *MxStartVec (void) ;
void    MxEVnonsymMx(double *U, double **_S);
void    MxIterate (double *p0, double *p8, double *S);
void    MxIterate_FULL( double *p0, double *p8, double *S, int lmins);
void    MxMemoryCleanUp (void);
void    MxExponent(double *p0, double *p8, double *U);

#endif 

/* End of file */
