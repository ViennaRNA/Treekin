/* calc.h */
/* Last changed Time-stamp: <2003-12-02 16:22:50 mtw> */
/*  static char rcsid[] = "$Id: calc.h,v 1.9 2005/06/21 10:08:31 mtw Exp $"; */

#ifndef _CALC_H_
#define _CALC_H_
#include "barparser.h"

void    MxInit (int d);
double *MxBar2Matrix (BarData *Data, double *);
double *MxEqDistr (BarData *Data);
double *MxEqDistrFULL (SubInfo *E);
double *MxMethodeFULL(InData *InD);
double *MxMethodeINPUT (BarData *Data, double *);
double *MxSymmetr (double *U, double *PI);
double *MxStartVec (void) ;
void    MxEVnonsymMx(double *U, double **_S);
void    MxIterate (double *p0, double *p8, double *S);
void    MxMemoryCleanUp (void);
void    MxExponent(double *p0, double *p8, double *U);
void    MxFPT(double *U, double *p8);

#endif 

/* End of file */
