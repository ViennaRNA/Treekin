/* calc.h */
/* Last changed Time-stamp: <2003-11-04 11:35:53 mtw> */
/*  static char rcsid[] = "$Id: calc.h,v 1.8 2003/11/18 17:27:59 mtw Exp $"; */

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

#endif 

/* End of file */
