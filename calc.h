/* calc.h */
/* Last changed Time-stamp: <2003-09-10 14:40:58 mtw> */
/*  static char rcsid[] = "$Id: calc.h,v 1.5 2003/09/10 13:53:11 mtw Exp $"; */

#ifndef _CALC_H_
#define _CALC_H_
#include "barparser.h"

void    MxInit (int d);
double *MxBar2Matrix (TypeBarData *Data, double *);
double *MxEqDistr (TypeBarData *Data);
double *MxEqDistrFULL (double *energi);
double *MxMethodeFULL(InData *InD);
double *MxMethodeINPUT (TypeBarData *Data, double *);
double *MxSymmetr (double *U, double *PI);
double *MxStartVec (void) ;
void    MxEVnonsymMx(double *U, double **_S);
void    MxIterate (double *p0, double *p8, double *S);
void    MxIteratehack(double *p0, double *S, double *T, double *EV_m);
void    MxIterate_FULL( double *p0, double *p8, double *S, int *assoc_gradbas, int lmin_nr);
void    MxMemoryCleanUp (void);
void    MxExponent(double *p0, double *p8, double *U);

#endif 

/* End of file */
