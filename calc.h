/* calc.h */
/* Last changed Time-stamp: <2003-07-19 11:51:46 mtw> */
/*  static char rcsid[] = "$Id: calc.h,v 1.2 2003/08/05 08:40:04 mtw Exp $"; */

#ifndef _CALC_H_
#define _CALC_H_
#include "barparser.h"

void    MxInit (int d);
double *MxBar2Matrix (TypeBarData *Data);
double *MxEqDistr (TypeBarData *Data);
double *MxEqDistrFULL (double *energi);
double *MxMethodeFULL(InData *InD);
double *MxMethodeINPUT (double *);
double *MxSymmetr (double *U, double *PI);
double *MxStartVec (void) ;
void    MxEVnonsymMx(double *U, double **_S);
void    MxIterate (double *p0, double *S);
void    MxIteratehack(double *p0, double *S, double *T, double *EV_m);
void    MxIterate_effective_lmins ( double *p0, double *p8, double *S, int *lmin_nr_so, int *assoc_gradbas);
void    MxMemoryCleanUp (void);

#endif 

/* End of file */
