/*=================================================================*/
/*=   mxccm.h                                                     =*/
/*=   header file fro matrix routines from ccmath library         =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 15:16:26 mtw>          =*/
/*=   $Id: mxccm.h,v 1.3 2006/03/15 14:18:20 mtw Exp $            =*/
/*=   ---------------------------------------------------------   =*/
/*=      (c) Daniel A. Atkinson, Michael Thomas Wolfinger         =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _MXCCM_H_
#define _MXCCM_H_

void eigen(double *a,double *ev,int n);
int qrevec(double *ev,double *v,double *d,int m);
void housev(double *a,double *d,double *dp,int n);
void mmul(double *c,double *a,double *b,int n);
void vmul(double *vp,double *mat,double *v,int n);
void trnm(double *a,int m);
void mcopy(double *a,double *b,int m);
int minv(double *a,int n);

#endif
