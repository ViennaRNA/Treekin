/* mxccm.h */
/* Last changed Time-stamp: <2001-10-11 17:17:44 mtw> */
/* static char rcsid[] = "$Id: mxccm.h,v 1.1 2003/07/14 07:42:20 mtw Exp $"; */

void eigen(double *a,double *ev,int n);
int qrevec(double *ev,double *v,double *d,int m);
void housev(double *a,double *d,double *dp,int n);
void mmul(double *c,double *a,double *b,int n);
void vmul(double *vp,double *mat,double *v,int n);
void trnm(double *a,int m);
void mcopy(double *a,double *b,int m);
int minv(double *a,int n);
