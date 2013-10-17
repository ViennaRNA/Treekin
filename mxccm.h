/*=================================================================*/
/*=   mxccm.h                                                     =*/
/*=   header file for matrix routines from ccmath library         =*/
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

/*extract all eigen values and vectors of a real symmetric matrix.*/
/*  input: a - symetric matrix, n - dimension
    output: a - eigenvectors, ev - eigenvalues*/
void eigen(double *a,double *ev,int n);

/*Perform a QR reduction of a real symmetric tridiagonal matrix
     to diagonal form and update an orthogonal transformation matrix.*/
int qrevec(double *ev,double *v,double *d,int m);

/*Transform a real symmetric matrix to tridiagonal form andŔ
     compute the orthogonal matrix of this transformation.*/
void housev(double *a,double *d,double *dp,int n);

/*Multiply two real square matrices C = A * B.*/
void mmul(double *c,double *a,double *b,int n);

/*Multiply a vector by a matrix vp = mat*v.*/
void vmul(double *vp,double *mat,double *v,int n);

/*Transpose a real square matrix in place A -> A~.*/
void trnm(double *a,int m);

/*Copy an array a = b.*/
void mcopy(double *a,double *b,int m);

/*Invert (in place) a general real matrix A -> Inv(A).*/
int minv(double *a,int n);

// new by Marcel
// multiply non-square matrices
void mmul_singular(double *c,double *a,double *b,int dim1, int dim2, int dim3, int verbose);

#endif
