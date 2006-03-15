/*=================================================================*/
/*=   exp_matrix.h                                                =*/
/*=   header file for matrix exponentials via pade approximation  =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 11:21:18 mtw>          =*/
/*=   $Id: exp_matrix.h,v 1.2 2006/03/15 11:08:15 mtw Exp $    =*/
/*=   ---------------------------------------------------------   =*/
/*=     (c) Michael Thomas Wolfinger, W. Andreas Svrcek-Seiler    =*/
/*=                  {mtw,svrci}@tbi.univie.ac.at                 =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _EXP_MATRIX_H_
#define _EXP_MATRIX_H_

void padexp(double *from,double *out,int n,int ord);

#endif
