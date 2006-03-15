/*=================================================================*/
/*=   exp_matrix.h                                                =*/
/*=   header file for matrix exponentials via pade approximation  =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 15:16:06 mtw>          =*/
/*=   $Id: exp_matrix.h,v 1.3 2006/03/15 14:18:20 mtw Exp $       =*/
/*=   ---------------------------------------------------------   =*/
/*=     (c) Michael Thomas Wolfinger, W. Andreas Svrcek-Seiler    =*/
/*=                  {mtw,svrci}@tbi.univie.ac.at                 =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _EXP_MATRIX_H_
#define _EXP_MATRIX_H_

void padexp(double *from,double *out,int n,int ord);

#endif
