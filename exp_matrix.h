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

#ifdef WITH_MPACK
#ifdef WITH_MPACK_GMP
   #include <gmpxx.h>
   #include <mpack/mblas_gmp.h>
   #include <mpack/mlapack_gmp.h>
#endif
#ifdef WITH_MPACK_QD
  #include <qd/qd_real.h>
   #include <mpack/mblas_qd.h>
   #include <mpack/mlapack_qd.h>
#endif
#ifdef WITH_MPACK_DD
    #include <qd/dd_real.h>
   #include <mpack/mblas_dd.h>
   #include <mpack/mlapack_dd.h>
#endif
#ifdef WITH_MPACK_MPFR
#include <mpack/mpreal.h>
   #include <mpack/mblas_mpfr.h>
   #include <mpack/mlapack_mpfr.h>
#endif
#ifdef WITH_MPACK___FLOAT28
#include <mpack/mutils___float128.h>
   #include <mpack/mblas___float128.h>
   #include <mpack/mlapack___float128.h>
#endif
#ifdef WITH_MPACK_DOUBLE
#include <mpack/mutils_double.h>
   #include <mpack/mblas_double.h>
   #include <mpack/mlapack_double.h>
#endif
#ifdef WITH_MPACK_LD
    #include <mpack/mblas_longdouble.h>
    #include <mpack/mlapack_longdouble.h>
#endif
#endif

#include <algorithm> //for copy
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mxccm.h"

#ifdef WITH_MPACK
# include "treekinCastableTypes.h"
#endif

using namespace std;

class ExpMatrix {

private:
  Mxccm mxccm;

  template<typename T>
  void ipmmul(T *a,T *b, T *c,int n);

public:
  ExpMatrix() : mxccm(){

  }
  ~ExpMatrix(){

  }

  template<typename T>
  void padexp(T *from,T *out,int n,int ord, double epsilon);
};



template<typename T>
void
ExpMatrix::ipmmul(T *a,T *b, T *c,int n)
{
  T *T1,*T2;
  size_t n_size = n*n;

  T1 = new T[n_size];
  T2 = new T[n_size];

  std::copy_n(a, n_size, T1);
  std::copy_n(b, n_size, T2);
  mxccm.mmul(T1,T2,c,n);
  std::copy_n(T1, n_size, a);
  delete[] T2;
  delete[] T1;
  return;
}

template<typename T>
void
ExpMatrix::padexp(T *from,T *out,int n,int ord, double epsilon)
{
  /* calculate exp(M) by a Pade approx */
  /* of order ord                      */
  /* exp(x) ~= p(x)/q(x)               */
  /* Sidje's trick makes               */
  /* exp(x) ~= 1 + P(x)/Q(x)           */

  int i,k;
  T *cp;
  T *mat_in,*MD,*MN,*X2;
  /*precon stuff */
  T maxnorm,ml;
  int ppow=0;
  T scfac;


  mat_in = new T[n*n];
  for (i=0; i<n*n; i++) mat_in[i] = from[i];

  /* get the infinity norm of the matrix for scaling/squaring*/
  maxnorm=-1;

  for (i=0; i<n*n; i++)
    maxnorm = maxnorm < (T)(abs(mat_in[i])) ? (T)(abs(mat_in[i])) : maxnorm ;

  ml =  (T)(log(maxnorm)/log(2.));

  /* yuck, heavy heuristics */
  /* calculate M/(2^n) to get exp (M) = (exp (M/2^n))^(2^n) afterwards */

  if (ml >=4.)  {
    ppow = (int) (T)(ml - 2.);

    scfac=1.;
    for (i=1; i<=ppow; i++) scfac *=2.;
    for (i=0; i<n*n; i++) mat_in[i]/=scfac;
  }
  /* cp holds the Pade coefficients ,
     X2 the squared input matrix    ,
     MN the numerator polynomial    ,
     MD the denominator polynomial */
  cp = new T[ord+1];
  X2 = new T[n*n];
  MN = new T[n*n];
  MD = new T[n*n];

  cp[0]=1.;
  for (i=1; i<=ord; i++)
    cp[i] = (T)(cp[i-1]*((T)(ord+1-i))/((T)((ord+ord+1-i)*i)));

  ipmmul(X2,mat_in,mat_in,n);

  for (i=0; i<n*n; i++) MN[i]=MD[i]=0.;

  for (i=0; i<n; i++) {
    MN[(n+1)*i]=cp[ord-1];
    MD[(n+1)*i]=cp[ord];
  }

  for (k=ord/2-1; k>0; k--) {
    ipmmul(MN,X2,MN,n);
    for (i=0; i<n; i++) MN[(n+1)*i]+=cp[2*k-1];
  }

  for (k=ord/2; k>0; k--) {
    ipmmul(MD,X2,MD,n);
    for (i=0; i<n; i++) MD[(n+1)*i] +=cp[2*k-2];
  }
  ipmmul(MN,mat_in,MN,n);

  for(i=0; i<n*n; i++) MD[i]=(T)(MD[i]-MN[i]);
  mxccm.minv(MD,n, epsilon);
  ipmmul(mat_in,MD,MN,n);
  for(i=0; i<n*n; i++) mat_in[i]*=2.0;
  for(i=0; i<n; i++) mat_in[(n+1)*i] +=1.0;

  if (ml >=6.) {
    for (i=1; i<=ppow; i++)
      ipmmul(mat_in,mat_in,mat_in,n);
  }

  std::copy_n(mat_in,n*n,out);

  delete[] MD;
  delete[] MN;
  delete[] X2;
  delete[] cp;
  delete[] mat_in;
}


#endif
