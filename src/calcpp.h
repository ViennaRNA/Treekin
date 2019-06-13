#ifndef _CALCPP_H_
#define _CALCPP_H_

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifdef WITH_MLAPACK
#ifdef WITH_MLAPACK_GMP
#include <gmpxx.h>
#include <mlapack/mblas_gmp.h>
#include <mlapack/mlapack_gmp.h>
#endif
#ifdef WITH_MLAPACK_QD
#include <qd/qd_real.h>
#include <mlapack/mblas_qd.h>
#include <mlapack/mlapack_qd.h>
#endif
#ifdef WITH_MLAPACK_DD
#include <qd/dd_real.h>
#include <mlapack/mblas_dd.h>
#include <mlapack/mlapack_dd.h>
#endif
#ifdef WITH_MLAPACK_MPFR
#include <mlapack/mpreal.h>
#include <mlapack/mblas_mpfr.h>
#include <mlapack/mlapack_mpfr.h>
#endif
#ifdef WITH_MLAPACK___FLOAT128
#include <mlapack/mutils___float128.h>
#include <mlapack/mblas___float128.h>
#include <mlapack/mlapack___float128.h>
#endif
#ifdef WITH_MLAPACK_DOUBLE
#include <mlapack/mblas_double.h>
#include <mlapack/mlapack_double.h>
#endif
#ifdef WITH_MLAPACK_LD
#include <mlapack/mblas_longdouble.h>
#include <mlapack/mlapack_longdouble.h>
#endif
#endif

#include <cmath>

#ifdef HAVE_LAPACKE_H
# include <lapacke.h>
#else
# ifdef HAVE_LAPACKE_LAPACKE_H
#   include <lapacke/lapacke.h>
# else
#   ifdef HAVE_OPENBLAS_LAPACKE_H
#     include <openblas/lapacke.h>
#   endif
# endif
#endif

#include "bardata.h"
#include "globals.h"

#ifdef WITH_MLAPACK
# include "treekinCastableTypes.h"
#endif

class Calccpp {
  treekin_options *opt;
  SubInfo *E;

public:
  Calccpp(treekin_options *opt,
          SubInfo         *E)
  {
    this->opt = opt;
    this->E   = E;
  }


  ~Calccpp()
  {}


  /* rescale is --minimal-rate hass been provided */
  template<class T>
  void MxRescale(T      *U,
                 int    dim,
                 double desired_rate);


  template<class T>
  void MxRescaleH(T       *U,
                  int     dim,
                  double  hard_rescale);


  template<class T>
  void MxTimes(T      *U,
               int    dim,
               double times);


  template<class T>
  int *MxErgoEigen(T    *U,
                   int  dim);


#if defined(WITH_MLAPACK_GMP)
  void MxEV_Mlapack_Sym_gmp(const mpf_class *U,
                          int             dim,
                          mpf_class       *evals,
                          mpf_class       *evecs,
                          int             precision);


#endif
#if defined(WITH_MLAPACK_QD)
  void MxEV_Mlapack_Sym_qd(const qd_real  *U,
                         int            dim,
                         qd_real        *evals,
                         qd_real        *evecs);


#endif
#if defined(WITH_MLAPACK_DD)
  void MxEV_Mlapack_Sym_dd(const dd_real  *U,
                         int            dim,
                         dd_real        *evals,
                         dd_real        *evecs);


#endif
#if defined(WITH_MLAPACK_MPFR)
  void MxEV_Mlapack_Sym_mpfr(const mpfr::mpreal *U,
                           int                dim,
                           mpfr::mpreal       *evals,
                           mpfr::mpreal       *evecs,
                           int                precision);


#endif
#if defined(WITH_MLAPACK___FLOAT128)
  void MxEV_Mlapack_Sym_float128(const __float128 *U,
                               int              dim,
                               __float128       *evals,
                               __float128       *evecs);


#endif
#if defined(WITH_MLAPACK_LD)
  void MxEV_Mlapack_Sym_longdouble(const long double  *U,
                                 int                dim,
                                 long double        *evals,
                                 long double        *evecs);


#endif
#if defined(WITH_MLAPACK_DOUBLE)
  void MxEV_Mlapack_Sym_double(const double *U,
                             int          dim,
                             double       *evals,
                             double       *evecs);


#endif

  template<typename T>
  void Mx_Dgesv(int *n,
                int *nrhs,
                T   *A,
                int *m,
                int *ipiv,
                T   *B,
                int *l,
                int *nfo);


  template<typename T>
  void Mx_Dgeevx(const char *balanc,
                 const char *jobvl,
                 const char *jobvr,
                 const char *sense,
                 int        *n,
                 T          *a,
                 int        *lda,
                 T          *wr,
                 T          *wi,
                 T          *vl,
                 int        *ldvl,
                 T          *vr,
                 int        *ldvr,
                 int        *ilo,
                 int        *ihi,
                 T          *scale,
                 T          *abnrm,
                 T          *rconde,
                 T          *rcondv,
                 T          *work,
                 int        *lwork,
                 int        *iwork,
                 int        *info);


  template<typename T>
  void Mx_Dsyevx(const char *jobz,
                 const char *range,
                 const char *uplo,
                 int        *n,
                 T          *a,
                 int        *lda,
                 T          *vl,
                 T          *vu,
                 int        *il,
                 int        *iu,
                 T          *abstol,
                 int        *m,
                 T          *w,
                 T          *z,
                 int        *ldz,
                 T          *work,
                 int        *lwork,
                 int        *iwork,
                 int        *ifail,
                 int        *info);
};

template<typename T>
inline
void
Calccpp::Mx_Dgesv(int *n,
                  int *nrhs,
                  T   *A,
                  int *m,
                  int *ipiv,
                  T   *B,
                  int *l,
                  int *nfo)
{
  /*DGESV computes the solution to a real system of linear equations A * X = B */
  switch (opt->mlapackMethod) {
#if defined(WITH_MLAPACK)
    case MLAPACK_GMP:
    case MLAPACK_QD:
    case MLAPACK_MPFR:
    case MLAPACK_FLOAT128:
    case MLAPACK_LD:
    case MLAPACK_DD:
    case MLAPACK_DOUBLE:
      Rgesv((mpackint) * n,
            (mpackint) * nrhs,
            A,
            (mpackint) * m,
            (mpackint *)ipiv,
            B,
            (mpackint) * l,
            (mpackint *)nfo);
      break;
#endif
    default:
      /*default standard lapack */
      double *AA = (double *)A;
      double *BB = (double *)B;

      ::dgesv_(n, nrhs, AA, m, ipiv, BB, l, nfo);
      break;
  }
}


template<typename T>
inline
void
Calccpp::Mx_Dsyevx(const char *jobz,
                   const char *range,
                   const char *uplo,
                   int        *n,
                   T          *a,
                   int        *lda,
                   T          *vl,
                   T          *vu,
                   int        *il,
                   int        *iu,
                   T          *abstol,
                   int        *m,
                   T          *w,
                   T          *z,
                   int        *ldz,
                   T          *work,
                   int        *lwork,
                   int        *iwork,
                   int        *ifail,
                   int        *info)
{
  /*DGESV computes the solution to a real system of linear equations A * X = B */
  switch (opt->mlapackMethod) {
#if defined(WITH_MLAPACK)
    case MLAPACK_GMP:
    case MLAPACK_QD:
    case MLAPACK_MPFR:
    case MLAPACK_FLOAT128:
    case MLAPACK_LD:
    case MLAPACK_DD:
    case MLAPACK_DOUBLE:
      Rsyevx(jobz, range, uplo, (mpackint) * n, a, (mpackint) * lda, *vl,
             *vu, (mpackint) * il, (mpackint) * iu, *abstol, (mpackint *)m, w,
             z, (mpackint) * ldz, work, (mpackint) * lwork, (mpackint *)iwork, (mpackint *)ifail,
             (mpackint *)info);
      break;
#endif
    default:
      /*default standard lapack */
      ::dsyevx_((char *)jobz, (char *)range, (char *)uplo, n, (double *)a, lda, (double *)vl,
              (double *)vu, il, iu, (double *)abstol, m, (double *)w,
              (double *)z, ldz, (double *)work, lwork, iwork, ifail, info);
      break;
  }
}


template<typename T>
inline
void
Calccpp::Mx_Dgeevx(const char *balanc,
                   const char *jobvl,
                   const char *jobvr,
                   const char *sense,
                   int        *n,
                   T          *a,
                   int        *lda,
                   T          *wr,
                   T          *wi,
                   T          *vl,
                   int        *ldvl,
                   T          *vr,
                   int        *ldvr,
                   int        *ilo,
                   int        *ihi,
                   T          *scale,
                   T          *abnrm,
                   T          *rconde,
                   T          *rcondv,
                   T          *work,
                   int        *lwork,
                   int        *iwork,
                   int        *info)
{
  /*DGESV computes the solution to a real system of linear equations A * X = B */
  switch (opt->mlapackMethod) {
#if defined(WITH_MLAPACK)
    case MLAPACK_GMP:
    case MLAPACK_QD:
    case MLAPACK_MPFR:
    case MLAPACK_FLOAT128:
    case MLAPACK_LD:
    case MLAPACK_DD:
    case MLAPACK_DOUBLE:
      Rgeevx(balanc, jobvl, jobvr, sense, (mpackint) * n, a, (mpackint) * lda, wr,
             wi, vl, *ldvl, vr, (mpackint) * ldvr,
             (mpackint *)ilo, (mpackint *)ihi, scale, abnrm, rconde,
             rcondv, work, (mpackint) * lwork,
             (mpackint *)iwork, (mpackint *)info);
      break;
#endif
    default:
      /*default standard lapack */
      ::dgeevx_((char *)balanc, (char *)jobvl, (char *)jobvr, (char *)sense, n, (double *)a, lda,
              (double *)wr, (double *)wi, (double *)vl, ldvl,
              (double *)vr, ldvr, ilo, ihi, (double *)scale, (double *)abnrm, (double *)rconde,
              (double *)rcondv, (double *)work,
              lwork, iwork, info);
      break;
  }
}


template<class T>
void
Calccpp::MxRescale(T      *U,
                   int    dim,
                   double desired_rate)
{
  /*first get the minimal rate */
  T minimal_rate = 1.0;

  for (int i = 0; i < dim * dim; i++)
    if (i % dim != i / dim && minimal_rate > U[i] && U[i] > 0.0)
      minimal_rate = U[i];

  /* calculate rescale factor: */
  T factor = (T)(log(desired_rate) / log(minimal_rate));
  fprintf(stderr,
          "rescale params: %20.10g %20.10g %20.10g\n",
          (double)minimal_rate,
          (double)desired_rate,
          (double)factor);

  if (desired_rate < minimal_rate)
    return;

  MxRescaleH(U, dim, factor);
}


template<class T>
void
Calccpp::MxRescaleH(T       *U,
                    int     dim,
                    double  hard_rescale)
{
  T factor = (T)hard_rescale;

  /* rescale! */
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      if (i != j)
        U[i * dim + j] = (T)std::pow(U[i * dim + j], factor);

  /* fix the diagonal: */
  for (int i = 0; i < dim; i++)
    U[dim * i + i] = 0.;
  for (int i = 0; i < dim; i++) {
    T tmp = 0.00;
    /* calculate column sum */
    for (int j = 0; j < dim; j++)
      tmp += U[dim * j + i];
    U[dim * i + i] = (T)(-tmp + (opt->useplusI ? 1.0 : 0.0));/* make U a stochastic matrix U = Q+I ?? */
  }
}


template<class T>
void
Calccpp::MxTimes(T      *U,
                 int    dim,
                 double times)
{
  /* multiply! */
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      if (i != j)
        U[i * dim + j] *= times;

  /* fix the diagonal: */
  for (int i = 0; i < dim; i++)
    U[dim * i + i] = 0.;
  for (int i = 0; i < dim; i++) {
    T tmp = 0.00;
    /* calculate column sum */
    for (int j = 0; j < dim; j++)
      tmp += U[dim * j + i];
    U[dim * i + i] = (T)(-tmp + (opt->useplusI ? 1.0 : 0.0));/* make U a stochastic matrix U = Q+I ?? */
  }
}


#endif
