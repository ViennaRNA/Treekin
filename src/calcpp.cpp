/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "calcpp.h"


extern "C" {
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
void dsyevx_( char* jobz, char* range, char* uplo, int* n, double* a,
                int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
                int* m, double* w, double* z, int* ldz, double* work, int* lwork,
                int* iwork, int* ifail, int* info );
void dgeevx_ (char *balanc, char *jobvl, char *jobvr,
              char *sense,
              int *n, double *A, int *lda,
              double *lambda_re, double *lambda_im,
              double *vl, int *ldvl,
              double *vr, int *ldvr,
              int *ilo, int *ihi,
              double *scale, double *abnrm,
              double *rconde, double *rcondv,
              double *work, int *lwork,
              int *iwork, int *info);
}

namespace treekin {

void default_dgesv(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info)
{
        dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}

void default_dsyevx( char* jobz, char* range, char* uplo, int* n, double* a,
                int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
                int* m, double* w, double* z, int* ldz, double* work, int* lwork,
                int* iwork, int* ifail, int* info )
{
        dsyevx_(jobz, range, uplo, n, a,
                lda, vl, vu, il, iu, abstol,
                m, w, z, ldz, work, lwork,
                iwork, ifail, info );
}

void default_dgeevx(char *balanc, char *jobvl, char *jobvr,
              char *sense,
              int *n, double *A, int *lda,
              double *lambda_re, double *lambda_im,
              double *vl, int *ldvl,
              double *vr, int *ldvr,
              int *ilo, int *ihi,
              double *scale, double *abnrm,
              double *rconde, double *rcondv,
              double *work, int *lwork,
              int *iwork, int *info)
{
        dgeevx_ (balanc, jobvl, jobvr,
              sense,
              n, A, lda,
              lambda_re, lambda_im,
              vl, ldvl,
              vr, ldvr,
              ilo, ihi,
              scale, abnrm,
              rconde, rcondv,
              work, lwork,
              iwork, info);
}

}



#ifdef WITH_MPACK_GMP
#include <gmp.h>

/**
 * Takes a symmetric matrix U as input and computes the eigenvalues and eigenvectors by using the MPACK library.
 * @param U - INPUT - the symmetric matrix
 * @param dim - INPUT - the number of rows or columns
 * @param evals - OUTPUT - the eigenvalues
 * @param evecs - OUTPUT - the eigenvectors
 * @param precision - INPUT - the precision is the number of bits for the values of the matrix
 */
void
Calccpp::MxEV_Mpack_Sym_gmp(const mpf_class *U,
                            int             dim,
                            mpf_class       *evals,
                            mpf_class       *evecs,
                            int             precision)
{
  mpackint  n = dim;
  mpackint  lwork, info;

  /*initialization of GMP */
  int       default_prec = precision;

  mpf_set_default_prec(default_prec);


  mpf_class *A = new mpf_class[n * n];
  /*setting A matrix */
  std::copy_n(U, n * n, A);

  /*work space query */
  lwork = -1;
  mpf_class *work = new mpf_class[1];

  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);
  lwork = (int)mpf_get_d(work[0].get_mpf_t());
  delete[]work;
  work = new mpf_class[std::max((mpackint)1, lwork)];
  /*inverse matrix */
  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);

  /*copy evecs */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      evecs[j + i * n] = A[i + j * n]; /*reverse order */

  delete[]work;
  delete[]A;
}


#endif

#ifdef WITH_MPACK_QD
void
Calccpp::MxEV_Mpack_Sym_qd(const qd_real  *U,
                           int            dim,
                           qd_real        *evals,
                           qd_real        *evecs)
{
  mpackint  n = dim;
  mpackint  lwork, info;


  qd_real   *A = new qd_real[n * n];

  /*setting A matrix */
  std::copy_n(U, n * n, A);

  /*work space query */
  lwork = -1;
  qd_real *work = new qd_real[1];

  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);
  lwork = (int)(work[0].x[0]);
  delete[]work;
  work = new qd_real[std::max((mpackint)1, lwork)];
  /*inverse matrix */
  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);

  /*copy evecs */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      evecs[j + i * n] = A[i + j * n]; /*reverse order */

  delete[]work;
  delete[]A;
}


#endif

#ifdef WITH_MPACK_DD
void
Calccpp::MxEV_Mpack_Sym_dd(const dd_real  *U,
                           int            dim,
                           dd_real        *evals,
                           dd_real        *evecs)
{
  mpackint  n = dim;
  mpackint  lwork, info;


  dd_real   *A = new dd_real[n * n];

  /*setting A matrix */
  std::copy_n(U, n * n, A);


  /*work space query */
  lwork = -1;
  dd_real *work = new dd_real[1];

  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);
  lwork = (int)(work[0].x[0]);
  delete[]work;
  work = new dd_real[std::max((mpackint)1, lwork)];
  /*inverse matrix */
  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);

  /*copy evecs */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      evecs[j + i * n] = A[i + j * n].x[0]; /*reverse order */

  delete[]work;
  delete[]A;
}


#endif

#ifdef WITH_MPACK_MPFR
void
Calccpp::MxEV_Mpack_Sym_mpfr(const mpreal *U,
                             int          dim,
                             mpreal       *evals,
                             mpreal       *evecs,
                             int          precision)
{
  mpackint  n = dim;
  mpackint  lwork, info;

  mpreal::set_default_prec(precision);


  mpreal    *A = new mpreal[n * n];
  /*setting A matrix */
  std::copy_n(U, n * n, A);


  /*work space query */
  lwork = -1;
  mpreal *work = new mpreal[1];

  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);
  lwork = (int)((double)(work[0]));
  delete[]work;
  work = new mpreal[std::max((mpackint)1, lwork)];
  /*inverse matrix */
  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);

  /*copy evecs */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      evecs[j + i * n] = A[i + j * n]; /*reverse order */

  delete[]work;
  delete[]A;
}


#endif

#ifdef WITH_MPACK___FLOAT128
void
Calccpp::MxEV_Mpack_Sym_float128(const __float128 *U,
                                 int              dim,
                                 __float128       *evals,
                                 __float128       *evecs)
{
  mpackint    n = dim;
  mpackint    lwork, info;


  __float128  *A = new __float128[n * n];

  /*setting A matrix */
  std::copy_n(U, n * n, A);

  /*work space query */
  lwork = -1;
  __float128 *work = new __float128[1];

  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);
  lwork = (int)(work[0]);
  delete[]work;
  work = new __float128[std::max((mpackint)1, lwork)];
  /*inverse matrix */
  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);

  /*copy evecs */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      evecs[j + i * n] = A[i + j * n]; /*reverse order */

  delete[]work;
  delete[]A;
}


#endif

#ifdef WITH_MPACK_LD
void
Calccpp::MxEV_Mpack_Sym_longdouble(const long double  *U,
                                   int                dim,
                                   long double        *evals,
                                   long double        *evecs)
{
  mpackint    n = dim;
  mpackint    lwork, info;


  long double *A = new long double[n * n];

  /*setting A matrix */
  std::copy_n(U, n * n, A);

  /*work space query */
  lwork = -1;
  long double *work = new long double[1];

  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);
  lwork = (int)(work[0]);
  delete[]work;
  work = new long double[std::max((mpackint)1, lwork)];
  /*inverse matrix */
  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);

  /*copy evecs */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      evecs[j + i * n] = A[i + j * n]; /*reverse order */

  delete[]work;
  delete[]A;
}


#endif

#ifdef WITH_MPACK_DOUBLE
void
Calccpp::MxEV_Mpack_Sym_double(const double *U,
                               int          dim,
                               double       *evals,
                               double       *evecs)
{
  mpackint  n = dim;
  mpackint  lwork, info;


  double    *A = new double[n * n];

  /*setting A matrix */
  std::copy_n(U, n * n, A);

  /*work space query */
  lwork = -1;
  double *work = new double[1];

  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);
  lwork = (int)(work[0]);
  delete[]work;
  work = new double[std::max((mpackint)1, lwork)];
  /*inverse matrix */
  Rsyev("V", "U", n, A, n, evals, work, lwork, &info);

  /*copy evecs */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      evecs[j + i * n] = A[i + j * n]; /*reverse order */

  delete[]work;
  delete[]A;
}


#endif
