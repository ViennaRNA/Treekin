

#include "calcpp.h"
#include <gmp.h>
/*
#ifdef WITH_MPACK
  #ifdef WITH_MPACK_GMP
      #include <mpack/mblas_gmp.h>
      #include <mpack/mlapack_gmp.h>
  #endif
  #ifdef WITH_MPACK_QD
      #include <mpack/mblas_qd.h>
      #include <mpack/mlapack_qd.h>
  #endif
  #ifdef WITH_MPACK_DD
      #include <mpack/mblas_dd.h>
      #include <mpack/mlapack_dd.h>
  #endif
  #ifdef WITH_MPACK_MPFR
      #include <mpack/mblas_mpfr.h>
      #include <mpack/mlapack_mpfr.h>
  #endif
  #ifdef WITH_MPACK___FLOAT128
      #include <mpack/mblas___float128.h>
      #include <mpack/mlapack___float128.h>
  #endif
  #ifdef WITH_MPACK_DOUBLE
      #include <mpack/mblas_double.h>
      #include <mpack/mlapack_double.h>
  #endif
  #ifdef WITH_MPACK_LD
      #include <mpack/mblas_longdouble.h>
      #include <mpack/mlapack_longdouble.h>
  #endif
#endif
*/


#ifdef WITH_MPACK_GMP
/**
 * Takes a symmetric matrix U as input and computes the eigenvalues and eigenvectors by using the MPACK library.
 * @param U - INPUT - the symmetric matrix
 * @param dim - INPUT - the number of rows or columns
 * @param evals - OUTPUT - the eigenvalues
 * @param evecs - OUTPUT - the eigenvectors
 * @param precision - INPUT - the precision is the number of bits for the values of the matrix
 */
void
Calccpp::MxEV_Mpack_Sym_gmp(const mpf_class *U, int dim, mpf_class *evals, mpf_class *evecs, int precision){
  mpackint n = dim;
  mpackint lwork, info;

 //initialization of GMP

  int default_prec = precision;
  mpf_set_default_prec(default_prec);

/*
  mpf_class eps;
      mpf_class one;
      unsigned long exp2;
      one = 1.0;
      exp2 = mpf_get_prec(one.get_mpf_t());
      mpf_div_2exp(eps.get_mpf_t(), one.get_mpf_t(), exp2);
  double blub = mpf_get_d(eps.get_mpf_t());
  printf("blub %10.4g",blub);
*/


  mpf_class *A = new mpf_class[n * n];
  mpf_class *w = new mpf_class[n];
    //setting A matrix
  for(int i =0; i < n; i++){
    for(int j=0; j < n; j++){
      A[i+j*n]=U[i+j*n];
    }
  }

  //work space query
  lwork = -1;
  mpf_class *work = new mpf_class[1];

  Rsyev("V", "U", n, A, n, w, work, lwork, &info);
  lwork = (int) mpf_get_d(work[0].get_mpf_t());
  delete[]work;
  work = new mpf_class[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, A, n, w, work, lwork, &info);


  //copy eigenvalues
  mpf_ptr tmp_ptr;
  double tmp_d;
  for(int i= 0; i < n; i++){
    //tmp_ptr= w[i].get_mpf_t();
    //tmp_d = mpf_get_d(tmp_ptr);
    evals[i] = (mpf_class)w[i];//tmp_d;
  }
  //copy evecs
  for(int i =0; i < n; i++){
      for(int j=0; j < n; j++){
        //tmp_ptr = A[i+j*n].get_mpf_t();
        //tmp_d = mpf_get_d(tmp_ptr);
        evecs[j+i*n]= (mpf_class)A[i+j*n] ;//tmp_d; //reverse order
      }
    }

  delete[]work;
  delete[]w;
  delete[]A;
}
#endif

#ifdef WITH_MPACK_QD
void
Calccpp::MxEV_Mpack_Sym_qd(const qd_real *U, int dim, qd_real *evals, qd_real *evecs){

  mpackint n = dim;
  mpackint lwork, info;


  qd_real *A = new qd_real[n * n];
  qd_real *w = new qd_real[n];
    //setting A matrix
  for(int i =0; i < n; i++){
    for(int j=0; j < n; j++){
      A[i+j*n]=U[i+j*n];
    }
  }

  //work space query
  lwork = -1;
  qd_real *work = new qd_real[1];

  Rsyev("V", "U", n, (qd_real*)A, n, w, work, lwork, &info);
  lwork = (int) (work[0].x[0]);
  delete[]work;
  work = new qd_real[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, (qd_real*)A, n, w, work, lwork, &info);


  //copy eigenvalues
  for(int i= 0; i < n; i++){
    evals[i] = w[i].x[0];
  }
  //copy evecs
  for(int i =0; i < n; i++){
      for(int j=0; j < n; j++){
        evecs[j+i*n]=A[i+j*n].x[0]; //reverse order
      }
    }

  delete[]work;
  delete[]w;
  delete[]A;
}
#endif

#ifdef WITH_MPACK_DD
void
Calccpp::MxEV_Mpack_Sym_dd(const dd_real *U, int dim, dd_real *evals, dd_real *evecs){

  mpackint n = dim;
  mpackint lwork, info;


  dd_real *A = new dd_real[n * n];
  dd_real *w = new dd_real[n];
    //setting A matrix
  for(int i =0; i < n; i++){
    for(int j=0; j < n; j++){
      A[i+j*n]=U[i+j*n];
    }
  }


  //work space query
  lwork = -1;
  dd_real *work = new dd_real[1];

  Rsyev("V", "U", n, (dd_real*)A, n, w, work, lwork, &info);
  lwork = (int) (work[0].x[0]);
  delete[]work;
  work = new dd_real[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n,(dd_real*) A, n, w, work, lwork, &info);


  //copy eigenvalues
  for(int i= 0; i < n; i++){
    evals[i] = w[i].x[0];
  }
  //copy evecs
  for(int i =0; i < n; i++){
      for(int j=0; j < n; j++){
        evecs[j+i*n]=A[i+j*n].x[0]; //reverse order
      }
    }


  delete[]work;
  delete[]w;
  delete[]A;

}
#endif

#ifdef WITH_MPACK_MPFR
void
Calccpp::MxEV_Mpack_Sym_mpfr(const mpreal *U, int dim, mpreal *evals, mpreal *evecs, int precision){

  mpackint n = dim;
  mpackint lwork, info;

  mpreal::set_default_prec(precision);


  mpreal *A = new mpreal[n * n];
  mpreal *w = new mpreal[n];
    //setting A matrix
  for(int i =0; i < n; i++){
    for(int j=0; j < n; j++){
      A[i+j*n]=U[i+j*n];
    }
  }


  //work space query
  lwork = -1;
  mpreal *work = new mpreal[1];

  Rsyev("V", "U", n, (mpreal*)A, n, w, work, lwork, &info);
  lwork = (int) ((double)(work[0]));
  delete[]work;
  work = new mpreal[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, (mpreal*)A, n, w, work, lwork, &info);


  //copy eigenvalues
  for(int i= 0; i < n; i++){
    evals[i] = (double)w[i];
  }
  //copy evecs
  for(int i =0; i < n; i++){
      for(int j=0; j < n; j++){
        evecs[j+i*n]=(double)A[i+j*n]; //reverse order
      }
    }


  delete[]work;
  delete[]w;
  delete[]A;

}
#endif

#ifdef WITH_MPACK___FLOAT128
void
Calccpp::MxEV_Mpack_Sym_float128(const __float128 *U, int dim, __float128 *evals, __float128 *evecs){

  mpackint n = dim;
  mpackint lwork, info;


  __float128 *A = new __float128[n * n];
  __float128 *w = new __float128[n];
    //setting A matrix
  for(int i =0; i < n; i++){
    for(int j=0; j < n; j++){
      A[i+j*n]=U[i+j*n];
    }
  }


  //work space query
  lwork = -1;
  __float128 *work = new __float128[1];

  Rsyev("V", "U", n, (__float128*)A, n, w, work, lwork, &info);
  lwork = (int) (work[0]);
  delete[]work;
  work = new __float128[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, (__float128*)A, n, w, work, lwork, &info);


  //copy eigenvalues
  for(int i= 0; i < n; i++){
    evals[i] = w[i];
  }
  //copy evecs
  for(int i =0; i < n; i++){
      for(int j=0; j < n; j++){
        evecs[j+i*n]=A[i+j*n]; //reverse order
      }
    }


  delete[]work;
  delete[]w;
  delete[]A;
}
#endif

#ifdef WITH_MPACK_LD
void
Calccpp::MxEV_Mpack_Sym_longdouble(const long double *U, int dim, long double *evals, long double *evecs){

  mpackint n = dim;
  mpackint lwork, info;


  long double *A = new long double[n * n];
  long double *w = new long double[n];
    //setting A matrix
  for(int i =0; i < n; i++){
    for(int j=0; j < n; j++){
      A[i+j*n]=U[i+j*n];
    }
  }

  //work space query
  lwork = -1;
  long double *work = new long double[1];

  Rsyev("V", "U", n, (long double*)A, n, w, work, lwork, &info);
  lwork = (int) (work[0]);
  delete[]work;
  work = new long double[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, (long double*)A, n, w, work, lwork, &info);


  //copy eigenvalues
  for(int i= 0; i < n; i++){
    evals[i] = w[i];
  }
  //copy evecs
  for(int i =0; i < n; i++){
      for(int j=0; j < n; j++){
        evecs[j+i*n]=A[i+j*n]; //reverse order
      }
    }


  delete[]work;
  delete[]w;
  delete[]A;
}
#endif

#ifdef WITH_MPACK_DOUBLE
void
Calccpp::MxEV_Mpack_Sym_double(const double *U, int dim, double *evals, double *evecs){

  mpackint n = dim;
  mpackint lwork, info;


  double *A = new double[n * n];
  double *w = new double[n];
    //setting A matrix
  for(int i =0; i < n; i++){
    for(int j=0; j < n; j++){
      A[i+j*n]=U[i+j*n];
    }
  }


  //work space query
  lwork = -1;
  double *work = new double[1];

  Rsyev("V", "U", n, (double*)A, n, w, work, lwork, &info);
  lwork = (int) (work[0]);
  delete[]work;
  work = new double[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n,(double*) A, n, w, work, lwork, &info);


  //copy eigenvalues
  for(int i= 0; i < n; i++){
    evals[i] = w[i];
  }

  //copy evecs
  for(int i =0; i < n; i++){
      for(int j=0; j < n; j++){
        evecs[j+i*n]=A[i+j*n]; //reverse order
      }
    }


  delete[]work;
  delete[]w;
  delete[]A;
}





#endif

