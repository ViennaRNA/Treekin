#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <set>
#include <queue>

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

extern "C" {
  #include "calc.h"
  //#include "expokit_wrappers.h"
}

//#include "expokit_wrappers.h"

using namespace std;

extern "C" void MxEgro(double **U, double **p0, int dim);
extern "C" double PrintProb(double *line, int dim, double time);
extern "C" double PrintProbNR(double *line, int dim, double time);
extern "C" double PrintProbFull(double *line, int dim, double time, int lmins);
extern "C" int ConvergenceReached(double *p8, double *pt, int dim, int full);
//extern "C" void TestExpokit(double *R, int dim, double *p0, double t_start, double t_end, double t_inc);
//extern "C" void TestEigen(double *R, int dim, double **evals, double **evecs);
extern "C" void MxRescale(double *U, int dim, double desired_rate);
extern "C" void MxRescaleH(double *U, int dim, double hard_rescale);
extern "C" void MxTimes(double *U, int dim, double times);
extern "C" int *MxErgoEigen(double *U, int dim);

//extern "C" void printmat_mpf(int N, int M, mpf_class * A, int LDA);
extern "C" void MxEV_Mpack_Sym_gmp(const double *U, int dim, double *evals, double *evecs, int precision);
extern "C" void MxEV_Mpack_Sym_qd(const double *U, int dim, double *evals, double *evecs);
extern "C" void MxEV_Mpack_Sym_dd(const double *U, int dim, double *evals, double *evecs);
extern "C" void MxEV_Mpack_Sym_mpfr(const double *U, int dim, double *evals, double *evecs, int precision);
extern "C" void MxEV_Mpack_Sym_float128(const double *U, int dim, double *evals, double *evecs);
extern "C" void MxEV_Mpack_Sym_longdouble(const double *U, int dim, double *evals, double *evecs);
extern "C" void MxEV_Mpack_Sym_double(const double *U, int dim, double *evals, double *evecs);


vector<int> reorganize; // reorganize array (so if LM 0 1 3 were reachable and 2 not, reorganize will contain r[0]=0 r[1]=1 r[2]=3), so r[x] = old position of x
int last_dim;

int *MxErgoEigen(double *U, int dim)
{
  int count = 1;
  vector<bool> reached(dim, false);
  reached[0] = true;

  queue<int> que_ergo;
  que_ergo.push(0);

  // fill the sets of non-empty states
  while(!que_ergo.empty() && count<dim) {
    int to_do = que_ergo.front();
    que_ergo.pop();

    // collect contingency
    for (int i=0; i<dim; i++) {
      if (i==to_do) continue;
      if (U[i*dim + to_do]!=0.0 && U[to_do*dim + i]!=0.0 && !reached[i]) {
        reached[i] = true;
        que_ergo.push(i);
        count++;
      }
    }
  }

  // return the array of unreached states:
  if (dim == count) return NULL;
  int *result = (int*) malloc(dim*sizeof(int));
  for (int i=0; i<dim; i++) result[i] = (reached[i]?1:0);
  return result;
}

void MxEgro(double **Up, double **p0p, int dim)
{
  // aliases
  double *U = *Up;
  double *p0 = *p0p;

  // lokalize first non-empty state
  int first = 0;
  while (p0[first]==0.0) first++;

  // make set of all non-empty states
  set<int> set_ergo;
  queue<int> que_ergo;

  set_ergo.insert(first);
  que_ergo.push(first);

  // fill the sets of non-empty states
  while(!que_ergo.empty() && (int)set_ergo.size()<dim) {
    int to_do = que_ergo.front();
    que_ergo.pop();

    // collect contingency
    for (int i=0; i<dim; i++) {
      if (i==to_do) continue;
      if (U[i*dim + to_do]>0.0 && (int)set_ergo.count(i)==0) {
        set_ergo.insert(i);
        que_ergo.push(i);
      }
    }
  }

  // check ergodicity
  if ((int)set_ergo.size()==dim) return; // all ok
  else {
    int i=first+1;
    while (i<dim) {
      if (p0[i]>0.0 && set_ergo.count(i)==0) { // disconnected :(
        fprintf(stderr, "ERROR: Matrix is non-ergodic and initial populations are disconected!! Exiting...\n");
        exit(-1);
      }
      i++;
    }

    // fill helper reorganize array
    reorganize.resize(set_ergo.size());
    i=0;
    for (set<int>::iterator it=set_ergo.begin(); it!=set_ergo.end(); it++) {
      reorganize[i++]=*it;
    }

    // reorganize matrix
    for (int i=0; i<(int)set_ergo.size(); i++) {
      for (int j=0; j<(int)set_ergo.size(); j++) {
        U[i*set_ergo.size()+j]=U[reorganize[i]*dim+reorganize[j]];
      }
    }
    last_dim = dim;
    dim = set_ergo.size();
    *Up = (double*)realloc(U, dim*dim*sizeof(double));

    // reorganize p0
    for (int i=0; i<(int)set_ergo.size(); i++) {
      p0[i]=p0[reorganize[i]];
    }
    *p0p = (double*)realloc(p0, dim*sizeof(double));

    // set dimension to global
    MxInit(dim);
    if (!opt.quiet) fprintf(stderr, "WARNING: Matrix is non-ergodic! Decreasing dimension to %d.\n", dim);
    //MxPrint(U, "Ergodic U", 'm');

    if (dim == 1) {
      PrintDummy(p0);
      exit(EXIT_SUCCESS);
    }
  }
}

double PrintProbFull(double *line, int dim, double time, int lmins)
{
  // for full process:
  vector<double> ptFULL (lmins, 0.0);
  if (reorganize.size() == 0) {
    for (int i=0; i<dim; i++) {
      ptFULL[E[i].ag] += line[i];
    }
  } else {
    for (int i=0; i<dim; i++) {
      ptFULL[E[reorganize[i]].ag] += line[i];
    }
  }

  // sum first
  double check = 0.0;
  for (int i=0; i<lmins; i++) {
    if(ptFULL[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, ptFULL[i]);
      if (opt.num_err == 'H') exit(EXIT_FAILURE);
      else if (opt.num_err == 'R') ptFULL[i] = 0.0;
    }
    /* map individual structure -> gradient basins */
    check += fabs(ptFULL[i]);
  }

  // check for overall propability
  if ( ((check-1) < -0.05) || ((check-1) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1.0 %s!\n", time, check, (opt.num_err == 'R'?"rescaling":"exiting") );
    if (opt.num_err == 'H' || check == 0.0) exit(EXIT_FAILURE);
  }

  // print:
  printf("%e ", time);
  for (int i=0; i<lmins; i++) {
    if (opt.num_err == 'R') printf("%e ", fabs(ptFULL[i])/check);
    else printf("%e ", fabs(ptFULL[i]));
  }
  printf("\n");

  return check;
}

double PrintProbNR(double *line, int dim, double time)
{
  double check = 0.0;

  // summ first
  for (int i=0; i<dim; i++) {
    if(line[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, line[i]);
      if (opt.num_err == 'H') exit(EXIT_FAILURE);
      else if (opt.num_err == 'R') line[i] = 0.0;
    }
    /* map individual structure -> gradient basins */
    check += fabs(line[i]);
  }

  // check for overall propability
  if ( ((check-1) < -0.05) || ((check-1) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1.0 %s!\n", time, check, (opt.num_err == 'R'?"rescaling":"exiting") );
    if (opt.num_err == 'H' || check == 0.0) exit(EXIT_FAILURE);
  }

  // print
  printf("%e ", time);
  for (int i=0; i<dim; i++) {
    if (opt.num_err == 'R') printf("%e ", fabs(line[i])/check);
    else printf("%e ", fabs(line[i]));
  }

  printf("\n");

  return check;
}

double PrintProb(double *line, int dim, double time)
{
  double check = 0.0;

  // sum it up
  for (int i=0; i<dim; i++) {
    if (line[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, line[i]);
      if (opt.num_err == 'H') exit(EXIT_FAILURE);
      else if (opt.num_err == 'R') line[i] = 0.0;
    }
    check += fabs(line[i]);
  }

  // check for overall probability
  if ( ((check-1) < -0.05) || ((check-1) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1.0 %s!\n", time, check, (opt.num_err == 'R' && check != 0.0?"rescaling":"exiting") );
    if (opt.num_err == 'H' || check == 0.0) exit(EXIT_FAILURE);
  }

  // print finally:
  printf("%e ", time);
  if (reorganize.size()==0) {
    for (int i=0; i<dim; i++) {
      /* map individual structure -> gradient basins */
      if (opt.num_err == 'R') printf("%e ", fabs(line[i])/check);
      else printf("%e ", fabs(line[i]));
    }
  } else {
    int j=0;
    for (int i=0; i<last_dim; i++) {
      if (j>(int)reorganize.size() || reorganize[j]!=i) {
        printf("%e ", 0.0);
      } else {

        if (opt.num_err == 'R') printf("%e ", fabs(line[j])/check);
        else printf("%e ", fabs(line[j]));
        j++;
      }
    }
  }
  printf("\n");

  return check;
}

/**
 * Takes a symmetric matrix U as input and computes the eigenvalues and eigenvectors by using the MPACK library.
 * @param U - INPUT - the symmetric matrix
 * @param dim - INPUT - the number of rows or columns
 * @param evals - OUTPUT - the eigenvalues
 * @param evecs - OUTPUT - the eigenvectors
 * @param precision - INPUT - the precision is the number of bits for the values of the matrix
 */
void
MxEV_Mpack_Sym_gmp(const double *U, int dim, double *evals, double *evecs, int precision){
#ifdef WITH_MPACK_GMP
  mpackint n = dim;
  mpackint lwork, info;

 //initialization of GMP
  int default_prec = precision;
  mpf_set_default_prec(default_prec);

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
  for(int i= 0; i < n; i++){
    evals[i] = mpf_get_d(w[i].get_mpf_t());
  }
  //copy evecs
  for(int i =0; i < n; i++){
      for(int j=0; j < n; j++){
        evecs[j+i*n]=mpf_get_d(A[i+j*n].get_mpf_t()); //reverse order
      }
    }

  delete[]work;
  delete[]w;
  delete[]A;
#endif
}

void
MxEV_Mpack_Sym_qd(const double *U, int dim, double *evals, double *evecs){
#ifdef WITH_MPACK_QD
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

  Rsyev("V", "U", n, A, n, w, work, lwork, &info);
  lwork = (int) (work[0].x[0]);
  delete[]work;
  work = new qd_real[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, A, n, w, work, lwork, &info);

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
#endif
}


void
MxEV_Mpack_Sym_dd(const double *U, int dim, double *evals, double *evecs){
#ifdef WITH_MPACK_DD
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

  Rsyev("V", "U", n, A, n, w, work, lwork, &info);
  lwork = (int) (work[0].x[0]);
  delete[]work;
  work = new dd_real[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, A, n, w, work, lwork, &info);

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
#endif
}


void
MxEV_Mpack_Sym_mpfr(const double *U, int dim, double *evals, double *evecs, int precision){
#ifdef WITH_MPACK_MPFR
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

  Rsyev("V", "U", n, A, n, w, work, lwork, &info);
  lwork = (int) ((double)(work[0]));
  delete[]work;
  work = new mpreal[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, A, n, w, work, lwork, &info);

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
#endif
}

void
MxEV_Mpack_Sym_float128(const double *U, int dim, double *evals, double *evecs){
#ifdef WITH_MPACK___FLOAT128
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

  Rsyev("V", "U", n, A, n, w, work, lwork, &info);
  lwork = (int) (work[0]);
  delete[]work;
  work = new __float128[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, A, n, w, work, lwork, &info);

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
#endif
}

void
MxEV_Mpack_Sym_longdouble(const double *U, int dim, double *evals, double *evecs){
#ifdef WITH_MPACK_LD
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

  Rsyev("V", "U", n, A, n, w, work, lwork, &info);
  lwork = (int) (work[0]);
  delete[]work;
  work = new long double[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, A, n, w, work, lwork, &info);

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
#endif
}

void
MxEV_Mpack_Sym_double(const double *U, int dim, double *evals, double *evecs){
#ifdef WITH_MPACK_DOUBLE
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

  Rsyev("V", "U", n, A, n, w, work, lwork, &info);
  lwork = (int) (work[0]);
  delete[]work;
  work = new double[std::max((mpackint) 1, lwork)];
  //inverse matrix
  Rsyev("V", "U", n, A, n, w, work, lwork, &info);

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
#endif
}


int ConvergenceReached(double *p8, double *pt, int dim, int full) {

  int pdiff_counter = 0;
  full = (full?1:0);
  /* now check if we have converged yet */
  for(int i = 0; i < dim; i++) {
    if (fabs(p8[i] - pt[i]) >= 0.000001) {
      pdiff_counter++;
      break;
    }
  }

  if (pdiff_counter < 1) /* all mins' pdiff lies within threshold */
    return 1;
  pdiff_counter = 0;
  /* end check of convergence */
  return false;
}

void MxRescale(double *U, int dim, double desired_rate)
{
  //first get the minimal rate
  double minimal_rate = 1.0;
  for (int i=0; i<dim*dim; i++) {
    if (i%dim != i/dim && minimal_rate > U[i] && U[i]>0.0) minimal_rate = U[i];
  }

  // calculate rescale factor:
  double factor = log(desired_rate)/log(minimal_rate);
  fprintf(stderr, "rescale params: %20.10g %20.10g %20.10g\n", minimal_rate, desired_rate, factor);

  if (desired_rate < minimal_rate) return ;

  MxRescaleH(U, dim, factor);
}

void MxRescaleH(double *U, int dim, double hard_rescale)
{

  double factor = hard_rescale;

  // rescale!
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      if (i!=j) U[i*dim+j] = std::pow(U[i*dim+j], factor);
    }
  }

  // fix the diagonal:
  for (int i = 0; i < dim; i++) U[dim*i+i] = 0;
  for (int i = 0; i < dim; i++) {
    double tmp = 0.00;
    // calculate column sum
    for(int j = 0; j < dim; j++)  tmp += U[dim*j+i];
    U[dim*i+i] = -tmp+(opt.useplusI?1.0:0.0);   // make U a stochastic matrix U = Q+I ??
  }
}

void MxTimes(double *U, int dim, double times)
{
  // multiply!
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      if (i!=j) U[i*dim+j] *= times;
    }
  }

  // fix the diagonal:
  for (int i = 0; i < dim; i++) U[dim*i+i] = 0;
  for (int i = 0; i < dim; i++) {
    double tmp = 0.00;
    // calculate column sum
    for(int j = 0; j < dim; j++)  tmp += U[dim*j+i];
    U[dim*i+i] = -tmp+(opt.useplusI?1.0:0.0);   // make U a stochastic matrix U = Q+I ??
  }
}

/*#include <iostream>
#include <Eigen/Eigenvalues>

using Eigen::MatrixXd;
using Eigen::EigenSolver;

void TestEigen(double *R, int dim, double **evals, double **evecs) {

  MatrixXd A = MatrixXd(dim, dim);

  for (int i=0; i<dim; i++) {
  for (int j=0; j<dim; j++) {
    A(i,j) = R[i*dim+j];
    fprintf(stderr, "%30.20g ", R[i*dim+j]);
  }
    fprintf(stderr, "\n");
  }

  EigenSolver<MatrixXd> es(A);

  MatrixXd D = es.pseudoEigenvalueMatrix();
  MatrixXd V = es.pseudoEigenvectors();
  cerr << "The input matrix: " << endl << A << endl;
  cerr << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
  cerr << "The pseudo-eigenvector matrix V is:" << endl << V << endl;
  cerr << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;

}*/

/*void TestExpokit(double *R, int dim, double *p0, double t_start, double t_end, double t_inc)
{

  int n = dim;
  int m = n-1;
  double w[n];
  double tol = 0.01; // change to something better maybe?
  double anorm[n*n]; //??
  int lwsp = n*(m+2)+5*(m+2)*(m+2)+7;
  int liwsp = max(m+2, 7);
  int iwsp[liwsp];
  int itrace = 1, iflag = 1;
  double wsp[lwsp];
  double res[n*n];

  // change the R into ia/ja matrix:

    // now get the non-zero ones
  vector<double> non_zero;
  vector<int> ia_vec;
  vector<int> ja_vec;
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      if (fabs(R[i*dim+j]-(i==j?1.0:0.0)) > 0.0) {
        non_zero.push_back(R[i*dim+j]-(i==j?1.0:0.0));
        ia_vec.push_back(i);
        ja_vec.push_back(j);
      }
    }
  }

    // convert them into C code
  int ia[non_zero.size()];
  int ja[non_zero.size()];
  double a[non_zero.size()];
  int nz = non_zero.size();
  for (unsigned int i=0; i<non_zero.size(); i++) {
    ia[i] = ia_vec[i];
    ja[i] = ja_vec[i];
    a[i] = non_zero[i];
  }


//  const int nzc = 9;
//  int ia[nzc] = {1,1,1,2,2,2,3,3,3};
//  int ja[nzc] = {1,2,3,1,2,3,1,2,3};
//  double a[nzc] = {-0.9, 0.5, 0.5, 0.5, -1.0, 0.5, 0.4, 0.5, -1.0};
//  int nz = nzc;

  // main loop
  for (double time=t_start; time<=t_end; time *= t_inc) {
    wrapsingledgexpv_(&n, &m, &time, p0, w, &tol, anorm, wsp, &lwsp, iwsp, &liwsp, &itrace, &iflag, ia, ja, a, &nz, res);
    PrintProb(w, dim, time);
  }

}*/
