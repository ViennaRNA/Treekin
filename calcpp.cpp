#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <set>
#include <queue>

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

vector<int> reorganize; // reorganize array (so if LM 0 1 3 were reachable and 2 not, reorganize will contain r[0]=0 r[1]=1 r[2]=3), so r[x] = old position of x
int last_dim;

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
  // do the printing:
  double check = 0.0;
  printf("%e ", time);
  for (int i=0; i<lmins; i++) {
    if(ptFULL[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, ptFULL[i]);
      if (opt.num_err == 'H') exit(EXIT_FAILURE);
      else if (opt.num_err == 'R') ptFULL[i] = 0.0;
    }
    /* map individual structure -> gradient basins */
    check += fabs(ptFULL[i]);
  }

  for (int i=0; i<lmins; i++) {
    if (opt.num_err == 'R') printf("%e ", fabs(ptFULL[i])/check);
    else printf("%e ", fabs(ptFULL[i]));
  }
  printf("\n");

  // check for overall propability
  if ( ((check-1) < -0.05) || ((check-1) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
    if (opt.num_err == 'H') exit(EXIT_FAILURE);
  }

  return check;
}

double PrintProbNR(double *line, int dim, double time)
{
  double check = 0.0;
  printf("%e ", time);

  for (int i=0; i<dim; i++) {
    if(line[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, line[i]);
      if (opt.num_err == 'H') exit(EXIT_FAILURE);
      else if (opt.num_err == 'R') line[i] = 0.0;
    }
    /* map individual structure -> gradient basins */
    check += fabs(line[i]);
  }

  for (int i=0; i<dim; i++) {
    if (opt.num_err == 'R') printf("%e ", fabs(line[i])/check);
    else printf("%e ", fabs(line[i]));
  }

  printf("\n");

  // check for overall propability
  if ( ((check-1) < -0.05) || ((check-1) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
    if (opt.num_err == 'H') exit(EXIT_FAILURE);
  }

  return check;
}

double PrintProb(double *line, int dim, double time)
{
  double check = 0.0;
  printf("%e ", time);

  for (int i=0; i<dim; i++) {
    if (line[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, line[i]);
      if (opt.num_err == 'H') exit(EXIT_FAILURE);
      else if (opt.num_err == 'R') line[i] = 0.0;
    }
    check += fabs(line[i]);
  }


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

  // check for overall probability
  if ( ((check-1) < -0.05) || ((check-1) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
    if (opt.num_err == 'H') exit(EXIT_FAILURE);
  }

  return check;
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
      if (i!=j) U[i*dim+j] = pow(U[i*dim+j], factor);
    }
  }

  // fix the diagonal:
  for (int i = 0; i < dim; i++) U[dim*i+i] = 0;
  for (int i = 0; i < dim; i++) {
    double tmp = 0.00;
    // calculate column sum
    for(int j = 0; j < dim; j++)  tmp += U[dim*j+i];
    U[dim*i+i] = -tmp+1.0;   // make U a stochastic matrix U = Q+I ??
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
