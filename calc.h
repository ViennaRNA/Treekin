/*=================================================================*/
/*=   calc.h                                                      =*/
/*=   header file for calculation routines for treekin            =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2017-07-27 18:30:26 ivo>          =*/
/*=   $Id: calc.h,v 1.16 2006/11/27 13:50:08 mtw Exp $            =*/
/*=   ---------------------------------------------------------   =*/
/*=     (c) Michael Thomas Wolfinger, W. Andreas Svrcek-Seiler    =*/
/*=                  {mtw,svrci}@tbi.univie.ac.at                 =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _CALC_H_
#define _CALC_H_

#include <cstdlib>
#include <cmath>
#include <algorithm> // fill_n
#include <set>
#include <queue>
#include <vector>
#include <errno.h>
#include <assert.h>

#include "mxccm.h"      /* functions for eigen-problems stolen from ccmath */
#include "globals.h"
#include "barparser.h"
#include "calcpp.h"
#include "bardata.h"
#include "exp_matrix.h"

template<typename T>
class Calc {
  int      dim;
  double   _kT;
  T  *evals;
  T  *evecs;
  T  *_sqrPI;   /* left constant array */
  T  *sqrPI_;   /* right constant array */
  T  *D;   /* matrix with degree of degeneracy */
  char     Aname[30];

  Globals *_globalParameters;
  treekin_options *_opt;
  int _lmins;
  SubInfo *_E;
  Mxccm mxccm;
  std::vector<int> reorganize; // reorganize array (so if LM 0 1 3 were reachable and 2 not, reorganize will contain r[0]=0 r[1]=1 r[2]=3), so r[x] = old position of x

  Calccpp *solver;

private:
  T *MxMethodeA (BarData *Data);
  T *MxMethodeFULL(T *);
  T *MxMethodeINPUT (BarData *Data, T *);
  T  max_saddle(int i, int j, BarData *Data);
  void    print_settings(void);
  char   *time_stamp(void);
  void    MxDoDegeneracyStuff(void);
  void    MxBinWrite(T *Mx, const char what[], char c);
  int     MxBinRead(T** Mx, const char what[], char c);
  void    MxASCIIWrite(T *Mx, const char *asciifile);
  void    MxASCIIWriteV(T *Mx, const char *asciifile);
  void    MxKotzOutMathematica(T *Mx);
  void    MxSortEig(T *evals, T *evecs);
  void    MxEVLapackSym(T *U);
  void    MxEVLapackNonSym(T *U);
  void    MxFixevecs(T *, T *);
  void    MxFixevecsAbsorb(T *, T *);
  void    MxDiagHelper(T *P8);
  void norm2(T *mx);
  int *MxErgoEigen(T *U, int dim);

public:
  Calc(Globals *globalParameters, SubInfo *E, int dim);
  ~Calc();

  /**
   * check for ergodicity
   */
  int MxEgro(T **U, T **p0, int dim);

  void PrintDummy(T *line);

  /**
   * print probabilities (adds zero (non-ergodic) columns)
   * returns sum of these probabilities
   */
  T PrintProb(T *line, int dim, double time);
  T PrintProbFull(T *line, int dim, double time, int lmins);
  T PrintProbNR(T *line, int dim, double time);

  void *MxNew (size_t size);
  void MxGetSpace(T **p8);
  T *MxBar2Matrix (BarData *Data, T *);
  void MxEqDistr (BarData *Data, T **p8);
  void MxEqDistrFULL (SubInfo *E, T *p8);
  void MxDiagonalize (T *U, T **_S, T *PI);
  void MxRecover(T **_S, T *P8);
  void MxStartVec (T **p0);
  int ConvergenceReached(T *p8, T *pt, int dim, int full);
  void MxIterate (T *p0, T *p8, T *S);
  void MxMemoryCleanUp (void);
  int MxExponent(T *p0, T *p8, T *U);

  void MxFPT(T *U, T *p8, FILE *out);

  void MxEqDistrFromLinSys(T *U, T **p8);
  void MxEqDistrFromDetailedBalance(T *U, T **p8);
  void MxFPTSimple(T *U);
  T *MxFPTOneState(T *U, int state);

  // not used! (just for past testing)
  T  MxFPTRandom(T *P, T *U, int src, int dst, int packets);
  void MxFPTRnd(T *U, int packets);

  void MxFPrint(T *mx, const char *name, char c, FILE *out, int pure);
  void MxPrint(T *mx, const char *name, char c);

};

template<typename T>
Calc<T>::Calc(Globals *globalParameters, SubInfo *E, int d){
   _globalParameters = globalParameters;
  _opt = &globalParameters->opt;
  _lmins = globalParameters->lmins;
  _E = E;
  _kT = 0.00198717*(273.15 + _opt->T);

  evals = NULL;
  evecs = NULL;
  _sqrPI = NULL;   /* left constant array */
  sqrPI_ = NULL;   /* right constant array */
  D = NULL;
  dim = d;
  solver = new Calccpp(_opt,_E);
}

template<typename T>
Calc<T>::~Calc(){
  delete solver;
}

template<typename T>
int Calc<T>::MxEgro(T **Up, T **p0p, int dim)
{
  // aliases
  T *U = *Up;
  T *p0 = *p0p;

  // lokalize first non-empty state
  int first = 0;
  while (p0[first]==(T)0.0) first++;

  // make set of all non-empty states
  std::set<int> set_ergo;
  std::queue<int> que_ergo;

  set_ergo.insert(first);
  que_ergo.push(first);

  // fill the sets of non-empty states
  while(!que_ergo.empty() && (int)set_ergo.size()<dim) {
    int to_do = que_ergo.front();
    que_ergo.pop();

    // collect contingency
    for (int i=0; i<dim; i++) {
      if (i==to_do) continue;
      if (U[i*dim + to_do]>(T)0.0 && (int)set_ergo.count(i)==0) {
        set_ergo.insert(i);
        que_ergo.push(i);
      }
    }
  }

  // check ergodicity
  if ((int)set_ergo.size()==dim) return dim; // all ok
  else {
    int i=first+1;
    while (i<dim) {
      if (p0[i]>(T)0.0 && set_ergo.count(i)==0) { // disconnected :(
        fprintf(stderr, "ERROR: Matrix is non-ergodic and initial populations are disconected!! Exiting...\n");
        exit(-1);
      }
      i++;
    }

    // fill helper reorganize array
    reorganize.resize(set_ergo.size());
    i=0;
    for (std::set<int>::iterator it=set_ergo.begin(); it!=set_ergo.end(); it++) {
      reorganize[i++]=*it;
    }

    // reorganize matrix
    for (int i=0; i<(int)set_ergo.size(); i++) {
      for (int j=0; j<(int)set_ergo.size(); j++) {
        U[i*set_ergo.size()+j]=U[reorganize[i]*dim+reorganize[j]];
      }
    }

    dim = set_ergo.size();
    //*Up = (T*)realloc(U, dim*dim*sizeof(T));
    *Up = new T[dim*dim];
    std::copy_n(U, dim*dim, *Up);
    delete[] U;
    U = *Up;

    // reorganize p0
    for (int i=0; i<(int)set_ergo.size(); i++) {
      p0[i]=p0[reorganize[i]];
    }
    //*p0p = (T*)realloc(p0, dim*sizeof(T));
    *p0p = new T[dim];
    std::copy_n(p0, dim, *p0p);
    delete[] p0;
    p0 = *p0p;

    if (!_opt->quiet) fprintf(stderr, "WARNING: Matrix is non-ergodic! Decreasing dimension to %d.\n", dim);
    //MxPrint(U, "Ergodic U", 'm');

    if (dim == 1) {
      PrintDummy(p0);
      return dim;
    }
  }
  // return new reduced dimension of the inputmatrix.
  return dim;
}

template<typename T>
void Calc<T>::PrintDummy(T *line)
{
  print_settings();
  PrintProb(line, 1, _opt->t0);
  PrintProb(line, 1, _opt->t8);
}

template<typename T>
T Calc<T>::PrintProbFull(T *line, int dim, double time, int lmins)
{
  // for full process:
  std::vector<T> ptFULL (lmins, 0.0);
  if (reorganize.size() == 0) {
    for (int i=0; i<dim; i++) {
      ptFULL[_E[i].ag] += line[i];
    }
  } else {
    for (int i=0; i<dim; i++) {
      ptFULL[_E[reorganize[i]].ag] += line[i];
    }
  }

  // sum first
  T check = 0.0;
  for (int i=0; i<lmins; i++) {
    if(ptFULL[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, (double)ptFULL[i]);
      if (_opt->num_err == 'H') exit(EXIT_FAILURE);
      else if (_opt->num_err == 'R') ptFULL[i] = 0.0;
    }
    /* map individual structure -> gradient basins */
    check += abs(ptFULL[i]);
  }

  // check for overall propability
  if ( ((check-1.) < -0.05) || ((check-1.) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1.0 %s!\n", time, (double)check, (_opt->num_err == 'R'?"rescaling":"exiting") );
    if (_opt->num_err == 'H' || check == 0.0) exit(EXIT_FAILURE);
  }

  // print:
  printf("%e ", time);
  for (int i=0; i<lmins; i++) {
    if (_opt->num_err == 'R') printf("%e ", (double)(T)(abs(ptFULL[i])/check));
    else printf("%e ", (double)(T)abs(ptFULL[i]));
  }
  printf("\n");

  return check;
}

template<typename T>
T Calc<T>::PrintProbNR(T *line, int dim, double time)
{
  T check = 0.0;

  // summ first
  for (int i=0; i<dim; i++) {
    if(line[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, line[i]);
      if (_opt->num_err == 'H') exit(EXIT_FAILURE);
      else if (_opt->num_err == 'R') line[i] = 0.0;
    }
    /* map individual structure -> gradient basins */
    check += abs(line[i]);
  }

  // check for overall propability
  if ( ((check-1.) < -0.05) || ((check-1.) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1.0 %s!\n", time, check, (_opt->num_err == 'R'?"rescaling":"exiting") );
    if (_opt->num_err == 'H' || check == 0.0) exit(EXIT_FAILURE);
  }

  // print
  printf("%e ", time);
  for (int i=0; i<dim; i++) {
    if (_opt->num_err == 'R') printf("%e ", abs(line[i])/check);
    else printf("%e ", abs(line[i]));
  }

  printf("\n");

  return check;
}

template<typename T>
T Calc<T>::PrintProb(T *line, int dim, double time)
{
  T check = (T)0.0;

  // sum it up
  for (int i=0; i<dim; i++) {
    if (line[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, (double)line[i]);
      if (_opt->num_err == 'H') exit(EXIT_FAILURE);
      else if (_opt->num_err == 'R') line[i] = 0.0;
    }
    check += abs(line[i]);
  }

  // check for overall probability
  if ( ((check-1.) < -0.05) || ((check-1.) > 0.05) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1.0 %s!\n", time, (double)check, (_opt->num_err == 'R' && (double)check != 0.0?"rescaling":"exiting") );
    if (_opt->num_err == 'H' || check == 0.0) exit(EXIT_FAILURE);
  }

  // print finally:
  printf("%e ", time);
  if (reorganize.size()==0) {
    for (int i=0; i<dim; i++) {
      /* map individual structure -> gradient basins */
      if (_opt->num_err == 'R') printf("%e ", (double)(T)(abs(line[i])/check));
      else printf("%e ", (double)(T)abs(line[i]));
    }
  } else {
    int j=0;
    for (int i=0; i<dim; i++) {
      if (j>(int)reorganize.size() || reorganize[j]!=i) {
        printf("%e ", 0.0);
      } else {

        if (_opt->num_err == 'R') printf("%e ", (double)(T)(abs(line[j])/check));
        else printf("%e ", (double)(T)abs(line[j]));
        j++;
      }
    }
  }
  printf("\n");

  return check;
}

template<typename T>
int * Calc<T>::MxErgoEigen(T *U, int dim)
{
  int count = 1;
  std::vector<bool> reached(dim, false);
  reached[0] = true;

  std::queue<int> que_ergo;
  que_ergo.push(0);

  // fill the sets of non-empty states
  while(!que_ergo.empty() && count<dim) {
    int to_do = que_ergo.front();
    que_ergo.pop();

    // collect contingency
    for (int i=0; i<dim; i++) {
      if (i==to_do) continue;
      if ((double)U[i*dim + to_do]!=0.0 && (double)U[to_do*dim + i]!=0.0 && !reached[i]) {
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


/*==*/
template<typename T>
T*
Calc<T>::MxBar2Matrix ( BarData *Data, T *R)
{
  T *U=NULL;
  if(_opt->want_degenerate) MxDoDegeneracyStuff();
  switch (_opt->method) {
  case 'A':
    U = MxMethodeA(Data);
    break;
  case 'F':
    U = MxMethodeFULL(R);
    break;
  case 'I':
    U = MxMethodeINPUT(Data, R);
    break;
  default:
    fprintf (stderr,
             "ERROR in MxBar2Matrix(): No handler 4 method %c\n", _opt->method);
    exit(EXIT_FAILURE);
  }
  if (_opt->dumpU) {
    MxASCIIWrite(U, "U.txt");
    MxBinWrite(U, "U", 'm');
  }
  return (U);
}
/*==*/
template<typename T>
void
Calc<T>::MxGetSpace (T **p8)
{
  *p8   = new T[dim];
  std::fill_n(*p8,dim,0.);
  evals = new T[dim];
  std::fill_n(evals,dim,0.);
  evecs = new T[dim*dim];
  std::fill_n(evecs,dim,0.);
  assert(evals!=NULL);assert(evecs!=NULL);assert(p8!=NULL);

if(!_opt->absrb) {
    _sqrPI = new T[dim*dim];
    std::fill_n(_sqrPI,dim*dim,0.);
    sqrPI_ = new T[dim*dim];
    std::fill_n(sqrPI_,dim*dim,0.);
  }
}

/*==*/
template<typename T>
void
Calc<T>::MxStartVec (T **p0)
{
  int i;
  T *pzero = new T[dim];
  std::fill_n(pzero,dim,0.0);
  if (_opt->pini) {
    for (i = 1; i < (int) *_opt->pini; i+=2)
      pzero[(int)_opt->pini[i]-1] = (T)_opt->pini[i+1];
    /* -1 because our lmins start with 1, not with 0 (as Data does ) */
  } else {
    // all into first state...
    pzero[0]=1.0;
  }

  if (_opt->want_verbose) MxPrint (pzero, "p0", 'v');
  *p0=pzero;
}

/*==*/
/* calculate equilibrium distribution */
template<typename T>
void
Calc<T>::MxEqDistr ( BarData *Data, T **p8 )
{
  int i;
  T Z = 0.;

  if(_opt->absrb) {
    for(i = 0; i < dim; i++)
      (*p8)[i] = 0.;
    (*p8)[dim-1] = 1.0; /* last entry is the 'new' absorbing state */
  }
  else {
    for(i = 0; i < dim; i++) Z += exp(-((Data[i].energy-Data[0].energy)/_kT));
    for(i = 0; i < dim; i++) (*p8)[i] = exp(-((Data[i].energy-Data[0].energy)/_kT));
  }

  /* now normalize the new p8 */
  for (i=0; i<dim; i++)
    *(*p8+i)  *= 1.0/Z; //sumsq;

  if(_opt->want_verbose) MxPrint (*p8, "p8", 'v');
  return;
}

/*==*/
template<typename T>
void
Calc<T>::MxEqDistrFULL (SubInfo *E, T *p8 )
{
  int i;
  T Z = 0.;

  if(_opt->absrb) {
    for(i = 0; i < dim; i++)
      p8[i] = 0.;
    p8[_opt->absrb-1] = 1.;
  }
  else {
    for(i = 0; i < dim; i++) Z += exp(-E[i].energy/_kT);
    for(i = 0; i < dim; i++) p8[i] = (T)(exp(-E[i].energy/_kT)/Z);
  }


  if(_opt->want_verbose) MxPrint (p8, "p8", 'v');
  return;
}

template<typename T>
void
Calc<T>::MxEqDistrFromLinSys( T *U, T **p8 )
{
  int i,j,n, nrhs, nfo, *ipiv=NULL;
  T *A=NULL, *B=NULL;

  if(_opt->absrb) {
    for(i = 0; i < dim; i++)
      *(*p8+i) = 0.;
    *(*p8+(dim-1)) = 1.0; /* last entry is the 'new' absorbing state */
  }
  else {
    n=dim-1;
    A    = new T[n*n];
    B    = new T[n];
    ipiv = (int *)malloc(n*sizeof(int));
    nrhs=1;

    if (_opt->want_verbose) MxPrint(U, "U", 'm' );

    for(i=1; i<=n; i++) /* all except first row */
      for(j=1; j<=n; j++)
        A[n*(i-1)+(j-1)]=U[dim*i+j] - (i==j && _opt->useplusI ? 1.0 : 0.0);  //U-I=Q
    for(n=0,i=1; i<dim; i++,n++)
      B[n]=-U[dim*i];
    dim=n;
    mxccm.trnm(A,n);

    if (_opt->want_verbose) MxPrint(A, "A in MxEqDistrFromLinSys", 'm' );
    if (_opt->want_verbose) MxPrint(B, "B in MxEqDistrFromLinSys", 'v' );

    //DGESV computes the solution to a real system of linear equations A * X = B
    solver->Mx_Dgesv(&n, &nrhs, A, &n, ipiv, B, &n, &nfo);

    if (nfo != 0) {
      fprintf(stderr, "dgesv exited with value %d (Cannot compute the equilibrium distribution) - check if the rates have good dimension(you have probably reached the numerical precision of treekin)\n", nfo);
      /// TODO switch to compute from detailed balance
      exit(EXIT_FAILURE);
    }
    if (_opt->want_verbose) MxPrint(B, "X in MxEqDistrFromLinSys", 'v' );
    dim=n+1;
    *p8[0]=1.;
    for(i=1; i<dim; i++) *(*p8+i)=B[i-1];

    if (_opt->want_verbose) MxPrint(*p8, "p8 in MxEqDistrFromLinSys before norm", 'v' );

    // now check if all are > 0.0
    for (i=dim-1; i!=0; i--) {
      if ((*p8)[i]<0.0) {
        if (i==0) exit(EXIT_FAILURE);
        fprintf(stderr, "Warning: p8[%5d] is negative (%e), setting it to a nearby value (%e)\n", i+1, (*p8)[i], (i==dim-1?(*p8)[i-1]:(*p8)[i+1]));
        (*p8)[i] = abs((i==dim-1?(*p8)[i-1]:(*p8)[i+1]));
      }
      if ((*p8)[i]==0.0) {
        fprintf(stderr, "Warning: p8[%5d] is zero (%e), setting it to a minimum value (%e)\n", i+1, (*p8)[i], 1e-20);
        (*p8)[i] = 1e-20;
      }
    }

    /* now make the vector stochastic (sum = 1.0) */
    T sum=0.0;
    for (i=0; i<dim; i++) {
      sum += (*p8)[i];
    }
    for (i=0; i<dim; i++) {
      (*p8)[i] /= sum;
    }

    delete[] A;
    delete[] B;
    free(ipiv);
  }
  if(_opt->want_verbose)
    MxPrint(*p8, "p8", 'v' );
}

/*==*/
template<typename T>
void
Calc<T>::MxDiagonalize ( T *U, T **_S, T *P8)
{
  int i,j;
  T *tmpMx=NULL;

  if(_opt->dumpMathematica == 1)  MxKotzOutMathematica(U);
  if(!_opt->absrb) {
    if (_opt->want_verbose) MxPrint (U, "input U", 'm');

    tmpMx = new T[dim*dim];
    //if (_opt->want_verbose) MxPrint (P8, "P8", 'v');
    MxDiagHelper(P8);
    mxccm.mmul(tmpMx, _sqrPI, U, dim);
    if (_opt->want_verbose) MxPrint (_sqrPI, "_sqrPI", 'm');
    if (_opt->want_verbose) MxPrint (tmpMx, "tmpMx = _sqrPI*U", 'm');
    if (_opt->want_verbose) MxPrint (sqrPI_, "sqrPI_", 'm');
    //memset(U,0,dim*dim*sizeof(T));
    fill_n(U,dim*dim,0);
    mxccm.mmul(U, tmpMx, sqrPI_, dim);
    if (_opt->want_verbose) MxPrint (U, "U = _sqrPI*U*sqrPI_ (this should be symmetric, but still uncorrected)", 'm');
    delete[] tmpMx;

    /* correct for numerical errors -- ensure U is symmetric */
    T err = 0.0;  // accumulated error
    for (i = 0; i < dim; i++) {
      for (j = i; j < dim; j++) {
        T err_inc = (T)(abs((U[dim*i+j]+U[dim*j+i])/2 - U[dim*i+j]));
        /*if (std::isnan(err_inc)) {
          fprintf(stderr, "%d %d\n", i,j);
        }*/
        err += err_inc;
        if (std::isnan((double)err_inc)) {
          fprintf(stderr, "Weird rates! r(%5d,%5d)=(%e, %e)\n", i, j, (double)U[dim*i+j], (double)U[dim*j+i]);
          //exit(-1);
        }
        U[dim*i+j] = (T)((U[dim*i+j]+U[dim*j+i])/2.);
        //U[dim*i+j] = std::sqrt(U[dim*i+j]*U[dim*j+i]);
        U[dim*j+i] = U[dim*i+j];
      }
    }
    if (std::isnan((double)err)
) {
      fprintf(stderr, "Rates are not good -- check your rates (maybe the equilibrium was not computed correctly? rerun with --equil-file or -v set)!!!... (%Lf)\n", (long double)err);
      exit(-1);
    }
    if (!_opt->quiet) fprintf(stderr, "Corrected numerical error: %e (%e per number)\n", (double)err, (double)(T)(err/(T)(dim*dim)));


    if (_opt->want_verbose) MxPrint (U, "force symmetrized U", 'm');
  }

    // get eigenv*
  if(_opt->absrb) {
    MxEVLapackNonSym(U);
  } else {
    switch(_opt->mpackMethod) {
#if defined(WITH_MPACK_GMP)
      case MPACK_GMP:
        solver->MxEV_Mpack_Sym_gmp((const mpf_class *)U, dim,(mpf_class *)evals,(mpf_class *)evecs, _opt->mpackMethod_Bits);
        break;
#endif
#if defined(WITH_MPACK_QD)
      case MPACK_QD:
        solver->MxEV_Mpack_Sym_qd((const qd_real *)U, dim,(qd_real *)evals,(qd_real *)evecs);
        break;
#endif
#if defined(WITH_MPACK_MPFR)
      case MPACK_MPFR:
        solver->MxEV_Mpack_Sym_mpfr((const mpfr::mpreal *)U, dim,(mpfr::mpreal *)evals,(mpfr::mpreal *)evecs, _opt->mpackMethod_Bits);
        break;
#endif
#if defined(WITH_MPACK___FLOAT128)
      case MPACK_FLOAT128:
        solver->MxEV_Mpack_Sym_float128((const __float128 *)U, dim,(__float128 *)evals,(__float128 *)evecs);
        break;
#endif
#if defined(WITH_MPACK_LD)
      case MPACK_LD:
        solver->MxEV_Mpack_Sym_longdouble((const long double *)U, dim,(long double *)evals,(long double *)evecs);
        break;
#endif
#if defined(WITH_MPACK_DD)
      case MPACK_DD:
        solver->MxEV_Mpack_Sym_dd((const dd_real *)U, dim,(dd_real *)evals,(dd_real *)evecs);
        break;
#endif
#if defined(WITH_MPACK_DOUBLE)
        case MPACK_DOUBLE:
          solver->MxEV_Mpack_Sym_double((const double *)U, dim,(double *)evals,(double *)evecs);
        break;
#endif
      default:
        //default standard lapack
        MxEVLapackSym(U);
        break;
    }
  }

  MxSortEig(evals, evecs);

  if (_opt->want_verbose) {
    MxPrint(evecs, "Eigenvectors of U (LAPACK)", 'm');
    MxPrint(evals, "Eigenvalues of U (LAPACK)", 'v');
  }

  // check for ergodicity of eigenvectors (if not ergodic, then not all the space can be reached)
  int *reachable = MxErgoEigen(evecs, dim);
  if (reachable) {
    fprintf(stderr, "WARNING: Eigenvector matrix is non-ergodic (unreachable states from 1st state):");
    for (i=0; i<dim; i++) if (!reachable[i]) fprintf(stderr, " %4d", i+1);
    fprintf(stderr, "\n         the overall probability can start to decrease if this state(s) is(are) populated!!!");

    free(reachable);
    //exit(EXIT_FAILURE);
  }

  if (!_opt->quiet) {
    if ((T)abs(evals[0] - ((_opt->useplusI)?1.:0.)) > (T)10*_opt->FEPS) {
      fprintf(stderr, "WARNING largest eigenvalue is %g, see evals.txt and evecs.txt\n", (double)evals[0]);
      MxASCIIWriteV(evals, "evals.txt");
      MxASCIIWrite(evecs, "evecs.txt");
    }
  }
  // fix evals[0] to be 0 (or 1.0)
  if (_opt->useplusI) evals[0] = 1.0;
  else evals[0] = 0.0;

  // check if we have only one eval = 0 / eval = 1
  i = 0;
  if (_opt->useplusI) while (evals[i] >= 1.0) evals[i++] = 1.0;
  else              while (evals[i] >= 0.0) evals[i++] = 0.0;

  if (i > 1) {
    fprintf(stderr, "WARNING: There are multiple evals=%d.0 (%4d) (multiple population traps -- maybe we compute crap)\n", _opt->useplusI?1:0, i);
  }

  if(_opt->absrb)
    MxFixevecsAbsorb(evecs,evals);
  //  else
  //  MxFixevecs(evecs,evals); /* so far it's not helping ... */

  *_S=evecs;
  if(_opt->wrecover) {
    MxBinWrite(evals, "evals", 'v');
    MxBinWrite(evecs, "evecs", 'm');
    MxASCIIWriteV(evals, "evals.txt");
    MxASCIIWrite(evecs, "evecs.txt");
    if (_opt->want_verbose) {
      T *CL; //, *tmpMx; ?
      int i,j;
      CL        = new T[dim*dim];
      mxccm.mmul (CL, sqrPI_, evecs, dim);
      MxASCIIWrite(CL, "evecR.txt");
      fprintf(stderr, "Sums of EVs of rate matrix R\n");
      for (i=0; i<dim; i++) {
        T sum, sum2, sum3;
        sum=sum2=sum3=0.;
        for (j=0; j<dim; j++) {
          sum += CL[dim*j+i];
        }
        fprintf(stderr, "%d %15.8g  ", i, (double)sum);
      }
      fprintf(stderr, "\n");
    delete[] CL;//free(CL);
    }
  }
  // compensate for the +I matrix
  if (_opt->useplusI) for(i = 0; i < dim; i++) evals[i] -= 1.0;  /* compensate 4 translation of matrix U */
}

/*==*/
template<typename T>
void
Calc<T>::MxDiagHelper(T *P8)
{
  int i,j;
  for(i = 0; i < dim; i++)
  {
    j=i;
    sqrPI_[dim*i+j] = (T)std::sqrt(P8[i]);          /* pos right */
    _sqrPI[dim*i+j] = (T)(1.0/(sqrPI_[dim*i+j]));  /* neg left */
  }
}

/*==*/
template<typename T>
void
Calc<T>::MxRecover (T **_S, T *P8)
{
  MxDiagHelper(P8);
  delete[] evecs;
  delete[] evals;

  if(MxBinRead(&evecs, "evecs", 'm') != dim ) {
    fprintf(stderr, "ERROR: MxBinRead() returns wrong dimension for evecs\n");
    exit(EXIT_FAILURE);
  }
  if(MxBinRead(&evals, "evals", 'v') != dim) {
    fprintf(stderr, "ERROR: MxBinRead() returns wrong dimension for evals\n");
    exit(EXIT_FAILURE);
  }
  *_S=evecs;
  if(_opt->want_verbose) {
    MxPrint(evals, "MxRecover: Eigenvalues", 'v');
    MxPrint(evecs, "MxRecover: Eigenvectors", 'm');
  }
  return;
}


template<class T>
int Calc<T>::ConvergenceReached(T *p8, T *pt, int dim, int full) {

  int pdiff_counter = 0;
  full = (full?1:0);
  /* now check if we have converged yet */
  for(int i = 0; i < dim; i++) {
    if (abs(p8[i] - pt[i]) >= 0.000001) {
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

/**
 * @param S which comes into this function contains the (right)
 * eigenvectors calculated by LAPACK
 */
template<typename T>
void
Calc<T>::MxIterate (T *p0, T *p8, T *S)
{
  /*  solve following equation 4 various times
    ***** NON-ABSORBING CASE: ******
    p(t)    = sqrPI_ * S * exp(time * EV) * St * _sqrPI * p(0)
    CL      = sqrPI_ * S
    CR      = St * _sqrPI
    tmpVec  = CR * p(0)
    tmpVec2 = exp(time * EV) * tmpVec
    p(t)    = CL * tmpVec2
    ******* ABSORBING CASE: *******
    p(t)    = S * exp(time * EV) * S_inv * p(0)
    tmpVec  = S_inv * p(0)
    tmpVec2 = exp(time * EV) *tmpVec
    p(t)    = S * tmpVec2
  */
  int i,  count = 0;
  T time, check = 0.;
  T *CL = NULL, *CR, *exptL, *tmpVec, *tmpVec2, *pt, *St, *pdiff;
  T *ptFULL = NULL;  /* prob dist 4 of the effective lmins of the tree at time t */
  T *p8FULL = NULL;  /* equ dist 4 gradient basins, full process */

  tmpVec    = new T[dim];
  std::fill_n(tmpVec,dim,0.);
  tmpVec2   = new T[dim];
  std::fill_n(tmpVec2,dim,0.);
  pt        = new T[dim];
  std::fill_n(pt,dim,0.);
  size_t dimSquared = dim*dim;
  exptL     = new T[dimSquared];
  std::fill_n(exptL,dimSquared,0.);
  if(_opt->method=='F') {
    size_t fullSize = _lmins+1;
    ptFULL    = new T[fullSize];
    std::fill_n(ptFULL,fullSize,0.);
    p8FULL    = new T[fullSize];
    std::fill_n(p8FULL,fullSize,0.);
    pdiff     = new T[fullSize];
    std::fill_n(pdiff,fullSize,0.);
  }
  else pdiff   = new T[dim];

  if(! _opt->absrb) { /* NON-absorbing case */
    CL        = new T[dim*dim];
    CR        = new T[dim*dim];
    St        = new T[dim*dim];
    mxccm.mcopy(St, S, dim*dim);
    mxccm.trnm(St, dim); /* transpose S (same as invert(S) ('cause of symmetry) minv(St, dim); )*/
    mxccm.mmul (CL, sqrPI_, S, dim);
    mxccm.mmul (CR, St, _sqrPI, dim);
    mxccm.vmul (tmpVec, CR, p0, dim);

    // test CL*CR = I
    if (_opt->want_verbose) {
      T *testMx, *tvec;
      testMx      = new T[dim*dim];
      tvec    = new T[dim];
      mxccm.mmul(testMx,CL,CR,dim);
      MxPrint(testMx, "CL*CR should be I", 'm');
      delete[] testMx;
      MxPrint(tmpVec, "CR*p0 should start with 1", 'v');
      mxccm.vmul(tvec, _sqrPI, p0, dim);
      mxccm.vmul(tmpVec2, St, tvec, dim);
      MxPrint(tmpVec2, "CR*p0 should start with 1", 'v');
      MxPrint(CL, "CL", 'm');
      MxPrint(CR, "CR", 'm');
      delete[] tvec;
    }
    delete[] St;
    delete[] CR;
  }
  else { /* absorbing case */
    T *S_inv = new T[dim*dim];
    mxccm.mcopy(S_inv, S, dim*dim);
    mxccm.minv(S_inv,dim, _opt->FEPS);
    if(_opt->want_verbose) MxPrint(evals, "evals in MxIterate", 'v');
    mxccm.vmul (tmpVec, S_inv, p0, dim);
    delete[] S_inv;
  }
  if(_opt->method=='F') { /* calculate equilibrium distribution once */
    for (i = 0; i < dim; i++)   p8FULL[_E[i].ag] += p8[i];
    for (i = 0; i < _lmins; i++) check += abs(p8FULL[i]);
    if ( ((check-1.) < -0.1) || ((check-1.) > 0.1) ) {
      fprintf(stderr, "overall equilibrium probability is %e != 1. !\n", (double)check);
      if (_opt->num_err == 'H' || check == 0.0) exit(EXIT_FAILURE);
      else if (_opt->num_err == 'R') {
        for (i=0; i<dim; i++) p8[i] /= check;
        check = 1.0;
      }
    }
  }
  check = 0.;

  /* solve fundamental equation */
  print_settings();
  if (_opt->t0 == 0.0) {
    if (_opt->method=='F')  PrintProbFull(p0, dim, 0.0, _lmins);
    else                  PrintProb(p0, dim, 0.0);
    _opt->t0 = _globalParameters->TZERO;
  }

  T underflow[dim], tt;
  for (i=0; i<dim; i++) underflow[i] = 0.0;

  // iterate
  for (tt = _opt->t0; tt < _opt->t8*_opt->tinc; tt *= _opt->tinc) {
    time = (tt<_opt->t8)?tt:(T)_opt->t8;
    for (i = 0; i < dim; i++) {
      errno = 0;
      exptL[dim*i+i] = (T)exp(time/_opt->times*evals[i]);
      if ((errno == ERANGE || std::isnan((double)exptL[dim*i+i])
) && underflow[i]==0.0) {
        //if (_opt->warnings) fprintf(stderr, "WARNING: underflow occured on %dth state at time %g -- exp(%g * %g) = %g\n", i+1, time/_opt->times, time/_opt->times, evals[i], exptL[dim*i+i]);
                //         the overall probability can start to decrease if this state is still populated!!!\n         p_%d(%g) = %g, so it seems this %s\n", i+1, time/_opt->times, time/_opt->times, evals[i], exptL[dim*i+i], i+1, time/_opt->times, pt[i], pt[i]>0.1?"is DEFINITELLY BAD":(pt[i]>0.001?"is POTENTIALLY BAD":"SHOULD BE OK"));
        underflow[i] = (T)(time/_opt->times);
        //exit(EXIT_FAILURE);
      }
    }
    mxccm.vmul(tmpVec2, exptL, tmpVec, dim);
    if(!_opt->absrb)  mxccm.vmul(pt, CL, tmpVec2, dim);
    else            mxccm.vmul(pt, S, tmpVec2, dim);

    /*if (count%100 == 0) {
      fprintf(stderr, "Time: %g\n", time/_opt->times);
      MxPrint(tmpVec2, "tmpVec2", 'v');
      MxPrint(pt, "pt", 'v');
    }*/

    count++;  /* # of iterations */

    if(_opt->method=='F') {
      //memset(ptFULL, 0, (_lmins+1)*sizeof(T)); //replaced by c++ style.
      fill_n(ptFULL, _lmins+1, 0);
      for (i = 0; i < dim; i++) {
        ptFULL[_E[i].ag] += pt[i];
      }
    }

    // print probabilities with respect to corrected ergodicity
    if (_opt->method=='F') check = PrintProbFull(pt, dim, time, _lmins);
    else                 check = PrintProb(pt, dim, time);

    //PrintProbNR(p8FULL, lmins, -1);
    int reached;
    if (_opt->method=='F') reached = ConvergenceReached(p8FULL, ptFULL, _lmins, 1);
    else                 reached = ConvergenceReached(p8, pt, dim, 0);
    fflush(stdout);
    if (reached) break;
  }

  // print underflow:
  if (_opt->warnings) {
    for (i=0; i<dim; i++) {
      if (underflow[i] > 0.0) fprintf(stderr, "underflow %5d at time %12g", i+1, (double)underflow[i]);
    }
  }


  if (time < _opt->t8) {
    if (_opt->method=='F') PrintProbFull(pt, dim, _opt->t8, _lmins);
    else                 PrintProb(pt, dim, _opt->t8);
  }
  printf("# of iterations: %d\n", count);

  /*** end solve fundamental equation ***/

  if(_opt->method=='F') {
    delete[] ptFULL;
    delete[] p8FULL;
    delete[] _E;
  }
  //free(evals);
  delete[] exptL;
  delete[] CL;
  delete[] tmpVec;
  delete[] tmpVec2;
  delete[] pdiff;
  delete[] pt;
}

/*==*/
template<typename T>
T*
Calc<T>::MxMethodeA (BarData *Data)
{
  /***************************************/
  /*           |       E_s           |   */
  /*           \   ____________      /   */
  /*            \  |   /\     |     /    */
  /*             \ |  /  \    |    /     */
  /*              \| /    \   |   /      */
  /*               \/      \  |  /       */
  /*               i        \ | /        */
  /*                         \/          */
  /*                         j           */
  /*                       -beta(E_s-E_j) */
  /*  rate:  j->i:  prop  e               */
  /*                       -beta(E_s-E_i) */
  /*         i->j   prop  e               */
  /****************************************/

  int i,j,real_abs = 0;
  T m_saddle, Zabs, abs_rate;

  T* U = new T[dim*dim];

  for( i = 0; i < dim; i++) {
    for( j = i+1; j < dim; j++) {
      m_saddle = max_saddle(i, j, Data);

      T tmp_ji = (T)(1.0*exp(-(m_saddle-Data[j].energy)/_kT));
      T tmp_ij = (T)(1.0*exp(-(m_saddle-Data[i].energy)/_kT));
      /* rate j -> i */
      U[dim*i+j] = tmp_ji;
      /* rate i -> j */
      U[dim*j+i] = tmp_ij;
    }
  }

  if(_opt->absrb) { /*==== absorbing  states ====*/
    dim++;
    fprintf(stderr, "dim increased to %i\n", dim);
    T * tmpU = new T[dim*dim];
    real_abs = _opt->absrb; /* the original absorbing lmin */
    real_abs--;
    _opt->absrb = dim; /* the 'new' abs state = last row/column of rate matrix */
    fprintf(stderr, "new absorbing lmin is: %i\n", _opt->absrb);
    Zabs = (T)(exp((-Data[real_abs].energy)/_kT));
    abs_rate = (T)(exp((-Data[real_abs].energy)/_kT)/Zabs);

    for(i = 0; i < (dim-1); i++) { /* all except the last row */
      for(j = 0; j < (dim-1); j++)
        tmpU[dim*i+j] = U[(dim-1)*i+j];
      tmpU[(dim-1)*j+(dim-1)] = 0.;
    }
    for(j = 0; j < dim; j++) /* last row */
      tmpU[dim*(dim-1)+j] = 0.;
    tmpU[dim*(dim-1)+real_abs] = abs_rate;
    delete[] U; //free(U);
    U = tmpU;
    if(_opt->want_verbose) MxPrint(U, "aufgeblasene Matrix", 'm');

  }  /*== end absorbing states ==*/

  /* set diagonal elements to 0 */
  for (i = 0; i < dim; i++) U[dim*i+i] = 0.;
  for (j = 0; j < dim; j++) {
    T tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++) tmp += U[dim*i+j];
    U[dim*j+j] = (T)(-tmp); /* make U a stochastic matrix */
    if (_opt->useplusI) U[dim*j+j] += 1.0;
  }
  if (_opt->want_verbose) MxPrint (U,"U with Methode A", 'm');
  return (U);
}

/*==*/
template<typename T>
T*
Calc<T>::MxMethodeFULL (T *R)
{
  int i, j;
  delete[] D; //free(D);

  if(_opt->absrb) { /*==== absorbing  states ====*/
    for(i = 0; i < dim; i++)
      R[dim*i+(_opt->absrb-1)] = 0. ;
  }              /*== end absorbing states ==*/

  /* set diagonal elements  to 0 */
  for (i = 0; i < dim; i++) R[dim*i+i] = 0.;
  for (j = 0; j < dim; j++) {
    T tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++) tmp += R[dim*i+j];
    R[dim*j+j] = (T)(-tmp); /* make U a stochastic matrix */
    if (_opt->useplusI) R[dim*j+j] += 1.0;
  }

  if (_opt->want_verbose) MxPrint (R, "R with Methode F", 'm');
  return R;
}

/*==*/
template<typename T>
T*
Calc<T>::MxMethodeINPUT (BarData *Data, T *Input)
{
  int i, j;
  T *U=NULL, Zabs, abs_rate;

  if (_opt->want_verbose) MxPrint(Input, "Input Matrix", 'm');

  if (_opt->absrb) {  /*==== absorbing  states ====*/
    dim++;
    fprintf(stderr, "dim increased to %i\n", dim);
    U = new T[dim*dim]; //(T *) MxNew(dim*dim*sizeof(T));
    _opt->real_abs = _opt->absrb; /* the original absorbing lmin */
    _opt->real_abs--;
    _opt->absrb = dim; /* the 'new' abs state = last row/column of rate matrix */
    fprintf(stderr, "new absorbing lmin is: %i\n", _opt->absrb);
    Zabs = exp((-Data[_opt->real_abs].FGr)/_kT);
    abs_rate = (T)(exp((-Data[_opt->real_abs].energy)/_kT)/Zabs);

    for(i = 0; i < (dim-1); i++) { /* all except the last row */
      for(j = 0; j < (dim-1); j++)
        U[dim*i+j] = Input[(dim-1)*i+j];
      U[(dim-1)*j+(dim-1)] = 0.;
    }
    for(j = 0; j < dim; j++) /* last row */
      U[dim*(dim-1)+j] = 0.;
    U[dim*(dim-1)+_opt->real_abs] = abs_rate;
    /*   if(_opt->want_verbose) MxPrint(U, "aufgeblasene Matrix", 'm'); */
  }      /*== end absorbing states ==*/
  else { /*== non-absorbing states ==*/
    U = new T[dim*dim]; //(T *) MxNew(dim*dim*sizeof(T));
    //memcpy(U, Input, dim*dim*sizeof(T));
    std::copy (Input, Input+dim*dim, U );
  }      /*== end non-absorbing states ==*/

  /* diagonal elements */
  for (i = 0; i < dim; i++) U[dim*i+i] = 0.0;
  //fprintf(stderr, "dim is %i\n", dim);
  for (i = 0; i < dim; i++) {
    T tmp = 0.00;
    // calculate column sum
    for(j = 0; j < dim; j++)  tmp += U[dim*j+i];
    U[dim*i+i] = (T)(-tmp);   // make U a stochastic matrix U = Q+I ??
    if (_opt->useplusI) U[dim*i+i] += 1.0;
  }

/*  // normalize each column to sum=1
  for (i = 0; i < dim; i++) {
    double tmp = 0.00;
    for(j = 0; j < dim; j++) tmp += U[dim*j+i];
    for(j = 0; j < dim; j++) U[dim*j+i] /= tmp;
  }*/


  if(_opt->want_verbose) MxPrint (U,"U with Methode I" , 'm');
  delete[] Input;
  return U;
}

/*==*/
template<typename T>
T
Calc<T>::max_saddle(int i, int j, BarData *Data)
{
  int tmp;

  if(Data[i].number > Data[j].father) { /* exchange i & j */
    tmp = i;
    i = j;
    j = tmp;
  }
  if(Data[i].number == Data[j].father) return (T)(Data[j].energy + Data[j].ediff);
  else {
    if((Data[j].energy + Data[j].ediff) > max_saddle((Data[j].father - 1), (Data[i].number - 1), Data)) return (Data[j].energy + Data[j].ediff);
    else return max_saddle((Data[j].father - 1), (Data[i].number - 1), Data);
  }
}


template<typename T>
void Calc<T>::MxFPrint(T *mx, const char *name, char c, FILE *out, int pure)
{
  int k, l;
  switch (c) {
  case 'm':    /* square matrix */
    if (!pure) fprintf(out,"%s:\n", name);
    for (k = 0; k < dim; k++) {
      for (l=0; l< dim; l++) {
        fprintf(out,"%15.7g ", (double)mx[dim*k+l]);
      }
      fprintf(out,"\n");
    }
    if (!pure) fprintf(out,"---\n");
    break;
  case 'v':
    if (!pure) fprintf(out,"%s:\n", name);
    for (k = 0; k < dim; k++) fprintf(out,"%15.7g ", (double)mx[k]);
    if (!pure) fprintf(out,"\n---\n");
    else fprintf(out, "\n");
    break;
  default:
    fprintf(out,"ERROR MxPrint(): no handler 4 type %c\n", c);
  }
}

/**
 * print matrix stored in ccmath-format
 */
template<typename T>
void
Calc<T>::MxPrint(T *mx, const char *name, char c)
{
  MxFPrint(mx, name, c, stderr, 0);
}


/*==*/
template<typename T>
void
Calc<T>::norm2(T *mx)
{
  /* normalize columns of matrix mx */
  /* (to euclidean norm 1)*/
  int i,j;
  T sumsq;
  T mxValue;
  for (j=0; j<dim; j++) {
    sumsq=0.0;
    for (i=0; i<dim; i++){
      mxValue=mx[dim*i+j];
      sumsq += mxValue*mxValue;
    }
    if(sumsq > 0.0)
      sumsq=(T)(1./std::sqrt(sumsq));
    for (i=0; i<dim; i++)
      mx[dim*i+j] *= sumsq;
  }
  return;
}

//#define abs(x) ((x)>0.0?(x):(-x))

/*==*/
template<typename T>
void
Calc<T>::MxFixevecsAbsorb(T *evecs, T *evals)
/* evecs: eigenvectors columns, dimension N   */
/* since the sum over each non-absorbing      */
/* column must vanish, replace the            */
/* problematic evecs with linear combinations */
/* satisfying that criterion                  */
/* rationale: with 1^T a row vector of 1s      */
/* Q=S*exp(L)*inv(S)                          */
/* 1^T*Q = 1^T  (probabilities sum up to 1)   */
/* 1^T*S*exp(L) = 1^T*S                       */
{
  int i,j,maxind,abscount;
  T maxent,csum;

  /* assume sorted evals        */
  /* and an absorbing state     */
  /* im the 1st evec=1st column */

  /* take care of possibly more than one abs. state*/
  /* all abs. states have been *assigned* eigenvalue one */
  /* for each of them, set the largest component of the */
  /* eigenvectors to one, all others to zero */
  abscount=0;
  for (i=0; i<dim; i++) {
    if (evals[i]==1.0)  {
      abscount++;
      maxent=0.;
      maxind=0;
      for (j=0; j<dim; j++) {
        if ((T)(abs(evecs[dim*j+i])) > maxent) {
          maxent=(T)(abs(evecs[dim*j+i]));
          maxind=j;
        }
      }
      evecs[dim*maxind+i]=1.0;
      for (j=0; j<dim; j++) {
        if (j==maxind) continue;
        evecs[dim*j+i]=0.0;
      }
    }
  }

  /* repair messed non abs. eigenvectors */
  /* using all abs. states equally */
  /* using all abs. states equally */
  for (j=abscount; j< dim; j++) {
    T mu=0.;
    for(i=0; i<dim; i++)
      mu += evecs[dim*i+j];
    for (i=0; i<abscount; i++)
      evecs[dim*i+j] -= mu/(T)abscount;
  }

  if (_opt->want_verbose) {
    MxPrint (evals, "evals_complete", 'v');
    MxPrint (evecs, "evecs_complete", 'm');
    fflush(stdout);
    fflush(stderr);
  }

  norm2(evecs);

  if (_opt->want_verbose) {
    MxPrint (evals, "evals_complete", 'v');
    MxPrint (evecs, "evecs_complete", 'm');
    fflush(stdout);
    fflush(stderr);
  }

  if (_opt->want_verbose) {
    fprintf(stderr,"colsums: ");
    for (i=0; i<dim; i++) {
      csum=0.0;
      for (j=0; j<dim; j++)
        csum += evecs[dim*j+i];
      fprintf(stderr,"%g ", (double)csum);
    }
    fprintf(stderr,"\n");
  }

  if (_opt->absrb && (abscount > 1)) {
    int i,j;
    if (!_opt->quiet) fprintf(stderr, "\nWARNING: found %d additional absorbing state(s): ", abscount-1);
    for(i=1; i<dim; i++) {
      if(evals[i] == 1.) {
        for(j=0; j<dim; j++) {
          if(evecs[dim*j+i]==1.)
            fprintf(stderr, " %5d", j+1);
        }
      }
    }
    fprintf(stderr, "\n");
  }
  fflush(stderr);
}


template<typename T>
void
Calc<T>::MxFixevecs(T *evecs, T *evals)
/* Component sum of eigenvectors of M has to be zero, except for largest eigenvalue
   Strategy to enforce this property:
   transform evecs (eigenvectors of the symmertic Matrix U)
   back into space of the Rate Matrix M. With S the Matrix of eigenvectors,
   we compute

       V = sqrPI_ * S

       V_0 should be the equilibrium vector, i.e.
       V[0,j]>0 and   Sum_j V[0,j] =1

       foreach V_i, i>0, compute d = Sum_j V[i,j]
       and V_i = V_i - d*V_0

   Now transform the vectors back into the space of the symmetric matrix U

       S = _sqrPI * V
*/
{
  int i,j;
  T *V, sum, sum0;

  V = new T[dim*dim]; //(T *) MxNew (dim*dim*sizeof(T));
  mxccm.mmul (V, sqrPI_, evecs, dim);

  if (_opt->want_verbose)
    MxPrint(V, "Eigenvectors of rate matrix M, before Fixevecs", 'm');

  fprintf(stderr, "Sums of EVs of rate matrix R before fixing\n");
  for (i=0; i<dim; i++) {
    T sum;
    sum=0.;
    for (j=0; j<dim; j++) {
      sum += V[dim*j+i];
    }
    fprintf(stderr, "%15.8g   ", i, sum);
  }
  fprintf(stderr, "\n");

  sum0=0.; i=0;
  for (j=0; j<dim; j++) {
    //    if (V[dim*j+i]<0) {
    //  if (!_opt->quiet)
    //  fprintf(stderr, "correcting negative equilib probability for state p[%d]=%g\n",j,V[dim*j+i]);
    //      V[dim*j+i]=0;
    //}
    sum0 += V[dim*j+i];
  }
  // for (j=0; j<dim; j++)
  //   V[dim*j+i] /= sum;

  for (i=1; i<dim; i++) {
    sum=0.;
    for (j=0; j<dim; j++)
      sum += V[dim*j+i];
    for (j=0; j<dim; j++)
      V[dim*j+i] -= sum /sum0 * V[dim*j+0];
  }

  fprintf(stderr, "Sums of EVs of rate matrix R after fixing\n");
  for (i=0; i<dim; i++) {
    T sum;
    sum=0.;
    for (j=0; j<dim; j++) {
      sum += V[dim*j+i];
    }
    fprintf(stderr, "%15.8g   ", i, (double)sum);
  }
  fprintf(stderr, "\n");

  mxccm.mmul(evecs, _sqrPI, V , dim);

  norm2(evecs);

  if (_opt->want_verbose)
    MxPrint(evecs, "Eigenvectors of symmetric matrix U, after Fixevecs", 'm');

  delete[] V; //free(V);
}


/**
 * sort evecs,eval
 */
template<typename T>
void
Calc<T>::MxSortEig(T *evals, T *evecs)
{
  int i,j,k;
  T p;

  for (i=0; i<dim; i++) {
    p=evals[k=i];
    for (j=i+1; j<dim; j++)
      if (evals[j] >= p) p=evals[k=j];
    if (k != i) {
      evals[k]=evals[i];
      evals[i]=p;
      for (j=0; j<dim; j++) {
        p=evecs[dim*j+i];
        evecs[dim*j+i]=evecs[dim*j+k];
        evecs[dim*j+k]=p;
      }
    }
  }
}

/*==*/
template<typename T>
void
Calc<T>::MxBinWrite(T *Mx, const char what[], char c)
{
  int i,j,len=-1;
  FILE *BINOUT=NULL;
  char *wosis=NULL, *binfile=NULL;
  const char *suffix = "bin";
  //size_t info;

  wosis=(char *)what;
  if (_opt->basename == NULL)
    len=strlen(suffix)+strlen(wosis)+4;
  else
    len=strlen(_opt->basename)+strlen(wosis)+strlen(suffix)+4;
  binfile = (char *)calloc(len, sizeof(char));
  assert(binfile != NULL);
if(_opt->basename != NULL) {
    strcpy(binfile, _opt->basename);
    strcat(binfile, ".");
    strcat(binfile, wosis);
    strcat(binfile, ".");
    strcat(binfile, suffix);
  }
  else { /* this should not happen */
    strcpy(binfile,wosis);
    strcat(binfile, ".");
    strcat(binfile,suffix);
  }
  if (!_opt->quiet) fprintf(stderr, "MxBinWrite: writing %s to %s\n", wosis, binfile);
  BINOUT = fopen(binfile, "w");
  if (!BINOUT) {
    fprintf(stderr, "ERROR: could not open file pointer for binary outfile %s\n", binfile);
    exit(EXIT_FAILURE);
  }
  /* first write dim to file */
  fwrite(&dim,sizeof(int),1,BINOUT);
  switch(c) {
  case 'm':  /* write matrix entries */
    for(i=0; i<dim; i++)
      for(j=0; j<dim; j++){
        double val = (double)Mx[dim*i+j];
        fwrite(&val,sizeof(double),1,BINOUT);
        //if (!_opt->quiet) fprintf(stderr, "\n");
      }
    break;
  case 'v': /* write vector entries */
    for(i=0; i<dim; i++){
      double val = (double)Mx[i];
      fwrite(&val,sizeof(double),1,BINOUT);
      //fprintf(stderr, "\n");
      }
    break;
  default:
    fprintf(stderr, "ERROR MxBinWrite(): no handler for type %c\n",c);
    exit(EXIT_FAILURE);
  }
  fclose(BINOUT);
  //  if(binfile) free(binfile);
}

/*==*/
template<typename T>
int
Calc<T>::MxBinRead(T **Mx, const char what[], char c)
{
  int dimension=0,len=-1;
  FILE *BININ=NULL;
  char *wosis=NULL, *binfile=NULL;
  const char *suffix="bin";
  double *data=NULL;
        int ref;
  //size_t info;

  wosis=(char*)what;
  if (_opt->basename == NULL)
    len=strlen(suffix)+strlen(wosis)+4;
  else
    len=strlen(_opt->basename)+strlen(wosis)+strlen(suffix)+4;
  binfile = (char *)calloc(len, sizeof(char));
  assert(binfile != NULL);
if(_opt->basename != NULL) {
  strcpy(binfile, _opt->basename);
  strcat(binfile, ".");
  strcat(binfile, wosis);
  strcat(binfile, ".");
  strcat(binfile, suffix);
}
else { /* this should not happen */
  strcpy(binfile,wosis);
  strcat(binfile, ".");
  strcat(binfile,suffix);
}
fprintf(stderr, "MxBinRead: reading %s from %s\n", wosis, binfile);
BININ = fopen(binfile, "r+");
if (!BININ) {
  fprintf(stderr, "ERROR: could not open file pointer for\
    binary infile %s\n", binfile);
  exit(EXIT_FAILURE);
}
/* read dimension from file */
ref = fread(&dimension,sizeof(int),1,BININ);
size_t numberOfValues;
switch(c) { /* read data */
  case 'm':
  numberOfValues = dimension*dimension;
  data = (double *)calloc(numberOfValues, sizeof(double));
  ref = fread((void*)data, sizeof(double), numberOfValues,BININ);
  break;
  case 'v':
  numberOfValues = dimension;
  data = (double *)calloc(numberOfValues, sizeof(double));
  ref = fread(data, sizeof(double), numberOfValues,BININ);
  break;
  default:
  fprintf(stderr, "ERROR MxBinRead(): no handler for type %c\n",c);
  exit(EXIT_FAILURE);
}

T* convertedValues= new T[numberOfValues]; //(T *)calloc(numberOfValues, sizeof(T));
for(size_t i = 0; i < numberOfValues; i++) {
  convertedValues[i] = (T)data[i];
}

*Mx = convertedValues;
fclose(BININ);
free(binfile);
return dimension;
}

/*==*/

template<typename T>
inline
void
Calc<T>::MxASCIIWrite(T *Mx, const char *asciifile)
{
int i, j;
FILE *ASCIIOUT;

ASCIIOUT = fopen(asciifile, "w");
if (!ASCIIOUT) {
  fprintf(stderr, "could not open file pointer 4 ASCII outfile\n");
  exit(EXIT_FAILURE);
}
for(i=0; i<dim; i++) {
  for(j=0; j<dim; j++) {
    fprintf(ASCIIOUT,"%15.10g ", (double)Mx[dim*j+i]);
  }
  fprintf(ASCIIOUT,"\n");
}
if (!_opt->quiet) fprintf(stderr, "matrix written to ASCII file\n");
fclose(ASCIIOUT);
}

template<typename T>
inline
void
Calc<T>::MxASCIIWriteV(T *Mx, const char *asciifile)
{
int i;
FILE *ASCIIOUT;

ASCIIOUT = fopen(asciifile, "w");
if (!ASCIIOUT) {
  fprintf(stderr, "could not open file pointer 4 ASCII outfile\n");
  exit(EXIT_FAILURE);
}
for(i=0; i<dim; i++) {
  fprintf(ASCIIOUT,"%15.10g ", (double)Mx[i]);
}
if (!_opt->quiet) fprintf(stderr, "vector written to ASCII file\n");
fclose(ASCIIOUT);
}

/*==*/

template<typename T>
int
Calc<T>::MxExponent(T *p0, T *p8, T *U)
{
ExpMatrix em;

int i,j, pdiff_counter = 0, count = 0;
T x, tt, time, *Uexp, *Umerk, *pt, *pdiff, check = 0.;

Umerk = new T[dim*dim]; //(T *) MxNew (dim*dim*sizeof(T));
Uexp = new T[dim*dim];//(T *) MxNew (dim*dim*sizeof(T));
pt = new T[dim];//(T *) MxNew (dim*sizeof(T));
pdiff = new T[dim];//(T *) MxNew (dim*sizeof(T));

//memcpy(Umerk, U, dim*dim*sizeof(T));
std::copy(U,U+dim*dim, Umerk);//U to Umerk

/* solve fundamental equation */
if (_opt->t0 == 0.0) {
  if (_opt->method=='F') PrintProbFull(p0, dim, 0.0, _lmins);
  else PrintProb(p0, dim, 0.0);
  _opt->t0 = _globalParameters->TZERO;
}

for (i=0; i<dim; i++) U[(dim+1)*i] -= 1;
print_settings();
for (tt = _opt->t0; tt < _opt->t8*_opt->tinc; tt *= _opt->tinc) {
  time = (tt<_opt->t8)? tt:(T)_opt->t8;
  //memcpy(U, Umerk, dim*dim*sizeof(T));
  std::copy(Umerk,Umerk+dim*dim, U);// Umerk to U;
  for (i=0; i<dim*dim; i++) U[i]*=time;
  em.padexp(U,Uexp,dim,30, _opt->FEPS);
  x = 0.;
  for(j=0; j<dim*dim; j++) x+=Uexp[j];
  for(j=0; j<dim*dim; j++) Uexp[j]*=(T)dim/x;
  mxccm.vmul(pt, Uexp, p0, dim);
  /* check convergence */
  for(i=0; i<dim; i++) {
    pdiff[i] = (T)(p8[i] - pt[i]);
    if (abs(pdiff[i]) >= 0.0001)
    pdiff_counter++;
  }
  if (pdiff_counter < 1)
  break;
  pdiff_counter = 0.;
  /* end check convergence */
  check = 0.;
  printf("%e ", (double)time); /* print p(t) to stdout */
  for (i = 0; i < dim; i++) {
    if(pt[i] < -0.00001) {
      fprintf(stderr, "prob of lmin %i has become negative!\n", i+1);
      exit(EXIT_FAILURE);
    }
    printf("%e ", (double)(T)abs(pt[i])); //abs for mpack
    check += abs(pt[i]);
  }
  printf("\n");

  count++; /* # of iterations */

  if ( ((check-1.) < -0.01) || ((check-1.) > 0.01) ) {
    fprintf(stderr, "overall probability at time %e is %e != 1.0 %s!\n", (double)time, (double)check, (_opt->num_err == 'R'?"rescaling":"exiting") );
    if (_opt->num_err == 'H') {
      //clean up before writing error.
      delete[] Uexp;
      delete[] pt;
      delete[] pdiff;
      delete[] Umerk;
      return 1;
      //exit(EXIT_FAILURE); //Exit in main.
    }
  }
  fill_n(pt, dim, 0);
  fill_n(pdiff, dim, 0);
  fill_n(Uexp, dim*dim, 0);
  //fill_n(U, dim*dim, 0);

}
printf("# of iterations: %d\n", count);
delete[] Uexp;
delete[] pt;
delete[] pdiff;
delete[] Umerk;

return 0;
}

/*==*/
template<typename T>
void
Calc<T>::MxFPT(T *U, T *p8, FILE *out)
{
int i,j;

//fprintf(stderr, "in MxFPT\n");
//if(_opt->want_verbose) MxPrint (U,"U" , 'm');
T *Z=NULL, *M=NULL;

//MxPrint(U, "U mxfpt", 'm');

if (_opt->absrb) {
  Z = new T[(dim-1)*(dim-1)];

  int i,j,nrhs,nfo,*ipiv=NULL;
  T *B = new T[(dim-1)];

  for (i=0; i<dim-1; i++) B[i] = 1.0;

  for(i = 0; i < dim-1; i++)
  for(j = 0; j < dim-1; j++) {
    Z[(dim-1)*i+j] = (T)((i==j && _opt->useplusI?1.0:0.0) - U[dim*i+j]); /* I-U = -Q  (U is transposed) without absorbing state*/
  }

  ipiv = (int *)malloc((dim-1)*sizeof(int));
  nrhs=1;

  int n=dim-1;
  //DGESV computes the solution to a real system of linear equations A * X = B (I think A must be transposed)
  solver->Mx_Dgesv(&n, &nrhs, Z, &n, ipiv, B, &n, &nfo);

  // fix non-zero in real absorbing state
  for(i=0; i<dim-1; i++) {
    B[i] -= B[_opt->real_abs];
  }
  dim--;
  MxFPrint(B, "Absorbing times: ", 'v', out, out!=stderr);
  dim++;

  delete[] B;
  delete[] Z;
  free(ipiv);

}
else { // non-absorbing case
  Z = new T[dim*dim];

  for(i = 0; i < dim; i++)
  for(j = 0; j < dim; j++) {
    Z[dim*i+j] = (T)((i==j&& _opt->useplusI?1.0:0.0) - U[dim*j+i] + p8[j]); /* I-U+W = W-Q  (U is transposed)*/
  }

  //if(_opt->want_verbose) MxPrint (Z,"I-U+W" , 'm');
  mxccm.minv(Z,dim, _opt->FEPS);
  //if(_opt->want_verbose)MxPrint (Z,"Fundamental matrix Z=inv(I-U+W)" , 'm');

  M = new T[dim*dim];
  for(i = 0; i < dim; i++) {
    for(j = 0; j < dim; j++) {
      M[dim*i+j] = (T)((Z[dim*j+j]-Z[dim*i+j])/p8[j]);
    }
  }

  MxFPrint(M, "First passage times (i-th row, j-th column represents fpt from i-th to j-th state)", 'm', out, out!=stderr);
  delete[] M;
  delete[] Z;
}
}

/*==*/

template<typename T>
void
Calc<T>::MxKotzOutMathematica(T *Mx)
{
int i,j;
FILE *MATHEMATICA_OUT=NULL;
const char *mathematica_file = "mxMat.txt";

MATHEMATICA_OUT = fopen(mathematica_file, "w");
if (!MATHEMATICA_OUT) {
  fprintf(stderr, "could not open file pointer 4 Matematica outfile\n");
  exit(EXIT_FAILURE);
}
fprintf(MATHEMATICA_OUT, "{");
for(i=0; i<dim; i++) {
  fprintf(MATHEMATICA_OUT, "{");
  for(j=0; j<dim; j++) {
    if (j != (dim-1))
    fprintf(MATHEMATICA_OUT, "%25.22f, ", (double)Mx[dim*j+i]);
    else
    fprintf(MATHEMATICA_OUT, "%25.22f}", (double)Mx[dim*j+i]);
  }
  if (i != (dim-1))
  fprintf(MATHEMATICA_OUT, ",\n");
  else
  fprintf(MATHEMATICA_OUT, "}\n");
}
fclose(MATHEMATICA_OUT);
}

template<typename T>
T * Calc<T>::MxFPTOneState(T *U, int state)
{
if (state>=dim) {
  fprintf(stderr, "State %d does not exist (--fpt %d)", state+1, state+1);
  return NULL;
}

int i,j,nrhs,nfo;
int *ipiv=NULL;
int n=dim-1;
T *A=NULL, *B=NULL;
A = new T[n*n];
B = new T[dim];

// A = Q(infinetisimal generator) = U^T-I with state-th column and row deleted (U is transposed!)
for(i=0; i<n; i++) {
  for (j=0; j<n; j++) {
    int ui = (i>=state?i+1:i);
    int uj = (j>=state?j+1:j);
    A[n*i+j] = (T)(U[dim*ui+uj] - (i==j && _opt->useplusI ?1.0:0.0));
  }
}

for (i=0; i<n; i++) B[i] = -1.0;

/*
 dim--;
 MxPrint(A, "A", 'm');
 MxPrint(B, "B", 'v');
 dim++;
 */

ipiv = (int *)malloc(n*sizeof(int));
nrhs=1;

//DGESV computes the solution to a real system of linear equations A * X = B (I think A must be transposed), solution goes to B
solver->Mx_Dgesv(&n, &nrhs, A, &n, ipiv, B, &n, &nfo);

for (i=dim-1; i>state; i--) {
  B[i] = B[i-1];
}
B[state] = 0.0;

delete[] A;
free(ipiv);
return B;
}

template<typename T>
inline
void Calc<T>::MxFPTSimple(T *U)
{
int i,j;
fprintf(stderr, "in MxFTPSimple\n");
T *M = new T[dim*dim];
T *B;
int p_done = 0;

for (i=0; i<dim; i++) {
  // calculate to one state and add it to global matrix
  B = MxFPTOneState(U, i);
  for (j=0; j<dim-1; j++) {
    M[dim*j+i] = B[j];
  }
  delete[] B; //free(B);
  // print progress
  if (i*100/dim >= p_done) {
    p_done = i*100/dim;
    fprintf(stderr, "%d%%\n", p_done);
  }
}

MxPrint(M, "FPT", 'm');

delete[] M;
}

template<typename T>
inline
void Calc<T>::MxFPTRnd(T *U, int packets)
{
srand(time(NULL));
T *P=NULL, *M=NULL;

P = (T *)malloc(dim*dim*sizeof(T));

int i, j;

for (i=0; i<dim; i++) {
  for (j=0; j<dim; j++) {
    //if (i==j) P[dim*i+j]
    P[dim*i+j] = U[dim*j+i]/(1.0-U[dim*i+i]);
  }
  P[dim*i+i] = 0.0;
}

//MxPrint(U, "U", 'm');
//#MxPrint(P, "P", 'm');

if (!_opt->absrb) {
  // ergodic option
  M = (T *)malloc(dim*dim*sizeof(T));
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      M[dim*i+j] = MxFPTRandom(P,U, i, j, packets);
      fprintf(stderr, "x");
    }
    fprintf(stderr, "\n");
  }

  MxPrint(M, "M (random)", 'm');
}
else {
  // absorbing option
  M = (T *)malloc(dim*sizeof(T));
  for (i=0; i<dim; i++) {
    M[i] = MxFPTRandom(P,U, i, dim-1, packets);
  }

  MxPrint(M, "M (random)", 'v');
}

free(M);
free(P);
}

template<typename T>
T Calc<T>::MxFPTRandom(T *P, T *U, int src, int dst, int packets)
{
if (src==dst) return 0.0;

T wtime = 0.0;
int i,j;
for (i=0; i<packets; i++) {
  int curr = src;
  T time = 0.0;
  while (curr != dst) {
    time += 1/(1.0-U[dim*curr+curr]);
    T next = random()/(T)RAND_MAX;
    j=0;
    while (next>0.0 && j<dim) {
      next -= P[dim*curr+j];
      j++;
    }
    j--;
    curr = j;
  }
  //fprintf(stderr, "time to consume %f\n", time);
  wtime += time;
}
return wtime/(T)packets;
}

template<typename T>
void
Calc<T>::MxEVLapackSym(T *U) {
/*   extern void dsyev_(char *jobz, char *uplo, int *n, double *A,int *lda, */
/*         double *w, double *work, int *lwork, int *info); */
T abstol;
abstol = _opt->FEPS;

if (_opt->want_verbose) fprintf(stderr, "FEPS precision %20.10g\n", (double)abstol);

T *work=NULL, vl, vu;   // unused or inputs
int il, iu, m, lwork, liwork, *iwork=NULL, *ifail=NULL, *isuppz=NULL, nfo;// unused or inputs
lwork = dim*(dim+26);
liwork = 10*dim;
work = (T *) malloc (lwork * sizeof(T));
iwork = (int *) malloc (liwork * sizeof(int));
ifail = (int *) malloc (dim * sizeof(int));
isuppz = (int *) malloc (2*dim * sizeof(int));

//if (!_opt->useplusI) for (il=0; il<dim; il++) U[il*dim+il]-=1.0;

// works only with: dim, U, evals, evecs
solver->Mx_Dsyevx("V", "A", "U",&dim, U, &dim, &vl, &vu, &il, &iu,
    &abstol, &m, evals, evecs, &dim, work, &lwork, iwork,
    ifail, &nfo);
/*
 #if 1
 dsyevx_("V", "A", "U",&dim, U, &dim, &vl, &vu, &il, &iu,
 &abstol, &m, evals, evecs, &dim, work, &lwork, iwork,
 ifail, &nfo);
 #else
 dsyevr_("V", "A", "U",&dim, U, &dim, &vl, &vu, &il, &iu,
 &abstol, &m, evals, evecs, &dim, isuppz, work, &lwork, iwork,
 &liwork, &nfo);
 #endif
 */

//dsyev_("V","U",&dim, S, &dim, evals, work, &lwork, &nfo);
//MxPrint(evecs, "Eigenvectors (bad)", 'm');
//MxPrint(evals, "Eigenvalues (bad)", 'v');
if(nfo != 0) {
  fprintf(stderr, "dsyevx exited with value %d (val=%20.10g) (Cannot compute eigenvalues) Try to:\n - lower --feps value (try --feps=-1.0)\n - try using --useplusI option (sometimes the eigenvalues are better computed like that)\n", nfo, (double)evals[nfo-1]);
  for (il=0; il<dim; il++) fprintf(stderr, "%20.10g ", (double)evecs[nfo*dim+il]);
  fprintf(stderr, "\n");

  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

//if (!_opt->useplusI) for (il=0; il<dim; il++) U[il*dim+il]+=1.0;

free(work);
free(iwork);
free(ifail);
free (isuppz);

// transpose output
mxccm.trnm(evecs, dim);
}

template<typename T>
void
Calc<T>::MxEVLapackNonSym(T *origU)
/* input: matrix origU, space for (right)evec-matrix S, */
/* since fortran sucks, we transpose the input matrix   */
/* and still want the right eigenvectors                */
/* matrix of right eigenvecs is transposed  to have */
/* the j-th evec in column c */
{
T *evals_re=NULL, *evals_im=NULL, *scale=NULL;
T *rconde=NULL, *rcondv=NULL, *work=NULL, abnrm;
int one, ilo, ihi, lwork, *iwork, nfo;
//int dimx2;
/* for sorting */
int i;

dim = dim+500;
one = 1;
//dimx2 = 2*dim;
lwork = dim*(dim+6);

evals_re = new T[dim];
evals_im = new T[dim];
scale = new T[dim];
rconde = new T[dim];
rcondv = new T[dim];
work = new T[lwork];
iwork = new int[2*(dim -2)];
dim=dim-500;
if ( (evals_re && evals_im && scale && rconde &&
        rcondv && work && iwork)==0) {
  fprintf(stderr,"no space for temporary lapack arrays!\n");
  exit(EXIT_FAILURE);
}

/* instead of more fiddling, we transpose the input */
mxccm.trnm(origU, dim);

T* tmpNull=NULL;
solver->Mx_Dgeevx("B", "N", "V", "V", &dim, origU, &dim, evals_re, evals_im, tmpNull, &one,
    evecs ,&dim, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work,
    &lwork, iwork, &nfo);

for (i=0; i<dim; i++) evals[i]=evals_re[i];
for (i=0; i<dim; i++)
if ((evals_re[i] != 0.0) && abs(evals_im[i]/evals_re[i])>_opt->FEPS)
fprintf(stderr,"eigenvalue %d is %g + i*%g, which is somewhat complex\n",
    i, (double)evals_re[i], (double)evals_im[i]);

/*transpose output*/
mxccm.trnm(evecs, dim);

delete[] evals_re;
delete[] evals_im;
delete[] scale;
delete[] rconde;
delete[] rcondv;
delete[] work;
delete[] iwork;
}

template<typename T>
inline
void Calc<T>::MxEqDistrFromDetailedBalance ( T *U, T **p8 )
{
T *res = *p8;
int *done;
int count = 1;  // num of solved states
int i=0;
int j=1;
int k;

done = new int[dim];
std::fill_n(done,dim,0);

for(k=1; k<dim; k++) {
  res[k] = 0.0;
}
res[0] = 1.0;

// while there are unsolved states
while (count != dim) {
  for (j=1; j<dim; j++) {
    if (res[j]>(T)0.) continue;
    if (U[dim*i+j]>(T)0.0) {
      if (U[dim*j+i]==(T)0.0) {
        fprintf(stderr, "Local balance is unsatisfiable at U[%d][%d]=%f\n", i, j, (double)U[dim*i+j]);
        exit(EXIT_FAILURE);
      }
      // local balance
      T tmp = U[dim*j+i];
      tmp *= res[i];
      tmp /= U[dim*i+j];
      res[j] = tmp;
      count++;
      if (res[j]>(T)1. || res[j]<=(T)0.) {
        fprintf(stderr, "Weird eqilibrium pop for state %d", j);
        MxPrint(res, "p8 (local balance)", 'v');
      }
    }
  }
  done[i] = 1;
  for (i=1; i<dim; i++) {
    if ((!done[i]) && (res[i]>(T)0.)) break;
  }

  if (i==dim && count<dim) {
    fprintf(stderr, "non-ergodic chain\n");
    MxPrint(res, "p8 (local balance)", 'v');
    exit(EXIT_FAILURE);
  }
}

delete[] done;

/* now make the vector stochastic (sum = 1.0) */
T sum=0.0;
for (i=0; i<dim; i++) {
  sum += res[i];
}
for (i=0; i<dim; i++) {
  res[i] /= sum;
}

MxPrint(res, "p8 (detailed balance)", 'v');
}

/*==*/
template<typename T>
void*
Calc<T>::MxNew ( size_t size )
{
void *mx = NULL;
if ( (mx = (void *) calloc (1, size)) == NULL )
fprintf (stderr, "ERROR: new_martix() allocation failed\n");

return mx;
}

/*==*/
template<typename T>
void
Calc<T>::MxMemoryCleanUp (void)
{
if(_opt->sequence) free(_opt->sequence);
if(_opt->basename) free(_opt->basename);
if(_opt->fpt_file) free(_opt->fpt_file);
delete[] _sqrPI;
delete[] sqrPI_;
delete[]evecs;
delete[]evals;
}

/*==*/
template<typename T>
void
Calc<T>::print_settings(void)
{
printf(
    "# Date: %s"
    "# Sequence: %s\n"
    "# Method: %c  Start Time: %.2f  Stop Time: %.2f Temperature: %.2f\n",
    time_stamp(),
    _opt->sequence,
    _opt->method,
    _opt->t0,
    _opt->t8,
    _opt->T
);
if(_opt->basename != NULL) printf("# basename: %s\n",_opt->basename);
else printf("# basename: <stdin>\n");
if (_opt->tinc) printf("# time increment: %.2f\n", _opt->tinc);
else printf("# time increment: %.2f \n", _opt->tinc);
if(_opt->want_degenerate == 1)printf("# degeneracy:  on\n");
else printf("# degeneracy:  off\n");
if (_opt->absrb < 1) printf("# absorbing lmin: none\n");
else printf("# absorbing lmin: %d\n", _opt->absrb);
if (_opt->n > 0) printf("# nlmins: %d\n", _opt->n);
else printf("# nlmins: till EOF\n");
}

/*==*/
template<typename T>
char*
Calc<T>::time_stamp(void)
{
time_t cal_time;
cal_time = time(NULL);
return ( ctime(&cal_time) );
}

/*==*/
template<typename T>
void
Calc<T>::MxDoDegeneracyStuff(void)
{
int i, j, b, nr, current, numsad = 1, count = 0;
TypeDegSaddle *saddle = NULL;
numsad = Barparser::ParseSaddleFile(&saddle);
/* loop over all elements of structure-array saddle: */
/* first we fill the upper triangle */
for (count = 0; count < numsad; count++) {
  current = 1;
  nr = saddle[count].list[0];
  /* only for saddles with a cc >= 1  AND those which connect at least 2 lmins */
  if(saddle[count].cc >= 1 && nr >= 2 && !(saddle[count].cc == 1 && nr == 2)) {
    while(current < nr) {
      for(b = current+1; b <= nr; b++) {
        /* skip in case a lmin which we don't see is connected by */
        /* the saddle (because we can only see --max x lmins) */
        if(saddle[count].list[current] > dim || saddle[count].list[b] > dim) {
          current++;
          continue;
        }
        /* FIRST: we consider the size of the cc the saddle belongs to */
        if(saddle[count].cc > 1) {
          D[dim * (saddle[count].list[current]-1) + (saddle[count].list[b]-1)] += (saddle[count].cc - 1);
          fprintf(stderr, "transition between %3d - %3d: adding %2d  cc\n",
              saddle[count].list[current], saddle[count].list[b], saddle[count].cc-1);
          /* -1 because if the size of cc == 1 there's nothing special */
          /* about it, i.e. there is just one saddle */
        }
        /* SECOND: we consider that the saddle connects several lmins */
        if(nr > 2) {
          D[dim * (saddle[count].list[current]-1) + (saddle[count].list[b]-1)] +=1;
          fprintf(stderr, "transition betweed %3d - %3d: adding  1 deg_saddle\n",
              saddle[count].list[current], saddle[count].list[b] );
          /* -1 because matrix starts with0, NOT with 1 */
        }
      }
      current++;
    }
  }
}

for(i = 0; i < dim; i++) /* make matrix symmetric */
{
  for(j = 0; j < dim; j++)
  if (i != j)
  D[dim*j+i] = D[dim*i+j];
}
if (_opt->want_verbose) {
  sprintf (Aname, "%s", "D (degeneracies)");
  MxPrint (D, Aname, 'm');
}
free(saddle);
}



#endif

/* End of file */
