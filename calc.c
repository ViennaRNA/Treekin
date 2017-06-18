/*=================================================================*/
/*=   calc.c                                                      =*/
/*=   main calculation and iteration routines for treekin         =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2017-06-02 16:53:08 ivo>          =*/
/*=   $Id: calc.c,v 1.41 2006/11/27 23:01:45 mtw Exp $            =*/
/*=   ---------------------------------------------------------   =*/
/*=     (c) Michael Thomas Wolfinger, W. Andreas Svrcek-Seiler    =*/
/*=                  {mtw,svrci}@tbi.univie.ac.at                 =*/
/*=                             treekin                           =*/
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include "exp_matrix.h" /* functions for matrix-exponent stuff */
#include "mxccm.h"      /* functions for eigen-problems stolen from ccmath */
#include "barparser.h"  /* functions for input */
#include "calc.h"       /* does all matrix stuff for markov process */
#include "globals.h"    /* contains getopt-stuff */

#include "calcpp.h"

#define SQ(X) ((X)*(X))

/* private function(s) */
static double *MxMethodeA (BarData *Data);
static double *MxMethodeFULL(double *);
static double *MxMethodeINPUT (BarData *Data, double *);
static double  max_saddle(int i, int j, BarData *Data);
static void    print_settings(void);
static char   *time_stamp(void);
static void    MxDoDegeneracyStuff(void);
static void    MxBinWrite(double *Mx, char what[], char T);
static int     MxBinRead(double** Mx, char what[], char T);
static void    MxASCIIWrite(double *Mx, char *asciifile);
static void    MxASCIIWriteV(double *Mx, char *asciifile);
static void    MxKotzOutMathematica(double *Mx);
static void    MxSortEig(double *evals, double *evecs);
static void    MxEVLapackSym(double *U);
static void    MxEVLapackNonSym(double *U);
static void    MxFixevecs(double *, double *);
static void    MxFixevecsAbsorb(double *, double *);
static void    MxDiagHelper(double *P8);

void MxFPrintD(double *mx, char *name, int dim1, int dim2, FILE *out);

/* private vars and arrays */
static int      dim = 0;
static double   _kT = 1.;
static double  *evals     = NULL;
static double  *evecs     = NULL;
static double  *_sqrPI    = NULL;   /* left constant array */
static double  *sqrPI_    = NULL;   /* right constant array */
static double  *D         = NULL;   /* matrix with degree of degeneracy */
static char     Aname[30];
static TypeDegSaddle *saddle = NULL;

/*==*/
void
MxInit (int d)
{
  _kT = 0.00198717*(273.15 + opt.T);
  if (d > 0 ) dim = d;
  else {
    fprintf(stderr, "dim <= 0\n");
    exit(EXIT_FAILURE);
  }
}

/*==*/
double*
MxBar2Matrix ( BarData *Data, double *R)
{
  double *U=NULL;
  if(opt.want_degenerate) MxDoDegeneracyStuff();
  switch (opt.method) {
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
             "ERROR in MxBar2Matrix(): No handler 4 method %c\n", opt.method);
    exit(EXIT_FAILURE);
  }
  if (opt.dumpU) {
    MxASCIIWrite(U, "U.txt");
    MxBinWrite(U, "U", 'm');
  }
  return (U);
}
/*==*/
void
MxGetSpace (double **p8)
{
  *p8   = (double *) MxNew (dim*sizeof(double));
  evals = (double *) MxNew (dim*sizeof(double));
  evecs = (double *) MxNew (dim*dim*sizeof(double));
  assert(evals!=NULL);
  assert(evecs!=NULL);
  assert(p8!=NULL);

  if(!opt.absrb) {
    _sqrPI = (double *) MxNew (dim*dim*sizeof(double));
    sqrPI_ = (double *) MxNew (dim*dim*sizeof(double));
  }
}

/*==*/
void
MxStartVec (double **p0)
{
  int i;
  double *pzero = NULL;

  pzero = (double *) MxNew(dim*sizeof(double));

  if (opt.pini) {
    for (i = 1; i < (int) *opt.pini; i+=2)
      pzero[(int)opt.pini[i]-1] = (double)opt.pini[i+1];
    /* -1 because our lmins start with 1, not with 0 (as Data does ) */
  } else {
    // all into first state...
    pzero[0]=1.0;
  }

  if (opt.want_verbose) MxPrint (pzero, "p0", 'v');
  *p0=pzero;
}

/*==*/
/* calculate equilibrium distribution */
void
MxEqDistr ( BarData *Data, double **p8 )
{
  int i;
  double Z = 0.;

  if(opt.absrb) {
    for(i = 0; i < dim; i++)
      (*p8)[i] = 0.;
    (*p8)[dim-1] = 1.0; /* last entry is the 'new' absorbing state */
  }
  else {
    for(i = 0; i < dim; i++) Z += exp(-((double)Data[i].energy/_kT));
    for(i = 0; i < dim; i++) (*p8)[i] = exp(-((double) Data[i].energy/_kT));
  }

  /* now normalize the new p8 */
  double sumsq=0.0;
  for (i=0; i<dim; i++)
    sumsq += SQ((*p8)[i]);
  if(sumsq > 0.0)
    sumsq=1./sqrtl(sumsq);
  for (i=0; i<dim; i++)
    *(*p8+i)  *= sumsq;

  if(opt.want_verbose) MxPrint (*p8, "p8", 'v');
  return;
}

/*==*/
void
MxEqDistrFULL (SubInfo *E, double *p8 )
{
  int i;
  double Z = 0.;

  if(opt.absrb) {
    for(i = 0; i < dim; i++)
      p8[i] = 0.;
    p8[opt.absrb-1] = 1.;
  }
  else {
    for(i = 0; i < dim; i++) Z += exp(-E[i].energy/_kT);
    for(i = 0; i < dim; i++) p8[i] = exp(-E[i].energy/_kT)/Z;
  }


  if(opt.want_verbose) MxPrint (p8, "p8", 'v');
  return;
}

void
MxEqDistrFromLinSys( double *U, double **p8 )
{
  extern void dgesv_(int *N,int *NRHS,double *A,int *LDA,int *IPIV,
                     double *B,int *LDB,int *INFO);

  int i,j,n, nrhs, nfo, *ipiv=NULL;
  double *A=NULL, *B=NULL;

  if(opt.absrb) {
    for(i = 0; i < dim; i++)
      *(*p8+i) = 0.;
    *(*p8+(dim-1)) = 1.0; /* last entry is the 'new' absorbing state */
  }
  else {
    n=dim-1;
    A    = (double *)malloc(n*n*sizeof(double));
    B    = (double *)malloc(n*sizeof(double));
    ipiv = (int *)malloc(n*sizeof(int));
    nrhs=1;

    if (opt.want_verbose) MxPrint(U, "U", 'm' );

    for(i=1; i<=n; i++) /* all except first row */
      for(j=1; j<=n; j++)
        A[n*(i-1)+(j-1)]=U[dim*i+j] - (i==j && opt.useplusI ? 1.0 : 0.0);  //U-I=Q
    for(n=0,i=1; i<dim; i++,n++)
      B[n]=-U[dim*i];
    dim=n;
    trnm(A,n);

    if (opt.want_verbose) MxPrint(A, "A in MxEqDistrFromLinSys", 'm' );
    if (opt.want_verbose) MxPrint(B, "B in MxEqDistrFromLinSys", 'v' );

    //DGESV computes the solution to a real system of linear equations A * X = B
    dgesv_(&n, &nrhs, A, &n, ipiv, B, &n, &nfo);
    if (nfo != 0) {
      fprintf(stderr, "dgesv exited with value %d (Cannot compute the equilibrium distribution) - check if the rates have good dimension(if you, you have probably reached the numerical precision of treekin)\n", nfo);
      /// TODO switch to compute from detailed balance
      exit(EXIT_FAILURE);
    }
    if (opt.want_verbose) MxPrint(B, "X in MxEqDistrFromLinSys", 'v' );
    dim=n+1;
    *p8[0]=1.;
    for(i=1; i<dim; i++) *(*p8+i)=B[i-1];

    if (opt.want_verbose) MxPrint(*p8, "p8 in MxEqDistrFromLinSys before norm", 'v' );

    // now check if all are > 0.0
    for (i=dim-1; i!=0; i--) {
      if ((*p8)[i]<0.0) {
        if (i==0) exit(EXIT_FAILURE);
        fprintf(stderr, "Warning: p8[%5d] is negative (%e), setting it to a nearby value (%e)\n", i+1, (*p8)[i], (i==dim-1?(*p8)[i-1]:(*p8)[i+1]));
        (*p8)[i] = fabs((i==dim-1?(*p8)[i-1]:(*p8)[i+1]));
      }
      if ((*p8)[i]==0.0) {
        fprintf(stderr, "Warning: p8[%5d] is zero (%e), setting it to a minimum value (%e)\n", i+1, (*p8)[i], 1e-20);
        (*p8)[i] = 1e-20;
      }
    }


    /* now make the vector stochastic (sum = 1.0) */
    long double sum=0.0;
    for (i=0; i<dim; i++) {
      sum += (*p8)[i];
    }
    for (i=0; i<dim; i++) {
      (*p8)[i] /= sum;
    }

    // now normalize the new p8
    /*long double sumsq=0.0;
    for (i=0; i<dim; i++)
      sumsq += SQ((*p8)[i]);
    if(sumsq > 0.0)
      sumsq=1./sqrtl(sumsq);
    for (i=0; i<dim; i++)
      *(*p8+i)  *= sumsq;
*/

    free(A);
    free(B);
    free(ipiv);
  }
  if(opt.want_verbose)
    MxPrint(*p8, "p8", 'v' );
}

/*==*/
void
MxDiagonalize ( double *U, double **_S, double *P8)
{
  int i,j;
  double *tmpMx=NULL;

  if(opt.dumpMathematica == 1)  MxKotzOutMathematica(U);
  if(!opt.absrb) {
    if (opt.want_verbose) MxPrint (U, "input U", 'm');

    tmpMx = (double *) MxNew (dim*dim*sizeof(double));
    //if (opt.want_verbose) MxPrint (P8, "P8", 'v');
    MxDiagHelper(P8);
    mmul(tmpMx, _sqrPI, U, dim);
    if (opt.want_verbose) MxPrint (_sqrPI, "_sqrPI", 'm');
    if (opt.want_verbose) MxPrint (tmpMx, "tmpMx = _sqrPI*U", 'm');
    if (opt.want_verbose) MxPrint (sqrPI_, "sqrPI_", 'm');
    memset(U,0,dim*dim*sizeof(double));
    mmul(U, tmpMx, sqrPI_, dim);
    if (opt.want_verbose) MxPrint (U, "U = _sqrPI*U*sqrPI_ (this should be symmetric, but still uncorrected)", 'm');
    free(tmpMx);

    /* correct for numerical errors -- ensure U is symmetirc */
    long double err = 0.0;  // acumulated error
    for (i = 0; i < dim; i++) {
      for (j = i; j < dim; j++) {
        long double err_inc = fabs((U[dim*i+j]+U[dim*j+i])/2 - U[dim*i+j]);
        /*if (isnan(err_inc)) {
          fprintf(stderr, "%d %d\n", i,j);
        }*/
        err += err_inc;
        if (isnan(err_inc)) {
          fprintf(stderr, "Weird rates! r(%5d,%5d)=(%e, %e)\n", i, j, U[dim*i+j], U[dim*j+i]);
          //exit(-1);
        }
        U[dim*i+j] = (U[dim*i+j]+U[dim*j+i])/2;
	//U[dim*i+j] = sqrt(U[dim*i+j]*U[dim*j+i]);
        U[dim*j+i] = U[dim*i+j];
      }
    }
    if (isnan(err)) {
      fprintf(stderr, "Rates are not good -- check your rates (maybe the equilibrium was not computed correctly? rerun with --equil-file or -v set)!!!... (%Lf)\n", err);
      exit(-1);
    }
    if (!opt.quiet) fprintf(stderr, "Corrected numerical error: %e (%e per number)\n", (double)err, (double)(err/(long double)(dim*dim)));


    if (opt.want_verbose) MxPrint (U, "force symmetrized U", 'm');
  }

    // get eigenv*
  if(opt.absrb) {
    MxEVLapackNonSym(U);
  } else {
    MxEVLapackSym(U);
  }

  MxSortEig(evals, evecs);

  /*// ### TEST EIGEN
  TestEigen(U, dim, NULL, NULL);
  MxPrint(evecs, "Eigenvectors of U (LAPACK)", 'm');
  MxPrint(evals, "Eigenvalues of U (LAPACK)", 'v');
  /// ### END test*/

  if (opt.want_verbose) {
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

  //for (i=0; i< 10; i++) fprintf(stderr, "%30.20g\n", evals[i]);

  if (!opt.quiet) {
    if (abs(evals[0] - ((opt.useplusI)?1:0)) > 10*opt.FEPS) {
      fprintf(stderr, "WARNING largest eigenvalue is %g, see evals.txt and evecs.txt\n", evals[0]);
      MxASCIIWriteV(evals, "evals.txt");
      MxASCIIWrite(evecs, "evecs.txt");
    }
  }
  // fix evals[0] to be 0 (or 1.0)
  if (opt.useplusI) evals[0] = 1.0;
  else evals[0] = 0.0;

  // check if we have only one eval = 0 / eval = 1
  i = 0;
  if (opt.useplusI) while (evals[i] >= 1.0) evals[i++] = 1.0;
  else              while (evals[i] >= 0.0) evals[i++] = 0.0;

  if (i > 1) {
    fprintf(stderr, "WARNING: There are multiple evals=%d.0 (%4d) (multiple population traps -- maybe we compute crap)\n", opt.useplusI?1:0, i);
  }

  /*
  for (i=0; i<dim; i++)
    if (evals[i]>(1.0-opt.FEPS)) evals[i]=1.0;
  */

  //for (i=0; i< 10; i++) fprintf(stderr, "%30.20g\n", evals[i]);

  if(opt.absrb)
    MxFixevecsAbsorb(evecs,evals);
  //  else
  //  MxFixevecs(evecs,evals); /* so far it's not helping ... */
  
  *_S=evecs;
  if(opt.wrecover) {
    MxBinWrite(evals, "evals", 'v');
    MxBinWrite(evecs, "evecs", 'm');
    MxASCIIWriteV(evals, "evals.txt");
    MxASCIIWrite(evecs, "evecs.txt");
    if (opt.want_verbose) {
      double *CL, *tmpMx;
      int i,j;
      CL        = (double *) MxNew (dim*dim*sizeof(double));
      mmul (CL, sqrPI_, evecs, dim);
      MxASCIIWrite(CL, "evecR.txt");
      fprintf(stderr, "Sums of EVs of rate matrix R\n");
      for (i=0; i<dim; i++) {
	double sum, sum2, sum3;
	sum=sum2=sum3=0.;
	for (j=0; j<dim; j++) {
	  sum += CL[dim*j+i];
	}
	fprintf(stderr, "%15.8g  ", i, sum);
      }
      fprintf(stderr, "\n");
    free(CL);
    }
  }
  // compensate for the +I matrix
  if (opt.useplusI) for(i = 0; i < dim; i++) evals[i] -= 1.0;  /* compensate 4 translation of matrix U */
}

/*==*/
static void
MxDiagHelper(double *P8)
{
  int i,j;
  for(i = 0; i < dim; i++)
    for(j = 0; j < dim; j++)
      if( i == j) {
        sqrPI_[dim*i+j] = sqrt(P8[i]);          /* pos right */
        _sqrPI[dim*i+j] = 1/(sqrPI_[dim*i+j]);  /* neg left */
      }
}

/*==*/
void
MxRecover (double **_S, double *P8)
{
  MxDiagHelper(P8);
  free(evecs);
  free(evals);

  if(MxBinRead(&evecs, "evecs", 'm') != dim ) {
    fprintf(stderr, "ERROR: MxBinRead() returns wrong dimension for evecs\n");
    exit(EXIT_FAILURE);
  }
  if(MxBinRead(&evals, "evals", 'v') != dim) {
    fprintf(stderr, "ERROR: MxBinRead() returns wrong dimension for evals\n");
    exit(EXIT_FAILURE);
  }
  *_S=evecs;
  if(opt.want_verbose) {
    MxPrint(evals, "MxRecover: Eigenvalues", 'v');
    MxPrint(evecs, "MxRecover: Eigenvectors", 'm');
  }
  return;
}

/*==*/
/* S which comes into this function contains the (right) eigenvectors
calculated by LAPACK */
void
MxIterate (double *p0, double *p8, double *S)
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
  double time, check = 0.;
  double *CL = NULL, *CR, *exptL, *tmpVec, *tmpVec2, *pt, *St, *pdiff;
  double *ptFULL = NULL;  /* prob dist 4 of the effective lmins of the tree at time t */
  double *p8FULL = NULL;  /* equ dist 4 gradient basins, full process */

  tmpVec    = (double *) MxNew (dim*sizeof(double));
  tmpVec2   = (double *) MxNew (dim*sizeof(double));
  pt        = (double *) MxNew (dim*sizeof(double));
  exptL     = (double *) MxNew (dim*dim*sizeof(double));
  if(opt.method=='F') {
    ptFULL    = (double *) MxNew ((lmins+1)*sizeof(double));
    p8FULL    = (double *) MxNew ((lmins+1)*sizeof(double));
    pdiff     = (double *) MxNew ((lmins+1)*sizeof(double));
  }
  else pdiff   = (double *) MxNew (dim*sizeof(double));

  if(! opt.absrb) { /* NON-absorbing case */
    CL        = (double *) MxNew (dim*dim*sizeof(double));
    CR        = (double *) MxNew (dim*dim*sizeof(double));
    St        = (double *) MxNew (dim*dim*sizeof(double));
    mcopy(St, S, dim*dim);
    trnm(St, dim); /* transpose S (same as invert(S) ('cause of symmetry) minv(St, dim); )*/
    mmul (CL, sqrPI_, S, dim);
    mmul (CR, St, _sqrPI, dim);
    vmul (tmpVec, CR, p0, dim);

    // test CL*CR = I
    if (opt.want_verbose) {
      double *testMx, *tvec;
      testMx      = (double *) MxNew (dim*dim*sizeof(double));
      tvec    = (double *) MxNew (dim*sizeof(double));
      mmul(testMx,CL,CR,dim);
      MxPrint(testMx, "CL*CR should be I", 'm');
      free(testMx);
      MxPrint(tmpVec, "CR*p0 should start with 1", 'v');
      vmul(tvec, _sqrPI, p0, dim);
      vmul(tmpVec2, St, tvec, dim);
      MxPrint(tmpVec2, "CR*p0 should start with 1", 'v');
      MxPrint(CL, "CL", 'm');
      free(tvec);
    }
    free(St);
    free(CR);

    //MxPrint(S, "S", 'm');
    //MxPrint(_sqrPI, "sqrPI", 'm');
    //MxPrint(CR, "CR", 'm');
    //MxPrint(tmpVec, "tmpVec", 'v');
  }
  else { /* absorbing case */
    double *S_inv;
    S_inv  = (double *) MxNew (dim*dim*sizeof(double));
    mcopy(S_inv, S, dim*dim);
    minv(S_inv,dim);
    if(opt.want_verbose) MxPrint(evals, "evals in MxIterate", 'v');
    vmul (tmpVec, S_inv, p0, dim);
    free(S_inv);
  }
  if(opt.method=='F') { /* calculate equilibrium distribution once */
    for (i = 0; i < dim; i++)   p8FULL[E[i].ag] += p8[i];
    for (i = 0; i < lmins; i++) check += fabs(p8FULL[i]);
    if ( ((check-1) < -0.1) || ((check-1) > 0.1) ) {
      fprintf(stderr, "overall equilibrium probability is %e != 1. !\n", check);
      if (opt.num_err == 'H' || check == 0.0) exit(EXIT_FAILURE);
      else if (opt.num_err == 'R') {
        for (i=0; i<dim; i++) p8[i] /= check;
        check = 1.0;
      }
    }
  }
  check = 0.;

  /* solve fundamental equation */
  print_settings();
  if (opt.t0 == 0.0) {
    if (opt.method=='F')  PrintProbFull(p0, dim, 0.0, lmins);
    else                  PrintProb(p0, dim, 0.0);
    opt.t0 = TZERO;
  }

  double underflow[dim], tt;
  for (i=0; i<dim; i++) underflow[i] = 0.0;

  // iterate
  for (tt = opt.t0; tt < opt.t8*opt.tinc; tt *= opt.tinc) {
    time = (tt<opt.t8)?tt:opt.t8;
    for (i = 0; i < dim; i++) {
      errno = 0;
      exptL[dim*i+i] = exp(time/opt.times*evals[i]);
      if ((errno == ERANGE || isnan(exptL[dim*i+i])) && underflow[i]==0.0) {
        //if (opt.warnings) fprintf(stderr, "WARNING: underflow occured on %dth state at time %g -- exp(%g * %g) = %g\n", i+1, time/opt.times, time/opt.times, evals[i], exptL[dim*i+i]);
                //         the overall probability can start to decrease if this state is still populated!!!\n         p_%d(%g) = %g, so it seems this %s\n", i+1, time/opt.times, time/opt.times, evals[i], exptL[dim*i+i], i+1, time/opt.times, pt[i], pt[i]>0.1?"is DEFINITELLY BAD":(pt[i]>0.001?"is POTENTIALLY BAD":"SHOULD BE OK"));
        underflow[i] = time/opt.times;
        //exit(EXIT_FAILURE);
      }
    }
    vmul(tmpVec2, exptL, tmpVec, dim);
    if(!opt.absrb)  vmul(pt, CL, tmpVec2, dim);
    else            vmul(pt, S, tmpVec2, dim);

    /*if (count%100 == 0) {
      fprintf(stderr, "Time: %g\n", time/opt.times);
      MxPrint(tmpVec2, "tmpVec2", 'v');
      MxPrint(pt, "pt", 'v');
    }*/

    count++;  /* # of iterations */

    if(opt.method=='F') {
      memset(ptFULL, 0, (lmins+1)*sizeof(double));
      for (i = 0; i < dim; i++) {
        ptFULL[E[i].ag] += pt[i];
      }
    }

    // print probabilities with respect to corrected ergodicity
    if (opt.method=='F') check = PrintProbFull(pt, dim, time, lmins);
    else                 check = PrintProb(pt, dim, time);

    //PrintProbNR(p8FULL, lmins, -1);
    int reached;
    if (opt.method=='F') reached = ConvergenceReached(p8FULL, ptFULL, lmins, 1);
    else                 reached = ConvergenceReached(p8, pt, dim, 0);
    fflush(stdout);
    if (reached) break;
  }

  // print underflow:
  if (opt.warnings) {
    for (i=0; i<dim; i++) {
      if (underflow[i] > 0.0) fprintf(stderr, "underflow %5d at time %12g", i+1, underflow[i]);
    }
  }


  if (time < opt.t8) {
    if (opt.method=='F') PrintProbFull(pt, dim, opt.t8, lmins);
    else                 PrintProb(pt, dim, opt.t8);
  }
  printf("# of iterations: %d\n", count);

  /*** end solve fundamental equation ***/

  if(opt.method=='F') {
    free(ptFULL);
    free(p8FULL);
    free(E);
  }
  free(evals);
  free(exptL);
  free(CL);
  free(tmpVec);
  free(tmpVec2);
  free(pdiff);
  if(pt != NULL) free(pt);
}

/*==*/
static double*
MxMethodeA (BarData *Data)
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
  double m_saddle, Zabs, abs_rate, *T, *U;

  U = (double *) MxNew (dim*dim*sizeof(double));

  for( i = 0; i < dim; i++) {
    for( j = i+1; j < dim; j++) {
      m_saddle = max_saddle(i, j, Data);
      /* rate j -> i */
      U[dim*i+j] = 1.0*exp(-(m_saddle-Data[j].energy)/_kT);
      /* rate i -> j */
      U[dim*j+i] = 1.0*exp(-(m_saddle-Data[i].energy)/_kT);
    }
  }

  if(opt.absrb) { /*==== absorbing  states ====*/
    dim++;
    fprintf(stderr, "dim increased to %i\n", dim);
    T = (double *) MxNew(dim*dim*sizeof(double));
    real_abs = opt.absrb; /* the original absorbing lmin */
    real_abs--;
    opt.absrb = dim; /* the 'new' abs state = last row/column of rate matrix */
    fprintf(stderr, "new absorbing lmin is: %i\n", opt.absrb);
    Zabs = exp((-Data[real_abs].energy)/_kT);
    abs_rate = exp((-Data[real_abs].energy)/_kT)/Zabs;

    for(i = 0; i < (dim-1); i++) { /* all except the last row */
      for(j = 0; j < (dim-1); j++)
        T[dim*i+j] = U[(dim-1)*i+j];
      T[(dim-1)*j+(dim-1)] = 0.;
    }
    for(j = 0; j < dim; j++) /* last row */
      T[dim*(dim-1)+j] = 0.;
    T[dim*(dim-1)+real_abs] = abs_rate;
    free(U);
    U = T;
    if(opt.want_verbose) MxPrint(U, "aufgeblasene Matrix", 'm');

  }  /*== end absorbing states ==*/

  /* set diagonal elements to 0 */
  for (i = 0; i < dim; i++) U[dim*i+i] = 0;
  for (j = 0; j < dim; j++) {
    double tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++) tmp += U[dim*i+j];
    U[dim*j+j] = -tmp; /* make U a stochastic matrix */
    if (opt.useplusI) U[dim*j+j] += 1.0;
  }
  if (opt.want_verbose) MxPrint (U,"U with Methode A", 'm');
  return (U);
}

/*==*/
static double*
MxMethodeFULL (double *R)
{
  int i, j;
  free(D);

  if(opt.absrb) { /*==== absorbing  states ====*/
    for(i = 0; i < dim; i++)
      R[dim*i+(opt.absrb-1)] = 0. ;
  }              /*== end absorbing states ==*/

  /* set diagonal elements  to 0 */
  for (i = 0; i < dim; i++) R[dim*i+i] = 0;
  for (j = 0; j < dim; j++) {
    double tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++) tmp += R[dim*i+j];
    R[dim*j+j] = -tmp; /* make U a stochastic matrix */
    if (opt.useplusI) R[dim*j+j] += 1.0;
  }

  if (opt.want_verbose) MxPrint (R, "R with Methode F", 'm');
  return R;
}

/*==*/
static double*
MxMethodeINPUT (BarData *Data, double *Input)
{
  int i, j;
  double *U=NULL, Zabs, abs_rate;

  if (opt.want_verbose) MxPrint(Input, "Input Matrix", 'm');

  if (opt.absrb) {  /*==== absorbing  states ====*/
    dim++;
    fprintf(stderr, "dim increased to %i\n", dim);
    U = (double *) MxNew(dim*dim*sizeof(double));
    opt.real_abs = opt.absrb; /* the original absorbing lmin */
    opt.real_abs--;
    opt.absrb = dim; /* the 'new' abs state = last row/column of rate matrix */
    fprintf(stderr, "new absorbing lmin is: %i\n", opt.absrb);
    Zabs = exp((-Data[opt.real_abs].FGr)/_kT);
    abs_rate = exp((-Data[opt.real_abs].energy)/_kT)/Zabs;

    for(i = 0; i < (dim-1); i++) { /* all except the last row */
      for(j = 0; j < (dim-1); j++)
        U[dim*i+j] = Input[(dim-1)*i+j];
      U[(dim-1)*j+(dim-1)] = 0.;
    }
    for(j = 0; j < dim; j++) /* last row */
      U[dim*(dim-1)+j] = 0.;
    U[dim*(dim-1)+opt.real_abs] = abs_rate;
    /*   if(opt.want_verbose) MxPrint(U, "aufgeblasene Matrix", 'm'); */
  }      /*== end absorbing states ==*/
  else { /*== non-absorbing states ==*/
    U = (double *) MxNew(dim*dim*sizeof(double));
    memcpy(U, Input, dim*dim*sizeof(double));
  }      /*== end non-absorbing states ==*/

  /* diagonal elements */
  for (i = 0; i < dim; i++) U[dim*i+i] = 0;
  //fprintf(stderr, "dim is %i\n", dim);
  for (i = 0; i < dim; i++) {
    double tmp = 0.00;
    // calculate column sum
    for(j = 0; j < dim; j++)  tmp += U[dim*j+i];
    U[dim*i+i] = -tmp;   // make U a stochastic matrix U = Q+I ??
    if (opt.useplusI) U[dim*i+i] += 1.0;
  }

/*  // normalize each column to sum=1
  for (i = 0; i < dim; i++) {
    double tmp = 0.00;
    for(j = 0; j < dim; j++) tmp += U[dim*j+i];
    for(j = 0; j < dim; j++) U[dim*j+i] /= tmp;
  }*/


  if(opt.want_verbose) MxPrint (U,"U with Methode I" , 'm');
  free(Input);
  return U;
}

/*==*/
static double
max_saddle(int i, int j, BarData *Data)
{
  int tmp;

  if(Data[i].number > Data[j].father) { /* exchange i & j */
    tmp = i;
    i = j;
    j = tmp;
  }
  if(Data[i].number == Data[j].father) return Data[j].energy + Data[j].ediff;
  else {
    if((Data[j].energy + Data[j].ediff) > max_saddle((Data[j].father - 1), (Data[i].number - 1), Data)) return (Data[j].energy + Data[j].ediff);
    else return max_saddle((Data[j].father - 1), (Data[i].number - 1), Data);
  }
}

/*==*/
void*
MxNew ( size_t size )
{
  void *mx = NULL;
  if ( (mx = (void *) calloc (1, size)) == NULL )
    fprintf (stderr, "ERROR: new_martix() allocation failed\n");

  return mx;
}

void MxFPrintD(double *mx, char *name, int dim1, int dim2, FILE *out)
{
  int k, l;
  fprintf(out,"%s:{\n", name);
  for (k = 0; k < dim1; k++) {
    if (k!=0) fprintf(out, ",");
    fprintf(out, "{");
    for (l=0; l< dim2; l++) {
      if (l!=0) fprintf(out, ",");
      fprintf(out,"%15.7g (%4d) ", mx[dim2*k+l], dim2*k+l);
    }
    fprintf(out,"}\n");
  }
  fprintf(out,"}-----------\n");
}

void MxFPrint(double *mx, char *name, char T, FILE *out, int pure)
{
  int k, l;
  switch (T) {
  case 'm':    /* square matrix */
    if (!pure) fprintf(out,"%s:\n", name);
    for (k = 0; k < dim; k++) {
      for (l=0; l< dim; l++) {
        fprintf(out,"%15.7g ", mx[dim*k+l]);
      }
      fprintf(out,"\n");
    }
    if (!pure) fprintf(out,"---\n");
    break;
  case 'v':
    if (!pure) fprintf(out,"%s:\n", name);
    for (k = 0; k < dim; k++) fprintf(out,"%15.7g ", mx[k]);
    if (!pure) fprintf(out,"\n---\n");
    else fprintf(out, "\n");
    break;
  default:
    fprintf(out,"ERROR MxPrint(): no handler 4 type %c\n", T);
  }
}

/*==*/
/* print matrix stored in ccmath-format */
void
MxPrint(double *mx, char *name, char T)
{
  MxFPrint(mx, name, T, stderr, 0);
}

void    PrintDummy(double *line)
{
  print_settings();
  PrintProb(line, 1, opt.t0);
  PrintProb(line, 1, opt.t8);
}

/*==*/
static void
print_settings(void)
{
  printf(
    "# Date: %s"
    "# Sequence: %s\n"
    "# Method: %c  Start Time: %.2f  Stop Time: %.2f Temperature: %.2f\n",
    time_stamp(),
    opt.sequence,
    opt.method,
    opt.t0,
    opt.t8,
    opt.T
  );
  if(opt.basename != NULL) printf("# basename: %s\n",opt.basename);
  else printf("# basename: <stdin>\n");
  if (opt.tinc) printf("# time increment: %.2f\n", opt.tinc);
  else printf("# time increment: %.2f \n", opt.tinc);
  if(opt.want_degenerate == 1)printf("# degeneracy:  on\n");
  else printf("# degeneracy:  off\n");
  if (opt.absrb < 1) printf("# absorbing lmin: none\n");
  else printf("# absorbing lmin: %d\n", opt.absrb);
  if (opt.n > 0) printf("# nlmins: %d\n", opt.n);
  else printf("# nlmins: till EOF\n");
}

/*==*/
static char*
time_stamp(void)
{
  time_t  cal_time;
  cal_time = time(NULL);
  return ( ctime(&cal_time) );
}

/*==*/
void
MxMemoryCleanUp (void)
{
  if(_sqrPI)       free(_sqrPI);
  if(sqrPI_)       free(sqrPI_);
  if(opt.sequence) free(opt.sequence);
  if(opt.basename) free(opt.basename);
  if(opt.fpt_file) free(opt.fpt_file);

  free_gengetopt();
}

/*==*/
static void
MxDoDegeneracyStuff(void)
{
  int i, j, b, nr, current, numsad = 1, count = 0;

  numsad = ParseSaddleFile(&saddle);
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
            D[dim * (saddle[count].list[current]-1) + (saddle[count].list[b]-1)]++;
            fprintf(stderr, "transition betweed %3d - %3d: adding  1 deg_saddle\n",
                    saddle[count].list[current], saddle[count].list[b] );
            /* -1 because matrix starts with0, NOT with 1 */
          }
        }
        current++;
      }
    }
  }

  for(i = 0; i < dim; i++)  /* make matrix symmetric */
    for(j = 0; j < dim; j++)
      if (i != j)
        D[dim*j+i] = D[dim*i+j];
  if (opt.want_verbose) {
    sprintf (Aname, "%s", "D (degeneracies)");
    MxPrint (D, Aname, 'm');
  }

  free(saddle);
}

/*==*/
static void
norm2(double *mx)
{
  /* normalize columns of matrix mx */
  /* (to euclidean norm 1)*/
  int i,j;
  long double sumsq;

  for (j=0; j<dim; j++) {
    sumsq=0.0;
    for (i=0; i<dim; i++)
      sumsq += SQ(mx[dim*i+j]);
    if(sumsq > 0.0)
      sumsq=1./sqrtl(sumsq);
    for (i=0; i<dim; i++)
      mx[dim*i+j] *= sumsq;
  }
  return;
}

#define abs(x) ((x)>0.0?(x):(-x))

/*==*/
static void
MxFixevecsAbsorb(double *evecs, double *evals)
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
  double maxent,csum;

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
        if (abs(evecs[dim*j+i]) > maxent) {
          maxent=abs(evecs[dim*j+i]);
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
    double mu=0.;
    for(i=0; i<dim; i++)
      mu += evecs[dim*i+j];
    for (i=0; i<abscount; i++)
      evecs[dim*i+j] -= mu/(double)abscount;
  }

  if (opt.want_verbose) {
    MxPrint (evals, "evals_complete", 'v');
    MxPrint (evecs, "evecs_complete", 'm');
    fflush(stdout);
    fflush(stderr);
  }

  norm2(evecs);

  if (opt.want_verbose) {
    MxPrint (evals, "evals_complete", 'v');
    MxPrint (evecs, "evecs_complete", 'm');
    fflush(stdout);
    fflush(stderr);
  }

  if (opt.want_verbose) {
    fprintf(stderr,"colsums: ");
    for (i=0; i<dim; i++) {
      csum=0.0;
      for (j=0; j<dim; j++)
        csum += evecs[dim*j+i];
      fprintf(stderr,"%g ", csum);
    }
    fprintf(stderr,"\n");
  }

  if (opt.absrb && (abscount > 1)) {
    int i,j;
    if (!opt.quiet) fprintf(stderr, "\nWARNING: found %d additional absorbing state(s): ", abscount-1);
    for(i=1; i<dim; i++) {
      if(evals[i] == 1) {
        for(j=0; j<dim; j++) {
          if(evecs[dim*j+i]==1)
            fprintf(stderr, " %5d", j+1);
        }
      }
    }
    fprintf(stderr, "\n");
  }
  fflush(stderr);
}

static void
MxFixevecs(double *evecs, double *evals)
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
  double *V, sum, sum0;

  V = (double *) MxNew (dim*dim*sizeof(double));
  mmul (V, sqrPI_, evecs, dim);

  if (opt.want_verbose)
    MxPrint(V, "Eigenvectors of rate matrix M, before Fixevecs", 'm');
  
  fprintf(stderr, "Sums of EVs of rate matrix R before fixing\n");
  for (i=0; i<dim; i++) {
    double sum;
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
    //  if (!opt.quiet)
    //	fprintf(stderr, "correcting negative equilib probability for state p[%d]=%g\n",j,V[dim*j+i]);
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
    double sum;
    sum=0.;
    for (j=0; j<dim; j++) {
      sum += V[dim*j+i];
    }
    fprintf(stderr, "%15.8g   ", i, sum);
  }
  fprintf(stderr, "\n");   

  mmul(evecs, _sqrPI, V , dim);

  norm2(evecs);
  
  if (opt.want_verbose)
    MxPrint(evecs, "Eigenvectors of symmetric matrix U, after Fixevecs", 'm');

  free(V); 
}
  

/*==*/
/* sort evecs,evals */
static void
MxSortEig(double *evals, double *evecs)
{
  int i,j,k;
  double p;

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
static void
MxBinWrite(double *Mx, char what[], char T)
{
  int i,j,len=-1;
  FILE *BINOUT=NULL;
  char *wosis=NULL, *binfile=NULL, *suffix = "bin";
  //size_t info;

  wosis=what;
  if (opt.basename == NULL)
    len=strlen(suffix)+strlen(wosis)+4;
  else
    len=strlen(opt.basename)+strlen(wosis)+strlen(suffix)+4;
  binfile = (char *)calloc(len, sizeof(char));
  assert(binfile != NULL);
  if(opt.basename != NULL) {
    strcpy(binfile, opt.basename);
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
  if (!opt.quiet) fprintf(stderr, "MxBinWrite: writing %s to %s\n", wosis, binfile);
  BINOUT = fopen(binfile, "w");
  if (!BINOUT) {
    fprintf(stderr, "ERROR: could not open file pointer for binary outfile %s\n", binfile);
    exit(EXIT_FAILURE);
  }
  /* first write dim to file */
  fwrite(&dim,sizeof(int),1,BINOUT);
  switch(T) {
  case 'm':  /* write matrix entries */
    for(i=0; i<dim; i++)
      for(j=0; j<dim; j++)
        fwrite(&Mx[dim*i+j],sizeof(double),1,BINOUT);
        //if (!opt.quiet) fprintf(stderr, "\n");
    break;
  case 'v': /* write vector entries */
    for(i=0; i<dim; i++)
      fwrite(&Mx[i],sizeof(double),1,BINOUT);
      //fprintf(stderr, "\n");
    break;
  default:
    fprintf(stderr, "ERROR MxBinWrite(): no handler for type %c\n",T);
    exit(EXIT_FAILURE);
  }
  fclose(BINOUT);
  //  if(binfile) free(binfile);
}

/*==*/
static int
MxBinRead(double **Mx, char what[], char T)
{
  int dimension=0,len=-1;
  FILE *BININ=NULL;
  char *wosis=NULL, *binfile=NULL, *suffix="bin";
  double *data=NULL;
	int ref;
  //size_t info;

  wosis=what;
  if (opt.basename == NULL)
    len=strlen(suffix)+strlen(wosis)+4;
  else
    len=strlen(opt.basename)+strlen(wosis)+strlen(suffix)+4;
  binfile = (char *)calloc(len, sizeof(char));
  assert(binfile != NULL);
  if(opt.basename != NULL) {
    strcpy(binfile, opt.basename);
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
  switch(T) { /* read data */
  case 'm':
    data = (double *)calloc(dimension*dimension, sizeof(double));
    ref = fread((void*)data, sizeof(double), dimension*dimension,BININ);
    break;
  case 'v':
    data = (double *)calloc(dimension, sizeof(double));
    ref = fread(data, sizeof(double), dimension,BININ);
    break;
  default:
    fprintf(stderr, "ERROR MxBinRead(): no handler for type %c\n",T);
    exit(EXIT_FAILURE);
  }

  *Mx = data;
  fclose(BININ);
  free(binfile);
  return dimension;
}

/*==*/
static void
MxASCIIWrite(double *Mx, char *asciifile)
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
      fprintf(ASCIIOUT,"%15.10g ", Mx[dim*j+i]);
    }
    fprintf(ASCIIOUT,"\n");
  }
  if (!opt.quiet) fprintf(stderr, "matrix written to ASCII file\n");
  fclose(ASCIIOUT);
}

static void
MxASCIIWriteV(double *Mx, char *asciifile)
{
  int i;
  FILE *ASCIIOUT;

  ASCIIOUT = fopen(asciifile, "w");
  if (!ASCIIOUT) {
    fprintf(stderr, "could not open file pointer 4 ASCII outfile\n");
    exit(EXIT_FAILURE);
  }
  for(i=0; i<dim; i++) {
      fprintf(ASCIIOUT,"%15.10g ", Mx[i]);
  }
  if (!opt.quiet) fprintf(stderr, "vector written to ASCII file\n");
  fclose(ASCIIOUT);
}

/*==*/
void
MxExponent(double *p0, double *p8, double *U)
{
  int i,j, pdiff_counter, count = 0;
  double x, tt, time, *Uexp, *Umerk, *pt, *pdiff, check = 0.;

  Umerk  = (double *) MxNew (dim*dim*sizeof(double));
  Uexp   = (double *) MxNew (dim*dim*sizeof(double));
  pt     = (double *) MxNew (dim*sizeof(double));
  pdiff  = (double *) MxNew (dim*sizeof(double));

  memcpy(Umerk, U, dim*dim*sizeof(double));

  /* solve fundamental equation */
  if (opt.t0 == 0.0) {
    if (opt.method=='F')  PrintProbFull(p0, dim, 0.0, lmins);
    else                  PrintProb(p0, dim, 0.0);
    opt.t0 = TZERO;
  }

  for (i=0; i<dim; i++) U[(dim+1)*i] -= 1;
  print_settings();
  for (tt = opt.t0; tt < opt.t8*opt.tinc; tt *= opt.tinc) {
    time = (tt<opt.t8)? tt:opt.t8;
    memcpy(U, Umerk, dim*dim*sizeof(double));
    for (i=0; i<dim*dim; i++) U[i]*=time;
    padexp(U,Uexp,dim,30);
    x = 0.;
    for(j=0; j<dim*dim; j++) x+=Uexp[j];
    for(j=0; j<dim*dim; j++) Uexp[j]*=(double)dim/x;
    vmul(pt, Uexp, p0, dim);
    /* check convergence */
    for(i=0; i<dim; i++) {
      pdiff[i] = p8[i] - pt[i];
      if (fabs(pdiff[i]) >= 0.0001)
        pdiff_counter++;
    }
    if (pdiff_counter < 1)
      break;
    pdiff_counter = 0.;
    /* end check convergence */
    check = 0.;
    printf("%e ", time);  /* print p(t) to stdout */
    for (i = 0; i < dim; i++) {
      if(pt[i] < -0.00001) {
        fprintf(stderr, "prob of lmin %i has become negative!\n", i+1);
        exit(EXIT_FAILURE);
      }
      printf("%e ", fabs(pt[i]));
      check += fabs(pt[i]);
    }
    printf("\n");

    count++;  /* # of iterations */

    if ( ((check-1) < -0.01) || ((check-1) > 0.01) ) {
      fprintf(stderr, "overall probability at time %e is %e != 1.0 %s!\n", time, check, (opt.num_err == 'R'?"rescaling":"exiting") );
      if (opt.num_err == 'H') exit(EXIT_FAILURE);
    }
    memset(pt,   0, dim*sizeof(double));
    memset(pdiff, 0, dim*sizeof(double));
    memset(Uexp, 0, dim*dim*sizeof(double));
    memset(U, 0, dim*dim*sizeof(double));

  }
  printf("# of iterations: %d\n", count);
  free(Uexp);
  free(pt);
}

/*==*/
void
MxFPT(double *U, double *p8, FILE *out)
{
  int i,j;

  //fprintf(stderr, "in MxFPT\n");
  //if(opt.want_verbose) MxPrint (U,"U" , 'm');
  double *Z=NULL, *M=NULL;

  //MxPrint(U, "U mxfpt", 'm');

  if (opt.absrb) {
    Z = (double *) MxNew((dim-1)*(dim-1)*sizeof(double));

    extern void dgesv_(int *N,int *NRHS,double *A,int *LDA,int *IPIV,
                     double *B,int *LDB,int *INFO);

    int i,j,nrhs,nfo,*ipiv=NULL;
    double *B=NULL;
    B    = (double *)malloc((dim-1)*sizeof(double));

    for (i=0; i<dim-1; i++) B[i] = 1.0;

    for(i = 0; i < dim-1; i++)
        for(j = 0; j < dim-1; j++) {
          Z[(dim-1)*i+j] = (i==j && opt.useplusI?1.0:0.0) - U[dim*i+j]; /* I-U = -Q  (U is transposed) without absorbing state*/
        }

    ipiv = (int *)malloc((dim-1)*sizeof(int));
    nrhs=1;

    int n=dim-1;
    //DGESV computes the solution to a real system of linear equations A * X = B (I think A must be transposed)
    dgesv_(&n, &nrhs, Z, &n, ipiv, B, &n, &nfo);

    // fix non-zero in real absorbing state
    for(i=0; i<dim-1; i++) {
      B[i] -= B[opt.real_abs];
    }
    dim--;
    MxFPrint(B, "Absorbing times: ", 'v', out, out!=stderr);
    dim++;

    free(B);
    free(Z);
    free(ipiv);

  } else { // non-absorbing case
    Z = (double *) MxNew(dim*dim*sizeof(double));

    for(i = 0; i < dim; i++)
      for(j = 0; j < dim; j++) {
        Z[dim*i+j] = (i==j&& opt.useplusI?1.0:0.0) - U[dim*j+i] + p8[j]; /* I-U+W = W-Q  (U is transposed)*/
      }

    //if(opt.want_verbose) MxPrint (Z,"I-U+W" , 'm');
    minv(Z,dim);
    //if(opt.want_verbose)MxPrint (Z,"Fundamental matrix Z=inv(I-U+W)" , 'm');

    M = (double *) MxNew(dim*dim*sizeof(double));
    for(i = 0; i < dim; i++) {
      for(j = 0; j < dim; j++) {
        M[dim*i+j] = (Z[dim*j+j]-Z[dim*i+j])/p8[j];
      }
    }

    MxFPrint(M, "First passage times (i-th row, j-th column represents fpt from i-th to j-th state)", 'm', out, out!=stderr);
    free(M);
    free(Z);
  }
}

/*==*/
static void
MxKotzOutMathematica(double *Mx)
{
  int i,j;
  FILE *MATHEMATICA_OUT=NULL;
  char *mathematica_file = "mxMat.txt";

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
        fprintf(MATHEMATICA_OUT, "%25.22f, ", Mx[dim*j+i]);
      else
        fprintf(MATHEMATICA_OUT, "%25.22f}", Mx[dim*j+i]);
    }
    if (i != (dim-1))
      fprintf(MATHEMATICA_OUT, ",\n");
    else
      fprintf(MATHEMATICA_OUT, "}\n");
  }
  fclose(MATHEMATICA_OUT);
}

double *MxFPTOneState(double *U, int state)
{
  extern void dgesv_(int *N,int *NRHS,double *A,int *LDA,int *IPIV,
                     double *B,int *LDB,int *INFO);

  if (state>=dim) {
    fprintf(stderr, "State %d does not exist (--fpt %d)", state+1, state+1);
    return NULL;
  }

  int i,j,nrhs,nfo;
  int *ipiv=NULL;
  int n=dim-1;
  double *A=NULL, *B=NULL;
  A    = (double *)malloc(n*n*sizeof(double));
  B    = (double *)malloc(dim*sizeof(double));

  // A = Q(infinetisimal generator) = U^T-I with state-th column and row deleted (U is transposed!)
  for(i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      int ui = (i>=state?i+1:i);
      int uj = (j>=state?j+1:j);
      A[n*i+j] = U[dim*ui+uj] - (i==j && opt.useplusI ?1.0:0.0);
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
  dgesv_(&n, &nrhs, A, &n, ipiv, B, &n, &nfo);


  for (i=dim-1; i>state; i--) {
    B[i] = B[i-1];
  }
  B[state] = 0.0;

  free(A);
  free(ipiv);
  return B;
}


void MxFPTSimple(double *U)
{
  int i,j;
  fprintf(stderr, "in MxFTPSimple\n");
  double *B=NULL, *M=NULL;
  M = (double *) MxNew(dim*dim*sizeof(double));

  int p_done = 0;

  for (i=0; i<dim; i++) {
    // calculate to one state and add it to global matrix
    B = MxFPTOneState(U, i);
    for (j=0; j<dim-1; j++) {
      M[dim*j+i] = B[j];
    }
    free(B);
    // print progress
    if (i*100/dim >= p_done) {
      p_done = i*100/dim;
      fprintf(stderr, "%d%%\n", p_done);
    }
  }

  MxPrint(M, "FPT", 'm');

  free(M);
}

void MxFPTRnd(double *U, int packets)
{
  srand(time(NULL));
  double *P=NULL, *M=NULL;

  P = (double *)malloc(dim*dim*sizeof(double));

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

  if (!opt.absrb) {
    // ergodic option
    M = (double *)malloc(dim*dim*sizeof(double));
    for (i=0; i<dim; i++) {
      for (j=0; j<dim; j++) {
        M[dim*i+j] = MxFPTRandom(P,U, i, j, packets);
        fprintf(stderr, "x");
      }
      fprintf(stderr, "\n");
    }

    MxPrint(M, "M (random)", 'm');
  } else {
    // absorbing option
    M = (double *)malloc(dim*sizeof(double));
    for (i=0; i<dim; i++) {
      M[i] = MxFPTRandom(P,U, i, dim-1, packets);
    }

    MxPrint(M, "M (random)", 'v');
  }

  free(M);
  free(P);
}

double MxFPTRandom(double *P, double *U, int src, int dst, int packets)
{
  if (src==dst) return 0.0;

  double wtime = 0.0;
  int i,j;
  for (i=0; i<packets; i++) {
      int curr = src;
      double time = 0.0;
      while (curr != dst) {
        time += 1/(1.0-U[dim*curr+curr]);
        double next = random()/(double)RAND_MAX;
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
  return wtime/(double)packets;
}

static void
MxEVLapackSym(double *U) {
  /*   extern void dsyev_(char *jobz, char *uplo, int *n, double *A,int *lda, */
  /*         double *w, double *work, int *lwork, int *info); */
  extern void dsyevx_(char *jobz,char *range,char *uplo,int *n,double *A,
                      int *lda,double *vl,double *vu,int *il,int *iu,
                      double *abstol,int *m,double *w,double *z,int *ldz,
                      double *work,int *lwork,int *iwork,int *ifail,int *info);

  extern void dsyevr_(char *jobz, char *range, char *uplo, int *n, double* a, int *lda, double *vl,
		     double *vu, int *il, int *iu,
		     double *abstol, int* m, double* w, double* z,
		     int *ldz, int* isuppz, double *work,int *lwork,int *iwork,int *liwork,int *info );
  
  double abstol;
  abstol = opt.FEPS;

  if (opt.want_verbose) fprintf(stderr, "FEPS precision %20.10g\n", abstol);


  double *work=NULL, vl, vu;   // unused or inputs
  int il, iu, m, lwork, liwork, *iwork=NULL, *ifail=NULL, *isuppz=NULL, nfo;  // unused or inputs
  lwork = dim*(dim+26);
  liwork = 10*dim;
  work   = (double *) malloc (lwork * sizeof(double));
  iwork  =    (int *) malloc (liwork * sizeof(int));
  ifail  =    (int *) malloc (dim * sizeof(int));
  isuppz =    (int *) malloc (2*dim * sizeof(int));
  

  //if (!opt.useplusI) for (il=0; il<dim; il++) U[il*dim+il]-=1.0;

  // works only with: dim, U, evals, evecs
#if 1
  dsyevx_("V", "A", "U",&dim, U, &dim, &vl, &vu, &il, &iu,
          &abstol, &m, evals, evecs, &dim, work, &lwork, iwork,
          ifail, &nfo);
#else
  dsyevr_("V", "A", "U",&dim, U, &dim, &vl, &vu, &il, &iu,
          &abstol, &m, evals, evecs, &dim, isuppz, work, &lwork, iwork,
          &liwork, &nfo);
#endif
  
  //dsyev_("V","U",&dim, S, &dim, evals, work, &lwork, &nfo);

  //MxPrint(evecs, "Eigenvectors (bad)", 'm');
  //MxPrint(evals, "Eigenvalues (bad)", 'v');

  if(nfo != 0) {
    fprintf(stderr, "dsyevx exited with value %d (val=%20.10g) (Cannot compute eigenvalues) Try to:\n - lower --feps value (try --feps=-1.0)\n - try using --useplusI option (sometimes the eigenvalues are better computed like that)\n", nfo, evals[nfo-1]);
    for (il=0; il<dim; il++) fprintf(stderr, "%20.10g ", evecs[nfo*dim+il]);
    fprintf(stderr, "\n");

    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
  }

  //if (!opt.useplusI) for (il=0; il<dim; il++) U[il*dim+il]+=1.0;


  free(work);
  free(iwork);
  free(ifail);
  free (isuppz);

  // transpose output
  trnm(evecs, dim);
}

static void
MxEVLapackNonSym(double *origU)
/* input: matrix origU, space for (right)evec-matrix S, */
/* since fortran sucks, we transpose the input matrix   */
/* and still want the right eigenvectors                */
/* matrix of right eigenvecs is transposed  to have */
/* the j-th evec in column c */
{
  extern void dgeevx_(char *B, char *JVL, char *JVR, char *SENS, int *N, double *A,
                      int *LDA,double *WR, double * WI, double *VL, int *LDVL,
                      double *VR, int *LDVR,int *ILO, int *IHI, double *scale,
                      double *ABNRM, double *RCONDE, double *RCONDV, double *WORK,
                      int *LWORK, int *IWORK, int *INFO);

  double *evals_re=NULL, *evals_im=NULL, *scale=NULL;
  double *rconde=NULL, *rcondv=NULL, *work=NULL, abnrm;
  int one, ilo, ihi, lwork, *iwork, nfo;
  //int dimx2;
  /* for sorting */
  int i;

  dim = dim+500;
  one = 1;
  //dimx2 = 2*dim;
  lwork = dim*(dim+6);

  evals_re = (double *) malloc (dim * sizeof(double));
  evals_im = (double *) malloc (dim * sizeof(double));
  scale    = (double *) malloc (dim * sizeof(double));
  rconde   = (double *) malloc (dim * sizeof(double));
  rcondv   = (double *) malloc (dim * sizeof(double));
  work     = (double *) malloc (lwork * sizeof(double));
  iwork    =    (int *) malloc (2*(dim -2) * sizeof(int));
  dim=dim-500;
  if ( (evals_re && evals_im && scale && rconde && \
        rcondv && work && iwork)==0) {
    fprintf(stderr,"no space for temporary lapack arrays!\n");
    exit(EXIT_FAILURE);
  }

  /* instead of more fiddling, we transpose the input */
  trnm(origU, dim);
  /*for (i=0; i<dim; i++)
    for (j=i+1; j<dim; j++) {
      tmp = origU[dim*i+j];
      origU[dim*i+j]=origU[dim*j+i];
      origU[dim*j+i]=tmp;
    }*/

  dgeevx_("B", "N", "V", "V", &dim, origU, &dim, evals_re, evals_im, NULL, &one,\
        evecs ,&dim, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work,\
          &lwork, iwork, &nfo);

  for (i=0; i<dim; i++) evals[i]=evals_re[i];
  for (i=0; i<dim; i++)
    if ((evals_re[i] != 0.0) && fabs(evals_im[i]/evals_re[i])>opt.FEPS)
      fprintf(stderr,"eigenvalue %d is %g + i*%g, which is somewhat complex\n",
              i,evals_re[i],evals_im[i]);

  /*transpose output*/
  trnm(evecs, dim);
  /*for (i=0; i<dim; i++)
    for (j=i+1; j<dim; j++) {
      tmp = evecs[dim*i+j];
      evecs[dim*i+j]=evecs[dim*j+i];
      evecs[dim*j+i]=tmp;
    }*/

  free(evals_re);
  free(evals_im);
  free(scale);
  free(rconde);
  free(rcondv);
  free(work);
  free(iwork);
}

void MxEqDistrFromLocalBalance ( double *U, double **p8 )
{
    double *res = *p8;
    int count = 1;  // num of solved states
    int i=0;
    int j=1;
    int k;

    for(k=1; k<dim; k++) {
      res[k] = 0.0;
    }
    res[0] = 1.0;

    // while there are unsolved states
    while (count != dim) {
      if (U[dim*i+j]!=0.0 && res[i]!=0.0 && res[j]==0.0) {
        if (U[dim*j+i]==0.0) {
          fprintf(stderr, "Local balance is unsatisfiable at U[%d][%d]=%f\n", i, j, U[dim*i+j]);
          return ;
        }
        // local balance
        long double tmp = U[dim*j+i];
        tmp *= res[i];
        tmp /= U[dim*i+j];
        res[j] = tmp;
        count ++;
      }
      j++;
      if (j==dim) {
        j=0;
        i++;
        if (i==dim && count != dim) {
          fprintf(stderr, "non-ergodic chain\n");
          return ;
        }
      }
    }

    /* now make the vector stochastic (sum = 1.0) */
    long double sum=0.0;
    for (i=0; i<dim; i++) {
      sum += res[i];
    }
    for (i=0; i<dim; i++) {
      res[i] /= sum;
    }

    MxPrint(res, "p8 (local balance)", 'v');

    /* make it unit vector */
    /* long double qsum=0.0; */
    /* for (i=0; i<dim; i++) { */
    /*   qsum += res[i]*res[i]; */
    /* } */
    /* qsum = sqrtl(qsum); */
    /* for (i=0; i<dim; i++) { */
    /*   res[i] /= qsum; */
    /* } */
}

int MxShorten(double **shorten, int nstates, int my_dim, int max) {
  // shorten the matrix from dimension my_dim to nstates:
  if (my_dim>nstates) {
    if (!opt.quiet) fprintf(stderr, "decreasing %d to %d\n", my_dim, nstates);

    // first we need to fix the diagonal entries tmp_rates[i][i] = sum_j tmp_rates[i][j]
    double *tmp_rates = *shorten;
    int i, j;
    for (i = 0; i < my_dim; i++) tmp_rates[my_dim*i+i] = 0.0;
    for (i = 0; i < my_dim; i++) {
      double tmp = 0.00;
      // calculate column sum
      for(j = 0; j < my_dim; j++)  tmp += tmp_rates[my_dim*j+i];
      tmp_rates[my_dim*i+i] = -tmp;
    }

    // shorten by some value max
    if (max!=1) {
      while (my_dim-max > nstates) {
        MxRShorten(shorten, my_dim, my_dim-max);
        if (!opt.quiet) fprintf(stderr, "%d done...\n", my_dim-max);
        my_dim -= max;
      }
      MxRShorten(shorten, my_dim, nstates);
      my_dim = nstates;
    } else {
      while (my_dim!=nstates) {
        MxOneShorten(shorten, my_dim);
        my_dim--;
        if (!opt.quiet && my_dim%100==0 && my_dim>0) fprintf(stderr, "%d done...\n", my_dim);
      }
    }
  }

  return my_dim;
}

void MxRShorten(double **shorten, int fulldim, int gdim)
{
  //does: shortened = GG - GB*BB^(-1)*BG, where matrix tmp_rates is split as:
  //tmp_rates = (GG | GB)
  //            (BG | BB)
  // GG has dimension gdim*gdim; tmp_rates fulldim*fulldim
  // create matrices:

  int bdim = fulldim - gdim;
  int i,j;
  double *tmp_rates = *shorten;

  double *gg = (double *)calloc(gdim*gdim,sizeof(double));
  double *bg = (double *)calloc(bdim*gdim,sizeof(double));
  double *bb = (double *)calloc(bdim*bdim,sizeof(double));
  double *gb = (double *)calloc(gdim*bdim,sizeof(double));

  // fill the matrices: (row = i; column = j)
  for (i=0; i<gdim; i++) {
    for (j=0; j<gdim; j++) {
      gg[gdim*i+j] = tmp_rates[fulldim*i+j];
    }
  }

  for (i=0; i<bdim; i++) {
    for (j=0; j<gdim; j++) {
      bg[gdim*i+j] = tmp_rates[fulldim*(i+gdim)+j];
    }
  }

  for (i=0; i<gdim; i++) {
    for (j=0; j<bdim; j++) {
      gb[bdim*i+j] = tmp_rates[fulldim*i+j+gdim];
    }
  }

  for (i=0; i<bdim; i++) {
    for (j=0; j<bdim; j++) {
      bb[bdim*i+j] = tmp_rates[fulldim*(i+gdim)+j+gdim];
    }
  }

  /*MxFPrintD(tmp_rates, "Q", my_dim, my_dim, stderr);
  MxFPrintD(gg, "GG", dim, dim, stderr);
  MxFPrintD(bg, "BG", bdim, dim, stderr);
  MxFPrintD(gb, "GB", dim, bdim, stderr);
  MxFPrintD(bb, "BB", bdim, bdim, stderr);
*/
  // result2 = gb*bb^(-1)*bg
  minv(bb, bdim);
  //MxFPrintD(bb, "BBinv", bdim, bdim, stderr);
  double *result = (double *)calloc(gdim*bdim,sizeof(double));
  mmul_singular(result, gb, bb, gdim, bdim, bdim, 0);
  //MxFPrintD(result, "gb*bb-1", dim, bdim, stderr);
  double *result2 = (double *)calloc(gdim*gdim,sizeof(double));
  mmul_singular(result2, result, bg, gdim, bdim, gdim, 1);

  if (opt.want_verbose) MxFPrintD(result2, "gb*bb-1*bg", gdim, gdim, stderr);


  // result2 = gg - result2
  for (i=0; i<gdim; i++) {
    for (j=0; j<gdim; j++) {
      result2[gdim*i+j] = gg[gdim*i+j] - result2[gdim*i+j];
    }
  }

  //MxFPrintD(result2, "matrix after shortening", dim ,dim, stderr);
  free(result);
  free(*shorten);
  free(gg);
  free(gb);
  free(bg);
  free(bb);
  *shorten = result2;
}

void MxOneShorten(double **shorten, int fulldim)
{
  //does: shortened = GG - GB*BB^(-1)*BG, where matrix tmp_rates is split as:
  //tmp_rates = (GG | GB)
  //            (BG | BB)
  // GG has dimension fulldim-1*fulldim-1; tmp_rates fulldim*fulldim
  // create matrices:

  int gdim = fulldim-1;
  double *result = (double *)calloc(gdim*gdim,sizeof(double));

  double *tmp_rates = *shorten;

  int i, j;
  double c = 1.0/tmp_rates[fulldim*gdim + gdim];
  for (i=0; i<gdim; i++) {
    for (j=0; j<gdim; j++) {
      // just x - a*c^-1*b
      result[gdim*i+j] = tmp_rates[fulldim*i+j] - c*tmp_rates[fulldim*gdim + j]*tmp_rates[fulldim*i + gdim];
    }
  }

  //MxFPrintD(tmp_rates, "Q", fulldim, fulldim, stderr);
  //MxFPrintD(result, "Q-1", gdim, gdim, stderr);

  free(*shorten);
  *shorten = result;
}

int MxReadBinRates(FILE *rate_file, double **rate_mx, int nstates, int max)
{
  int dimension = 0;
	int ref;
  /* read dimension from file */
  ref = fread(&dimension,sizeof(int),1,rate_file);
  double *data = (double *)calloc(dimension*dimension, sizeof(double));
  ref = fread((void*)data, sizeof(double), dimension*dimension,rate_file);

  *rate_mx = data;
  fclose(rate_file);

  // decrease dimension
  dimension = MxShorten(rate_mx, nstates, dimension, max);

  return dimension;
}


/* End of file */

