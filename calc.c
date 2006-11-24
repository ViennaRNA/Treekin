/*=================================================================*/
/*=   calc.c                                                      =*/
/*=   main calculation and iteration routines for treekin         =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-11-24 15:43:38 mtw>          =*/
/*=   $Id: calc.c,v 1.37 2006/11/24 15:00:39 mtw Exp $            =*/
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
#include "exp_matrix.h" /* functions for matrix-exponent stuff */
#include "mxccm.h"      /* functions for eigen-problems stolen from ccmath */
#include "barparser.h"  /* functions for input */
#include "calc.h"       /* does all matrix stuff for markov process */
#include "globals.h"    /* contains getopt-stuff */

#define SQ(X) ((X)*(X))
#define FEPS 1.e-15

/* private function(s) */
static void   *MxNew (size_t size);
static double *MxMethodeA (BarData *Data);
static double *MxMethodeFULL(double *);
static double *MxMethodeINPUT (BarData *Data, double *);
static double  max_saddle(int i, int j, BarData *Data);
static void    print_settings(void);
static char   *time_stamp(void);
static void    MxDoDegeneracyStuff(void);
static void    MxBinWrite (double *matrix);
static void    MxASCIIWrite(double *matrix);
static void    MxKotzOutMathematica(double *matrix);
static void    MxSortEig(double *evals, double *evecs);
static void    MxEVLapackSym(double *U);
static void    MxEVLapackNonSym(double *U);
static void    MxFixevecs(double *, double *);

/* private vars and arrays */
static int      dim = 0;
static double   _kT = 1.;
static double  *evals     = NULL;
static double  *evecs     = NULL;
static double  *_sqrPI    = NULL;   /* left constant array */
static double  *sqrPI_    = NULL;   /* right constant array */
static double  *D         = NULL;   /* matrix with degree of degeneacy */
static char     Aname[30];
static TypeDegSaddle *saddle = NULL;

/*==*/
void
MxInit (int d)
{
  _kT = 0.00198717*(273.15 + opt.T);
  if (d > 0 ) dim = d;
  else { fprintf(stderr, "dim <= 0\n"); exit(EXIT_FAILURE); }
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
  if (opt.dumpU)
   /*    MxBinWrite(U); */
    MxASCIIWrite(U);
  return (U);
}
/*==*/
void
MxGetSpace (double **p8)
{
  *p8   = (double *) MxNew (dim*sizeof(double));
  evals = (double *) MxNew (dim*sizeof(double));
  evecs = (double *) MxNew (dim*dim*sizeof(double));
  assert(evals!=NULL); assert(evecs!=NULL); assert(p8!=NULL);
  
  if(!opt.absrb){
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
  for (i = 1; i < (int) *opt.pini; i+=2)
    pzero[(int)opt.pini[i]-1] = (double)opt.pini[i+1];
  /* -1 because our lmins start with 1, not with 0 (as Data does ) */
  
  if (opt.want_verbose) MxPrint (pzero, "p0", 'v');
  *p0=pzero;
}

/*==*/
/* calculate equilibrium distribution */
void
MxEqDistr ( BarData *Data, double *p8 )
{
  int i;
  double Z = 0.;
  
  if(opt.absrb){
    for(i = 0; i < dim; i++)
      p8[i] = 0.;
    p8[dim-1] = 1.0; /* last entry is the 'new' absorbing state */
  }
  else{
    for(i = 0; i < dim; i++) Z += exp(-((double)Data[i].FGr/_kT));
    for(i = 0; i < dim; i++) p8[i] = exp(-((double) Data[i].FGr/_kT))/Z;
  }
  if(opt.want_verbose) MxPrint (p8, "p8", 'v');
  return;
}

/*==*/
void
MxEqDistrFULL (SubInfo *E, double *p8 )
{
  int i;
  double Z = 0.;

  if(opt.absrb){
    for(i = 0; i < dim; i++)
      p8[i] = 0.;
    p8[opt.absrb-1] = 1.;
  }
  else{
    for(i = 0; i < dim; i++) Z += exp(-E[i].energy/_kT);
    for(i = 0; i < dim; i++) p8[i] = exp(-E[i].energy/_kT)/Z;
  }
  if(opt.want_verbose) MxPrint (p8, "p8", 'v');
  return;
}

/*==*/
void
MxEqDistrFromLinSys ( double *U, double **p8 )
{
  extern void dgesv_(int *N,int *NRHS,double *A,int *LDA,int *IPIV,
		     double *B,int *LDB,int *INFO);
  int i,j,n, nrhs, nfo, *ipiv=NULL;
  double *A=NULL, *B=NULL;
  long double sumsq;
  
  n=dim-1;
  A    = (double *)malloc(n*n*sizeof(double));
  B    = (double *)malloc(n*sizeof(double));
  ipiv = (int *)malloc(n*sizeof(int));
  nrhs=1;
  
  for(i=1;i<=n;i++)   /* all except first row */
    for(j=0;j<n;j++)
      A[n*(i-1)+j]=U[dim*i+j];
  for(n=0,i=1;i<dim;i++,n++) 
    B[n]=U[dim*i+(dim-1)];
  dim=n;
  if(opt.want_verbose){
    MxPrint(A, "A in MxEqDistrFromLinSys", 'm' );
    MxPrint(B, "B in MxEqDistrFromLinSys", 'v' );
  }

  dgesv_(&n, &nrhs, A, &n, ipiv, B, &n, &nfo);
  if(nfo != 0){
    fprintf(stderr, "dgesv exited with value %d\n", nfo);
    exit(EXIT_FAILURE);
  }
  if(opt.want_verbose){
    MxPrint(B, "B in MxEqDistrFromLinSys", 'v' );
  }
  dim=n+1;
  *p8[0]=1.;
  for(i=1;i<dim;i++)
    *(*p8+i)=B[i-1];
  if(opt.want_verbose){
    MxPrint(*p8, "p8 in MxEqDistrFromLinSys before norm", 'v' );
  }
  /* now normalize the new p8 */
   sumsq=0.0;
  for (i=0;i<dim;i++)
    sumsq += SQ(*(*p8+i));
  if(sumsq > 0.0)
    sumsq=1./sqrtl(sumsq);
  for (i=0;i<dim;i++)
    *(*p8+i)  *= sumsq;
  if(opt.want_verbose){
    MxPrint(*p8, "p8 in MxEqDistrFromLinSys after  norm", 'v' );
   }
  free(A);
  free(B);
  free(ipiv);
  exit(EXIT_SUCCESS);
}

/*==*/
void
MxDiagonalize ( double *U, double **_S, double *P8 )
{
  int i,j;
  double csum, *tmpMx=NULL;

  if(opt.dumpMathematica == 1)  MxKotzOutMathematica(U);
  if(!opt.absrb){
    tmpMx = (double *) MxNew (dim*dim*sizeof(double));
    
    for(i = 0; i < dim; i++) {
      for(j = 0; j < dim; j++) {
	if( i == j) {
	  sqrPI_[dim*i+j] = sqrt(P8[i]);            /* pos right */
	  _sqrPI[dim*i+j] = 1/(sqrPI_[dim*i+j]);}}} /* neg left */
    mmul(tmpMx, _sqrPI, U, dim);
    memset(U,0,dim*dim*sizeof(double));
    mmul(U, tmpMx, sqrPI_, dim);
    free(tmpMx);
    /* correct for numerical errors */
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
	U[dim*i+j] = (U[dim*i+j]+U[dim*j+i])/2;
	U[dim*j+i] = U[dim*i+j];
      }
    }
    if (opt.want_verbose) MxPrint (U, "force symmetrized U", 'm');
  }
  
  if (opt.dumpU){
    fprintf(stderr, "in MxDiagonalize: writing U to mx.bin\n");
    MxBinWrite(U);
  }

  if(opt.absrb)
    MxEVLapackNonSym(U);
  else    
    MxEVLapackSym(U);
  
  MxSortEig(evals, evecs);
  
  if (opt.want_verbose){
    MxPrint(evecs, "Eigenvectors of U (LAPACK)", 'm');
    MxPrint(evals, "Eigenvalues of  U (LAPACK)", 'v');
  }
  for (i=0;i<dim;i++)
    if (evals[i]>(1.0-FEPS)) evals[i]=1.0;

  if(opt.absrb)
    MxFixevecs(evecs,evals);
  *_S=evecs;
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
  int i,  count = 0, pdiff_counter = 0, tzero_flag=0;
  double time, check = 0.;
  double *CL = NULL, *CR, *exptL, *tmpVec, *tmpVec2, *pt, *St, *pdiff;
  double *ptFULL = NULL;  /* prob dist 4 of the effective lmins of the tree at time t */
  double *p8FULL = NULL;  /* equ dist 4 gradient basins, full process */
  
  tmpVec    = (double *) MxNew (dim*sizeof(double));
  tmpVec2   = (double *) MxNew (dim*sizeof(double));
  pt        = (double *) MxNew (dim*sizeof(double));
  exptL     = (double *) MxNew (dim*dim*sizeof(double));
  if(opt.method=='F'){
    ptFULL    = (double *) MxNew ((lmins+1)*sizeof(double));
    p8FULL    = (double *) MxNew ((lmins+1)*sizeof(double));
    pdiff     = (double *) MxNew ((lmins+1)*sizeof(double));
  }
  else pdiff   = (double *) MxNew (dim*sizeof(double));
  
  if(! opt.absrb){  /* NON-absorbing case */
    CL        = (double *) MxNew (dim*dim*sizeof(double));
    CR        = (double *) MxNew (dim*dim*sizeof(double));
    St        = (double *) MxNew (dim*dim*sizeof(double));
    mcopy(St, S, dim*dim);  trnm(St, dim); /* transpose S */
    mmul (CL, sqrPI_, S, dim);
    mmul (CR, St, _sqrPI, dim);
    vmul (tmpVec, CR, p0, dim);
    free(St);
    free(CR);
  }
  else{  /* absorbing case */
    double *S_inv;
    S_inv  = (double *) MxNew (dim*dim*sizeof(double));
    mcopy(S_inv, S, dim*dim);
    minv(S_inv,dim);
    if(opt.want_verbose) MxPrint(evals, "evals in MxIterate", 'v');
    vmul (tmpVec, S_inv, p0, dim);
    free(S_inv);
  }
  if(opt.method=='F'){  /* calculate equilibrium distribution once */
    for (i = 0; i < dim; i++)   p8FULL[E[i].ag] += p8[i];
    for (i = 0; i < lmins; i++) check += fabs(p8FULL[i]);
    if ( ((check-1) < -0.1) || ((check-1) > 0.1) ){
      fprintf(stderr, "overall equilibrium probability is %e != 1. ! exiting\n", check);
      exit(EXIT_FAILURE);
    }
  }
  check = 0.;
  for(i = 0; i < dim; i++) evals[i] -= 1;  /* compensate 4 translation of matrix U */
  
  
  /* solve fundamental equation */
  print_settings();
  if (opt.t0 == 0.0) {
    printf(" %e ", 0.0);
    for (i = 0; i < dim; i++) printf("%e ", fabs(p0[i]));
    printf("\n");
    opt.t0 = TZERO;
    tzero_flag=1;
  }
      
  for (time = opt.t0; time <= opt.t8; time *= opt.tinc) {
    for (i = 0; i < dim; i++)
      exptL[dim*i+i] = exp(time*evals[i]);
    vmul(tmpVec2, exptL, tmpVec, dim);
    if(!opt.absrb)  vmul(pt, CL, tmpVec2, dim);
    else            vmul(pt, S, tmpVec2, dim);
    
    count++;  /* # of iterations */
    
    printf(" %e ", time);  /* print p(t) to stdout */
    for (i = 0; i < dim; i++){
      if(pt[i] < -0.01){
	fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n",
		i+1, time, pt[i]);
	exit(EXIT_FAILURE);
      }
      /* map individual structure -> gradient basins */
      if(opt.method=='F') ptFULL[E[i].ag] += pt[i]; 
      else   printf("%e ", fabs(pt[i]));
      check += fabs(pt[i]); 
    }
    if(opt.method=='F') 
      for(i = 1; i <= lmins; i++) printf("%e ", fabs(ptFULL[i]));
    printf("\n");
    
    if ( ((check-1) < -0.05) || ((check-1) > 0.05) ){
      fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
      exit(EXIT_FAILURE);
    }
    check = 0.;
    /* now check if we have converged yet */
    if(opt.method=='F'){
      for(i = 1; i <= lmins; i++){
	pdiff[i] = p8FULL[i] - ptFULL[i];
	if (fabs(pdiff[i]) >= 0.001)
	  pdiff_counter++; /* # of lmins whose pdiff is > the threshold */
      }
    }
    else{
      for(i = 0; i < dim; i++){
	pdiff[i] = p8[i] - pt[i];
	if (fabs(pdiff[i]) >= 0.001)
	  pdiff_counter++;
      }
    }
    if (pdiff_counter < 1) /* all mins' pdiff lies within threshold */
      break;
    pdiff_counter = 0;
    memset(pdiff, 0, (lmins+1)*sizeof(double));
    /* end check of convergence */

    if(opt.method=='F') memset(ptFULL, 0, (lmins+1)*sizeof(double));
    fflush(stdout);
  }
  if ( tzero_flag && (time < opt.t8) ){
    printf(" %e ", opt.t8);
    for (i = 0; i < dim; i++)
       printf("%e ", fabs(pt[i]));
     printf("\n");
  }
  printf("# of iterations: %d\n", count);

  /*** end solve fundamental equation ***/

  if(opt.method=='F'){
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
  free(p0);
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
  /*  	                   -beta(E_s-E_j) */
  /*  rate:  j->i:  prop  e               */
  /*                       -beta(E_s-E_i) */
  /*         i->j   prop  e               */
  /****************************************/
  
  int i,j,real_abs = 0;
  double m_saddle, Zabs, abs_rate, *T, *U;

  U = (double *) MxNew (dim*dim*sizeof(double));
  D  = (double *) MxNew ((dim+1)*(dim+1)*sizeof(double)); /* BRAUCH'MA DES IMMA ? */
  for(i=0;i<((dim+1)*(dim+1));i++) D[i] = 1.;
  
  for( i = 0; i < dim; i++)
    for( j = i+1; j < dim; j++){
      m_saddle = max_saddle(i, j, Data);
      /* rate j -> i */
      U[dim*i+j] = D[dim*i+j]*exp(-(m_saddle-Data[j].FGr)/_kT); 
      /* rate i -> j */
      U[dim*j+i] = D[dim*i+j]*exp(-(m_saddle-Data[i].FGr)/_kT);
    }
  if(D != NULL) free(D);

  if(opt.absrb){ /*==== absorbing  states ====*/
    dim++;
    fprintf(stderr, "dim increased to %i\n", dim);
    T = (double *) MxNew(dim*dim*sizeof(double));
    real_abs = opt.absrb; /* the original absorbing lmin */
    real_abs--;
    opt.absrb = dim; /* the 'new' abs state = last row/column of rate matrix */
    fprintf(stderr, "new absorbing lmin is: %i\n", opt.absrb);
    Zabs = exp((-Data[real_abs].FGr)/_kT);
    abs_rate = exp((-Data[real_abs].energy)/_kT)/Zabs;

    for(i = 0; i < (dim-1); i++){ /* all except the last row */ 
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
    /* calculate colum sum */
    for(i = 0; i < dim; i++) tmp += U[dim*i+j];
    U[dim*j+j] = -tmp+1; /* make U a stochastic matrix */
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
  
  if(opt.absrb){ /*==== absorbing  states ====*/
    for(i = 0; i < dim; i++)
      R[dim*i+(opt.absrb-1)] = 0. ;
  }              /*== end absorbing states ==*/
  
  /* set diagonal elements  to 0 */
  for (i = 0; i < dim; i++) R[dim*i+i] = 0;
  for (j = 0; j < dim; j++) {
    double tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++) tmp += R[dim*i+j];
     R[dim*j+j] = -tmp+1; /* make U a stochastic matrix */
  }

  if (opt.want_verbose) MxPrint (R, "R with Methode F", 'm');
  return R;
}

/*==*/
static double*
MxMethodeINPUT (BarData *Data, double *Input)
{  
  int i, j, real_abs = 0;
  double *U=NULL, Zabs, abs_rate;

  if (opt.want_verbose) MxPrint(Input, "Input Matrix", 'm');
    
  if (opt.absrb) {  /*==== absorbing  states ====*/
    dim++;
    fprintf(stderr, "dim increased to %i\n", dim);
    U = (double *) MxNew(dim*dim*sizeof(double));
    real_abs = opt.absrb; /* the original absorbing lmin */
    real_abs--;
    opt.absrb = dim; /* the 'new' abs state = last row/column of rate matrix */
    fprintf(stderr, "new absorbing lmin is: %i\n", opt.absrb);
    Zabs = exp((-Data[real_abs].FGr)/_kT);
    abs_rate = exp((-Data[real_abs].energy)/_kT)/Zabs;

    for(i = 0; i < (dim-1); i++){ /* all except the last row */ 
      for(j = 0; j < (dim-1); j++)
	U[dim*i+j] = Input[(dim-1)*i+j];
      U[(dim-1)*j+(dim-1)] = 0.;
    }
    for(j = 0; j < dim; j++) /* last row */
      U[dim*(dim-1)+j] = 0.;
    U[dim*(dim-1)+real_abs] = abs_rate;
  /*   if(opt.want_verbose) MxPrint(U, "aufgeblasene Matrix", 'm'); */
  }      /*== end absorbing states ==*/
  else{  /*== non-absorbing states ==*/
    U = (double *) MxNew(dim*dim*sizeof(double));
    memcpy(U, Input, dim*dim*sizeof(double));
  }      /*== end non-absorbing states ==*/

  /* diagonal elements */
  for (i = 0; i < dim; i++) U[dim*i+i] = 0;
  fprintf(stderr, "dim is %i\n", dim);
  for (i = 0; i < dim; i++) {
    double tmp = 0.00;
    /* calculate column sum */
     for(j = 0; j < dim; j++)  tmp += U[dim*j+i];
    U[dim*i+i] = -tmp+1.;   /* make U a stochastic matrix */
  }
  if(opt.want_verbose) MxPrint (U,"U with Methode I" , 'm');
  if (opt.dumpU) MxBinWrite(U);

  free(Input);
  return U;
}

/*==*/
static double
max_saddle(int i, int j, BarData *Data)
{
  int tmp;

  if(Data[i].number > Data[j].father){ /* exchange i & j */
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
static void*
MxNew ( size_t size )
{
  void *mx = NULL;
  if ( (mx = (void *) calloc (1, size)) == NULL )
    fprintf (stderr, "ERROR: new_martix() allocation failed\n");

  return mx;
}

/*==*/
/* print matrix stored in ccmath-fromat */
void
MxPrint(double *mx, char *name, char T)
{
  int k, l;
  switch (T) {
  case 'm':    /* square matrix */
    fprintf(stderr,"%s:\n", name);
    for (k = 0; k < dim; k++) {
      for (l=0; l< dim; l++) {
	fprintf(stderr,"%10.4g ", mx[dim*k+l]);
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"---\n");
    break;
  case 'v':
    fprintf(stderr,"%s:\n", name);    
    for (k = 0; k < dim; k++) fprintf(stderr,"%15.10f ", mx[k]);
    fprintf(stderr,"\n---\n");
    break;
  default:
    fprintf(stderr,"ERROR MxPrint(): no handler 4 type %c\n", T);
  }
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
  fclose(opt.INFILE);
}

/*==*/
static void
MxDoDegeneracyStuff(void)
{
  int i, j, b, nr, current, numsad = 1, count = 0;

  numsad = ParseSaddleFile(&saddle);
  /* loop over all elements of structure-array saddle: */
  /* first we fill the upper triangle */
  for (count = 0;count < numsad; count++) {
    current = 1;
    nr = saddle[count].list[0];
    /* only for saddles with a cc >= 1  AND those which connect at least 2 lmins */
    if(saddle[count].cc >= 1 && nr >= 2 && !(saddle[count].cc == 1 && nr == 2)) {
      while(current < nr){
	for(b = current+1; b <= nr; b++){
	  /* skip in case a lmin which we don't see is connected by */
	  /* the saddle (because we can only see --max x lmins) */
	  if(saddle[count].list[current] > dim || saddle[count].list[b] > dim){
	    current++;
	    continue;
	  }
	  /* FIRST: we consider the size of the cc the saddle belongs to */
	  if(saddle[count].cc > 1){
	    D[dim * (saddle[count].list[current]-1) + (saddle[count].list[b]-1)] += (saddle[count].cc - 1);
	    fprintf(stderr, "transition between %3d - %3d: adding %2d  cc\n",
		    saddle[count].list[current], saddle[count].list[b], saddle[count].cc-1);
	    /* -1 because if the size of cc == 1 there's nothing special */
	    /* about it, i.e. there is just one saddle */
	  }
	  /* SECOND: we consider that the saddle connects several lmins */
	  if(nr > 2){
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
  
  for (j=0;j<dim;j++){
    sumsq=0.0;
    for (i=0;i<dim;i++)
      sumsq += SQ(mx[dim*i+j]);
    if(sumsq > 0.0)
      sumsq=1./sqrtl(sumsq);
    for (i=0;i<dim;i++)
      mx[dim*i+j] *= sumsq;
  }
  return;
}
	  
/*==*/
static void
MxFixevecs(double *evecs, double *evals)
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
  for (i=0;i<dim;i++) {
    if (evals[i]==1.0) 	{
      abscount++;
      maxent=0.;
      maxind=0;
      for (j=0;j<dim;j++) {
	if (evecs[dim*j+i] > maxent) {
	  maxent=evecs[dim*j+i];
	  maxind=j;
	}
      }
      evecs[dim*maxind+i]=1.0;
      for (j=0;j<dim;j++) {
	if (j==maxind) continue;
	evecs[dim*j+i]=0.0;
      }
    }
  }

  /* repair messed non abs. eigenvectors */
  /* using all abs. states equally */
  for (j=abscount;j< dim;j++){
    double mu=0.;
    for(i=0;i<dim;i++)
      mu += evecs[dim*i+j];
    for (i=0;i<abscount;i++)
      evecs[dim*i+j] -= mu/(double)abscount;
  }
  norm2(evecs);
  
  if (opt.want_verbose){
    MxPrint (evals, "evals_complete", 'v');
    MxPrint (evecs, "evecs_complete", 'm');
    fflush(stdout); fflush(stderr);
  }
  
  fprintf(stderr,"colsums: ");
  for (i=0;i<dim;i++) {
    csum=0.0;
    for (j=0;j<dim;j++)	
      csum += evecs[dim*j+i];
    fprintf(stderr,"%g ", csum);
  }
  fprintf(stderr,"\n");

  if (opt.absrb && (abscount > 1)){
    int i,j;
    fprintf(stderr, "\nWARNING: found %d additional absorbing state(s): ",
	    abscount-1);
    for(i=1;i<dim;i++){
      if(evals[i] == 1){
	for(j=0;j<dim;j++){
	  if(evecs[dim*j+i]==1)
	    fprintf(stderr, " %5d", j+1);
	}
      }
    }
    fprintf(stderr, "\n");
  }
  fflush(stderr);
}

/*==*/
/* sort evecs,evals */
static void
MxSortEig(double *evals, double *evecs)
{
  int i,j,k;
  double p;
  
  for (i=0;i<dim;i++) {
    p=evals[k=i];
    for (j=i+1;j<dim;j++)
      if (evals[j] >= p) p=evals[k=j];
    if (k != i) {
      evals[k]=evals[i];
      evals[i]=p;
      for (j=0;j<dim;j++) {
	p=evecs[dim*j+i];
	evecs[dim*j+i]=evecs[dim*j+k];
	evecs[dim*j+k]=p;
      }
    }
  }
}

/*==*/
static void
MxBinWrite(double *matrix)
{
  int i, j;
  FILE *BINOUT;
  char *binfile = "mx.bin";

  BINOUT = fopen(binfile, "w");
  if (!BINOUT){
    fprintf(stderr, "could not open file pointer 4 binary outfile\n");
    exit(EXIT_FAILURE);
  }
  /* first write dim to file */
  fwrite(&dim,sizeof(int),1,BINOUT);
  /* then write matrix entries */
  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      fwrite(&matrix[dim*i+j],sizeof(double),1,BINOUT);

  fprintf(stderr, "matrix written to binfile\n");
  fclose(BINOUT);
}

/*==*/
static void
MxASCIIWrite(double *matrix)
{
  int i, j;
  FILE *ASCIIOUT;
  char *asciifile = "mx.txt";

  ASCIIOUT = fopen(asciifile, "w");
  if (!ASCIIOUT){
    fprintf(stderr, "could not open file pointer 4 ASCII outfile\n");
    exit(EXIT_FAILURE);
  }
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      fprintf(ASCIIOUT,"%15.10g ", matrix[dim*i+j]);
    }
    fprintf(ASCIIOUT,"\n");
  }
  fprintf(stderr, "matrix written to ASCII file\n");
  fclose(ASCIIOUT);
}

/*==*/
void
MxExponent(double *p0, double *p8, double *U)
{
  int i,j, pdiff_counter = 0;
  double x, time, *Uexp, *Umerk, *pt, *pdiff, check = 0.;

  Umerk  = (double *) MxNew (dim*dim*sizeof(double));
  Uexp   = (double *) MxNew (dim*dim*sizeof(double));
  pt     = (double *) MxNew (dim*sizeof(double));
  pdiff  = (double *) MxNew (dim*sizeof(double));

  memcpy(Umerk, U, dim*dim*sizeof(double));
  
  for (i=0; i<dim; i++) U[(dim+1)*i] -= 1;
  print_settings();
  for (time = opt.t0; time <= opt.t8; time *= opt.tinc) {
    memcpy(U, Umerk, dim*dim*sizeof(double));
    for (i=0; i<dim*dim; i++) U[i]*=time;
    padexp(U,Uexp,dim,30);
    x = 0.;
    for(j=0;j<dim*dim;j++) x+=Uexp[j];
    for(j=0;j<dim*dim;j++) Uexp[j]*=(double)dim/x;
    vmul(pt, Uexp, p0, dim);
    /* check convergence */
    for(i=0; i<dim; i++){
      pdiff[i] = p8[i] - pt[i];
      if (fabs(pdiff[i]) >= 0.0001)
	pdiff_counter++;
      }
    if (pdiff_counter < 1)
      break;
    pdiff_counter = 0.;
    /* end check convergence */
    check = 0.;
    printf(" %e ", time);  /* print p(t) to stdout */
    for (i = 0; i < dim; i++){
      if(pt[i] < -0.00001){
	fprintf(stderr, "prob of lmin %i has become negative!\n", i+1);
	exit(EXIT_FAILURE);
      }
      printf("%e ", fabs(pt[i]));
      check += fabs(pt[i]); 
    }
    printf("\n");
  
    if ( ((check-1) < -0.01) || ((check-1) > 0.01) ){
      fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
      exit(EXIT_FAILURE);
    }
    memset(pt,   0, dim*sizeof(double));
    memset(pdiff, 0, dim*sizeof(double));
    memset(Uexp, 0, dim*dim*sizeof(double));
    memset(U, 0, dim*dim*sizeof(double));

  }
  free(Uexp);
  free(pt);
}

/*==*/
void
MxFPT(double *U, double *p8)
{
  int i,j, val;
  double *M, *Q, *W, *Z, *c,*d;
  
  M = (double *) MxNew (dim*dim*sizeof(double));
  Q = (double *) MxNew (dim*dim*sizeof(double));
  W = (double *) MxNew (dim*dim*sizeof(double));
  Z = (double *) MxNew (dim*dim*sizeof(double));
  c = (double *) MxNew (dim*sizeof(double));
  d  = (double *) MxNew (dim*sizeof(double));

  memcpy(Q,U, dim*dim*sizeof(double));
  for(i=0; i < dim; i++)
    Q[dim*i+i] -= 1;
  
  for (j = 0; j < dim; j++)
    for(i = 0; i < dim; i++)
      W[dim*i+j] = p8[i];
  MxPrint(Q, "Q (MxFPT)", 'm'); MxPrint(W, "W (MxFPT)", 'm'); 

  for(i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      Z[dim*i+j] = W[dim*i+j] - Q[dim*i+j];
  MxPrint(Z, "W-Q == I-P+W (MxFPT)", 'm');

  val = minv(Z, dim);
  if (val != 0){
    fprintf(stderr, "Z is singular, can't invert it\n");
    exit(EXIT_FAILURE);
  }
  MxPrint(Z, "Z (MxFPT)", 'm');
  free(Q);  free(W);
  for(i = 0; i < dim; i++) c[i]=1.;
  MxPrint(c, "c", 'v');
  vmul(d, Z, c,dim);
  MxPrint(d, "d", 'v');
  for(i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      M[dim*i+j] = (Z[dim*j+j] - Z[dim*i+j])/p8[j]; /* !~!!! */
  MxPrint(M, "M (MxFPT)", 'm');
  free(M); free(Z); free(c); free(d);
}

/*==*/
static void
MxKotzOutMathematica(double *matrix)
{
  int i,j;
  FILE *MATHEMATICA_OUT=NULL;
  char *mathematica_file = "mxMat.txt";
  
  MATHEMATICA_OUT = fopen(mathematica_file, "w");
  if (!MATHEMATICA_OUT){
    fprintf(stderr, "could not open file pointer 4 Matematica outfile\n");
    exit(EXIT_FAILURE);
  }
  fprintf(MATHEMATICA_OUT, "{");
  for(i=0;i<dim;i++){
    fprintf(MATHEMATICA_OUT, "{");
    for(j=0;j<dim;j++){
      if (j != (dim-1))
	fprintf(MATHEMATICA_OUT, "%25.22f, ", matrix[dim*i+j]);
      else
	fprintf(MATHEMATICA_OUT, "%25.22f}", matrix[dim*i+j]);
    }
    if (i != (dim-1))
      fprintf(MATHEMATICA_OUT, ",\n");
    else
      fprintf(MATHEMATICA_OUT, "}\n");
  }  
  fclose(MATHEMATICA_OUT);
}

/*==*/
void
MxFirstPassageTime(double *U, double *p8)
{
#define MFPT
#ifdef MFPT
  int i,j;
  { /* calculate MFPT */
    /* absorbing case: Z = inv(I-P+W) */
    /*   double *t=NULL, *c=NULL, *Q=NULL, *N=NULL; */
    /*     Q = (double *) MxNew(dim*dim*sizeof(double)); */
    /*     N = (double *) MxNew(dim*dim*sizeof(double)); */
    /*     c = (double *) MxNew (dim*sizeof(double)); */
    /*     t = (double *) MxNew (dim*sizeof(double)); */
    /*     for(i = 0; i < dim; i++){ */
    /*       for(j = 0; j < dim; j++){ */
    /* 	Q[dim*i+j] = Input[dim*i+j]; */
    /* 	if(i==j) Q[dim*i+j] = 1- Q[dim*i+j]; */
    /*       } */
    /*     } */
    /*     if(opt.want_verbose) MxPrint (Q,"Q" , 'm'); */
    /*     mcopy(N, Q, dim*dim); */
    /*     minv(N,dim); */
    /*     if(opt.want_verbose) MxPrint (N,"N" , 'm'); */
    /*     for(i = 0; i < dim; i++) c[i]=1; */
    /*     if(opt.want_verbose) MxPrint(c, "c", 'v'); */
    /*     vmul (t, N, c, dim); */
    /*     if(opt.want_verbose) MxPrint(t, "t", 'v'); */
    /*     free(Q);free(N);free(c);free(t); */
    /*==*/
    
    /* ergodic case: Z = inv(I-P+W) */
    fprintf(stderr, "in MxFirstPassageTime\n");
    if(opt.want_verbose) MxPrint (U,"U" , 'm');
    double *Z=NULL, *M=NULL, *I=NULL;
    Z = (double *) MxNew(dim*dim*sizeof(double));
    M = (double *) MxNew(dim*dim*sizeof(double));
    I = (double *) MxNew(dim*dim*sizeof(double));

    for(i=0; i < dim; i++)
      I[dim*i+i] = 1;
    
    for(i = 0; i < dim; i++)
      for(j = 0; j < dim; j++){
	  Z[dim*i+j] = I[dim*i+j]-U[dim*i+j]; /* I-U */
      }
    if(opt.want_verbose) MxPrint (Z,"I-U" , 'm');
    for(i = 0; i < dim; i++)
      for(j = 0; j < dim; j++)
	Z[dim*i+j] += p8[i];  /* I-U+W */
    
    if(opt.want_verbose) MxPrint (Z,"I-U+W" , 'm');
    minv(Z,dim);
    if(opt.want_verbose) MxPrint (Z,"Fundamental matrix Z=inv(I-U+W)" , 'm');
    
    for(i = 0; i < dim; i++){
      for(j = 0; j < dim; j++){
	M[dim*i+j] = (Z[dim*i+i]-Z[dim*i+j])/p8[j];
      }
    }
    if(opt.want_verbose) MxPrint (M,"M" , 'm');
    free(Z);free(M);free(I);
  }
#endif
  return;
}

static void
MxEVLapackSym(double *U){
/*   extern void dsyev_(char *jobz, char *uplo, int *n, double *A,int *lda, */
/* 		     double *w, double *work, int *lwork, int *info); */
  extern void dsyevx_(char *jobz,char *range,char *uplo,int *n,double *A,
		      int *lda,double *vl,double *vu,int *il,int *iu,
		      double *abstol,int *m,double *w,double *z,int *ldz,
		      double *work,int *lwork,int *iwork,int *ifail,int *info);
  double *work=NULL;
  double abstol=FEPS ,tmp, vl,vu;
  int i,j,lwork, il,iu,nfo,m, *iwork=NULL, *ifail=NULL;

  lwork = dim*(dim+6);
  work   = (double *) malloc (lwork * sizeof(double));
  iwork  =    (int *) malloc (5*(dim) * sizeof(int));
  ifail  =    (int *) malloc (dim * sizeof(int));
  
  dsyevx_("V","A","U",&dim,U,&dim,&vl,&vu,&il,&iu,
	 &abstol,&m, evals,evecs,&dim,work, &lwork, iwork,
	 ifail, &nfo );
  //dsyev_("V","U",&dim, S, &dim, evals, work, &lwork, &nfo);
  if(nfo != 0){
    fprintf(stderr, "dsyevx exited with value %d\n", nfo);
    exit(EXIT_FAILURE);
  }

  free(work);
  free(iwork);
  free(ifail);
  for (i=0;i<dim;i++)
    for (j=i+1;j<dim;j++){
      tmp = evecs[dim*i+j];evecs[dim*i+j]=evecs[dim*j+i];evecs[dim*j+i]=tmp;}
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
  
  double *evals_re=NULL, *evals_im=NULL,*scale=NULL;
  double *rconde=NULL,*rcondv=NULL,*work=NULL, abnrm;
  int one,dimx2,ilo,ihi,lwork, *iwork,nfo;
  /* for sorting */
  int i,j,k;
  float p;
  double tmp;
  
  dim=dim+500;
  one=1;
  dimx2 = 2*dim;
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
	rcondv && work && iwork)==0){
    fprintf(stderr,"no space for temporary lapack arrays!\n");
    exit(EXIT_FAILURE);
  }
 
  /* instead of more fiddling, we transpose the input */
  for (i=0;i<dim;i++)
    for (j=i+1;j<dim;j++){
      tmp = origU[dim*i+j];origU[dim*i+j]=origU[dim*j+i];origU[dim*j+i]=tmp;}
  
  dgeevx_("B","N","V","V",&dim,origU,&dim,evals_re,evals_im,NULL,&one,\
	  evecs ,&dim, &ilo, &ihi,  scale,&abnrm ,rconde, rcondv, work, \
	  &lwork  , iwork, &nfo);
    
  for (i=0;i<dim;i++) evals[i]=evals_re[i];
  for (i=0;i<dim;i++) 
    if ((evals_re[i] != 0.0) && fabs(evals_im[i]/evals_re[i])>1.e-15)
      fprintf(stderr,"eigenvalue %d is %g + i*%g, which is somewhat complex\n",
	      i,evals_re[i],evals_im[i]);

  /*transpose output*/
  for (i=0;i<dim;i++)
    for (j=i+1;j<dim;j++){
      tmp = evecs[dim*i+j];evecs[dim*i+j]=evecs[dim*j+i];evecs[dim*j+i]=tmp;}
  
  free(evals_re);
  free(evals_im);
  free(scale);
  free(rconde);
  free(rcondv);
  free(work);
  free(iwork);
}


/* End of file */
