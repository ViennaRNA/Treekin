/* calc.c */
/* Last changed Time-stamp: <2003-10-03 19:17:12 mtw> */
/* static char rcsid[] = "$Id: calc.c,v 1.18 2003/10/06 12:19:34 mtw Exp $"; */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "exp_matrix.h" /* functions for matrix-exponent stuff */
#include "mxccm.h"      /* functions for eigen-problems stolen from ccmath */
#include "barparser.h"  /* functions for input */
#include "matrix.h"     /* basic matrix functions from meschach */
#include "matrix2.h"    /* advances matrix functions from meschach */
#include "calc.h"       /* does all matrix stuff for markov process */
#include "globals.h"    /* contains getopt-stuff */

#define TOL 0.000000000000001
#define ABS_VAL 0.00001

/* private function(s) */
static void   *MxNew (size_t size);
static double *MxMethodeA (TypeBarData *Data);
static void    MxPrint(double *mx, char *name, char T);
static void    MxPrintMeschachMat(MAT *matrix, char *name);
static void    MxPrintMeschachVec(VEC* vector, char *name);
static void    MxMeschach2ccmath(MAT *meschach_matrix, double **origM);
static void    MxMeschach2ccmathVec(VEC *meschach_vector, double **origV);
static MAT    *Mxccmath2Meschach(double* ccmath_matrix);
static double  max_saddle(int i, int j, TypeBarData *Data);
static void    print_settings(void);
static char   *time_stamp(void);
static void    MxDoDegeneracyStuff(void);
static int     Mxempty(MAT *matrix);
static void    MxBinWrite (double *matrix);

/* private vars and arrays */
static int      dim = 0;
static double   _kT = 1.;
static double  *EV        = NULL;   /* array 4 eigenvalues */
static double  *EV_mesch  = NULL;   /* ev from meschach-routine, ie in other order */
static double  *_sqrPI    = NULL;   /* left constant array */
static double  *sqrPI_    = NULL;   /* right constant array */
static double  *D         = NULL;   /* matrix with degree of degeneacy */
static char     Aname[30];
static TypeDegSaddle *saddle;

/*==*/
void MxInit (int d) {
  int i;
  _kT = 0.00198717*(273.15 + opt.T);
  if (d > 0 ) dim = d;
  else { fprintf(stderr, "dim <= 0\n"); exit(1); }
  
  EV     = (double *) MxNew (dim*sizeof(double));
  _sqrPI = (double *) MxNew (dim*dim*sizeof(double));
  sqrPI_ = (double *) MxNew (dim*dim*sizeof(double));
  D      = (double *) MxNew (dim*dim*sizeof(double));
  for(i=0;i<dim*dim;i++) D[i] = 1.;
}

/*==*/
double *MxBar2Matrix ( TypeBarData *Data, double *R) {

  double *U;

  if(opt.want_degenerate) MxDoDegeneracyStuff();
  
  switch (opt.method) {
  case 'A':
    U = MxMethodeA(Data);
    break;
  case 'I':
    U = MxMethodeINPUT(Data, R);
    break;
  default:
    fprintf (stderr,
	     "ERROR in MxBar2Matrix(): No handler 4 method %c\n", opt.method);
    exit(1);
  }
  if (opt.dumpU)
    MxBinWrite(U);
  return (U);
}

/*==*/
double *MxStartVec (void) {

  int i;
  double *p0;

  p0 = (double *) MxNew(dim*sizeof(double));
  for (i = 1; i < (int) *opt.pini; i+=2)
    p0[(int)opt.pini[i]-1] = (double)opt.pini[i+1];
  /* -1 because our lmins start with 1, not with 0 (as Data does ) */

  MxPrint (p0, "p0", 'v');
  return (p0);
}

/*==*/
/* calculate equilibrium distribution */
double *MxEqDistr ( TypeBarData *Data ) {
  int i;
  double *p8, Z = 0.;
  
  p8 = (double *) MxNew (dim*sizeof(double));

  if(opt.absrb) dim--; /* because it was inceased before and we */
  /* only have dim lmins in bar-file, not dim+1 */
  for(i = 0; i < dim; i++) Z += exp(-((double)Data[i].FGr/_kT));
  for(i = 0; i < dim; i++) p8[i] = exp(-((double) Data[i].FGr/_kT))/Z;
  if(opt.absrb) dim++;
  
  if(opt.absrb){
    double tmp = 0.;
    for(i = 0; i < dim; i++){
      p8[i] = ABS_VAL;
      tmp += p8[i];
    }
    p8[opt.absrb-1] = 1.0-tmp;
  }
  if(opt.want_verbose) MxPrint(p8, "p8", 'v');
  return (p8);
}

/*==*/
double *MxEqDistrFULL (SubInfo *E) {
  int i;
  double *p8, Z = 0.;
  
  p8 = (double *) MxNew (dim*sizeof(double));

  /*  for(i = 0; i < dim; i++) fprintf (stderr, "E[%i].energy: %lf\n", i, E[i].energy); */
  for(i = 0; i < dim; i++) Z += exp(-E[i].energy/_kT);
  for(i = 0; i < dim; i++) p8[i] = exp(-E[i].energy/_kT)/Z;

  if(opt.absrb){
    double tmp = 0.;
    for(i = 0; i < dim; i++){
      p8[i] = ABS_VAL;
      tmp += p8[i];
    }
    p8[opt.absrb-1] = 1.0-tmp;
    /* fprintf(stderr, "set p8[%i] to %e\n", opt.absrb-1, p8[opt.absrb-1]); */
  }
  MxPrint (p8, "p8", 'v');
  return (p8);
}

/*==*/
double *MxSymmetr ( double *U, double *P8 ) {

  int i, j;
  double *S, *tmpMx;
  
  S     = (double *) MxNew (dim*dim*sizeof(double));
  tmpMx = (double *) MxNew (dim*dim*sizeof(double));
    
  for(i = 0; i < dim; i++) {
    for(j = 0; j < dim; j++) {
      if( i == j) {
	sqrPI_[dim*i+j] = sqrt(P8[i]);            /* pos right */
	_sqrPI[dim*i+j] = 1/(sqrPI_[dim*i+j]);}}} /* neg left */
  fprintf(stderr,"dim is %d\n", dim);
  mmul(tmpMx, _sqrPI, U, dim); mmul(S, tmpMx, sqrPI_, dim);
    if (opt.want_verbose) {
          MxPrint (_sqrPI, "_sqrPI (= negative Wurzl von pi)", 'm');
          MxPrint (sqrPI_, "sqrPI_ (= Wurzl von pi)", 'm');
          MxPrint (tmpMx, "tmpMx (= negative Wurzl von pi * U)", 'm');
          MxPrint (S, "S before (= tmpMx * Wurzl von pi)", 'm'); }
  free(tmpMx);
  /* correct for numerical errors */
  for (i = 0; i < dim; i++) {
    for (j = 0; j < dim; j++) {
      S[dim*i+j] = (S[dim*i+j]+S[dim*j+i])/2;
      S[dim*j+i] = S[dim*i+j];
    }
  }
  
  if (opt.want_verbose) MxPrint (S, "force symmetrized S", 'm');
  
  eigen(S, EV, dim);  /* S is overwritten with its eigenvectors */ 
  if (opt.want_verbose){
    MxPrint (S, "Eigenvectors of S", 'm');MxPrint(EV, "Eigenvalues of S", 'v'); }
  for (i=0;i<dim;i++){
    if (EV[i] > 1.){
      fprintf(stderr, "\nEV[%i] is > 1: %20.18g, now setting to 1.\n", i, EV[i]);
      EV[i] = 1.;
    }
  }
  
  /* compensate 4 translation of matrix U */
  for(i = 0; i < dim; i++) EV[i] = EV[i] - 1;
  return (S);
}

/*==*/
/* S which comes into this function contains the (right) eigenvectors  */
/* calculated either by ccmath (non-abs) or meschach (absorbing case) */
void MxIterate ( double *p0, double *p8, double *S) {
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
  int i, count = 0, pdiff_counter = 0;
  double time, check = 0.;
  double *CL, *CR, *exptL, *tmpVec, *tmpVec2, *St;
  double *pt, *pdiff;  /* probability distribution/difference 4 time t */
  
  CL        = (double *) MxNew (dim*dim*sizeof(double));
  exptL     = (double *) MxNew (dim*dim*sizeof(double));
  tmpVec    = (double *) MxNew (dim*sizeof(double));
  tmpVec2   = (double *) MxNew (dim*sizeof(double));
  pt        = (double *) MxNew (dim*sizeof(double));
  pdiff     = (double *) MxNew (dim*sizeof(double));

  if(! opt.absrb){  /* NON-absorbing case */
    St     = (double *) MxNew (dim*dim*sizeof(double));
    CR     = (double *) MxNew (dim*dim*sizeof(double));
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
    for(i = 0; i < dim; i++) EV_mesch[i] -= 1;  /* compensate 4 translation of matrix U */
    free(EV);         /* was allocated in MxInit */
    EV = EV_mesch;    /* let EV point at EV_mesch */
    if(opt.want_verbose) MxPrint(EV_mesch, "EV_mesch in MxIterate", 'v');
    vmul (tmpVec, S_inv, p0, dim);
    free(S_inv);
  }
  
  /*** solve fundamental equation ***/  /* using logarithmic time-scale */
  print_settings();
  for (time = opt.t0; time <= opt.t8; time *= opt.tinc) {
    for (i = 0; i < dim; i++) 
      exptL[dim*i+i] = exp(time*EV[i]);

    vmul(tmpVec2, exptL, tmpVec, dim);
    if(!opt.absrb)  vmul(pt, CL, tmpVec2, dim);
    else            vmul(pt, S, tmpVec2, dim);
 
    count++;  /* # of iterations */
    
    printf(" %e ", time);  /* print p(t) to stdout */
    for (i = 0; i < dim; i++){
      if(pt[i] < -0.01){
	fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, pt[i]);
	exit(866);
      }
      printf("%e ", fabs(pt[i]));
      check += fabs(pt[i]); 
    }
    printf("\n");

    if ( ((check-1) < -0.05) || ((check-1) > 0.05) ){
      fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
      exit(888);
    }
    check = 0.;
    
    /* now check if we have converged yet */
    for(i=0; i<dim; i++){
      pdiff[i] = p8[i] - pt[i];
      if (fabs(pdiff[i]) >= 0.001)
	pdiff_counter++; /* # of lmins whose pdiff is > the threshold */
    }
    if (pdiff_counter < 1) /* all mins' pdiff lies within threshold */
      break;
    pdiff_counter = 0.;
    memset(pdiff, 0, dim*sizeof(double));
    /* end check of convergence */
  }
  printf("# of iterations: %d\n", count);
  /*** end solve fundamental equation ***/ 

  free(EV);
  free(exptL);
  free(CL);
  free(tmpVec);
  free(tmpVec2);
  free(pt);
  free(pdiff);
  free(p0);
/*    free(EV_mesch); */
}


/*==*/
void MxIterate_FULL (double *p0, double *p8, double *S, int lmins) {
  /*
    solve following equation 4 various times
    p(t) = sqrPI_ * S * exp(time * EV) * St * _sqrPI * p(0)
    CL = sqrPI_ * S
    CR = St * _sqrPI
    tmpVec = CR * p(0)
  */
  int i,  pdiff_counter = 0;
  double time, *CL, *CR, *exptL, *tmpMx, *tmpVec, *tmpVec2, *pt, *St;
  double *ptFULL;    /* prob dist 4 of the effective lmins of the tree at time t */
  double *p8FULL;    /* equ dist 4 gradient basins, full process */
  double *pdiffFULL; /* population prob difference between p8FULL and ptFULL */ 
  double check = 0. , checkp8 = 0.;
  
  St        = (double *) MxNew (dim*dim*sizeof(double));
  CL        = (double *) MxNew (dim*dim*sizeof(double));
  CR        = (double *) MxNew (dim*dim*sizeof(double));
  tmpVec    = (double *) MxNew (dim*sizeof(double));
  tmpVec2   = (double *) MxNew (dim*sizeof(double));
  pt        = (double *) MxNew (dim*sizeof(double));
  exptL     = (double *) MxNew (dim*dim*sizeof(double));
  tmpMx     = (double *) MxNew (dim*dim*sizeof(double));
  ptFULL    = (double *) MxNew ((lmins+1)*sizeof(double));
  p8FULL    = (double *) MxNew ((lmins+1)*sizeof(double));
  pdiffFULL = (double *) MxNew ((lmins+1)*sizeof(double));
  
  mcopy(St, S, dim*dim);  trnm(St, dim); /* transpose S */
  mmul (CL, sqrPI_, S, dim);
  mmul (CR, St, _sqrPI, dim);
  free(St);
  
  vmul (tmpVec, CR, p0, dim);
  free(CR);
  /* calculate equilibrium distribution once */
  for (i = 0; i < dim; i++)
    p8FULL[E[i].ag] += p8[i]; /* eq distr of the gradient basins */
  for (i = 0; i < lmins; i++) checkp8 += fabs(p8FULL[i]);
  if ( ((checkp8-1) < -0.1) || ((checkp8-1) > 0.1) ){
    fprintf(stderr, "overall equilibrium probability is %e != 1. ! exiting\n", checkp8);
    exit(888);
  }
  
   /* solve fundamental equation */
  for (time = opt.t0; time <= opt.t8; time *= opt.tinc) {
    for (i = 0; i < dim; i++) 
      exptL[dim*i+i] = exp(time*EV[i]);

    vmul(tmpVec2, exptL, tmpVec, dim);
    vmul(pt, CL, tmpVec2, dim);
    
    /*    mmul (tmpMx, CL, exptL, dim); vmul (pt, tmpMx, tmpVec, dim); */
        
    for (i = 0; i < dim; i++){
      if(pt[i] < -0.01){
	fprintf(stderr, "prob of lmin %i has become negative: %6.4f\n", i+1,pt[i]);exit(866);
      }
      ptFULL[E[i].ag] += pt[i]; /* map individual structure -> gradient basins */
      check += fabs(pt[i]); 
    }
 
    printf(" %e ", time);
    for(i = 1; i <= lmins; i++) printf("%e ", fabs(ptFULL[i]));
    printf("\n");
    
    if ( ((check-1) < -0.01) || ((check-1) > 0.01) ){
      fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
      exit(888);
    }
    
    /* now check if we have converged yet */
    /*  printf("#---------------"); */
    for(i = 1; i <= lmins; i++){
      pdiffFULL[i] = p8FULL[i] - ptFULL[i];
      /*  printf("%7.4f ", pdiffFULL[i]); */
      if (fabs(pdiffFULL[i]) >= 0.001)
	pdiff_counter++; /* # of lmins whose pdiff is > the threshold */
    }
    /*   printf("  pdiff_counter: %i", pdiff_counter); */
    /*     printf("\n"); */
    if (pdiff_counter < 1) /* all mins' pdiff lies within threshold */
      break;
    pdiff_counter = 0;
    memset(pdiffFULL, 0, (lmins+1)*sizeof(double));
    /* end check of convergence */
    
    memset(ptFULL, 0, (lmins+1)*sizeof(double));
    check = 0.;
    fflush(stdout);
  }
  /* end solve fundamental equation */ 
  free(E);
  free(tmpVec2);
  free(EV);
  free(exptL);
  free(CL);
  free(tmpVec);
  free(tmpMx);
  free(pt);
  free(ptFULL);
  free(p8FULL);
  free(pdiffFULL);
  free(p0);
}

/*==*/
static double *MxMethodeA (TypeBarData *Data) {
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
  
  int i,j,real_abs = 0;;
  double m_saddle, Zabs, abs_rate, *T, *U;

  U = (double *) MxNew (dim*dim*sizeof(double));

  for( i = 0; i < dim; i++)
    for( j = i+1; j < dim; j++){
      m_saddle = max_saddle(i, j, Data);
      /* rate j -> i */
      U[dim*i+j] = D[dim*i+j]*exp(-(m_saddle-Data[j].FGr)/_kT); 
      /* rate i -> j */
      U[dim*j+i] = D[dim*i+j]*exp(-(m_saddle-Data[i].FGr)/_kT);
    }
  if(D != NULL) free(D);

  MxPrint(U, "original rate matrix", 'm');
  
  if(opt.absrb){ /*==== absorbing  states ====*/
    dim++;
    fprintf(stderr, "dim inceased to %i\n", dim);
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
    MxPrint(U, "aufgeblasene Matrix", 'm');
    
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
extern double *MxMethodeFULL (InData *InData){

  int a, i, j;
  double *U, biggest;
  
  U = (double *) MxNew (dim*dim*sizeof(double));
  free(D);
  biggest = (double)InData[0].rate;
  
  for(a = 0; a <= in_nr; a++){
    i = InData[a].i;
    j = InData[a].j;
    U[dim*i+j] = InData[a].rate;
    if(InData[a].rate >= biggest) biggest = InData[a].rate;
  }

  for(a = 0; a < dim*dim; a++) U[a] /= (biggest);

  if(opt.absrb){ /*==== absorbing  states ====*/
    for(i = 0; i < dim; i++)
      U[dim*i+(opt.absrb-1)] = 0. ;
  } /*== end absorbing states ==*/
  
  /* set diagonal elements  to 0 */
  for (i = 0; i < dim; i++) U[dim*i+i] = 0;
  for (j = 0; j < dim; j++) {
    double tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++) tmp += U[dim*i+j];
     U[dim*j+j] = -tmp+1; /* make U a stochastic matrix */
  }

  MxPrint (U, "U with Methode F", 'm');
  return U;
}

/*==*/
double *MxMethodeINPUT (TypeBarData *Data, double *Input){
  
  int i, j, real_abs = 0;
  double *U, Zabs, abs_rate;

   MxPrint(Input, "Input Matrix", 'm');

  if (opt.want_verbose) MxPrint(Input, "Input Matrix", 'm');

  if (opt.absrb) {  /*==== absorbing  states ====*/
    dim++;
    fprintf(stderr, "dim inceased to %i\n", dim);
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
    MxPrint(U, "aufgeblasene Matrix", 'm');
  }      /*== end absorbing states ==*/
  else{  /*== non-absorbing states ==*/
    U = (double *) MxNew(dim*dim*sizeof(double));
    for(i = 0; i < dim; i++)
      for(j = 0; j < dim; j++)
	U[dim*i+j] = Input[dim*i+j];
    MxPrint(U, "input Matrix before Verschoenerung", 'm');
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
  
  MxPrint (U,"U with Methode I" , 'm');
  if (opt.dumpU) MxBinWrite(U);

  free(Input);
  return U;
}

/*==*/
static double max_saddle(int i, int j,TypeBarData *Data){

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
static void *MxNew ( size_t size ) {

  void *mx;
  if ( (mx = (void *) calloc (1, size)) == NULL )
    fprintf (stderr, "ERROR: new_martix() allocation failed\n");

  return mx;
}

/*==*/
/* print matrix stored in ccmath-fromat */
static void MxPrint(double *mx, char *name, char T) {
  int k, l;
  switch (T) {
  case 'm':    /* square matrix */
    fprintf(stderr,"%s:\n", name);
    for (k = 0; k < dim; k++) {
      for (l=0; l< dim; l++) {
	fprintf(stderr,"%10.5f ", mx[dim*k+l]);
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"---\n");
    break;
  case 'v':
    fprintf(stderr,"%s:\n", name);    
    for (k = 0; k < dim; k++) fprintf(stderr,"%10.5f ", mx[k]);
    fprintf(stderr,"\n---\n");
    break;
  default:
    fprintf(stderr,"ERROR MxPrint(): no handler 4 type %c\n", T);
  }
}

/*==*/
static void print_settings(void) {

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
static char *time_stamp(void) {
  time_t  cal_time;
  
  cal_time = time(NULL);
  return ( ctime(&cal_time) );
}

/*==*/
void MxMemoryCleanUp (void) {
  if(_sqrPI)       free(_sqrPI);
  if(sqrPI_)       free(sqrPI_);
  if(opt.sequence) free(opt.sequence);
  fclose(opt.INFILE);
}

/*==*/
void MxDoDegeneracyStuff(void){

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
	    fprintf(stderr, "transition between %3d - %3d: adding %2d  cc\n",saddle[count].list[current], saddle[count].list[b], saddle[count].cc-1);
	    /* -1 because if the size of cc == 1 there's nothing special */
	    /* about it, i.e. there is just one saddle */
	  }
	  /* SECOND: we consider that the saddle connects several lmins */
	  if(nr > 2){
	    D[dim * (saddle[count].list[current]-1) + (saddle[count].list[b]-1)]++;
	    fprintf(stderr, "transition betweed %3d - %3d: adding  1 deg_saddle\n",saddle[count].list[current], saddle[count].list[b] );
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
static void MxMeschach2ccmath(MAT *meschach_matrix, double **origS){
  int  i, j;

  if((dim != meschach_matrix->m) != (dim != meschach_matrix->n)) {
    fprintf(stderr, "meschach-matrix is not square...\n");
    exit(1);}

  for (i = 0; i < dim; i++){
    for (j = 0; j < dim; j++){
      *(*origS+(dim*i+j)) = meschach_matrix->me[i][j];
    }
  }
 /*  M_FREE(meschach_matrix); */
}

/*==*/
static void MxMeschach2ccmathVec(VEC *meschach_vector, double **origEV){
  int i;
  for(i = 0; i < dim; i++) *(*origEV+i) = meschach_vector->ve[i];
 /*  V_FREE(meschach_vector); */
}

/*==*/
static MAT *Mxccmath2Meschach(double* ccmath_matrix){
  int i,j;
  MAT *tmp;
  tmp = m_get(dim,dim);
  for (i = 0; i < dim; i++){
    for (j = 0; j < dim; j++){
      tmp->me[i][j] = ccmath_matrix[dim*i+j];
    }
  }
  return tmp;
}

/*==*/
/* print matrix stored in meschach-format */
static void MxPrintMeschachMat(MAT *matrix, char *name) {
  int i,j;  
    fprintf(stderr, "%s (meschach):\n", name);
    for (i = 0; i < dim; i++){
      for (j = 0; j < dim; j++)
	fprintf(stderr,"%10.5f ", matrix->me[i][j]);
      fprintf(stderr, "\n");
    }
     fprintf(stderr,"---\n");
}

/*==*/
/* print vector stored in meschach-format */
static void MxPrintMeschachVec(VEC* vector, char *name){
  int i;
  fprintf(stderr, "%s (meschach):\n", name);
  for(i = 0; i < dim; i++) fprintf(stderr,"%10.5f ", vector->ve[i]);
  fprintf(stderr,"\n---\n");
}

/*==*/
void MxEVnonsymMx(double *origU, double **_S){

  int check = 1;
  double *tmp, *tmp_vec;
  MAT *A, *T, *Q, *X_re, *X_im;
  VEC *evals_re, *evals_im;
 
  tmp =     (double *) MxNew (dim*dim*sizeof(double));
  tmp_vec = (double *) MxNew (dim * sizeof(double));

  A = Mxccmath2Meschach(origU); /* transform U into meschach-format */
  /*    MxPrint(origU, "we'll calculate Eigenvectors and values of this matrix", 'm'); */
  Q = m_get(A->m, A->n);
  T = m_copy(A, MNULL);
  schur(T, Q);   /* compute Schur form: A = Q*T*Q^T */
  /*extract eigenvalues */
  evals_re = v_get(A->m);
  evals_im = v_get(A->m);
  schur_evals(T, evals_re, evals_im);
  /* Q not needed for eigenvalues */
  X_re = m_get(A->m, A->n);
  X_im = m_get(A->m, A->n);
  schur_vecs(T, Q, X_re, X_im);
  /* k'th eigenvector is k'th column of (X_re + i*X_im) */
  M_FREE(T);
  M_FREE(Q);
  M_FREE(A);
  if((check = Mxempty(X_im)) != 1)
    MxPrintMeschachMat(X_im, "imag Eigenvectors NOT EMPTY !!!");
  if(opt.want_verbose){
    MxPrintMeschachMat(X_re, "real Eigenvectors");
    MxPrintMeschachVec(evals_re, "real Eigenvalues");
  }
  MxMeschach2ccmath(X_re, &tmp);
  *_S = tmp;
  M_FREE(X_re);
  M_FREE(X_im);
  MxMeschach2ccmathVec(evals_re, &tmp_vec);
  EV_mesch = tmp_vec;
  V_FREE(evals_re);
  V_FREE(evals_im);
}

/*==*/
static int Mxempty(MAT *matrix){
  int i,j;

  for(i = 0; i < dim; i++)
    for(j = 0; j < dim; j++)
      if(matrix->me[i][j] != 0.) return 0;
  return 1;

}

/*==*/
static void MxBinWrite(double *matrix){
  int i, j;
  FILE *BINOUT;
  char *binfile = "matrixU.bin";

  BINOUT = fopen(binfile, "w");
  if (!BINOUT){
    fprintf(stderr, "could not open file pointer 4 binary outfile\n");
    exit(101);
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
void  MxExponent(double *p0, double *p8, double *U){
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
	exit(866);
      }
      printf("%e ", fabs(pt[i]));
      check += fabs(pt[i]); 
    }
    printf("\n");
  
    if ( ((check-1) < -0.01) || ((check-1) > 0.01) ){
      fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
      exit(888);
    }
    memset(pt,   0, dim*sizeof(double));
    memset(pdiff, 0, dim*sizeof(double));
    memset(Uexp, 0, dim*dim*sizeof(double));
    memset(U, 0, dim*dim*sizeof(double));

  }
  free(Uexp);
  free(pt);
}


/* End of file */
