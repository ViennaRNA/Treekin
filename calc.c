/* calc.c */
/* Last changed Time-stamp: <2003-08-27 18:42:32 mtw> */
/* static char rcsid[] = "$Id: calc.c,v 1.7 2003/08/27 17:08:26 mtw Exp $"; */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h> 
#include "mxccm.h"     /* functions for eigen-problems stolen from ccmath */
#include "barparser.h" /* functions for input */
#include "matrix.h"    /* basic matrix functions from meschach */
#include "matrix2.h"   /* advances matrix functions from meschach */
#include "calc.h"      /* does all matrix stuff for markov process */
#include "globals.h"   /* contains getopt-stuff */

#define TOL 0.000000000000001
#define ABS_VAL 0.00001

/* private function(s) */
static void   *MxNew (size_t size);
static double *MxMethodeA (TypeBarData *Data);
static double *MxMethodeB (TypeBarData *Data);
static double *MxMethodeC (TypeBarData *Data);
static void    MxPrint(double *mx, char *name, char T);
static void    MxPrintMeschachMat(MAT *matrix, char *name);
static void    MxPrintMeschachVec(VEC* vector, char *name);
static void    MxMeschach2ccmath(MAT *meschach_matrix, double **origM);
static void    MxMeschach2ccmathVec(VEC *meschach_vector, double **origV);
static MAT    *Mxccmath2Meschach(double* ccmath_matrix);
static double  max_saddle(int i, int j, TypeBarData *Data);
static int   **build_subtree_list(TypeBarData *Data);
static int     yet_seen(int *start, int number);
static void    calc_effective_bsize(TypeBarData *Data, int **P);
static void    recalculate_free_energy(TypeBarData *Data, int **lmin, double kT);
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

  _kT = 0.00198717*(273.15 + opt.T);
  if (d > 0 ) dim = d;
  else { fprintf(stderr, "dim <= 0\n"); exit(1); }

  EV = (double *) MxNew (dim*sizeof(double));
  _sqrPI = (double *) MxNew (dim*dim*sizeof(double));
  sqrPI_ = (double *) MxNew (dim*dim*sizeof(double));
}

/*==*/
double *MxBar2Matrix ( TypeBarData *Data) {

  double *U;

  if(opt.want_degenerate) MxDoDegeneracyStuff();
  
  switch (opt.method) {
  case 'A':
    U = MxMethodeA(Data);
    break;
  case 'B':
    U = MxMethodeB(Data);
    break;
  case 'C':
    U = MxMethodeC(Data);
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

  if (opt.want_verbose) {
    sprintf(Aname, "%s", "p0");
    MxPrint (p0, Aname, 'v');
  }
  return (p0);
}

/*==*/
/* calculate equilibrium distribution */
/* originally, we used Data[i].energy here, but now we use FGr to */
/* get the same equilibrium distribution as in the full proess */
double *MxEqDistr ( TypeBarData *Data ) {

  int i;
  double *p8, Z = 0.;
  
  p8 = (double *) MxNew (dim*sizeof(double));

  for(i = 0; i < dim; i++) Z += exp(-((double)Data[i].FGr/_kT));
  for(i = 0; i < dim; i++) p8[i] = exp(-((double) Data[i].FGr/_kT))/Z;

  if (opt.want_verbose) MxPrint (p8, "p8", 'v');

  if(opt.absrb > 0){
    double tmp = 0.;
    for(i = 0; i < dim; i++){
      p8[i] = ABS_VAL;
      tmp += p8[i];
    }
    p8[opt.absrb-1] = 1.0-tmp;
    if(opt.want_verbose) MxPrint (p8, "p8 for absorbing state", 'v'); 
  }
  MxPrint(p8, "p8", 'v');
  return (p8);
}

/*==*/
/* calculate equilibrium distribution from the energi-array */ 
double *MxEqDistrFULL (double *energi) {
  int i;
  double *p8;
  double Z = 0., test = 0.;
  
  p8 = (double *) MxNew (dim*sizeof(double));

  for(i = 0; i < dim; i++) Z += exp(-((double)energi[i]/_kT));
  for(i = 0; i < dim; i++) p8[i] = exp(-((double) energi[i]/_kT))/Z;

  /* this one is for check reasons only */
  for(i = 0; i < dim; i++) test += p8[i];
  /* end check */

  if (opt.want_verbose) {
    sprintf (Aname, "%s", "p8");
    MxPrint (p8, Aname, 'v');
  }

  return (p8);
}

/*==*/
double *MxSymmetr ( double *U, double *P8 ) {

  int i, j;
  double *S;
  double *tmpMx;
  
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
  if (opt.want_verbose)
    MxPrint (S, "Eigenvectors of S", 'm');MxPrint(EV, "Eigenvalues of S", 'v'); 
  
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
    p(t)   = sqrPI_ * S * exp(time * EV) * St * _sqrPI * p(0)
    CL     = sqrPI_ * S
    CR     = St * _sqrPI
    tmpVec = CR * p(0)
    ******* ABSORBING CASE: *******
    p(t)   = S * exp(time * EV) * S_inv * p(0) 
    tmpMx  = S * exp(time * EV) 
    tmpVec = S_inv * p(0)
    p(t)   = tmpMx * tmpVec 
  */
  int i, j, count = 0, pdiff_counter = 0;
  double time, check = 0.;
  double *CL, *CR, *exptL, *tmpMx, *tmpVec, *St, *S_inv; /* transposed of S */
  double *pt, *pdiff;  /* probability distribution/difference 4 time t */
  
  St     = (double *) MxNew (dim*dim*sizeof(double));
  CL     = (double *) MxNew (dim*dim*sizeof(double));
  CR     = (double *) MxNew (dim*dim*sizeof(double));
  exptL  = (double *) MxNew (dim*dim*sizeof(double));
  tmpMx  = (double *) MxNew (dim*dim*sizeof(double));
  S_inv  = (double *) MxNew (dim*dim*sizeof(double));
  tmpVec = (double *) MxNew (dim*sizeof(double));
  pt     = (double *) MxNew (dim*sizeof(double));
  pdiff  = (double *) MxNew (dim*sizeof(double));
  
  if(! opt.absrb){  /* NON-absorbing case */
    mcopy(St, S, dim*dim);  trnm(St, dim); /* transpose S */
    mmul (CL, sqrPI_, S, dim);
    mmul (CR, St, _sqrPI, dim);
    vmul (tmpVec, CR, p0, dim);
    free(St);
    free(CR);
    free(S_inv);
  }
  else{  /* absorbing case */
    mcopy(S_inv, S, dim*dim);
    minv(S_inv,dim);
    for(i = 0; i < dim; i++) EV_mesch[i] -= 1;  /* compensate 4 translation of matrix U */
    EV = EV_mesch;    /* let EV point at EV_mesch */
    MxPrint(EV_mesch, "EV_mesch in MxIterate", 'v');
    vmul (tmpVec, S_inv, p0, dim);
    free(S_inv);
  }
  
  /*** solve fundamental equation ***/  /* using logarithmic time-scale */
  print_settings();
  for (time = opt.t0; time <= opt.t8; time *= opt.tinc) {
    for (i = 0; i < dim; i++) {
      for (j = 0 ; j < dim; j++) {
	if ( i == j) exptL[dim*i+j] = exp(time*EV[i]);
	else exptL[dim*i+j] = 0.;
      }
    }
    if(!opt.absrb) mmul (tmpMx, CL, exptL, dim);
    else mmul (tmpMx, S, exptL, dim);
    vmul (pt, tmpMx, tmpVec, dim);
    count++;  /* # of iterations */
    
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

    if ( ((check-1) < -0.05) || ((check-1) > 0.05) ){
      fprintf(stderr, "overall probability at time %e is %e != 1. ! exiting\n", time,check );
      exit(888);
    }
    check = 0.;
    
    /* now check if we have converged yet */
    /*  printf("#---------------"); */
    for(i=0; i<dim; i++){
      pdiff[i] = p8[i] - pt[i];
      /*   printf("%7.4f ", pdiff[i]); */
      if (fabs(pdiff[i]) >= 0.0001)
	pdiff_counter++; /* # of lmins whose pdiff is > the threshold */
    }
    /*  printf(" # pdiff_counter: %i", pdiff_counter); */
    /*   printf("\n"); */
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
  free(tmpMx);
  free(pt);
  free(pdiff);
  free(p0);
/*    free(EV_mesch); */
}


/*==*/
void MxIterate_FULL (double *p0, double *p8, double *S,  int *assoc_gradbas, int lmin_nr) {
  /*
    solve following equation 4 various times
    p(t) = sqrPI_ * S * exp(time * EV) * St * _sqrPI * p(0)
    CL = sqrPI_ * S
    CR = St * _sqrPI
    tmpVec = CR * p(0)
  */
  int i, j, lmins, pdiff_counter = 0;
  double time, *CL, *CR, *exptL, *tmpMx, *tmpVec, *pt, *St;
  double *ptFULL;    /* prob dist 4 of the effective lmins of the tree at time t */
  double *p8FULL;    /* euqilibrium distr 4 gradient basins, calculated from full process */
  double *pdiffFULL; /* population prob difference between p8FULL and ptFULL */ 
  double check = 0.;
  double checkp8 = 0.;
  
  lmins = lmin_nr; /* # of gradient basins == # of lmins in tree */

  St        = (double *) MxNew (dim*dim*sizeof(double));
  CL        = (double *) MxNew (dim*dim*sizeof(double));
  CR        = (double *) MxNew (dim*dim*sizeof(double));
  tmpVec    = (double *) MxNew (dim*sizeof(double));
  exptL     = (double *) MxNew (dim*dim*sizeof(double));
  tmpMx     = (double *) MxNew (dim*dim*sizeof(double));
  pt        = (double *) MxNew (dim*sizeof(double));
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
    p8FULL[assoc_gradbas[i]] += p8[i]; /* eq distr of the gradient basins */
  for (i = 0; i < lmins; i++) checkp8 += fabs(p8FULL[i]);
  if ( ((checkp8-1) < -0.1) || ((checkp8-1) > 0.1) ){
    fprintf(stderr, "overall equilibrium probability is %e != 1. ! exiting\n", checkp8);
    exit(888);
  }
  
   /* solve fundamental equation */
  for (time = opt.t0; time <= opt.t8; time *= opt.tinc) {
    for (i = 0; i < dim; i++) {
      for (j = 0 ; j < dim; j++) {
	if ( i == j) exptL[dim*i+j] = exp(time*EV[i]);
	else exptL[dim*i+j] = 0.;
      }
    }
    
    mmul (tmpMx, CL, exptL, dim);
    vmul (pt, tmpMx, tmpVec, dim);
        
    for (i = 0; i < dim; i++){
      if(pt[i] < -0.001){
	fprintf(stderr, "prob of lmin %i has become negative: %6.4f\n", i+1,pt[i]);exit(866);
      }
      ptFULL[assoc_gradbas[i]] += pt[i]; /* map individual structure -> gradient basins */
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
      if (fabs(pdiffFULL[i]) >= 0.0001)
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

  int i, j;
  double *U;

  U = (double *) MxNew (dim*dim*sizeof(double));
  
  for (i = 0; i < dim; i++) {
    if (Data[i].father > 0) {  /* merged local minima */
      j = Data[i].father - 1;
      if(i != j) {
	/* rate i -> j */
	U[dim*j+i] = exp(-Data[i].ediff/_kT);                         
        /* rate j -> i */
        U[dim*i+j] = exp(-(Data[i].ediff+Data[i].energy-Data[j].energy)/_kT);
        continue;
      }
    }
    if(Data[i].father == 0) {  /* unmerged local minima */
      j = 0;
      if(i == j) continue;
      /* rate i -> j */
      U[dim*j+i] = exp(-Data[i].ediff/_kT);
      /* rate j -> i */
      U[dim*i+j] = exp(-(Data[i].ediff+Data[i].energy-Data[j].energy)/_kT);
    } 
  }
  
  /* diagonal elements */
  for (j = 0; j < dim; j++) {
    double tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++){ 
      tmp += U[dim*i+j];
    }
    /* make U a stochastic matrix */
    U[dim*j+j] = -tmp+1;
  }
  
  if (opt.want_verbose) {
    sprintf (Aname, "%s", "U with Methode A");
    MxPrint (U, Aname, 'm');
  }
  
  return (U);
}

/*==*/
static double *MxMethodeB (TypeBarData *Data) {
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
  
  int i,j;
  double m_saddle;
  double *U;

  U = (double *) MxNew (dim*dim*sizeof(double));

  for(i = 0; i < (dim*dim); i++) U[i] = -1;
  
  if(opt.want_degenerate) {
    for( i = 0; i < dim; i++)
      for( j = i+1; j < dim; j++){
	m_saddle = max_saddle(i, j, Data);
	/* rate j -> i */
	U[dim*i+j] = D[dim*i+j]*exp(-(m_saddle-Data[j].FGr)/_kT); 
	/* rate i -> j */
	U[dim*j+i] = D[dim*i+j]*exp(-(m_saddle-Data[i].FGr)/_kT);
      }
    if(D != NULL) free(D);
  }
  else {
    for( i = 0; i < dim; i++)
      for( j = i+1; j < dim; j++){
	m_saddle = max_saddle(i, j, Data);
	/* rate j -> i */
	U[dim*i+j] = exp(-(m_saddle-Data[j].FGr)/_kT);
	/* rate i -> j */
	U[dim*j+i] = exp(-(m_saddle-Data[i].FGr)/_kT);
      }
  }
  /*===========================*/
  /*==== absorbing  states ====*/
  /*===========================*/
  if(opt.absrb > 0)
    for(i = 0; i < dim; i++)  /* make column opt.absrb absorbing by */
      U[dim*i+(opt.absrb-1)] *= exp(-10/_kT); 
  /*==========================*/
  /*== end absorbing states ==*/
  /*==========================*/
  
  for(i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      if( i == j) U[dim*i+j] = 0;  /* set diagonal elements to 0 */
  
  /* diagonal elements */
  for (j = 0; j < dim; j++) {
    double tmp = 0.00;
    /* calculate colum sum */
    for(i = 0; i < dim; i++) tmp += U[dim*i+j];
    U[dim*j+j] = -tmp+1; /* make U a stochastic matrix */
  }
  
  if (opt.want_verbose) MxPrint (U,"U with Methode B", 'm');
  
  return (U);
}

/*==*/
static double *MxMethodeC (TypeBarData *Data){

 int i,j;
 int **P; 
 double m_saddle; 
 double *U;

 U = (double *) MxNew (dim*dim*sizeof(double));
 
 P =  build_subtree_list(Data);   /* info: **P shows the number of subtrees :-) */

 calc_effective_bsize(Data, P);
 recalculate_free_energy(Data, P, _kT);

 for(i = 0; i < (dim*dim); i++) U[i] = -1;
 
 for( i = 0; i < dim; i++)
   for( j = i+1; j < dim; j++){
      m_saddle = max_saddle(i, j, Data);
      /* rate j -> i */
      U[dim*i+j] = (Data[i].eff_bsize*Data[j].eff_bsize)*exp(-(m_saddle-Data[j].F_eff)/_kT);
      /* rate i -> j */
      U[dim*j+i] = U[dim*i+j] * exp((Data[i].energy-Data[j].energy)/_kT); /* back-rate from detailed balance */ 
     /*  U[dim*j+i] = exp(-(m_saddle-Data[i].energy)/_kT); */ /* exactly the same */ 
   }
 
 for(i = 0; i < dim; i++)
   for (j = 0; j < dim; j++)
     if( i == j) U[dim*i+j] = 0;  /* set diagonal elements to 0 */
 
 /* diagonal elements */
 for (j = 0; j < dim; j++) {
   double tmp = 0.00;
   /* calculate column sum */
   for(i = 0; i < dim; i++){ 
     tmp += U[dim*i+j];
   }
   /* make U a stochastic matrix */
   U[dim*j+j] = -tmp+1;
 }
 
 if (opt.want_verbose) {
   sprintf (Aname, "%s", "U with Methode C");
   MxPrint (U, Aname, 'm');
 }

 free(P);
 return (U);
}

/*==*/ 
extern double *MxMethodeFULL (InData *InData){

  int a, i, j;
  extern int in_nr;
  double *Q;
  float biggest;
  
  Q = (double *) MxNew (dim*dim*sizeof(double));

  biggest = InData[0].rate;
  
  for(a = 0; a <= in_nr; a++){
    i = InData[a].i;
    j = InData[a].j;
    Q[dim*i+j] = InData[a].rate;
    if(InData[a].rate >= biggest) biggest = InData[a].rate;
  }

  for(a = 0; a < dim*dim; a++)
    Q[a] /= biggest;
  
  /* diagonal elements */
  for (j = 0; j < dim; j++) {
    double tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++){ 
      tmp += Q[dim*i+j];
    }
    /* make Q a stochastic matrix */
    Q[dim*j+j] = -tmp+1;
  }

  if (opt.want_verbose) {
    sprintf (Aname, "%s", "U with Methode F");
    MxPrint (Q, Aname, 'm');
  }
  return Q;
}

/*==*/
extern double *MxMethodeINPUT (double *Input){
  int i,j;
  double *test;

  test = (double *)calloc(dim,sizeof(double));
  
  if (opt.want_verbose)
    MxPrint(Input, "Input Matrix", 'm'); 

  for(i=0; i<dim*dim; i++)
    Input[i] /= 100;

  /* diagonal elements */
  for (i=0; i<dim; i++) Input[dim*i+i] = 0;
  for (j = 0; j < dim; j++) {
    double tmp = 0.00;
    /* calculate column sum */
    for(i = 0; i < dim; i++)
      tmp += Input[dim*i+j];
    /* make Q a stochastic matrix */
    Input[dim*j+j] = -tmp+1.;
  }
  
  if (opt.want_verbose)
    MxPrint (Input,"U with Methode I" , 'm');

  free(test);
  return Input;
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
static int **build_subtree_list(TypeBarData *Data){
  
  int i,j, test = -1, a,b = 0;
  int **deepest;
  int *seen;
  
  deepest = (int **)calloc((dim+1),sizeof(int *));
  seen = (int *)calloc(dim, sizeof(int));
  
  for(i = 0; i <= dim; i++){        /* dim+1 rows */
    *(deepest+i) = (int *)calloc((dim), sizeof(int));
    for(j = 0; j < dim; j++)        /* dim colums */ 
      *(*(deepest+i)+j) = test;
  }

  **deepest = 0;

  i = j = 0;
  for(i = 1; i <= dim; i++){   /* the one to merge with */ 
    j = 0;                     /* # of lmin in this subtree */
    for(a = 1; a < dim; a++){
      if(yet_seen(seen, Data[a].number)) continue;
      if(Data[a].father == i){
	*(*(deepest+i)+j) = Data[a].number;
	seen[b++] = Data[a].number;
	j++;
      }
    }
    if( **(deepest+i) > 0   )(**deepest)++;
  }
  
  free(seen);
  return deepest;  
}

/*==*/
/* checks if number is already in the array and returns 1 if yes, 0 if not */
static int yet_seen(int *start, int number){
  
  int i;
  
  for(i = 0; i < dim; i++)
    if(start[i] == number) return 1;

  return 0;
}

/*==*/
/* calculates the effective bsize of each lmin without its children */ 
static void calc_effective_bsize(TypeBarData *Data, int **lmin){

  int i,j, test = 0;

  for(i = 0; i < dim; i++)
    Data[i].eff_bsize = Data[i].bsize;
    
  for(i = 1; i <= dim; i++){ /* look over all lmins */
    j = 0;
    while(*(*(lmin+i)+j) != -1){
      Data[i].eff_bsize -= Data[(*(*(lmin+i)+j))-1].bsize;
      j++;
    }
  }

  /* check if we still have the same # of structures in sum */
  for(i = 0; i < dim; i++)
    test += Data[i].eff_bsize;
  if ( test != Data[0].bsize){
    fprintf (stderr, "bsize has changed in calc_effective_bsize \n");
    exit(1);
  }
}

/*==*/
static void recalculate_free_energy(TypeBarData *Data, int **lmin, double kT){

  int i,j;
  float mfe;

   mfe = Data[0].energy;
   /*    printf("mfe: %7.4f\n", mfe); */

   for(i = 0; i < dim; i++){
     Data[i].F_eff = Data[i].F;
     Data[i].Z = exp((-Data[i].F)/kT);  /* calculate partition function Z of lmin */
     Data[i].Z_eff = Data[i].Z;
     /*printf("lmin %2d  F: %7.4f Z: %7.4f\n", Data[i].number, Data[i].F, Data[i].Z);*/
   }
  
   for(i = 1; i <= dim; i++){ /* look over all lmins */
     j = 0;
     while(*(*(lmin+i)+j) != -1){
       Data[i-1].Z_eff -= Data[(*(*(lmin+i)+j))-1].Z;
       j++;
     }
     Data[i-1].F_eff = -kT*log(Data[i-1].Z_eff);
     /*   printf("F_eff: %7.4f, Z_eff: %7.4f \n",Data[i-1].F_eff, Data[i-1].Z_eff); */
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

  int i, j,  numsad, count = 0, b, nr, current;
  numsad = 1;

  D = (double *) MxNew (dim*dim*sizeof(double));
  for(i = 0; i < (dim*dim); i++) D[i] = 1; /* initialize D */
  
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
  MxPrintMeschachMat(X_re, "real Eigenvectors");
  if((check = Mxempty(X_im)) != 1)
    MxPrintMeschachMat(X_im, "imag Eigenvectors NOT EMPTY !!!");
  MxPrintMeschachVec(evals_re, "real Eigenvalues");
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

/* End of file */

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

  fclose(BINOUT);
}


