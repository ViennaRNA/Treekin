/*=================================================================*/
/*=   exp_matrix.c                                                =*/
/*=   routines for calculating matrix exponents via pade approx.  =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 12:52:09 mtw>          =*/
/*=   $Id: exp_matrix.c,v 1.4 2006/03/15 11:52:30 mtw Exp $       =*/
/*=   ---------------------------------------------------------   =*/
/*=     (c) W. Andreas Svrcek-Seiler, Michael Thomas Wolfinger    =*/
/*=                  {svrci,mtw}@tbi.univie.ac.at                 =*/
/*=                             treekin                           =*/
/*=================================================================*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "exp_matrix.h"

/* private function(s) with long double */

static void
trnm(long double *a,int n)
{
  long double s,*p,*q;
  int i,j,e;
  for(i=0,e=n-1; i<n-1 ;++i,--e,a+=n+1){
    for(p=a+1,q=a+n,j=0; j<e ;++j){
      s= *p; *p++ = *q; *q=s; q+=n;
     }
   }
}

static void
mmul(long double *c,long double *a,long double *b,int n)
{
  long double *p,*q,s; int i,j,k;
  trnm(b,n);
  for(i=0; i<n ;++i,a+=n){
    for(j=0,q=b; j<n ;++j){
      for(k=0,p=a,s=0.; k<n ;++k) s+= *p++ * *q++;
      *c++ =s;
     }
   }
  trnm(b,n);
}


static int
minv(long double *a,int n)
{
  int lc,*le; long double s,t,tq=0.,zr=1.e-50;
  long double *pa,*pd,*ps,*p,*q,*q0;
  int i,j,k,m;
  le=(int *)malloc(n*sizeof(int));
  q0=(long double *)malloc(n*sizeof(long double));
  for(j=0,pa=pd=a; j<n ;++j,++pa,pd+=n+1){
    if(j>0){
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
      for(i=1; i<n ;++i){ lc=i<j?i:j;
        for(k=0,p=pa+i*n-j,q=q0,t=0.; k<lc ;++k) t+= *p++ * *q++;
      	q0[i]-=t;
       }
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
     }
    s=fabs(*pd); lc=j;
    for(k=j+1,ps=pd; k<n ;++k){
      if((t=fabs(*(ps+=n)))>s){ s=t; lc=k;}
     }
    tq=tq>s?tq:s; if(s<zr*tq){ free(le-j); free(q0); return -1;}
    *le++ =lc;
    if(lc!=j){
      for(k=0,p=a+n*j,q=a+n*lc; k<n ;++k){
        t= *p; *p++ = *q; *q++ =t;
       }
     }
    for(k=j+1,ps=pd,t=1./ *pd; k<n ;++k) *(ps+=n)*=t;
    *pd=t;
   }
  for(j=1,pd=ps=a; j<n ;++j){
    for(k=0,pd+=n+1,q= ++ps; k<j ;++k,q+=n) *q*= *pd;
   }
  for(j=1,pa=a; j<n ;++j){ ++pa;
    for(i=0,q=q0,p=pa; i<j ;++i,p+=n) *q++ = *p;
    for(k=0; k<j ;++k){ t=0.;
      for(i=k,p=pa+k*n+k-j,q=q0+k; i<j ;++i) t-= *p++ * *q++;
      q0[k]=t;
     }
    for(i=0,q=q0,p=pa; i<j ;++i,p+=n) *p= *q++;
   }
  for(j=n-2,pd=pa=a+n*n-1; j>=0 ;--j){ --pa; pd-=n+1;
    for(i=0,m=n-j-1,q=q0,p=pd+n; i<m ;++i,p+=n) *q++ = *p;
    for(k=n-1,ps=pa; k>j ;--k,ps-=n){ t= -(*ps);
      for(i=j+1,p=ps,q=q0; i<k ;++i) t-= *++p * *q++;
      q0[--m]=t;
     }
    for(i=0,m=n-j-1,q=q0,p=pd+n; i<m ;++i,p+=n) *p= *q++;
   }
  for(k=0,pa=a; k<n-1 ;++k,++pa){
    for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
    for(j=0,ps=a; j<n ;++j,ps+=n){
      if(j>k){ t=0.; p=ps+j; i=j;}
      else{ t=q0[j]; p=ps+k+1; i=k+1;}
      for(; i<n ;) t+= *p++ *q0[i++];
      q0[j]=t;
     }
    for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
   }
  for(j=n-2,le--; j>=0 ;--j){
    for(k=0,p=a+j,q=a+ *(--le); k<n ;++k,p+=n,q+=n){
      t=*p; *p=*q; *q=t;
     }
   }
  free(le); free(q0);
  return 0;
}


static void
ipmmul(long double *a,long double *b, long double *c,int n)
{
  long double *T1,*T2;
  
  T1 = (long double *) malloc(n*n*sizeof(long double));
  T2 =  (long double *) malloc(n*n*sizeof(long double));
  
  memcpy(T1,a,sizeof(long double)*n*n);
  memcpy(T2,b,sizeof(long double)*n*n);
  
  mmul(T1,T2,c,n);
  memcpy(a,T1,sizeof(long double)*n*n);
  free(T2);
  free(T1);
  return;
}

void
padexp(double *from,double *out,int n,int ord)
{
  /* calculate exp(M) by a Pade approx */
  /* of order ord                      */
  /* exp(x) ~= p(x)/q(x)               */
  /* Sidje's trick makes               */
  /* exp(x) ~= 1 + P(x)/Q(x)           */
  
  int i,k;
  long double *cp;
  long double *mat_in,*MD,*MN,*X2;
  /*precon stuff */
  long double maxnorm,ml;
  int ppow=0;
  long  double scfac;
  
  
  mat_in = (long double *) malloc(n*n*sizeof(long double));
  for (i=0;i<n*n;i++) mat_in[i] = (long double) from[i];
  
  /* get the infinity norm of the matrix for scaling/squaring*/
  maxnorm=-1;
  
  for (i=0;i<n*n;i++)
    maxnorm = maxnorm < fabs(mat_in[i]) ? fabs(mat_in[i]) : maxnorm ;
  
  ml =  log(maxnorm)/log(2.);
  
  /* yuck, heavy heuristics */
  /* calculate M/(2^n) to get exp (M) = (exp (M/2^n))^(2^n) afterwards */
  
  if (ml >=4.)  {
    ppow = (int) (ml - 2.);
    
    scfac=1.;
    for (i=1;i<=ppow;i++) scfac *=2.;
    for (i=0;i<n*n;i++) mat_in[i]/=scfac;
  }
  /* cp holds the Pade coefficients ,
     X2 the squared input matrix    ,
     MN the numerator polynomial    ,
     MD the denominator polyniomial */
  cp = (long double *) malloc((ord+1)*sizeof(long double));
  X2 =  (long double *) malloc(n*n*sizeof(long double));
  MN =  (long double *) malloc(n*n*sizeof(long double));
  MD = (long double *) malloc(n*n*sizeof(long double));
  
  cp[0]=1.;
  for (i=1;i<=ord;i++)
    cp[i] = cp[i-1]*((long double)(ord+1-i))/((long double)((ord+ord+1-i)*i));
  
  ipmmul(X2,mat_in,mat_in,n);
  
  for (i=0;i<n*n;i++) MN[i]=MD[i]=0.;
  
  for (i=0;i<n;i++){
    MN[(n+1)*i]=cp[ord-1];
    MD[(n+1)*i]=cp[ord];
  }
  
  for (k=ord/2-1;k>0;k--){
    ipmmul(MN,X2,MN,n);
    for (i=0;i<n;i++) MN[(n+1)*i]+=cp[2*k-1];
  }
  
  for (k=ord/2;k>0;k--){
    ipmmul(MD,X2,MD,n);
    for (i=0;i<n;i++) MD[(n+1)*i] +=cp[2*k-2];
  }
  ipmmul(MN,mat_in,MN,n);
  
  for(i=0;i<n*n;i++) MD[i]=MD[i]-MN[i];
  minv(MD,n);
  ipmmul(mat_in,MD,MN,n);
  for(i=0;i<n*n;i++) mat_in[i]*=2.0;
  for(i=0;i<n;i++) mat_in[(n+1)*i] +=1.0;
  
  if (ml >=6.){
    for (i=1;i<=ppow; i++)
      ipmmul(mat_in,mat_in,mat_in,n);
  }
  
  for (i=0;i<n*n;i++) out[i] = (double) mat_in[i];
  free(MD);
  free(MN);
  free(X2);
  free(cp);
  free(mat_in);
}  
