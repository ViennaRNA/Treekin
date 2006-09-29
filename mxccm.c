/*=================================================================*/
/*=   mxccm.c                                                     =*/
/*=   matrix routines from ccmath library                         =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-09-29 13:46:17 mtw>          =*/
/*=   $Id: mxccm.c,v 1.5 2006/09/29 16:35:20 mtw Exp $            =*/
/*=   ---------------------------------------------------------   =*/
/*=      (c) Daniel A. Atkinson, Michael Thomas Wolfinger         =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#include <stdlib.h>
#include <math.h>
#include "mxccm.h"

void
eigen(double *a,double *ev,int n)
{
  double *dp;
  dp=(double *)calloc(n,sizeof(double));
  housev(a,ev,dp,n);
  qrevec(ev,a,dp,n); trnm(a,n);
  free(dp);
}

int
qrevec(double *ev,double *evec,double *dp,int n)
{
  double cc,sc=0.,d,x,y,h,tzr=1.e-15;
  int i,j,k,m,mqr=8*n;
  double *p;
  for(j=0,m=n-1;;++j){
    while(1){ if(m<1) return 0; k=m-1;
      if(fabs(dp[k])<=fabs(ev[m])*tzr) --m;
      else{ x=(ev[k]-ev[m])/2.; h=sqrt(x*x+dp[k]*dp[k]);
        if(m>1 && fabs(dp[m-2])>fabs(ev[k])*tzr) break;
        if((cc=sqrt((1.+x/h)/2.))!=0.) sc=dp[k]/(2.*cc*h); else sc=1.;
        x+=ev[m]; ev[m--]=x-h; ev[m--]=x+h;
        for(i=0,p=evec+n*(m+1); i<n ;++i,++p){
          h=p[0]; p[0]=cc*h+sc*p[n]; p[n]=cc*p[n]-sc*h;
         }
       }
     }
    if(j>mqr) return -1;
    if(x>0.) d=ev[m]+x-h; else d=ev[m]+x+h;
    cc=1.; y=0.; ev[0]-=d;
    for(k=0; k<m ;++k){
      x=ev[k]*cc-y; y=dp[k]*cc; h=sqrt(x*x+dp[k]*dp[k]);
      if(k>0) dp[k-1]=sc*h;
      ev[k]=cc*h; cc=x/h; sc=dp[k]/h; ev[k+1]-=d; y*=sc;
      ev[k]=cc*(ev[k]+y)+ev[k+1]*sc*sc+d;
      for(i=0,p=evec+n*k; i<n ;++i,++p){
        h=p[0]; p[0]=cc*h+sc*p[n]; p[n]=cc*p[n]-sc*h;
       }
     }
    ev[k]=ev[k]*cc-y; dp[k-1]=ev[k]*sc; ev[k]=ev[k]*cc+d;
   }
  return 0;
}

void
housev(double *a,double *d,double *dp,int n)
{
  double sc,x,y,h;
  int i,j,k,m,e;
  double *qw,*qs,*pc,*p;
  qs=(double *)calloc(n,sizeof(double));
  for(j=0,pc=a; j<n-2 ;++j,pc+=n+1){
    m=n-j-1;
    for(i=1,sc=0.; i<=m ;++i) sc+=pc[i]*pc[i];
    if(sc>0.){ sc=sqrt(sc);
      if((x= *(pc+1))<0.){ y=x-sc; h=1./sqrt(-2.*sc*y);}
      else{ y=x+sc; h=1./sqrt(2.*sc*y); sc= -sc;}
      for(i=0,qw=pc+1; i<m ;++i){
        qs[i]=0.; if(i) qw[i]*=h; else qw[i]=y*h;
       }
      for(i=0,e=j+2,p=pc+n+1,h=0.; i<m ;++i,p+=e++){
        qs[i]+=(y=qw[i])* *p++;
        for(k=i+1; k<m ;++k){
          qs[i]+=qw[k]* *p; qs[k]+=y* *p++;
         }
        h+=y*qs[i];
       }
      for(i=0; i<m ;++i){
        qs[i]-=h*qw[i]; qs[i]+=qs[i];
       }
      for(i=0,e=j+2,p=pc+n+1; i<m ;++i,p+=e++){
        for(k=i; k<m ;++k) *p++ -=qw[i]*qs[k]+qs[i]*qw[k];
       }
     }
    d[j]= *pc; dp[j]=sc;
   }
  d[j]= *pc; dp[j]= *(pc+1); d[j+1]= *(pc+=n+1);
  free(qs);
  for(i=0,m=n+n,p=pc; i<m ;++i) *p-- =0.;
  *pc=1.; *(pc-=n+1)=1.; qw=pc-n;
  for(m=2; m<n ;++m,qw-=n+1){
    for(j=0,p=pc,*pc=1.; j<m ;++j,p+=n){
      for(i=0,qs=p,h=0.; i<m ;) h+=qw[i++]* *qs++;
      for(i=0,qs=p,h+=h; i<m ;) *qs++ -=h*qw[i++];
     }
    for(i=0,p=qw+m; i<n ;++i) *(--p)=0.;
    *(pc-=n+1)=1.;
   }
}

void
mmul(double *c,double *a,double *b,int n)
{
  double *p,*q,s; int i,j,k;
  trnm(b,n);
  for(i=0; i<n ;++i,a+=n){
    for(j=0,q=b; j<n ;++j){
      for(k=0,p=a,s=0.; k<n ;++k) s+= *p++ * *q++;
      *c++ =s;
     }
   }
  trnm(b,n);
}

void
trnm(double *a,int n)
{
  double s,*p,*q;
  int i,j,e;
  for(i=0,e=n-1; i<n-1 ;++i,--e,a+=n+1){
    for(p=a+1,q=a+n,j=0; j<e ;++j){
      s= *p;
      *p = *q;
      p++;
      *q=s;
      q+=n;
     }
   }
}

void
vmul(double *vp,double *mat,double *v,int n)
{
  double s,*q; int k,i;
  for(k=0; k<n ;++k){
    for(i=0,q=v,s=0.; i<n ;++i) s+= *mat++ * *q++;
    *vp++ =s;
   }
}

void
mcopy(double *a,double *b,int m)
{
  double *p,*q; int k;
  for(p=a,q=b,k=0; k<m ;++k) *p++ = *q++;
}

int
minv(double *a,int n)
{
  int lc,*le; double s,t,tq=0.,zr=1.e-15;
  double *pa,*pd,*ps,*p,*q,*q0;
  int i,j,k,m;
  le=(int *)malloc(n*sizeof(int));
  q0=(double *)malloc(n*sizeof(double));
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

int
solv(double *a,double *b,int n)
{
  int i,j,k,lc; double *ps,*p,*q,*pa,*pd;
  double *q0,s,t,tq=0.,zr=1.e-15;
  q0=(double *)calloc(n,sizeof(double));
  for(j=0,pa=a,pd=a; j<n ;++j,++pa,pd+=n+1){
    if(j){
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
    tq=tq>s?tq:s; if(s<zr*tq){ free(q0); return -1;}
    if(lc!=j){ t=b[j]; b[j]=b[lc]; b[lc]=t;
      for(k=0,p=a+n*j,q=a+n*lc; k<n ;++k){
        t= *p; *p++ = *q; *q++ =t;
      }
    }
    for(k=j+1,ps=pd,t=1./ *pd; k<n ;++k) *(ps+=n)*=t;
  }
  for(j=1,ps=b+1; j<n ;++j){
    for(k=0,p=a+n*j,q=b,t=0.; k<j ;++k) t+= *p++ * *q++;
    *ps++ -=t;
  }
  for(j=n-1,--ps,pd=a+n*n-1; j>=0 ;--j,pd-=n+1){
    for(k=j+1,p=pd,q=b+j,t=0.; k<n ;++k) t+= *++p * *++q;
    *ps-=t; *ps-- /= *pd;
  }
  free(q0); return 0;
}
