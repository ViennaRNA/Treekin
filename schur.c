/*=================================================================*/
/*=   schur.c                                                     =*/
/*=   routines for computing the Schur decomposition of a real    =*/
/*=   non-symmetric matrix from meschach library                  =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 11:00:52 mtw>          =*/
/*=   $Id: schur.c,v 1.2 2006/03/15 11:08:15 mtw Exp $    =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

/**************************************************************************
**
** Copyright (C) 1993 David E. Stewart & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

#include	<stdio.h>
#include	"matrix.h"
#include	<math.h>


static void
hhldr3(double x, double y, double z, double *nu1, double *beta, double *newval)
{
  double	alpha;  
  if ( x >= 0.0 )
    alpha = sqrt(x*x+y*y+z*z);
  else
    alpha = -sqrt(x*x+y*y+z*z);
  *nu1 = x + alpha;
  *beta = 1.0/(alpha*(*nu1));
  *newval = alpha;
}

static void
hhldr3cols(MAT *A, int k, int j0, double beta,double nu1, double nu2, double nu3)
{
  double	**A_me, ip, prod;
  int	j, n;
  
  if ( k < 0 || k+3 > A->m || j0 < 0 )
    error(E_BOUNDS,"hhldr3cols");
  A_me = A->me;
  n = A->n;
  
  for ( j = j0; j < n; j++ ) {
    ip = nu1*m_entry(A,k,j)+nu2*m_entry(A,k+1,j)+nu3*m_entry(A,k+2,j);
    prod = ip*beta;
    m_add_val(A,k  ,j,-prod*nu1);
    m_add_val(A,k+1,j,-prod*nu2);
    m_add_val(A,k+2,j,-prod*nu3);
  }
}

static void
hhldr3rows(MAT *A, int k, int i0, double beta,double nu1, double nu2, double nu3)
{
  double	**A_me, ip, prod;
  int	i, m;
  
  if ( k < 0 || k+3 > A->n )
    error(E_BOUNDS,"hhldr3rows");
  A_me = A->me;		m = A->m;
  i0 = min(i0,m-1);
  
  for ( i = 0; i <= i0; i++ ){
    ip = nu1*m_entry(A,i,k)+nu2*m_entry(A,i,k+1)+nu3*m_entry(A,i,k+2);
    prod = ip*beta;
    m_add_val(A,i,k  , - prod*nu1);
    m_add_val(A,i,k+1, - prod*nu2);
    m_add_val(A,i,k+2, - prod*nu3);
  }
}

/* schur -- computes the Schur decomposition of the matrix A in situ
	-- optionally, gives Q matrix such that Q^T.A.Q is upper triangular
	-- returns upper triangular Schur matrix */
MAT*
schur(MAT *A, MAT *Q)
{
  int		i, j, iter, k, k_min, k_max, k_tmp, n, split;
  double	beta2, c, discrim, dummy, nu1, s, tmp, x, y, z;
  double	**A_me;
  double	sqrt_macheps;
  static	VEC	*diag=VNULL, *beta=VNULL;
  
  if ( ! A )
    error(E_NULL,"schur");
  if ( A->m != A->n || ( Q && Q->m != Q->n ) )
    error(E_SQUARE,"schur");
  if ( Q != MNULL && Q->m != A->m )
    error(E_SIZES,"schur");
  n = A->n;
  diag = v_resize(diag,A->n);
  beta = v_resize(beta,A->n);
  /*MEM_STAT_REG(diag,TYPE_VEC);
    MEM_STAT_REG(beta,TYPE_VEC); */
  /* compute Hessenberg form */
  Hfactor(A,diag,beta);
  
  /* save Q if necessary */
  if ( Q )
    Q = makeHQ(A,diag,beta,Q);
  V_FREE(diag); /* added by mtw */
  V_FREE(beta); /* added by mtw */
  makeH(A,A);
  
  sqrt_macheps = sqrt(MACHEPS);
  
  k_min = 0;	A_me = A->me;
  
  while ( k_min < n ){
    double    a00, a01, a10, a11;
    double  scale, t, numer, denom;
    
    /* find k_max to suit:
       submatrix k_min..k_max should be irreducible */
    k_max = n-1;
    for ( k = k_min; k < k_max; k++ )
      if ( m_entry(A,k+1,k) == 0.0 ){
	k_max = k;
	break;
      }
    
    if ( k_max <= k_min ) {
      k_min = k_max + 1;
      continue;	/* outer loop */
    }
    
    /* check to see if we have a 2 x 2 block
       with complex eigenvalues */
    if ( k_max == k_min + 1 )
      {
	a00 = m_entry(A,k_min,k_min);
	a01 = m_entry(A,k_min,k_max);
	a10 = m_entry(A,k_max,k_min);
	a11 = m_entry(A,k_max,k_max);
	tmp = a00 - a11;
	discrim = tmp*tmp + 4*a01*a10;
	if ( discrim < 0.0 )
	  {	/* yes -- e-vals are complex
		   -- put 2 x 2 block in form [a b; c a];
		   then eigenvalues have real part a & imag part sqrt(|bc|) */
	    numer = - tmp;
	    denom = ( a01+a10 >= 0.0 ) ?
	      (a01+a10) + sqrt((a01+a10)*(a01+a10)+tmp*tmp) :
	      (a01+a10) - sqrt((a01+a10)*(a01+a10)+tmp*tmp);
	    if ( denom != 0.0 ) {   /* t = s/c = numer/denom */
	      t = numer/denom;
	      scale = c = 1.0/sqrt(1+t*t);
	      s = c*t;
	    }
	    else {
	      c = 1.0;
	      s = 0.0;
	    }
	    rot_cols(A,k_min,k_max,c,s,A);
	    rot_rows(A,k_min,k_max,c,s,A);
	    if ( Q != MNULL )
	      rot_cols(Q,k_min,k_max,c,s,Q);
	    k_min = k_max + 1;
	    continue;
	  }
	else /* discrim >= 0; i.e. block has two real eigenvalues */
	  {	/* no -- e-vals are not complex;
		   split 2 x 2 block and continue */
		/* s/c = numer/denom */
	    numer = ( tmp >= 0.0 ) ?
	      - tmp - sqrt(discrim) : - tmp + sqrt(discrim);
	    denom = 2*a01;
	    if ( fabs(numer) < fabs(denom) ) {   /* t = s/c = numer/denom */
	      t = numer/denom;
	      scale = c = 1.0/sqrt(1+t*t);
	      s = c*t;
	    }
	    else if ( numer != 0.0 ) {   /* t = c/s = denom/numer */
	      t = denom/numer;
	      scale = 1.0/sqrt(1+t*t);
	      c = fabs(t)*scale;
	      s = ( t >= 0.0 ) ? scale : -scale;
	    }
	    else { /* numer == denom == 0 */
	      c = 0.0;
	      s = 1.0;
	    }
	    rot_cols(A,k_min,k_max,c,s,A);
	    rot_rows(A,k_min,k_max,c,s,A);
	    if ( Q != MNULL )
	      rot_cols(Q,k_min,k_max,c,s,Q);
	    k_min = k_max + 1;	/* go to next block */
	    continue;
	  }
      }
    
    /* now have r x r block with r >= 2:
       apply Francis QR step until block splits */
    split = FALSE;
    iter = 0;
    while ( ! split ) {
      iter++;
      
      /* set up Wilkinson/Francis complex shift */
      k_tmp = k_max - 1;
      
      a00 = m_entry(A,k_tmp,k_tmp);
      a01 = m_entry(A,k_tmp,k_max);
      a10 = m_entry(A,k_max,k_tmp);
      a11 = m_entry(A,k_max,k_max);
      
      /* treat degenerate cases differently
	 -- if there are still no splits after five iterations
	 and the bottom 2 x 2 looks degenerate, force it to
	 split */
      if ( iter >= 5 &&
	   fabs(a00-a11) < sqrt_macheps*(fabs(a00)+fabs(a11)) &&
	   (fabs(a01) < sqrt_macheps*(fabs(a00)+fabs(a11)) ||
	    fabs(a10) < sqrt_macheps*(fabs(a00)+fabs(a11))) )
	{
	  if ( fabs(a01) < sqrt_macheps*(fabs(a00)+fabs(a11)) )
	    m_set_val(A,k_tmp,k_max,0.0);
	  if ( fabs(a10) < sqrt_macheps*(fabs(a00)+fabs(a11)) )
	    {
	      m_set_val(A,k_max,k_tmp,0.0);
	      split = TRUE;
	      continue;
	    }
	}
      
      s = a00 + a11;
      t = a00*a11 - a01*a10;
      
      /* break loop if a 2 x 2 complex block */
      if ( k_max == k_min + 1 && s*s < 4.0*t ) {
	split = TRUE;
	continue;
      }
      
      /* perturb shift if convergence is slow */
      if ( (iter % 10) == 0 ){
	s += iter*0.02;
	t += iter*0.02;
      }
      
      /* set up Householder transformations */
      k_tmp = k_min + 1;
      a00 = m_entry(A,k_min,k_min);
      a01 = m_entry(A,k_min,k_tmp);
      a10 = m_entry(A,k_tmp,k_min);
      a11 = m_entry(A,k_tmp,k_tmp);
      x = a00*a00 + a01*a10 - s*a00 + t;
      y = a10*(a00+a11-s);
      if ( k_min + 2 <= k_max )
	z = a10* /* m_entry(A,k_min+2,k_tmp) */ A->me[k_min+2][k_tmp];
      else
	z = 0.0;
      
      for ( k = k_min; k <= k_max-1; k++ ){
	if ( k < k_max - 1 ){
	  hhldr3(x,y,z,&nu1,&beta2,&dummy);
	  tracecatch(hhldr3cols(A,k,max(k-1,0),  beta2,nu1,y,z),"schur");
	  tracecatch(hhldr3rows(A,k,min(n-1,k+3),beta2,nu1,y,z),"schur");
	  if ( Q != MNULL )
	    hhldr3rows(Q,k,n-1,beta2,nu1,y,z);
	}
	else
	  {
	    givens(x,y,&c,&s);
	    rot_cols(A,k,k+1,c,s,A);
	    rot_rows(A,k,k+1,c,s,A);
	    if ( Q )
	      rot_cols(Q,k,k+1,c,s,Q);
	  }
	x = m_entry(A,k+1,k);
	if ( k <= k_max - 2 )
	  y = m_entry(A,k+2,k);
	else
	  y = 0.0;
	if ( k <= k_max - 3 )
	  z = m_entry(A,k+3,k);
	else
	  z = 0.0;
      }
      for ( k = k_min; k <= k_max-2; k++ ) {
	/* zero appropriate sub-diagonals */
	m_set_val(A,k+2,k,0.0);
	if ( k < k_max-2 )
	  m_set_val(A,k+3,k,0.0);
      }
      
      /* test to see if matrix should split */
      for ( k = k_min; k < k_max; k++ )
	if ( fabs(A_me[k+1][k]) < MACHEPS*
	     (fabs(A_me[k][k])+fabs(A_me[k+1][k+1])) )
	  {	A_me[k+1][k] = 0.0;	split = TRUE;	}
    }
  }
  
  /* polish up A by zeroing strictly lower triangular elements
     and small sub-diagonal elements */
  for ( i = 0; i < A->m; i++ )
    for ( j = 0; j < i-1; j++ )
      A_me[i][j] = 0.0;
  for ( i = 0; i < A->m - 1; i++ )
    if ( fabs(A_me[i+1][i]) < MACHEPS*
	 (fabs(A_me[i][i])+fabs(A_me[i+1][i+1])) )
      A_me[i+1][i] = 0.0;
  
  return A;
}

/* schur_vals -- compute real & imaginary parts of eigenvalues
   -- assumes T contains a block upper triangular matrix
   as produced by schur()
   -- real parts stored in real_pt, imaginary parts in imag_pt */
void
schur_evals(MAT *T, VEC *real_pt, VEC *imag_pt)
{
  int	i, n;
  double	discrim, **T_me;
  double	diff, sum, tmp;
  
  if ( ! T || ! real_pt || ! imag_pt )
    error(E_NULL,"schur_evals");
  if ( T->m != T->n )
    error(E_SQUARE,"schur_evals");
  n = T->n;	T_me = T->me;
  real_pt = v_resize(real_pt,(u_int)n);
  imag_pt = v_resize(imag_pt,(u_int)n);
  
  i = 0;
  while ( i < n ) {
    if ( i < n-1 && T_me[i+1][i] != 0.0 ) {   /* should be a complex eigenvalue */
      sum  = 0.5*(T_me[i][i]+T_me[i+1][i+1]);
      diff = 0.5*(T_me[i][i]-T_me[i+1][i+1]);
      discrim = diff*diff + T_me[i][i+1]*T_me[i+1][i];
      if ( discrim < 0.0 ) {	/* yes -- complex e-vals */
	real_pt->ve[i] = real_pt->ve[i+1] = sum;
	imag_pt->ve[i] = sqrt(-discrim);
	imag_pt->ve[i+1] = - imag_pt->ve[i];
      }
      else {	/* no -- actually both real */
	tmp = sqrt(discrim);
	real_pt->ve[i]   = sum + tmp;
	real_pt->ve[i+1] = sum - tmp;
	imag_pt->ve[i]   = imag_pt->ve[i+1] = 0.0;
      }
      i += 2;
    }
    else {   /* real eigenvalue */
      real_pt->ve[i] = T_me[i][i];
      imag_pt->ve[i] = 0.0;
      i++;
    }
  }
}

/* schur_vecs -- returns eigenvectors computed from the real Schur
		decomposition of a matrix
	-- T is the block upper triangular Schur matrix
	-- Q is the orthognal matrix where A = Q.T.Q^T
	-- if Q is null, the eigenvectors of T are returned
	-- X_re is the real part of the matrix of eigenvectors,
		and X_im is the imaginary part of the matrix.
	-- X_re is returned */
MAT*
schur_vecs(MAT *T, MAT *Q, MAT *X_re, MAT *X_im)
{
  int	i, j, limit;
  double	t11_re, t11_im, t12, t21, t22_re, t22_im;
  double	l_re, l_im, det_re, det_im, invdet_re, invdet_im,
    val1_re, val1_im, val2_re, val2_im,
    tmp_val1_re, tmp_val1_im, tmp_val2_re, tmp_val2_im, **T_me;
  double	sum, diff, discrim, magdet, norm, scale;
  static VEC	*tmp1_re=VNULL, *tmp1_im=VNULL, *tmp2_re=VNULL, *tmp2_im=VNULL;

  if ( ! T || ! X_re )
    error(E_NULL,"schur_vecs");
  if ( T->m != T->n || X_re->m != X_re->n ||
       ( Q != MNULL && Q->m != Q->n ) ||
       ( X_im != MNULL && X_im->m != X_im->n ) )
    error(E_SQUARE,"schur_vecs");
  if ( T->m != X_re->m ||
       ( Q != MNULL && T->m != Q->m ) ||
       ( X_im != MNULL && T->m != X_im->m ) )
    error(E_SIZES,"schur_vecs");
  
  tmp1_re = v_resize(tmp1_re,T->m);
  tmp1_im = v_resize(tmp1_im,T->m);
  tmp2_re = v_resize(tmp2_re,T->m);
  tmp2_im = v_resize(tmp2_im,T->m);
  /* MEM_STAT_REG(tmp1_re,TYPE_VEC);
  MEM_STAT_REG(tmp1_im,TYPE_VEC);
  MEM_STAT_REG(tmp2_re,TYPE_VEC);
  MEM_STAT_REG(tmp2_im,TYPE_VEC);*/
  
  T_me = T->me;
  i = 0;
  while ( i < T->m ) {
    if ( i+1 < T->m && T->me[i+1][i] != 0.0 ) {	/* complex eigenvalue */
      sum  = 0.5*(T_me[i][i]+T_me[i+1][i+1]);
      diff = 0.5*(T_me[i][i]-T_me[i+1][i+1]);
      discrim = diff*diff + T_me[i][i+1]*T_me[i+1][i];
      l_re = l_im = 0.0;
      if ( discrim < 0.0 ) {	/* yes -- complex e-vals */
	l_re = sum;
	l_im = sqrt(-discrim);
      }
      else /* not correct Real Schur form */
	error(E_RANGE,"schur_vecs");
    }
    else {
      l_re = T_me[i][i];
      l_im = 0.0;
    }
    
    v_zero(tmp1_im);
    v_rand(tmp1_re);
    sv_mlt(MACHEPS,tmp1_re,tmp1_re);

    /* solve (T-l.I)x = tmp1 */
    limit = ( l_im != 0.0 ) ? i+1 : i;
    /* printf("limit = %d\n",limit); */
    for ( j = limit+1; j < T->m; j++ )
      tmp1_re->ve[j] = 0.0;
    j = limit;
    while ( j >= 0 ) {
      if ( j > 0 && T->me[j][j-1] != 0.0 ) {   /* 2 x 2 diagonal block */
	val1_re = tmp1_re->ve[j-1] -
	  __ip__(&(tmp1_re->ve[j+1]),&(T->me[j-1][j+1]),limit-j);
	val1_im = tmp1_im->ve[j-1] -
	  __ip__(&(tmp1_im->ve[j+1]),&(T->me[j-1][j+1]),limit-j);
	val2_re = tmp1_re->ve[j] -
	  __ip__(&(tmp1_re->ve[j+1]),&(T->me[j][j+1]),limit-j);
	val2_im = tmp1_im->ve[j] -
	  __ip__(&(tmp1_im->ve[j+1]),&(T->me[j][j+1]),limit-j);
	
	t11_re = T_me[j-1][j-1] - l_re;
	t11_im = - l_im;
	t22_re = T_me[j][j] - l_re;
	t22_im = - l_im;
	t12 = T_me[j-1][j];
	t21 = T_me[j][j-1];
	
	scale =  fabs(T_me[j-1][j-1]) + fabs(T_me[j][j]) +
	  fabs(t12) + fabs(t21) + fabs(l_re) + fabs(l_im);
	
	det_re = t11_re*t22_re - t11_im*t22_im - t12*t21;
	det_im = t11_re*t22_im + t11_im*t22_re;
	magdet = det_re*det_re+det_im*det_im;
	if ( sqrt(magdet) < MACHEPS*scale ){
	  det_re = MACHEPS*scale;
	  magdet = det_re*det_re+det_im*det_im;
	}
	invdet_re =   det_re/magdet;
	invdet_im = - det_im/magdet;
	tmp_val1_re = t22_re*val1_re-t22_im*val1_im-t12*val2_re;
	tmp_val1_im = t22_im*val1_re+t22_re*val1_im-t12*val2_im;
	tmp_val2_re = t11_re*val2_re-t11_im*val2_im-t21*val1_re;
	tmp_val2_im = t11_im*val2_re+t11_re*val2_im-t21*val1_im;
	tmp1_re->ve[j-1] = invdet_re*tmp_val1_re -
	  invdet_im*tmp_val1_im;
	tmp1_im->ve[j-1] = invdet_im*tmp_val1_re +
	  invdet_re*tmp_val1_im;
	tmp1_re->ve[j]   = invdet_re*tmp_val2_re -
	  invdet_im*tmp_val2_im;
	tmp1_im->ve[j]   = invdet_im*tmp_val2_re +
	  invdet_re*tmp_val2_im;
	j -= 2;
      }
      else {
	t11_re = T_me[j][j] - l_re;
	t11_im = - l_im;
	magdet = t11_re*t11_re + t11_im*t11_im;
	scale = fabs(T_me[j][j]) + fabs(l_re);
	if ( sqrt(magdet) < MACHEPS*scale ){
	  t11_re = MACHEPS*scale;
	  magdet = t11_re*t11_re + t11_im*t11_im;
	}
	invdet_re =   t11_re/magdet;
	invdet_im = - t11_im/magdet;
	val1_re = tmp1_re->ve[j] -
	  __ip__(&(tmp1_re->ve[j+1]),&(T->me[j][j+1]),limit-j);
	val1_im = tmp1_im->ve[j] -
	  __ip__(&(tmp1_im->ve[j+1]),&(T->me[j][j+1]),limit-j);
	tmp1_re->ve[j] = invdet_re*val1_re - invdet_im*val1_im;
	tmp1_im->ve[j] = invdet_im*val1_re + invdet_re*val1_im;
	j -= 1;
      }
    }
    
    norm = v_norm_inf(tmp1_re) + v_norm_inf(tmp1_im);
    sv_mlt(1/norm,tmp1_re,tmp1_re);
    if ( l_im != 0.0 )
      sv_mlt(1/norm,tmp1_im,tmp1_im);
    mv_mlt(Q,tmp1_re,tmp2_re);
    if ( l_im != 0.0 )
      mv_mlt(Q,tmp1_im,tmp2_im);
    if ( l_im != 0.0 )
      norm = sqrt(in_prod(tmp2_re,tmp2_re)+in_prod(tmp2_im,tmp2_im));
    else
      norm = v_norm2(tmp2_re);
    sv_mlt(1/norm,tmp2_re,tmp2_re);
    if ( l_im != 0.0 )
      sv_mlt(1/norm,tmp2_im,tmp2_im);
    
    if ( l_im != 0.0 ) {
      if ( ! X_im )
	error(E_NULL,"schur_vecs");
      set_col(X_re,i,tmp2_re);
      set_col(X_im,i,tmp2_im);
      sv_mlt(-1.0,tmp2_im,tmp2_im);
      set_col(X_re,i+1,tmp2_re);
      set_col(X_im,i+1,tmp2_im);
      i += 2;
    }
    else {
      set_col(X_re,i,tmp2_re);
      if ( X_im != MNULL )
	set_col(X_im,i,tmp1_im);	/* zero vector */
      i += 1;
    }
  }
  V_FREE(tmp1_re); /* added by mtw */
  V_FREE(tmp1_im); /* added by mtw */
  V_FREE(tmp2_re); /* added by mtw */
  V_FREE(tmp2_im); /* added by mtw */
  
  return X_re;
}



/* hhvec -- calulates Householder vector to eliminate all entries after the
	i0 entry of the vector vec. It is returned as out. May be in-situ */
VEC*
hhvec(VEC *vec, u_int i0, Real *beta, VEC *out, Real *newval)
{
  Real	norm;
  
  out = _v_copy(vec,out,i0);
  norm = sqrt(_in_prod(out,out,i0));
  if ( norm <= 0.0 ) {
    *beta = 0.0;
    return (out);
  }
  *beta = 1.0/(norm * (norm+fabs(out->ve[i0])));
  if ( out->ve[i0] > 0.0 )
    *newval = -norm;
  else
    *newval = norm;
  out->ve[i0] -= *newval;
  
  return (out);
}

/* hhtrvec -- apply Householder transformation to vector -- may be in-situ */
VEC*
hhtrvec(VEC *hh, double beta, u_int i0, VEC *in, VEC *out)
{
  /* hh = Householder vector */
  Real	scale;
  if ( hh==(VEC *)NULL || in==(VEC *)NULL )
    error(E_NULL,"hhtrvec");
  if ( in->dim != hh->dim )
    error(E_SIZES,"hhtrvec");
  if ( i0 > in->dim )
    error(E_BOUNDS,"hhtrvec");
  
  scale = beta*_in_prod(hh,in,i0);
  out = v_copy(in,out);
  __mltadd__(&(out->ve[i0]),&(hh->ve[i0]),-scale,(int)(in->dim-i0));
  return (out);
}

/* hhtrrows -- transform a matrix by a Householder vector by rows
	starting at row i0 from column j0 -- in-situ */
MAT*
hhtrrows(MAT *M, u_int i0, u_int j0, VEC *hh, double beta)
{
  Real	ip, scale;
  int	i /*, j */;
  
  if ( M==(MAT *)NULL || hh==(VEC *)NULL )
    error(E_NULL,"hhtrrows");
  if ( M->n != hh->dim )
    error(E_RANGE,"hhtrrows");
  if ( i0 > M->m || j0 > M->n )
    error(E_BOUNDS,"hhtrrows");
  
  if ( beta == 0.0 )	return (M);
  
  /* for each row ... */
  for ( i = i0; i < M->m; i++ )	{	/* compute inner product */
    ip = __ip__(&(M->me[i][j0]),&(hh->ve[j0]),(int)(M->n-j0));
    scale = beta*ip;
    if ( scale == 0.0 )
      continue;
    
    /* do operation */
    __mltadd__(&(M->me[i][j0]),&(hh->ve[j0]),-scale, (int)(M->n-j0));
  }
  
  return (M);
}


/* hhtrcols -- transform a matrix by a Householder vector by columns
	starting at row i0 from column j0 -- in-situ */
MAT*
hhtrcols(MAT *M, u_int i0, u_int j0, VEC *hh, double beta)
{
  int	i;
  static	VEC	*w = VNULL;
  
  if ( M==(MAT *)NULL || hh==(VEC *)NULL )
    error(E_NULL,"hhtrcols");
  if ( M->m != hh->dim )
    error(E_SIZES,"hhtrcols");
  if ( i0 > M->m || j0 > M->n )
    error(E_BOUNDS,"hhtrcols");
  
  if ( beta == 0.0 )	return (M);
  
  w = v_resize(w,M->n);
  /*  MEM_STAT_REG(w,TYPE_VEC);*/
  v_zero(w);
  
  for ( i = i0; i < M->m; i++ )
    if ( hh->ve[i] != 0.0 )
      __mltadd__(&(w->ve[j0]),&(M->me[i][j0]),hh->ve[i],
		 (int)(M->n-j0));
  for ( i = i0; i < M->m; i++ )
    if ( hh->ve[i] != 0.0 )
      __mltadd__(&(M->me[i][j0]),&(w->ve[j0]),-beta*hh->ve[i],
		 (int)(M->n-j0));
  V_FREE(w);  /* added by mtw */
  return (M);
}


/* givens -- returns c,s parameters for Givens rotation to
		eliminate y in the vector [ x y ]' */
void
givens(x,y,c,s)
double  x,y;
Real	*c,*s;
{
	Real	norm;

	norm = sqrt(x*x+y*y);
	if ( norm == 0.0 )
	{	*c = 1.0;	*s = 0.0;	}	/* identity */
	else
	{	*c = x/norm;	*s = y/norm;	}
}

/* rot_rows -- premultiply mat by givens rotation described by c,s */
MAT	*rot_rows(mat,i,k,c,s,out)
MAT	*mat,*out;
u_int	i,k;
double	c,s;
{
	u_int	j;
	Real	temp;

	if ( mat==(MAT *)NULL )
		error(E_NULL,"rot_rows");
	if ( i >= mat->m || k >= mat->m )
		error(E_RANGE,"rot_rows");
	if ( mat != out )
		out = m_copy(mat,m_resize(out,mat->m,mat->n));

	for ( j=0; j<mat->n; j++ )
	{
		/* temp = c*out->me[i][j] + s*out->me[k][j]; */
		temp = c*m_entry(out,i,j) + s*m_entry(out,k,j);
		/* out->me[k][j] = -s*out->me[i][j] + c*out->me[k][j]; */
		m_set_val(out,k,j, -s*m_entry(out,i,j) + c*m_entry(out,k,j));
		/* out->me[i][j] = temp; */
		m_set_val(out,i,j, temp);
	}

	return (out);
}

/* rot_cols -- postmultiply mat by givens rotation described by c,s */
MAT	*rot_cols(mat,i,k,c,s,out)
MAT	*mat,*out;
u_int	i,k;
double	c,s;
{
	u_int	j;
	Real	temp;

	if ( mat==(MAT *)NULL )
		error(E_NULL,"rot_cols");
	if ( i >= mat->n || k >= mat->n )
		error(E_RANGE,"rot_cols");
	if ( mat != out )
		out = m_copy(mat,m_resize(out,mat->m,mat->n));

	for ( j=0; j<mat->m; j++ )
	{
		/* temp = c*out->me[j][i] + s*out->me[j][k]; */
		temp = c*m_entry(out,j,i) + s*m_entry(out,j,k);
		/* out->me[j][k] = -s*out->me[j][i] + c*out->me[j][k]; */
		m_set_val(out,j,k, -s*m_entry(out,j,i) + c*m_entry(out,j,k));
		/* out->me[j][i] = temp; */
		m_set_val(out,j,i,temp);
	}

	return (out);
}


/* Hfactor -- compute Hessenberg factorisation in compact form.
	-- factorisation performed in situ
	-- for details of the compact form see QRfactor.c and matrix2.doc */
MAT*
Hfactor(MAT *A, VEC *diag, VEC *beta)
{
  static	VEC	*tmp1 = VNULL;
  int	k, limit;  
  if ( ! A || ! diag || ! beta )
    error(E_NULL,"Hfactor");
  if ( diag->dim < A->m - 1 || beta->dim < A->m - 1 )
    error(E_SIZES,"Hfactor");
  if ( A->m != A->n )
    error(E_SQUARE,"Hfactor");
  limit = A->m - 1;
  
  tmp1 = v_resize(tmp1,A->m);
  /*MEM_STAT_REG(tmp1,TYPE_VEC);*/
  
  for ( k = 0; k < limit; k++ ) {
    get_col(A,(u_int)k,tmp1);
    hhvec(tmp1,k+1,&beta->ve[k],tmp1,&A->me[k+1][k]);
    v_set_val(diag,k,v_entry(tmp1,k+1));
    hhtrcols(A,k+1,k+1,tmp1,v_entry(beta,k));
    hhtrrows(A,0  ,k+1,tmp1,v_entry(beta,k));
  }
  V_FREE(tmp1); /* added by mtw */
  return (A);
}

/* makeHQ -- construct the Hessenberg orthogonalising matrix Q;
	-- i.e. Hess M = Q.M.Q'	*/
MAT	*makeHQ(H, diag, beta, Qout)
MAT	*H, *Qout;
VEC	*diag, *beta;
{
	int	i, j, limit;
	static	VEC	*tmp1 = VNULL, *tmp2 = VNULL;

	if ( H==(MAT *)NULL || diag==(VEC *)NULL || beta==(VEC *)NULL )
		error(E_NULL,"makeHQ");
	limit = H->m - 1;
	if ( diag->dim < limit || beta->dim < limit )
		error(E_SIZES,"makeHQ");
	if ( H->m != H->n )
		error(E_SQUARE,"makeHQ");
	Qout = m_resize(Qout,H->m,H->m);

	tmp1 = v_resize(tmp1,H->m);
	tmp2 = v_resize(tmp2,H->m);
	/*	MEM_STAT_REG(tmp1,TYPE_VEC);
		MEM_STAT_REG(tmp2,TYPE_VEC);*/

	for ( i = 0; i < H->m; i++ )
	{
		/* tmp1 = i'th basis vector */
		for ( j = 0; j < H->m; j++ )
			/* tmp1->ve[j] = 0.0; */
		    v_set_val(tmp1,j,0.0);
		/* tmp1->ve[i] = 1.0; */
		v_set_val(tmp1,i,1.0);

		/* apply H/h transforms in reverse order */
		for ( j = limit-1; j >= 0; j-- )
		{
			get_col(H,(u_int)j,tmp2);
			/* tmp2->ve[j+1] = diag->ve[j]; */
			v_set_val(tmp2,j+1,v_entry(diag,j));
			hhtrvec(tmp2,beta->ve[j],j+1,tmp1,tmp1);
		}

		/* insert into Qout */
		set_col(Qout,(u_int)i,tmp1);
	}
	V_FREE(tmp1);  /* added by mtw */
	V_FREE(tmp2);  /* added by mtw */
	return (Qout);
}

/* makeH -- construct actual Hessenberg matrix */
MAT	*makeH(H,Hout)
MAT	*H, *Hout;
{
	int	i, j, limit;

	if ( H==(MAT *)NULL )
		error(E_NULL,"makeH");
	if ( H->m != H->n )
		error(E_SQUARE,"makeH");
	Hout = m_resize(Hout,H->m,H->m);
	Hout = m_copy(H,Hout);

	limit = H->m;
	for ( i = 1; i < limit; i++ )
		for ( j = 0; j < i-1; j++ )
			/* Hout->me[i][j] = 0.0;*/
		    m_set_val(Hout,i,j,0.0);

	return (Hout);
}



#define SQ(X) ((X)*(X))

/* _v_norm2 -- computes (scaled) 2-norm (Euclidean norm) of vectors */
double	_v_norm2(x,scale)
VEC	*x, *scale;
{
	int	i, dim;
	Real	s, sum;

	if ( x == (VEC *)NULL )
		error(E_NULL,"_v_norm2");
	dim = x->dim;

	sum = 0.0;
	if ( scale == (VEC *)NULL )
		for ( i = 0; i < dim; i++ )
			sum += SQ(x->ve[i]);
	else if ( scale->dim < dim )
		error(E_SIZES,"_v_norm2");
	else
		for ( i = 0; i < dim; i++ )
		{	s = scale->ve[i];
			sum += ( s== 0.0 ) ? SQ(x->ve[i]) :
							SQ(x->ve[i]/s);
		}

	return sqrt(sum);
}

#define	max(a,b)	((a) > (b) ? (a) : (b))

/* _v_norm_inf -- computes (scaled) infinity-norm (supremum norm) of vectors */
double	_v_norm_inf(x,scale)
VEC	*x, *scale;
{
	int	i, dim;
	Real	s, maxval, tmp;

	if ( x == (VEC *)NULL )
		error(E_NULL,"_v_norm_inf");
	dim = x->dim;

	maxval = 0.0;
	if ( scale == (VEC *)NULL )
		for ( i = 0; i < dim; i++ )
		{	tmp = fabs(x->ve[i]);
			maxval = max(maxval,tmp);
		}
	else if ( scale->dim < dim )
		error(E_SIZES,"_v_norm_inf");
	else
		for ( i = 0; i < dim; i++ )
		{	s = scale->ve[i];
			tmp = ( s== 0.0 ) ? fabs(x->ve[i]) : fabs(x->ve[i]/s);
			maxval = max(maxval,tmp);
		}

	return maxval;
}
