
#include	<limits.h>
#include 	"matrix.h"

/* m_get -- gets an mxn matrix (in MAT form) by dynamic memory allocation */
MAT	*m_get(int m, int n) {
  MAT	*matrix;
  int	i;
  
  if (m < 0 || n < 0)
    error(E_NEG,"m_get");
  
  if ((matrix=NEW(MAT)) == (MAT *)NULL )
    error(E_MEM,"m_get");
   
  matrix->m = m;
  matrix->n = matrix->max_n = n;
  matrix->max_m = m;
  matrix->max_size = m*n;

  if ((matrix->base = NEW_A(m*n,Real)) == (Real *)NULL ) {
    free(matrix);
    error(E_MEM,"m_get");
  }
  if ((matrix->me = (Real **)calloc(m,sizeof(Real *))) == (Real **)NULL ) {
    free(matrix->base);
    free(matrix);
    error(E_MEM,"m_get");
  }
  
  /* set up pointers */
  for ( i=0; i<m; i++ )
    matrix->me[i] = &(matrix->base[i*n]);
  return (matrix);
}


/* v_get -- gets a VEC of dimension 'dim'
   -- Note: initialized to zero */
VEC	*v_get(int size) {
  VEC	*vector;
  
  if (size < 0)
    error(E_NEG,"v_get");
  
  if ((vector=NEW(VEC)) == (VEC *)NULL )
    error(E_MEM,"v_get");
  
  vector->dim = vector->max_dim = size;
  if ((vector->ve=NEW_A(size,Real)) == (Real *)NULL ) {
    free(vector);
    error(E_MEM,"v_get");
  }
  
  return (vector);
}

/* m_free -- returns MAT & asoociated memory back to memory heap */
int	m_free(MAT *mat) {
   
  if ( mat==(MAT *)NULL || (int)(mat->m) < 0 ||
       (int)(mat->n) < 0 )
    /* don't trust it */
    return (-1);
  
  if ( mat->base != (Real *)NULL ) {
    free((char *)(mat->base));
  }
  if ( mat->me != (Real **)NULL ) {
    free((char *)(mat->me));
  }
  
  free((char *)mat);
  
  return (0);
}

/* v_free -- returns VEC & asoociated memory back to memory heap */
int	v_free(VEC *vec){
  if ( vec==(VEC *)NULL || (int)(vec->dim) < 0 )
    /* don't trust it */
    return (-1);
  
  if ( vec->ve == (Real *)NULL ) {
    free((char *)vec);
  }
  else {
    free((char *)vec->ve);
    free((char *)vec);
  }
  
  return (0);
}



/* m_resize -- returns the matrix A of size new_m x new_n; A is zeroed
   -- if A == NULL on entry then the effect is equivalent to m_get() */
MAT	*m_resize(MAT *A, int new_m, int new_n) {
  int	i;
  int	new_max_m, new_max_n, new_size, old_m, old_n;
  
  if (new_m < 0 || new_n < 0)
    error(E_NEG,"m_resize");
  
  if ( ! A )
    return m_get(new_m,new_n);
  
  /* nothing was changed */
  if (new_m == A->m && new_n == A->n)
    return A;
  
  old_m = A->m;	old_n = A->n;
  if ( new_m > A->max_m ) {	/* re-allocate A->me */
    
    A->me = RENEW(A->me,new_m,Real *);
    if ( ! A->me )
      error(E_MEM,"m_resize");
  }
  new_max_m = max(new_m,A->max_m);
  new_max_n = max(new_n,A->max_n);
  
   new_size = new_max_m*new_max_n;
   if ( new_size > A->max_size ){	/* re-allocate A->base */
     
     A->base = RENEW(A->base,new_size,Real);
     if ( ! A->base )
       error(E_MEM,"m_resize");
     A->max_size = new_size;
   }
   
   /* now set up A->me[i] */
   for ( i = 0; i < new_m; i++ )
     A->me[i] = &(A->base[i*new_n]);
   
   /* now shift data in matrix */
   if ( old_n > new_n ) {
     for ( i = 1; i < min(old_m,new_m); i++ )
       MEM_COPY((char *)&(A->base[i*old_n]),
		(char *)&(A->base[i*new_n]),
		sizeof(Real)*new_n);
   }
   else if ( old_n < new_n ) {
     for ( i = (int)(min(old_m,new_m))-1; i > 0; i-- )
       {   /* copy & then zero extra space */
	 MEM_COPY((char *)&(A->base[i*old_n]),
		  (char *)&(A->base[i*new_n]),
		  sizeof(Real)*old_n);
	 __zero__(&(A->base[i*new_n+old_n]),(new_n-old_n));
       }
     __zero__(&(A->base[old_n]),(new_n-old_n));
     A->max_n = new_n;
   }
   /* zero out the new rows.. */
   for ( i = old_m; i < new_m; i++ )
     __zero__(&(A->base[i*new_n]),new_n);
   
   A->max_m = new_max_m;
   A->max_n = new_max_n;
   A->max_size = A->max_m*A->max_n;
   A->m = new_m;
   A->n = new_n;
   
   return A;
}

/* v_resize -- returns the vector x with dim new_dim
   -- x is set to the zero vector */
VEC	*v_resize(VEC *x, int new_dim) {
  
  if (new_dim < 0)
    error(E_NEG,"v_resize");
  
  if ( ! x )
    return v_get(new_dim);
  
  /* nothing is changed */
  if (new_dim == x->dim)
    return x;
  
  if ( x->max_dim == 0 )	/* assume that it's from sub_vec */
    return v_get(new_dim);
  
  if ( new_dim > x->max_dim ) {
    
    x->ve = RENEW(x->ve,new_dim,Real);
    if ( ! x->ve )
      error(E_MEM,"v_resize");
    x->max_dim = new_dim;
  }
  
  if ( new_dim > x->dim )
    __zero__(&(x->ve[x->dim]),new_dim - x->dim);
  x->dim = new_dim;
  
  return x;
}


/* _in_prod -- inner product of two vectors from i0 downwards */
double	_in_prod(VEC *a, VEC *b, unsigned int i0) {
  u_int	limit;
  
  if ( a==(VEC *)NULL || b==(VEC *)NULL )
    error(E_NULL,"_in_prod");
  limit = min(a->dim,b->dim);
  if ( i0 > limit )
    error(E_BOUNDS,"_in_prod");
  
  return __ip__(&(a->ve[i0]),&(b->ve[i0]),(int)(limit-i0));
  /*****************************************
	a_v = &(a->ve[i0]);		b_v = &(b->ve[i0]);
	for ( i=i0; i<limit; i++ )
		sum += a_v[i]*b_v[i];
		sum += (*a_v++)*(*b_v++);

	return (double)sum;
  ******************************************/
}

/* sv_mlt -- scalar-vector multiply -- may be in-situ */
VEC	*sv_mlt(double scalar, VEC *vector, VEC *out) {
  if ( vector==(VEC *)NULL )
    error(E_NULL,"sv_mlt");
  if ( out==(VEC *)NULL || out->dim != vector->dim )
    out = v_resize(out,vector->dim);
  if ( scalar == 0.0 )
    return v_zero(out);
  if ( scalar == 1.0 )
    return v_copy(vector,out);
  
  __smlt__(vector->ve,(double)scalar,out->ve,(int)(vector->dim));
  /**************************************************
	dim = vector->dim;
	out_ve = out->ve;	vec_ve = vector->ve;
	for ( i=0; i<dim; i++ )
		out->ve[i] = scalar*vector->ve[i];
		(*out_ve++) = scalar*(*vec_ve++);
  **************************************************/
  return (out);
}


/* m_add -- matrix addition -- may be in-situ */
MAT	*m_add(MAT *mat1, MAT *mat2, MAT *out) {
  u_int	m,n,i;
  
  if ( mat1==(MAT *)NULL || mat2==(MAT *)NULL )
    error(E_NULL,"m_add");
  if ( mat1->m != mat2->m || mat1->n != mat2->n )
    error(E_SIZES,"m_add");
  if ( out==(MAT *)NULL || out->m != mat1->m || out->n != mat1->n )
    out = m_resize(out,mat1->m,mat1->n);
  m = mat1->m;	n = mat1->n;
  for ( i=0; i<m; i++ ) {
    __add__(mat1->me[i],mat2->me[i],out->me[i],(int)n);
    /**************************************************
		for ( j=0; j<n; j++ )
			out->me[i][j] = mat1->me[i][j]+mat2->me[i][j];
    **************************************************/
  }
  
  return (out);
}

/* mv_mlt -- matrix-vector multiplication 
   -- Note: b is treated as a column vector */
VEC	*mv_mlt(MAT *A, VEC *b, VEC *out){
  u_int	i, m, n;
  Real	**A_v, *b_v /*, *A_row */;
  
  if ( A==(MAT *)NULL || b==(VEC *)NULL )
    error(E_NULL,"mv_mlt");
  if ( A->n != b->dim )
    error(E_SIZES,"mv_mlt");
  if ( b == out )
    error(E_INSITU,"mv_mlt");
  if ( out == (VEC *)NULL || out->dim != A->m )
    out = v_resize(out,A->m);
  
  m = A->m;		n = A->n;
  A_v = A->me;		b_v = b->ve;
  for ( i=0; i<m; i++ ) {
    out->ve[i] = __ip__(A_v[i],b_v,(int)n);
    /**************************************************
		A_row = A_v[i];		b_v = b->ve;
		for ( j=0; j<n; j++ )
			sum += (*A_row++)*(*b_v++);
		out->ve[i] = sum;
    **************************************************/
  }
  
  return out;
}


/* get_col -- gets a specified column of a matrix and retruns it as a vector */
VEC	*get_col(MAT *mat, unsigned int col,VEC *vec) {
  unsigned int	i;
  
  if ( mat==(MAT *)NULL )
    error(E_NULL,"get_col");
  if ( col >= mat->n )
    error(E_RANGE,"get_col");
  if ( vec==(VEC *)NULL || vec->dim<mat->m )
    vec = v_resize(vec,mat->m);
  
  for ( i=0; i<mat->m; i++ )
    vec->ve[i] = mat->me[i][col];
  
  return (vec);
}

/* _set_col -- sets column of matrix to values given in vec (in situ) */
MAT	*_set_col(MAT *mat, unsigned int col, VEC *vec, unsigned int i0) {
  unsigned int	i,lim;
  
  if ( mat==(MAT *)NULL || vec==(VEC *)NULL )
    error(E_NULL,"_set_col");
  if ( col >= mat->n )
    error(E_RANGE,"_set_col");
  lim = min(mat->m,vec->dim);
  for ( i=i0; i<lim; i++ )
    mat->me[i][col] = vec->ve[i];
  
  return (mat);
}


#define MODULUS	LONG_MAX
#define MZ	0L

/* v_zero -- zero the vector x */
VEC	*v_zero(VEC *x){
  if ( x == VNULL )
    error(E_NULL,"v_zero");
  
  __zero__(x->ve,x->dim);
  return x;
}

static long mrand_list[56];
static int  started = FALSE;
static int  inext = 0, inextp = 31;

/* mrand -- pseudo-random number generator */
double mrand(void) {
  long	lval;
  static Real  factor = 1.0/((Real)MODULUS);
  
  if ( ! started )
    smrand(3127);
  
  inext = (inext >= 54) ? 0 : inext+1;
  inextp = (inextp >= 54) ? 0 : inextp+1;
  
  lval = mrand_list[inext]-mrand_list[inextp];
  if ( lval < 0L )
    lval += MODULUS;
  mrand_list[inext] = lval;
  
  return (double)lval*factor;
}

/* mrandlist -- fills the array a[] with len random numbers */
void	mrandlist(Real *a, int len) {
  int		i;
  long	lval;
  static Real  factor = 1.0/((Real)MODULUS);
  
  if ( ! started )
    smrand(3127);
  
  for ( i = 0; i < len; i++ )
    {
      inext = (inext >= 54) ? 0 : inext+1;
      inextp = (inextp >= 54) ? 0 : inextp+1;
      
      lval = mrand_list[inext]-mrand_list[inextp];
      if ( lval < 0L )
	lval += MODULUS;
      mrand_list[inext] = lval;
      
      a[i] = (Real)lval*factor;
    }
}

/* smrand -- set seed for mrand() */
void smrand( int seed) {
  int		i;
  
  mrand_list[0] = (123413*seed) % MODULUS;
  for ( i = 1; i < 55; i++ )
    mrand_list[i] = (123413*mrand_list[i-1]) % MODULUS;
  
  started = TRUE;
  
  for ( i = 0; i < 55*55; i++ )
    mrand();
}
#undef MODULUS
#undef MZ
#undef FAC

/* v_rand -- initialises x to be a random vector, components
   independently & uniformly ditributed between 0 and 1 */
VEC	*v_rand( VEC *x){
  
  if ( ! x )
    error(E_NULL,"v_rand");
  
  mrandlist(x->ve,x->dim);
  
  return x;
}


/* __ip__ -- inner product */
double	__ip__(double *dp1, double *dp2, int len) {
  int	i;
  double     sum;
  sum = 0.0;
  for ( i = 0; i < len; i++ )
    sum  += dp1[i]*dp2[i];
  
  return sum;
}

/* __mltadd__ -- scalar multiply and add c.f. v_mltadd() */
void	__mltadd__(double *dp1, double *dp2, double s, int len) {
  int	i;
  for ( i = 0; i < len; i++ )
    dp1[i] += s*dp2[i];
}

/* __smlt__ scalar multiply array c.f. sv_mlt() */
void	__smlt__(double *dp, double s, double *out, int len) {
  int	i;
  for ( i = 0; i < len; i++ )
    out[i] = s*dp[i];
}

/* __add__ -- add arrays c.f. v_add() */
void	__add__(double *dp1, double *dp2, double *out, int len) {
  int	i;
  for ( i = 0; i < len; i++ )
    out[i] = dp1[i] + dp2[i];
}

/* __zero__ -- zeros an array of floating point numbers */
void	__zero__(double *dp, int len){
    memset(((char *)dp),'\0',(len*sizeof(double)));
}



/* _m_copy -- copies matrix into new area */
MAT	*_m_copy(MAT *in, MAT *out, u_int i0, u_int j0) {
  u_int	i /* ,j */;
  
  if ( in==MNULL )
    error(E_NULL,"_m_copy");
  if ( in==out )
    return (out);
  if ( out==MNULL || out->m < in->m || out->n < in->n )
    out = m_resize(out,in->m,in->n);
  
  for ( i=i0; i < in->m; i++ )
    MEM_COPY(&(in->me[i][j0]),&(out->me[i][j0]), (in->n - j0)*sizeof(Real));
  return (out);
}

/* _v_copy -- copies vector into new area */
VEC	*_v_copy(VEC *in, VEC *out, u_int i0) {
  if ( in==VNULL )
    error(E_NULL,"_v_copy");
  if ( in==out )
    return (out);
  if ( out==VNULL || out->dim < in->dim )
    out = v_resize(out,in->dim);
  
  MEM_COPY(&(in->ve[i0]),&(out->ve[i0]),(in->dim - i0)*sizeof(Real));
  return (out);
}
