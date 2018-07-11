/*=================================================================*/
/*=   matrix.h                                                    =*/
/*=   header file for matrix-calculation stuff from treekin       =*/
/*=   mainly type definitions from meschach library               =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 15:16:18 mtw>          =*/
/*=   $Id: matrix.h,v 1.3 2006/03/15 14:18:20 mtw Exp $           =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

/**************************************************************************
**
** Copyright (C) 1993 David E. Stewart & Zbigniew Leyk, all rights reserved.
**
**           Meschach Library
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

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include        <stdio.h>
#include        <malloc.h>
#include        <stdlib.h>
#include        <stddef.h>
#include        <string.h>
#include        <float.h>

#include        "err.h"

/* standard copy & zero functions */
#define MEM_COPY(from, to, size)  memmove((to), (from), (size))

#define Real double

#define MACHEPS DBL_EPSILON

extern int     isatty(int);


/* vector definition */
typedef struct  {
  unsigned int  dim, max_dim;
  double        *ve;
} VEC;

/* matrix definition */
typedef struct  {
  unsigned int  m, n;
  unsigned int  max_m, max_n, max_size;
  double        **me, *base; /* base is base of alloc'd mem */
} MAT;

/* allocate one object of given type */
#define NEW(type) ((type *)calloc((size_t)1, (size_t)sizeof(type)))

/* allocate num objects of given type */
#define NEW_A(num, type) ((type *)calloc((size_t)(num), (size_t)sizeof(type)))

/* re-allocate arry to have num objects of the given type */
#define RENEW(var, num, type) \
  ((var) = (type *)((var) ? \
                    realloc((char *)(var), (size_t)((num) *sizeof(type))) : \
                    calloc((size_t)(num), (size_t)sizeof(type))))

/* type independent min and max operations */
#define max(a, b)  ((a) > (b) ? (a) : (b))
#define min(a, b)  ((a) > (b) ? (b) : (a))

#undef TRUE
#define TRUE  1
#undef FALSE
#define FALSE 0

/* Dynamic memory allocation */

/* Should use M_FREE/V_FREE in programs instead of m/v_free()
 * as this is considerably safer -- also provides a simple type check ! */

/* get/resize vector to given dimension */
extern VEC *v_get(int), *v_resize(VEC *, int);


/* get/resize matrix to be m x n */
extern MAT *m_get(int,
                  int), *m_resize(MAT *, int, int);


/* free (de-allocate)  matrices and vectors */
extern int m_free(MAT *), v_free(VEC *);


/* MACROS */

/* macros that also check types and sets pointers to NULL */
#define M_FREE(mat) (m_free(mat), (mat) = (MAT *)NULL)
#define V_FREE(vec) (v_free(vec), (vec) = (VEC *)NULL)

/* Entry level access to data structures */
/* returns x[i] */
#define v_entry(x, i)    ((x)->ve[i])

/* x[i] <- val */
#define v_set_val(x, i, val)  ((x)->ve[i] = (val))

/* returns A[i][j] */
#define m_entry(A, i, j)    ((A)->me[i][j])

/* A[i][j] <- val */
#define m_set_val(A, i, j, val)  ((A)->me[i][j] = (val))

/* A[i][j] <- A[i][j] + val */
#define m_add_val(A, i, j, val)  ((A)->me[i][j] += (val))

/* MACROS */

/* Copying routines */
/* copy in to out starting at out[i0][j0] */
extern MAT *_m_copy(MAT           *in,
                    MAT           *out,
                    unsigned int  i0,
                    unsigned int  j0),
*m_move(MAT * in, int, int, int, int, MAT * out, int, int),
*vm_move(VEC * in, int, MAT * out, int, int, int, int);


/* copy in to out starting at out[i0] */
extern VEC *_v_copy(VEC           *in,
                    VEC           *out,
                    unsigned int  i0),
*v_move(VEC * in, int, int, VEC * out, int),
*mv_move(MAT * in, int, int, int, int, VEC * out, int);


/* MACROS */
#define m_copy(in, out)  _m_copy(in, out, 0, 0)
#define v_copy(in, out)  _v_copy(in, out, 0)


/* Initialisation routines -- to be zero, ones, random or identity */
extern VEC *v_zero(VEC *), *v_rand(VEC *);


/* Basic vector operations */

extern VEC *sv_mlt(double,
                   VEC *,
                   VEC *),    /* out <- s.x */
*mv_mlt(MAT *, VEC *, VEC *), /* out <- A.x */
*vm_mlt(MAT *, VEC *, VEC *), /* out^T <- x^T.A */
*v_add(VEC *, VEC *, VEC *),  /* out <- x + y */
*v_sub(VEC *, VEC *, VEC *);  /* out <- x - y */


/* returns inner product starting at component i0 */
extern double  _in_prod(VEC           *x,
                        VEC           *y,
                        unsigned int  i0),
__ip__(double *, double *, int);


/* see v_mltadd(), v_add(), v_sub() and v_zero() */
extern void  __mltadd__(double *,
                        double *,
                        double,
                        int),
__add__(double *, double *, double *, int),
__smlt__(double *, double, double *, int),
__zero__(double *, int);


/* MACRO */
/* usual way of computing the inner product */
#define in_prod(a, b)  _in_prod(a, b, 0)

/* Norms */
/* scaled vector norms -- scale == NULL implies unscaled */
/* returns (scaled) Euclidean norm */
extern double _v_norm2(VEC  *x,
                       VEC  *scale);


/* returns max_i |x[i]/scale[i]| */
extern double _v_norm_inf(VEC *x,
                          VEC *scale);


/* MACROS */
#define v_norm2(x)  _v_norm2(x, VNULL)
#define v_norm_inf(x) _v_norm_inf(x, VNULL)


extern MAT *m_add(MAT *A,
                  MAT *B,
                  MAT *out);                /* out <- A + B */


extern MAT *_set_col(MAT          *A,
                     unsigned int i,
                     VEC          *out,
                     unsigned int j0);


extern VEC *get_row(MAT *,
                    unsigned int,
                    VEC *);


extern VEC *get_col(MAT *,
                    unsigned int,
                    VEC *);


/* MACROS */
/* col j of A <- vec */
#define set_col(mat, col, vec)  _set_col(mat, col, vec, 0)

/* miscellaneous functions */

double  mrand(void);                   /* returns random # in [0,1) */


void  smrand(int seed);              /* seeds mrand() */


void    mrandlist(double  *x,
                  int     len);        /* generates len random numbers */


/* miscellaneous constants */
#define VNULL ((VEC *)NULL)
#define MNULL ((MAT *)NULL)


/* Hessenberg factorisation of A -- for schur() */
extern MAT *Hfactor(MAT *A,
                    VEC *diag1,
                    VEC *diag2);


extern MAT *makeHQ(MAT  *A,
                   VEC  *diag1,
                   VEC  *diag2,
                   MAT  *Qout),
/* Hout is the Hessenberg matrix in Hessenberg factorisation */
*makeH(MAT * A, MAT * Hout);


extern VEC  *hhvec(VEC *, u_int, Real *, VEC *, Real *);
extern VEC  *hhtrvec(VEC *, double, u_int, VEC *, VEC *);
extern MAT  *hhtrrows(MAT *, u_int, u_int, VEC *, double);
extern MAT  *hhtrcols(MAT *, u_int, u_int, VEC *, double);

extern void  givens(double,
                    double,
                    Real *,
                    Real *);


extern VEC  *rot_vec(VEC *, u_int, u_int, double, double, VEC *);   /* in situ */
extern MAT  *rot_rows(MAT *, u_int, u_int, double, double, MAT *);  /* in situ */
extern MAT  *rot_cols(MAT *, u_int, u_int, double, double, MAT *);  /* in situ */


/* eigenvalue routines */

/* computes real Schur form = Q^T.A.Q */
extern MAT *schur(MAT *A,
                  MAT *Q);


/* computes real and imaginary parts of the eigenvalues
 *             of A after schur() */
extern void  schur_evals(MAT  *A,
                         VEC  *re_part,
                         VEC  *im_part);


/* computes real and imaginary parts of the eigenvectors
 *            of A after schur() */
extern MAT *schur_vecs(MAT  *T,
                       MAT  *Q,
                       MAT  *X_re,
                       MAT  *X_im);


#endif
