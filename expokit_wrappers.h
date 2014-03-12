
// not needed apparently
#ifndef EXPOKIT_WRAPPERS_H
#define EXPOKIT_WRAPPERS_H

/*
#######################################################
# expokit_wrappers.h
# Wrapper functions for expokit functions.
#######################################################
# Assembled, not really authored, by
# Nick Matzke, matzke@berkeley.edu
# June 2012
#######################################################
*/

/*
Acknowledgements/sources:

1. Copied in part from a file in Lagrange:
http://code.google.com/p/lagrange/
https://github.com/blackrim/lagrange

Specifically:
 * RateMatrix.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: smitty
 *

...and the my_*.f wrappers for the EXPOKIT *.f code files.

*/

/*
2. Also copied in part (to get the .h file) from:

Python package "Pyprop":
http://code.google.com/p/pyprop/
http://pyprop.googlecode.com/svn/trunk/core/krylov/expokit/expokitpropagator.cpp
http://www.koders.com/python/fidCA95B5A4B2FB77455A72B8A361CF684FFE48F4DC.aspx?s=fourier+transform

Specifically:
pyprop/core/krylov/expokit/f2c/expokit.h
*/

/*
3. EXPOKIT package is available at:
http://www.maths.uq.edu.au/expokit/

Copyright:
http://www.maths.uq.edu.au/expokit/copyright
...or...
expokit_copyright.txt in this install
*/

/*
Note: Because these functions are defined with
extern"C", they can be sourced from R with e.g.
these commands:

################################
# Begin R code
################################
# Specify the *.so file (a 'shared object' file, in this case
# rexpokit.so), which is produced by the compiler linking
# the individual *.o files
wrapper_so_fn = paste(wd, "rexpokit.so", sep="")

# Unload the previous *.so file (if needed), then load the new one
dyn.unload(wrapper_so_fn)
dyn.load(wrapper_so_fn)

# Alledgely, this will list the loaded routines, but
# it doesn't seem to work
(dlls = getLoadedDLLs())
getDLLRegisteredRoutines(dll="rexpokit", addNames = TRUE)

# However, they are available, as you can see:
is.loaded("wrapalldmexpv_")
is.loaded("wrapsingledmexpv_")
is.loaded("wrapdgpadm_")
###############################
*/


// not needed apparently
//#include <core/common.h>

// not needed apparently
//namespace expokit_wrappers {

extern"C" {
// These functions are defined in my_matexp.f, which calls the following functions/

// DMEXPV contains an additional check ensuring sums to 1, for Markov-chain applications

// This exponentiates a large, sparse matrix (sparse = lots of zeros)
// Before input, the matrix should be transposed and
// put into coordinate list (COO) format:
// http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29
//
// i.e.:
// ia = row number
// ja = col number
// a  = corresponding value of that cell (zeros are excluded)
void wrapalldmexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
		int *ia, int *ja, double *a, int *nz, double * res);

// This returns just one row (?) of the transition matrix, useful in Lagrange;
// same inputs as wrapalldmexpv_.
void wrapsingledmexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
			int *ia, int *ja, double *a, int *nz, double * res);

// DGEXPV contains an additional check ensuring sums to 1, for Markov-chain applications
void wrapalldgexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
		int *ia, int *ja, double *a, int *nz, double * res);

// This returns just one row (?) of the transition matrix, useful in Lagrange;
// same inputs as wrapalldmexpv_.
void wrapsingledgexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
			int *ia, int *ja, double *a, int *nz, double * res);


// The myDMEXPV etc. functions provide direct access to EXPOKIT functions;
// This should be faster, especially for sparse matrices
// Here, you input v (starting probabilities) and it fills in w, which are the
// output probabilities (in output list item #5)
void myDMEXPV_(int* n, int* m, double* t, double* v, double* w, double* tol,
	double* anorm, double* wsp, int* lwsp, int* iwsp, int* liwsp, int* itrace, int* iflag,
		int* ia, int* ja, double* a, int* nz );

void myDGEXPV_(int* n, int* m, double* t, double* v, double* w, double* tol,
	double* anorm, double* wsp, int* lwsp, int* iwsp, int* liwsp, int* itrace, int* iflag,
		int* ia, int* ja, double* a, int* nz );


// This exponentiates a small, dense (no zeros) matrix with the padm
// approximation.  Here, the input matrix is just the
// matrix, tranposed and then put into a list of numbers, e.g.:
//
/*
# R code
# transpose the matrix in R
tmat = t(mat)
# convert into list of numbers (as numeric or double)
tmat_nums = as.double(tmat)
*/
void wrapdgpadm_(int * ideg,int * m,double * t,double * H,int * ldh,
	double * wsp,int * lwsp,int * ipiv,int * iexph,int *ns,int *iflag );

} // end extern "C"

// not needed apparently
//};
#endif
