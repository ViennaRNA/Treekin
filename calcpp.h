#ifndef _CALCPP_H_
#define _CALCPP_H_


// check for ergodicity
void MxEgro(double **U, double **p0, int dim);

// print probabilities (adds zero (non-ergodic) columns)
// returns sum of these probabilities
double PrintProb(double *line, int dim, double time);
double PrintProbFull(double *line, int dim, double time, int lmins);
double PrintProbNR(double *line, int dim, double time);

int ConvergenceReached(double *p8, double *pt, int dim, int full);

// rescale is --minimal-rate hass been provided
void MxRescale(double *U, int dim, double desired_rate);
void MxRescaleH(double *U, int dim, double hard_rescale);
void MxTimes(double *U, int dim, double times);
int *MxErgoEigen(double *U, int dim);
//void TestExpokit(double *R, int dim, double *p0, double t_start, double t_end, double t_inc);

//void TestEigen(double *R, int dim, double **evals, double **evecs);

//void printmat_mpf(int N, int M, mpf_class * A, int LDA);
void MxEV_Mpack_Sym_gmp(const double *U, int dim, double *evals, double *evecs, int precision);
void MxEV_Mpack_Sym_qd(const double *U, int dim, double *evals, double *evecs);
void MxEV_Mpack_Sym_dd(const double *U, int dim, double *evals, double *evecs);
void MxEV_Mpack_Sym_mpfr(const double *U, int dim, double *evals, double *evecs, int precision);
void MxEV_Mpack_Sym_float128(const double *U, int dim, double *evals, double *evecs);
void MxEV_Mpack_Sym_longdouble(const double *U, int dim, double *evals, double *evecs);
void MxEV_Mpack_Sym_double(const double *U, int dim, double *evals, double *evecs);

#endif
