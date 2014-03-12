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


void TestExpokit(double *R, int dim, double *p0, double t_start, double t_end, double t_inc);


#endif
