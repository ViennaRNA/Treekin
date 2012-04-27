#ifndef _CALCPP_H_
#define _CALCPP_H_


// check for ergodicity
void MxEgro(double *U, double *p0, int dim);

// print probabilities (adds zero (non-ergodic) columns)
// returns sum of these probabilities
double PrintProb(double *line, int dim, double time);

#endif
