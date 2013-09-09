#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <set>
#include <queue>

extern "C" {
  #include "calc.h"
}

using namespace std;

extern "C" void MxEgro(double **U, double **p0, int dim);
extern "C" double PrintProb(double *line, int dim, double time);
extern "C" double PrintProbNR(double *line, int dim, double time);
extern "C" double PrintProbFull(double *line, int dim, double time, int lmins);

vector<int> reorganize; // reorganize array (so if LM 0 1 3 were reachable and 2 not, reorganize will contain r[0]=0 r[1]=1 r[2]=3), so r[x] = old position of x
int last_dim;

void MxEgro(double **Up, double **p0p, int dim)
{
  // aliases
  double *U = *Up;
  double *p0 = *p0p;

  // lokalize first non-empty state
  int first = 0;
  while (p0[first]==0.0) first++;

  // make set of all non-empty states
  set<int> set_ergo;
  queue<int> que_ergo;

  set_ergo.insert(first);
  que_ergo.push(first);

  // fill the sets of non-empty states
  while(!que_ergo.empty() && (int)set_ergo.size()<dim) {
    int to_do = que_ergo.front();
    que_ergo.pop();

    // collect contingency
    for (int i=0; i<dim; i++) {
      if (i==to_do) continue;
      if (U[i*dim + to_do]>0.0 && (int)set_ergo.count(i)==0) {
        set_ergo.insert(i);
        que_ergo.push(i);
      }
    }
  }

  // check ergodicity
  if ((int)set_ergo.size()==dim) return; // all ok
  else {
    int i=first+1;
    while (i<dim) {
      if (p0[i]>0.0 && set_ergo.count(i)==0) { // disconnected :(
        fprintf(stderr, "ERROR: Matrix is non-ergodic and initial populations are disconected!! Exiting...\n");
        exit(-1);
      }
      i++;
    }

    // fill helper reorganize array
    reorganize.resize(set_ergo.size());
    i=0;
    for (set<int>::iterator it=set_ergo.begin(); it!=set_ergo.end(); it++) {
      reorganize[i++]=*it;
    }

    // reorganize matrix
    for (int i=0; i<(int)set_ergo.size(); i++) {
      for (int j=0; j<(int)set_ergo.size(); j++) {
        U[i*set_ergo.size()+j]=U[reorganize[i]*dim+reorganize[j]];
      }
    }
    last_dim = dim;
    dim = set_ergo.size();
    *Up = (double*)realloc(U, dim*dim*sizeof(double));

    // reorganize p0
    for (int i=0; i<(int)set_ergo.size(); i++) {
      p0[i]=p0[reorganize[i]];
    }
    *p0p = (double*)realloc(p0, dim*sizeof(double));

    // set dimension to global
    MxInit(dim);
    if (!opt.quiet) fprintf(stderr, "WARNING: Matrix is non-ergodic! Decreasing dimension to %d.\n", dim);
    //MxPrint(U, "Ergodic U", 'm');

    if (dim == 1) {
      PrintDummy(p0);
      exit(EXIT_SUCCESS);
    }
  }

}

double PrintProbFull(double *line, int dim, double time, int lmins)
{
  // for full process:
  vector<double> ptFULL (lmins, 0.0);
  if (reorganize.size() == 0) {
    for (int i=0; i<dim; i++) {
      ptFULL[E[i].ag] += line[i];
    }
  } else {
    for (int i=0; i<dim; i++) {
      ptFULL[E[reorganize[i]].ag] += line[i];
    }
  }
  // do the printing:
  double check = 0.0;
  printf("%e ", time);
  for (int i=0; i<lmins; i++) {
    if(ptFULL[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, ptFULL[i]);
      exit(EXIT_FAILURE);
    }
    /* map individual structure -> gradient basins */
    else   printf("%e ", fabs(ptFULL[i]));
    check += fabs(ptFULL[i]);
  }
  printf("\n");
  return check;
}

double PrintProbNR(double *line, int dim, double time)
{
  double check = 0.0;
  printf("%e ", time);

  for (int i=0; i<dim; i++) {
    if(line[i] < -0.01) {
      fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, line[i]);
      exit(EXIT_FAILURE);
    }
    /* map individual structure -> gradient basins */
    else printf("%e ", fabs(line[i]));
    check += fabs(line[i]);
  }

  printf("\n");
  return check;
}

double PrintProb(double *line, int dim, double time)
{
  double check = 0.0;
  printf("%e ", time);

  if (reorganize.size()==0) {
    for (int i=0; i<dim; i++) {
      if(line[i] < -0.01) {
        fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, line[i]);
        exit(EXIT_FAILURE);
      }
      /* map individual structure -> gradient basins */
      else printf("%e ", fabs(line[i]));
      check += fabs(line[i]);
    }
  } else {
    int j=0;
    for (int i=0; i<last_dim; i++) {
      if (j>(int)reorganize.size() || reorganize[j]!=i) {
        printf("%e ", 0.0);
      } else {
        if(line[j] < -0.01) {
          fprintf(stderr, "prob of lmin %i at time %e has become negative: %e \n", i+1, time, line[i]);
          exit(EXIT_FAILURE);
        }
        printf("%e ", fabs(line[j]));
        check += fabs(line[j]);
        j++;
      }
    }
  }
  printf("\n");
  return check;
}
