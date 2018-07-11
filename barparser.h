/*=================================================================*/
/*=   barparser.h                                                 =*/
/*=   header file for parsing bar-files etc.                      =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-11-25 22:47:42 mtw>          =*/
/*=   $Id: barparser.h,v 1.12 2006/11/27 13:46:29 mtw Exp $       =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _BARPARSER_H_
#define _BARPARSER_H_

#include <cstdlib>
#include <cstring>
#include "globals.h"
#include "bardata.h"
#include "mxccm.h"

/*#include "calc.h" */

class Barparser {
  static const int LMINBASE;
  static const int ARRAYSIZE;

  treekin_options *_opt;
  int _lmins;
  SubInfo *_E;

private:
  static char *my_getline(FILE *fp);


public:
  Barparser(treekin_options *opt,
            int             lmins,
            SubInfo         *E)
  {
    _opt    = opt;
    _lmins  = lmins;
    _E      = E;
  }


  ~Barparser()
  {}


  /* functions */
  int ParseBarfile(FILE     *fp,
                   BarData  **lmin);


  static int ParseSaddleFile(TypeDegSaddle **my_saddle);


  template<typename T>
  void MxFPrintD(T          *mx,
                 const char *name,
                 int        dim1,
                 int        dim2,
                 FILE       *out);


  template<typename T>
  void MxOneShorten(T   **shorten,
                    int fulldim);


  template<typename T>
  void MxRShorten(T   **shorten,
                  int fulldim,
                  int gdim);


  /* matrix shortening routine (when nstates specified) */
  template<typename T>
  int MxShorten(T   **shorten,
                int nstates,
                int my_dim,
                int max);


  /* read rates from binary file */
  int MxReadBinRates(FILE   *rate_file,
                     double **rate_mx,
                     int    nstates,
                     int    max);


  int ParseInfile(FILE    *infile_fp,
                  FILE    *mr_FP,
                  double  **microrates);


  int ParseRatesFile(FILE   *rates_FP,
                     double **Raten);


  /* visualisation */
  template<typename T>
  void VisualizeRates(char    *filename,
                      T       *R,
                      BarData *Data,
                      int     dim);
};


template<typename T>
inline
void
Barparser::MxFPrintD(T          *mx,
                     const char *name,
                     int        dim1,
                     int        dim2,
                     FILE       *out)
{
  int k, l;

  fprintf(out, "%s:{\n", name);
  for (k = 0; k < dim1; k++) {
    if (k != 0)
      fprintf(out, ",");

    fprintf(out, "{");
    for (l = 0; l < dim2; l++) {
      if (l != 0)
        fprintf(out, ",");

      fprintf(out, "%15.7g (%4d) ", mx[dim2 * k + l], dim2 * k + l);
    }
    fprintf(out, "}\n");
  }
  fprintf(out, "}-----------\n");
}


template<typename T>
inline
void
Barparser::MxRShorten(T   **shorten,
                      int fulldim,
                      int gdim)
{
  /*does: shortened = GG - GB*BB^(-1)*BG, where matrix tmp_rates is split as: */
  /*tmp_rates = (GG | GB) */
  /*            (BG | BB) */
  /* GG has dimension gdim*gdim; tmp_rates fulldim*fulldim */
  /* create matrices: */

  int bdim = fulldim - gdim;
  int i, j;
  T   *tmp_rates = *shorten;

  T   *gg = new T[gdim * gdim]; /*(T *)calloc(gdim*gdim,sizeof(T)); */
  T   *bg = new T[bdim * gdim]; /*(T *)calloc(bdim*gdim,sizeof(T)); */
  T   *bb = new T[bdim * bdim]; /*(T *)calloc(bdim*bdim,sizeof(T)); */
  T   *gb = new T[gdim * bdim]; /*(T *)calloc(gdim*bdim,sizeof(T)); */

  /* fill the matrices: (row = i; column = j) */
  for (i = 0; i < gdim; i++)
    for (j = 0; j < gdim; j++)
      gg[gdim * i + j] = tmp_rates[fulldim * i + j];

  for (i = 0; i < bdim; i++)
    for (j = 0; j < gdim; j++)
      bg[gdim * i + j] = tmp_rates[fulldim * (i + gdim) + j];

  for (i = 0; i < gdim; i++)
    for (j = 0; j < bdim; j++)
      gb[bdim * i + j] = tmp_rates[fulldim * i + j + gdim];

  for (i = 0; i < bdim; i++)
    for (j = 0; j < bdim; j++)
      bb[bdim * i + j] = tmp_rates[fulldim * (i + gdim) + j + gdim];

  /*MxFPrintD(tmp_rates, "Q", my_dim, my_dim, stderr);
   * MxFPrintD(gg, "GG", dim, dim, stderr);
   * MxFPrintD(bg, "BG", bdim, dim, stderr);
   * MxFPrintD(gb, "GB", dim, bdim, stderr);
   * MxFPrintD(bb, "BB", bdim, bdim, stderr);
   */
  /* result2 = gb*bb^(-1)*bg */
  Mxccm::minv(bb, bdim, _opt->FEPS);
  /*MxFPrintD(bb, "BBinv", bdim, bdim, stderr); */
  T *result = new T[gdim * bdim]; /*(T *)calloc(gdim*bdim,sizeof(T)); */
  Mxccm::mmul_singular(result, gb, bb, gdim, bdim, bdim, 0);
  /*MxFPrintD(result, "gb*bb-1", dim, bdim, stderr); */
  T *result2 = new T[gdim * gdim]; /*(T *)calloc(gdim*gdim,sizeof(T)); */
  Mxccm::mmul_singular(result2, result, bg, gdim, bdim, gdim, 1);

  if (_opt->want_verbose)
    MxFPrintD(result2, "gb*bb-1*bg", gdim, gdim, stderr);

  /* result2 = gg - result2 */
  for (i = 0; i < gdim; i++)
    for (j = 0; j < gdim; j++)
      result2[gdim * i + j] = gg[gdim * i + j] - result2[gdim * i + j];

  /*MxFPrintD(result2, "matrix after shortening", dim ,dim, stderr); */
  delete[] result;    /*free(result); */
  delete[] *shorten;  /*free(*shorten); */
  delete[] gg;        /*free(gg); */
  delete[] gb;        /*free(gb); */
  delete[] bg;        /*free(bg); */
  delete[] bb;        /*free(bb); */
  *shorten = result2;
}


template<typename T>
inline
void
Barparser::MxOneShorten(T   **shorten,
                        int fulldim)
{
  /*does: shortened = GG - GB*BB^(-1)*BG, where matrix tmp_rates is split as: */
  /*tmp_rates = (GG | GB) */
  /*            (BG | BB) */
  /* GG has dimension fulldim-1*fulldim-1; tmp_rates fulldim*fulldim */
  /* create matrices: */

  int gdim    = fulldim - 1;
  T   *result = new T[gdim * gdim]; /*(T *)calloc(gdim*gdim,sizeof(T)); */

  T   *tmp_rates = *shorten;

  int i, j;
  T   c = 1.0 / tmp_rates[fulldim * gdim + gdim];

  for (i = 0; i < gdim; i++) {
    for (j = 0; j < gdim; j++)
      /* just x - a*c^-1*b */
      result[gdim * i + j] = tmp_rates[fulldim * i + j] - c * tmp_rates[fulldim * gdim + j] *
                             tmp_rates[fulldim * i + gdim];
  }

  /*MxFPrintD(tmp_rates, "Q", fulldim, fulldim, stderr); */
  /*MxFPrintD(result, "Q-1", gdim, gdim, stderr); */

  delete[] *shorten; /*free(*shorten); */
  *shorten = result;
}


template<typename T>
inline
int
Barparser::MxShorten(T    **shorten,
                     int  nstates,
                     int  my_dim,
                     int  max)
{
  /* shorten the matrix from dimension my_dim to nstates: */
  if (my_dim > nstates) {
    if (!_opt->quiet)
      fprintf(stderr, "decreasing %d to %d\n", my_dim, nstates);

    /* first we need to fix the diagonal entries tmp_rates[i][i] = sum_j tmp_rates[i][j] */
    T   *tmp_rates = *shorten;
    int i, j;
    for (i = 0; i < my_dim; i++)
      tmp_rates[my_dim * i + i] = 0.0;
    for (i = 0; i < my_dim; i++) {
      T tmp = 0.00;
      /* calculate column sum */
      for (j = 0; j < my_dim; j++)
        tmp += tmp_rates[my_dim * j + i];
      tmp_rates[my_dim * i + i] = -tmp;
    }

    /* shorten by some value max */
    if (max != 1) {
      while (my_dim - max > nstates) {
        MxRShorten(shorten, my_dim, my_dim - max);
        if (!_opt->quiet)
          fprintf(stderr, "%d done...\n", my_dim - max);

        my_dim -= max;
      }
      MxRShorten(shorten, my_dim, nstates);
      my_dim = nstates;
    } else {
      while (my_dim != nstates) {
        MxOneShorten(shorten, my_dim);
        my_dim--;
        if (!_opt->quiet && my_dim % 100 == 0 && my_dim > 0)
          fprintf(stderr, "%d done...\n", my_dim);
      }
    }
  }

  return my_dim;
}


/* visualisation */
template<typename T>
inline
void
Barparser::VisualizeRates(char    *filename,
                          T       *R,
                          BarData *Data,
                          int     dim)
{
  FILE  *dotf;

  char  fname[100];
  int   i, j;

  strcpy(fname, filename);
  strcat(fname, ".dot");
  dotf = fopen(fname, "w");

  if (dotf) {
    fprintf(dotf, "digraph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
    /* nodes: */
    for (i = 0; i < dim; i++) {
      char energy[20] = "";
      if (Data != NULL)
        sprintf(energy, "\\n%.2f", Data[i].energy);

      fprintf(dotf, "\"%d\" [label=\"%d%s\"]\n", i + 1, i + 1, energy);
    }
    fprintf(dotf, "\n");

    /* edges: */
    int count = 1;
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        if (i != j && R[dim * j + i] > 0.) {
          fprintf(dotf,
                  "\"%d\" -> \"%d\" [label=\"%10.4g\"]\n",
                  i + 1,
                  j + 1,
                  (double)R[dim * j + i]);
          count++;
        }
      }
    }
    fprintf(dotf, "\n}\n");

    fclose(dotf);
  }

  char  syst[200];
  int   ref;
  sprintf(syst, "dot -Tps < %s > %s.eps", fname, filename);
  ref = system(syst);
}


#endif

/* End of file */
