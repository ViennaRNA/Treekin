/*=================================================================*/
/*=   globals.h                                                   =*/
/*=   header file for global routines from treekin                =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-11-24 19:17:45 mtw>          =*/
/*=   $Id: globals.h,v 1.15 2006/11/24 18:24:11 mtw Exp $         =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include <cstdio>
#include <getopt.h>
#include "treekin_cmdline.h"

typedef enum {MPACK_GMP = 2, MPACK_MPFR = 4, MPACK_DD = 8, MPACK_QD = 16, MPACK_FLOAT128 = 32, MPACK_LD = 64, MPACK_DOUBLE = 128 } MPackMethod;

typedef struct {         /* command-line options */
  char *basename;        /* base name of processed file */
  int absrb;             /* make one lmin absorbing */
  int real_abs;          /* remember absorbing state */
  int want_linear;       /* logarithmic time-scale */
  int want_verbose;      /* verbose output */
  int want_degenerate;   /* consider degeneracy */
  int dumpU;             /* dump U to a binary file */
  int dumpMathematica;   /* dump U to a Mathematica-readable file */
  int matexp;            /* use matrix-exponent routines */
  int binrates;          /* assume binary rates file */
  int rrecover;          /* recover from previous diagonalization */
  int wrecover;          /* write recovery file */
  int n;                 /* read only n lmins */
  double T;              /* Temperature */
  double t0;             /* start time */
  double t8;             /* stop time */
  double tinc;           /* time increment */
  double *pini;          /* start population at lmin xxx */
  char method;           /* method to build transition matrix */
  FILE *INFILE;          /* input file (usually bar-file) */
  FILE *RATFILE;         /* input file (rates) */
  char *sequence;        /* sequence */
  int  fpt;              /* switch for fpt-related calculations */
  int  fpt_num;          /* state to count ftp's to (-1 for all states)*/
  char *fpt_file;        /* output file for first passage times */
  char *vis_file;        /* output file for visualisation */
  int quiet;             // be quiet?
  int just_sh;           // just shorten
  int max_decrease;      // how many states to decrease at once
  char num_err;          // method for numerical error handling
  double FEPS;           // precision
  int useplusI;          // wheter to use the plus I matrix or no (eigenvalues are all +1)
  double minimal_rate;   // rescale to minimal rate?
  double hard_rescale;   // hard rescale rates?
  char *equil_file;      // file for equilibrium distribution
  double times;          // multiply the rates matrix?
  int mpackMethod_Bits;  // use mpack library for eigenvalue computation with the given number of bits
  int mpackMethod;       // MPackMethod (compromise between precision and speed)
  int warnings;          // all warnings on?
} treekin_options;

class Globals {

private:
  static Globals* instance;
  Globals() : lmins(0),opt(){}
  ~Globals(){}
  void ini_globs(void);
  void set_parameters(void);
  void display_settings(void);
  void to_basename(char *arg);
  gengetopt_args_info args_info;

public:
  treekin_options opt;
  int lmins; /* # of lmins in the barrier tree, needed in FULL process */

  static double TZERO;

  static Globals* initGlobals(){
    if(instance == NULL){
      instance = new Globals();
    }
    return instance;
  }
  static void destroy();

  void parse_commandline(int argc, char **argv);
  void free_gengetopt();
};



#endif
/* End of file */
