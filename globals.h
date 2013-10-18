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

#define TZERO 0.1

//typedef enum { false = 0, true = 1 } bool;

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
} treekin_options;

treekin_options opt;
int lmins; /* # of lmins in the barrier tree, needed in FULL process */

void parse_commandline(int argc, char **argv);

void free_gengetopt();

void MxFPrint(double *mx, char *name, char T, FILE *out, int pure);
void MxPrint(double *mx, char *name, char T);
#endif
/* End of file */
