/*=================================================================*/
/*=   globals.h                                                   =*/
/*=   header file for global routines from treekin                =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-09-29 13:59:55 mtw>          =*/
/*=   $Id: globals.h,v 1.11 2006/09/29 16:33:32 mtw Exp $          =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _GLOBALS_H_
#define _GLOBALS_H_

typedef struct {         /* command-line options */
  int absrb;             /* make one lmin absorbing */
  int want_linear;       /* logarithmic time-scale */
  int want_verbose;      /* verbose output */
  int want_degenerate;   /* consider degeneracy */
  int dumpU;             /* dump U to a binary file */
  int dumpMathematica;   /* dump U to a Mathematica-readable file */
  int matexp;            /* use matrix-exponent routines */
  int binrates;          /* assume binary rates file */
  int n;                 /* read only n lmins */
  double T;               /* Temperature */
  double t0;              /* start time */
  double t8;              /* stop time */
  double tinc;            /* time increment */
  double *pini;           /* start population at lmin xxx */
  char method;           /* method to build transition matrix */
  FILE *INFILE;          /* input file (usually bar-file) */
  FILE *RATENFILE;       /* input file containing rates from barriers */
  char *sequence;        /* sequence */
} treekin_options;

treekin_options opt;
int lmins; /* # of lmins in the barrier tree, needed in FULL process */

void parse_commandline(int argc, char **argv);

void MxPrint(double *mx, char *name, char T);
#endif
/* End of file */
