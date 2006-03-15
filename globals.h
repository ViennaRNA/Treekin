/*=================================================================*/
/*=   globals.h                                                   =*/
/*=   header file for global routines from treekin                =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 18:04:13 mtw>          =*/
/*=   $Id: globals.h,v 1.8 2006/03/15 18:04:29 mtw Exp $          =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include <glib.h>

typedef struct {         /* command-line options */
  short absrb;           /* make one lmin absorbing */
  short want_linear;     /* logarithmic time-scale */
  short want_verbose;    /* verbose output */
  short want_degenerate; /* consider degeneracy */
  short dumpU;           /* dump U to a binary file */
  short matexp;          /* use matrix-exponent routines */
  int n;                 /* read only n lmins */
  float T;               /* Temperature */
  float t0;              /* start time */
  float t8;              /* stop time */
  float tinc;            /* time increment */
  float *pini;           /* start population at lmin xxx */
  char method;           /* method to build transition matrix */
  FILE *INFILE;          /* input file (usually bar-file) */
  FILE *RATENFILE;       /* input file containing rates from barriers */
  char *sequence;        /* sequence */
  char *program_name;    /* name of executable file */ 
} treekin_options;

treekin_options opt;
int lmins; /* # of lmins in the barrier tree, needed in FULL process */

void parse_commandline(int argc, char **argv);

#endif
/* End of file */
