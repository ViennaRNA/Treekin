/* globals.h */
/* Last changed Time-stamp: <2003-07-12 18:29:29 mtw> */
/*  static char rcsid[] = "$Id: globals.h,v 1.1 2003/07/14 07:42:20 mtw Exp $"; */

typedef struct {         /* command-line options */
  short absrb;           /* make one lmin absorbing */
  short want_linear;     /* logarithmic time-scale */
  short want_verbose;    /* verbose output */
  short want_degenerate; /* consider degeneracy */
  short dumpU;           /* dump U to a binary file */
  int n;                 /* read only n lmins */
  float T;               /* Temperature */
  float t0;              /* start time */
  float t8;              /* stop time */
  float tinc;            /* time increment */
  float *pini;           /* start population at lmin xxx */
  char method;           /* method to build transition matrix */
  FILE *INFILE;          /* input file (usually bar-file) */
  char *sequence;        /* sequence */
  char *program_name;    /* name of executable file */ 
} markov_options;

markov_options opt;

void parse_commandline(int argc, char **argv);

/* End of file */
