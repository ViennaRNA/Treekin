/* globals.h */
/* Last changed Time-stamp: <2003-10-09 11:15:39 mtw> */
/*  static char rcsid[] = "$Id: globals.h,v 1.5 2003/10/09 17:01:35 mtw Exp $"; */

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
} markov_options;

markov_options opt;
int in_nr; /* counts # of elements read from FULL process input */
int lmins; /* # of lmins in the barrier tree, needed in FULL process */

void parse_commandline(int argc, char **argv);

/* End of file */
