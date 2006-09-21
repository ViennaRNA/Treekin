/*=================================================================*/
/*=   globals.c                                                   =*/
/*=   global routines for treekin                                 =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-09-21 15:52:56 mtw>          =*/
/*=   $Id: globals.c,v 1.10 2006/09/21 13:59:57 mtw Exp $          =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include "globals.h"
#include "treekin_cmdline.h"

static void ini_globs(void);
static void set_parameters(void);
static void display_settings(void);
static int  check_pini_prob(float*);

static struct gengetopt_args_info args_info;

/*==============================*/
void
parse_commandline(int argc, char **argv)
{
  ini_globs();
  if (cmdline_parser (argc, argv, &args_info) != 0){
    fprintf(stderr, "error while parsing command-line options\n");
    exit(EXIT_FAILURE);
  }
  set_parameters();
  if (args_info.inputs_num)
    opt.INFILE = fopen(args_info.inputs[0], "r");
  else
    opt.INFILE = stdin;
} 

/*==============================*/
static void
set_parameters(void)
{
  if(strncmp(args_info.method_arg, "F", 1)==0)
    opt.method = 'F';
  else if (strncmp(args_info.method_arg, "I", 1)==0)
    opt.method = 'I';

  if (args_info.p0_given){
    int i, j=1, lmintmp;
    float poptmp = 0., probtmp = 0., *pinitmp=NULL;
    pinitmp = (float *)calloc(2*args_info.p0_given+1, sizeof(float));
    *pinitmp = 1;
    for (i=0; i<args_info.p0_given; i++,j+=2){
      if (sscanf(args_info.p0_arg[i], "%d=%f",&lmintmp, &poptmp) == 0)
	exit(EXIT_FAILURE);
      if(lmintmp <1){
	fprintf(stderr, "States in --p0 must be >=1\n");
	exit (EXIT_FAILURE);	
      }
      else if (poptmp >0. && poptmp < 1.01){
	*(pinitmp + j)     = (float) lmintmp;
	*(pinitmp + j + 1) = poptmp;
	*pinitmp += 2;
      }  
    }
    opt.pini = pinitmp;
    for (i = 2; i < *pinitmp; i +=2)
      probtmp += pinitmp[i];
    if (probtmp <= 0.99 || probtmp >= 1.01){
      fprintf(stderr, "Values of --p0 must sum up to 1 (currently %5.2f)\n",
	      probtmp);
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.absorb_given){
    if( (opt.absrb = args_info.absorb_arg) <= 0 ){
      fprintf(stderr, "Value of --absorb must be >= 0\n");
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.t0_given){
    if( (opt.t0 = args_info.t0_arg) <= 0. ){
      fprintf(stderr, "Value of --t0 must be >= 0.\n");
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.t8_given){
    if( (opt.t8 = args_info.t8_arg) <= 0. ){
      fprintf(stderr, "Value of --t8 must be >= 0.\n");
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.tinc_given){
    if( (opt.tinc = args_info.tinc_arg) <= 0. ){
      fprintf(stderr, "Value of --tinc must be >= 0.\n");
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.Temp_given){
    if( (opt.T = args_info.Temp_arg) < -273.15 ){
      fprintf(stderr, "Value of --Temp must be > -273.15\n");
      exit (EXIT_FAILURE);
    }
  }
  
  if (args_info.nstates_given){
    if( (opt.n = args_info.nstates_arg) <= 1 ){
      fprintf(stderr, "Value of --nstates must be >= 1\n");
      exit (EXIT_FAILURE);
    }
  }
 
  if (args_info.degeneracy_given) opt.want_degenerate = 1;
  if (args_info.verbose_given) opt.want_verbose = 1;
  if (args_info.umatrix_given) opt.dumpU = 1;
  if (args_info.mathematicamatrix_given) opt.dumpMathematica = 1;
  if (args_info.exponent_given) opt.matexp = 1;
  if (args_info.info_given){
    display_settings();
    exit(EXIT_SUCCESS);
  }
}

/*==============================*/
static void
ini_globs(void)
{
  opt.absrb           =          0;
  opt.T               =         37.;
  opt.t0              =          0.1;
  opt.t8              = 1000000000.;
  opt.want_degenerate =          0; 
  opt.tinc            =          1.02;
  opt.method          =          'A';
  opt.dumpU           =          0;
  opt.dumpMathematica =          0;
  opt.matexp          =          0;
}

/*==============================*/
static void
display_settings(void)
{
  int i,j;
  fprintf(stderr,
          "Settings:\n"
	  "--absorb   = %3i\n"
	  "--t0       = %8.4f\n"
	  "--t8       = %7.2g\n"
	  "--tinc     = %6.2f\n"
	  "--Temp     = %6.2f\n"
	  "--method   = %c\n"
	  "--nstates  = %d\n"
	  "-d         = %d\n"
	  "-e         = %d\n"
	  "-u         = %d\n"
	  "-x         = %d\n"
	  "-v         = %d\n",
	  opt.absrb,
	  opt.t0,
	  opt.t8,
	  opt.tinc,
	  opt.T,
	  opt.method,
	  opt.n,
	  opt.want_degenerate,
	  opt.matexp,
	  opt.dumpU,
	  opt.dumpMathematica,
	  opt.want_verbose);
  for (i=0,j=1;i<args_info.p0_given;i++,j+=2){
    fprintf(stderr,
	    "--p0      %4i=%.2f\n",
	    (int)opt.pini[j],opt.pini[j+1]);
  }
}

/* End of file */


