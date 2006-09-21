/* treekin_cmdline.h */

/* File autogenerated by gengetopt version 2.17  */

#ifndef TREEKIN_CMDLINE_H
#define TREEKIN_CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
#define CMDLINE_PARSER_PACKAGE PACKAGE
#endif

#ifndef CMDLINE_PARSER_VERSION
#define CMDLINE_PARSER_VERSION VERSION
#endif

struct gengetopt_args_info
{
  const char *help_help; /* Print help and exit help description.  */
  const char *version_help; /* Print version and exit help description.  */
  int absorb_arg;	/* Make a state absorbing.  */
  char * absorb_orig;	/* Make a state absorbing original value given at command line.  */
  const char *absorb_help; /* Make a state absorbing help description.  */
  char * method_arg;	/* Select method to build transition matrix:\nA ==>Arrhenius-like kinetics\nF ==> Full process kinetics (whole subopt)\nI==>use rates from barriers (default='A').  */
  char * method_orig;	/* Select method to build transition matrix:\nA ==>Arrhenius-like kinetics\nF ==> Full process kinetics (whole subopt)\nI==>use rates from barriers original value given at command line.  */
  const char *method_help; /* Select method to build transition matrix:\nA ==>Arrhenius-like kinetics\nF ==> Full process kinetics (whole subopt)\nI==>use rates from barriers help description.  */
  float t0_arg;	/* Start time.  */
  char * t0_orig;	/* Start time original value given at command line.  */
  const char *t0_help; /* Start time help description.  */
  float t8_arg;	/* Stop time.  */
  char * t8_orig;	/* Stop time original value given at command line.  */
  const char *t8_help; /* Stop time help description.  */
  float Temp_arg;	/* Temperatur in Celsius.  */
  char * Temp_orig;	/* Temperatur in Celsius original value given at command line.  */
  const char *Temp_help; /* Temperatur in Celsius help description.  */
  int nstates_arg;	/* Read <int> states.  */
  char * nstates_orig;	/* Read <int> states original value given at command line.  */
  const char *nstates_help; /* Read <int> states help description.  */
  char ** p0_arg;	/* Set initial population of state <int> to <float>\nCan be given multiple times\n(NOTE: sum of <float> must equal 1).  */
  char ** p0_orig;	/* Set initial population of state <int> to <float>\nCan be given multiple times\n(NOTE: sum of <float> must equal 1) original value given at command line.  */
  int p0_min; /* Set initial population of state <int> to <float>\nCan be given multiple times\n(NOTE: sum of <float> must equal 1)'s minimum occurreces */
  int p0_max; /* Set initial population of state <int> to <float>\nCan be given multiple times\n(NOTE: sum of <float> must equal 1)'s maximum occurreces */
  const char *p0_help; /* Set initial population of state <int> to <float>\nCan be given multiple times\n(NOTE: sum of <float> must equal 1) help description.  */
  float tinc_arg;	/* Time scaling factor (for log time-scale).  */
  char * tinc_orig;	/* Time scaling factor (for log time-scale) original value given at command line.  */
  const char *tinc_help; /* Time scaling factor (for log time-scale) help description.  */
  int degeneracy_flag;	/* Consider degeracy in transition rates (default=off).  */
  const char *degeneracy_help; /* Consider degeracy in transition rates help description.  */
  int exponent_flag;	/* Use matrix-expontent routines, NO diagonalization (default=off).  */
  const char *exponent_help; /* Use matrix-expontent routines, NO diagonalization help description.  */
  int umatrix_flag;	/* Dump transition matrix U to a binary file mx.bin (default=off).  */
  const char *umatrix_help; /* Dump transition matrix U to a binary file mx.bin help description.  */
  int mathematicamatrix_flag;	/* Dump transition matrix U to Mathematica-readable file mxMat.txt (default=off).  */
  const char *mathematicamatrix_help; /* Dump transition matrix U to Mathematica-readable file mxMat.txt help description.  */
  int info_flag;	/* show settings (default=off).  */
  const char *info_help; /* show settings help description.  */
  int verbose_flag;	/* verbose output (default=off).  */
  const char *verbose_help; /* verbose output help description.  */
  
  int help_given ;	/* Whether help was given.  */
  int version_given ;	/* Whether version was given.  */
  int absorb_given ;	/* Whether absorb was given.  */
  int method_given ;	/* Whether method was given.  */
  int t0_given ;	/* Whether t0 was given.  */
  int t8_given ;	/* Whether t8 was given.  */
  int Temp_given ;	/* Whether Temp was given.  */
  int nstates_given ;	/* Whether nstates was given.  */
  unsigned int p0_given ;	/* Whether p0 was given.  */
  int tinc_given ;	/* Whether tinc was given.  */
  int degeneracy_given ;	/* Whether degeneracy was given.  */
  int exponent_given ;	/* Whether exponent was given.  */
  int umatrix_given ;	/* Whether umatrix was given.  */
  int mathematicamatrix_given ;	/* Whether mathematicamatrix was given.  */
  int info_given ;	/* Whether info was given.  */
  int verbose_given ;	/* Whether verbose was given.  */

  char **inputs ; /* unamed options */
  unsigned inputs_num ; /* unamed options number */
} ;

extern const char *gengetopt_args_info_purpose;
extern const char *gengetopt_args_info_usage;
extern const char *gengetopt_args_info_help[];

int cmdline_parser (int argc, char * const *argv,
  struct gengetopt_args_info *args_info);
int cmdline_parser2 (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

void cmdline_parser_print_help(void);
void cmdline_parser_print_version(void);

void cmdline_parser_init (struct gengetopt_args_info *args_info);
void cmdline_parser_free (struct gengetopt_args_info *args_info);

int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);

extern char *cmdline_parser_method_values[] ;	/* Possible values for method.  */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* TREEKIN_CMDLINE_H */
