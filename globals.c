/* globals.c */
/* Last changed Time-stamp: <2003-08-25 15:25:28 mtw> */
/* static char rcsid[] = "$Id: globals.c,v 1.3 2003/08/27 14:59:08 mtw Exp $"; */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include "globals.h"

static int decode_switches (int argc, char **argv);
static int check_pini_prob(float*);
static void ini_globs(void);
static void display_settings(void);
static void usage(int status);

static struct option const long_options[] =
{ {"verbose",     no_argument, 0, 'v'},
  {"help",        no_argument, 0, 'h'},
  {"absorb",      required_argument, 0, 0},
  {"t0",          required_argument, 0, 0},
  {"t8",          required_argument, 0, 0},
  {"tinc",        required_argument, 0, 0},
  {"nolog",       no_argument, 0, 0},
  {"Temp",        required_argument, 0, 0}, 
  {"method",      required_argument, 0, 0},
  {"p0",          required_argument, 0, 0},
  {"nlmins",      required_argument, 0, 0},
  {NULL, 0, NULL, 0}
};


/*==============================*/
static int decode_switches (int argc, char **argv) {
  int c, atmp, ntmp, p0flag = 0, i = 1, tinc_set = 0;
  float t0tmp, t8tmp, tinctmp, Ttmp, *pinitmp=NULL;
  char mtmp;
  int option_index = 0;

  ini_globs();
  
  while ((c = getopt_long (argc, argv,
			   "v"    /* verbose output */
			   "h"    /* help */
			   "d"     /* degeneracy */
			   "u",    /* dump U to bin-file */
			   long_options, &option_index)) != EOF)
    {
      switch (c)
	{
	case 0:
	  if (strcmp(long_options[option_index].name,"absorb")==0) {
	    atmp = -1;
	    if (sscanf(optarg, "%d", &atmp) == 0)
	      usage(EXIT_FAILURE);
	    else if (atmp >= 0)
	      opt.absrb = atmp;
	    else {
	      fprintf(stderr, "Value of --absorb must be >= 0, and not >%d<\n", atmp);
	      usage (EXIT_FAILURE);
	    }
	  }
	  
	  if (strcmp(long_options[option_index].name,"t0")==0) {
	    t0tmp = -1.0;
	    if (sscanf(optarg, "%f", &t0tmp) == 0)
	      usage(EXIT_FAILURE);
	    else if (t0tmp >= 0.)
	      opt.t0 = t0tmp;
	     else {
	      fprintf(stderr, "Value of --t0 must be >= 0. , and not >%f<\n", t0tmp);
	      usage (EXIT_FAILURE);
	    }
	  }
	  
	  if (strcmp(long_options[option_index].name,"t8")==0) {
	    t8tmp = -1.0;
	    if (sscanf(optarg, "%f", &t8tmp) == 0)
	      usage(EXIT_FAILURE);
	    else if (t8tmp >= 0.)
	      opt.t8 = t8tmp;
	     else {
	      fprintf(stderr, "Value of --t8 must be >= 0. , and not >%f<\n", t8tmp);
	      usage (EXIT_FAILURE);
	    }
	  }

	  if (strcmp(long_options[option_index].name,"tinc")==0) {
	    tinctmp = -1.0;
	    if (sscanf(optarg, "%f", &tinctmp) == 0)
	      usage(EXIT_FAILURE);
	    else if (tinctmp >= 0.){
	      tinc_set = 1;
	      opt.tinc = tinctmp;
	    }
	    else {
	      fprintf(stderr, "Value of --tinc must be >= 0. , and not >%f<\n", tinctmp);
	      usage (EXIT_FAILURE);
	    }
	  }

	  if (strcmp(long_options[option_index].name,"Temp")==0) {
	    Ttmp = -1.0;
	    if (sscanf(optarg, "%f", &Ttmp) == 0)
	      usage(EXIT_FAILURE);
	    else if ( Ttmp >= 0 )
	      opt.T = Ttmp;
	    else {
	      fprintf(stderr, "Value of --Temp must be >= 0, and not >%.2f<\n", Ttmp);
	      usage (EXIT_FAILURE);
	    }
	  }

	  if(strcmp(long_options[option_index].name,"method")==0) {
	    mtmp = 'Z';
	    if (sscanf(optarg, "%c", &mtmp) == 0)
	      usage(EXIT_FAILURE);
	    else if ( mtmp == 'A' || mtmp == 'B' || mtmp == 'C' || mtmp == 'F' || mtmp == 'I' )
	      opt.method = mtmp;
	    else {
	      fprintf(stderr, "Value of --method must be A|B|C|F, and not  >%d<\n", mtmp);
	      usage (EXIT_FAILURE);
	    }
	  }
	  
	  if (strcmp(long_options[option_index].name,"nlmins")==0) {
	    ntmp = -1;
	    if (sscanf(optarg, "%d", &ntmp) == 0)
	      usage(EXIT_FAILURE);
	    else if (ntmp >= 0)
	      opt.n = ntmp;
	    else {
	      fprintf(stderr, "Value of --nlmins must be >= 0, and not >%d<\n", ntmp);
	      usage (EXIT_FAILURE);
	    }
	  }
	  
	  if (strcmp(long_options[option_index].name,"p0")==0) {
            int lmintmp;
	    float poptmp = 0.;
	    if (!p0flag) {
	      pinitmp = (float *)calloc(3, sizeof(float));
	      if (pinitmp == NULL) {
		fprintf(stderr, "Could not allocate pinitmp!\n");
		exit(EXIT_FAILURE);
	      }
	      *pinitmp = 1;
	    }
	    else if (p0flag > 0) {
	      pinitmp = (float *)realloc(pinitmp, (*pinitmp+2)*sizeof(float));
	      if (pinitmp == NULL) {
		fprintf(stderr, "Could not re-allocate pinitmp!\n");
		exit(EXIT_FAILURE);
	      }
	    }
	    else {
	      fprintf(stderr, "p0flag is < 0, this should't happen >%d<\n", p0flag);
	      usage (EXIT_FAILURE);
	    }
	    if (sscanf(optarg, "%d=%f",&lmintmp, &poptmp) == 0)
	      usage(EXIT_FAILURE);
	    else if (poptmp >0. && poptmp < 1.01){
	      *(pinitmp + i)     = (float) lmintmp;
	      *(pinitmp + i + 1) = poptmp;
	      *pinitmp += 2;
	      i += 2;
	    };
	    
	    p0flag += 2;
	    opt.pini = pinitmp;
	  } 
	  break;
	  
	case 'd':
	  opt.want_degenerate = 1;
	  break;
	  
	case 'v':
	  opt.want_verbose = 1;
	  break;

	case 'u':
	  opt.dumpU = 1;
	  break;

	case 'h':
	  usage(0);
	  
	default:
	  usage (EXIT_FAILURE);
	}
    }
  if (i > 1){  /* just if --p0 was used */
    if (check_pini_prob(pinitmp)) {
      fprintf(stderr, "population prob must be == 1 in --p0 \n");
      exit(EXIT_FAILURE);
    }
  }
  return optind;
}

/*==============================*/
static void ini_globs(void) {
  opt.absrb           =      0;
  opt.T               =    37.;
  opt.t0              =      1;
  opt.t8              = 10000.;
   opt.want_degenerate =      0; 
  opt.tinc            =   1.02;
  opt.method          =    'B';
  opt.dumpU           =      0;
}

/*==============================*/
static void usage(int status) {
  printf("\n%s - simulate kinetic folding of RNA with a Markov process\n",
	 opt.program_name);
  
  printf("Usage: %s [OPTION]... [FILE]  , where [FILE] must have\n", opt.program_name);
  printf("barriers output-format, i.e. barriers --saddle --bsize --ssize\n"); 
  printf(
         "\nOptions:\n"
	 "--absorb <int>          Make lmin <int> an absorbing lmin (default none)\n"
	 "--t0 <float>            Start Time\n"
	 "--t8 <float>            Stop Time\n"
	 "--Temp <float>          Temperature in Celsius\n"
	 "--method <char>         Select method to build transition matrix\n"
	 "--nlmins <int>          Read <int> local minima (default till EOF)\n"
	 "--p0 (<int>=<float>)    Populate local minimum <int> with <float>\n"
	 "                        (NOTE: sum of <float> must equal 1)\n"
	 "--tinc <float>          Time scaling-factor (for log time-scale)\n"
	 "-d                      Consider degeracy in transition rates\n"
	 "-u                      Dump transition matrix U to a binary file matrix.bin\n"
	 "-v, --verbose           Verbose output (for debugging)\n"
	 "\n"
	 );
  display_settings();
  exit (status);
}

/*==============================*/
static void display_settings(void) {
  fprintf(stderr,
          "Settings:\n"
	  "--absorb   = %d\n"
	  "--t0       = %.2f\n"
	  "--t8       = %.2f\n"
	  "--tinc     = %.2f\n"
	  "--Temp     = %.2f\n"
	  "--method   = %c\n"
	  "--nlmins   = %d\n"
	  "-d         = %d\n"
	  "-u         = %d\n"
	  "--verbose  = %d\n",
	  opt.absrb,
	  opt.t0,
	  opt.t8,
	  opt.tinc,
	  opt.T,
	  opt.method,
	  opt.n,
	  opt.want_degenerate,
	  opt.dumpU,
	  opt.want_verbose);
}

/*==============================*/
static int check_pini_prob(float *pini_tmp){

  int i;
  float probtmp = 0.;

  for (i = 2; i < *pini_tmp; i +=2)
    probtmp += pini_tmp[i];
  if (probtmp <= 1.01 && probtmp >= 0.99) return 0;
  else return 1;
}

/*==============================*/
void parse_commandline(int argc, char **argv){

  int i;
  
  opt.program_name = argv[0];
  
  i = decode_switches (argc, argv);
  if (i==argc-1) {
    opt.INFILE = fopen(argv[i], "r");
    if (opt.INFILE==NULL) fprintf(stderr, "can't open file \n");
  } else {
    if (i < argc) usage(EXIT_FAILURE);
    opt.INFILE = stdin;
  }
}

/* End of file */


