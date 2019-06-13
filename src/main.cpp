/*=================================================================*/
/*=   main.c                                                      =*/
/*=   main file for treekin                                       =*/
/*=                                                               =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mxccm.h"
#include "calc.h"
using namespace std;

#ifdef WITH_MLAPACK
# include "treekinCastableTypes.h"
using namespace treekinCastableTypes;
#endif


template<class T>
T *
convertRates(double *rates,
             size_t dim)
{
  size_t  length  = dim * dim;
  T       *data   = new T[length];

  for (size_t i = 0; i < length; i++)
    data[i] = (T)(rates[i]);
  free(rates);
  return data;
}


template<class T>
int
treekin_main_precision(Globals *globalParameters)
{
  clock_t         clck1 = clock();

  treekin_options *opt  = &globalParameters->opt;
  int             lmins = globalParameters->lmins;
  SubInfo         *E; /*will be filled in Barparser::ParseInfile */

  BarData         *Data = NULL;
  /*  U - matrix (Q+I)^T, where Q is infinitesimal generator (^T - transposed)
   *  S - eigenvectors of U
   *  p0 - distribution in the beginning
   *  p8 - stable (end) distribution
   *  R - tmp rates matrix (transposed) */
  T         *U = NULL, *S = NULL, *p0 = NULL, *p8 = NULL, *R = NULL;
  double    *tmpReadRates = NULL;
  int       dim;

  Barparser *bp = new Barparser(opt, dim, E);

  switch (opt->method) {
    case 'I':
      if (!opt->quiet)
        fprintf(stderr,
                "Using rates matrix from STDIN for constructing transition matrix\n");

      if (opt->binrates)
        dim = bp->MxReadBinRates(opt->RATES, &tmpReadRates, opt->n, opt->max_decrease);
      else
        dim = bp->ParseRatesFile(opt->RATES, &tmpReadRates);

      if (opt->absrb) {
        int dimb;
        if (!opt->BARFILE) {
          fprintf(stderr,
                  "ERROR: bar file must be provided via --bar option for computing proper free energy/partition function for the absorbing state\n");
          exit(EXIT_FAILURE);
        } else {
          dimb = bp->ParseBarfile(opt->BARFILE, &Data); /* NB: opt.BARFILE; not opt.RATES ! */
          if (dim != dimb) {
            fprintf(stderr,
                    "ERROR: dimension mismatch among input rates file and bar file %d != %d\n",
                    dim, dimb);
            exit(EXIT_FAILURE);
          }
        }
      }

      if (dim == 0) {
        fprintf(stderr, "ERROR: Rate file empty!\n");
        globalParameters->free_gengetopt();
        exit(EXIT_FAILURE);
      }

      R = convertRates<T>(tmpReadRates, (size_t)dim);
      /* graph visualization */
      if (opt->vis_file)
        bp->VisualizeRates(opt->vis_file, R, Data, dim);

      break;
    case 'A':
      if (!opt->quiet)
        fprintf(stderr, "Using bar file from STDIN for constructing transition matrix\n");

      dim = bp->ParseBarfile(opt->BARFILE, &Data);
      break;
  }

  delete bp;

  Calc<T> *calc = new Calc<T>(globalParameters, E, dim);

  /* create rate matrix U which is used throughout the program */
  U = (T *)calc->MxBar2Matrix(Data, R);

  Calccpp *calccpp = new Calccpp(opt, E);

  /* rescale in case minimal_rate has been set */
  if (opt->hard_rescale != 1.0)
    calccpp->MxRescaleH(U, dim, opt->hard_rescale);
  else if (opt->minimal_rate > 0.0)
    calccpp->MxRescale(U, dim, opt->minimal_rate);

  /* multiply in case times was set */
  if (opt->times != 1.0)
    calccpp->MxTimes(U, dim, opt->times);

  delete calccpp;

  /* create initial population probability vector */
  calc->MxStartVec(&p0);

  /* check for ergodicity and adjust accordingly */
  dim = calc->MxEgro(&U, &p0, dim);
  if (opt->want_verbose)
    calc->MxPrint(U, "Ergodic U", 'm');

  /* allocate space for other matrices */
  calc->MxGetSpace(&p8);

  if (!opt->quiet)
    fprintf(stderr, "Time to initialize: %.2f secs.\n", (clock() - clck1) / (double)CLOCKS_PER_SEC);

  clck1 = clock();

  /* compute equilibrium distribution */
  if (opt->method == 'F') {
    calc->MxEqDistrFULL(E, p8);
  } else {
    if (opt->method == 'A')
      calc->MxEqDistr(Data, &p8);
    else
      calc->MxEqDistrFromDetailedBalance(U, &p8);

    /*MxEqDistrFromLinSys(U, &p8); */
  }

  /* dump equilibrium */
  if (opt->equil_file) {
    FILE *equil = fopen(opt->equil_file, "w");
    if (equil) {
      calc->MxFPrint(p8, "Equilibrium distribution", 'v', equil, 1);
      fclose(equil);
    }
  }

  /* compute first passage times */
  if (opt->fpt) {
    /* output file? */
    FILE *fpt = stderr;
    if (opt->fpt_file != NULL) {
      fpt = fopen(opt->fpt_file, "w");
      if (fpt == NULL) {
        if (!opt->quiet)
          fprintf(stderr,
                  "Cannot open file \"%s\" for FPT output, using \"stderr\" instead\n",
                  opt->fpt_file);

        fpt = stderr;
      }
    }

    /* all or just one state? */
    if ((opt->fpt_num == -1) || opt->absrb) {
      /* all FPTs */
      calc->MxFPT(U, p8, fpt);
    } else {
      /* just to state opt.ftp_num */
      T *res = calc->MxFPTOneState(U, opt->fpt_num - 1);  /* starting from 0 */
      if (res != NULL) {
        if (fpt == stderr)
          fprintf(fpt, "First passage times to state number %d:\n", opt->fpt_num);

        int k;
        for (k = 0; k < dim; k++)
          fprintf(fpt, "%15.7f ", (double)res[k]);
        if (fpt == stderr)
          fprintf(fpt, "\n---\n");
        else
          fprintf(fpt, "\n");

        delete[] res; /*free(res); */
      }
    }

    fclose(fpt);
    if (!opt->quiet)
      fprintf(stderr,
              "Time to compute fpt: %.2f secs.\n",
              (clock() - clck1) / (double)CLOCKS_PER_SEC);

    clck1 = clock();
  }

  int exponentError = 0;
  if (!opt->just_sh) {
    /* diagonalization + iteration */
    if (opt->matexp) {
      exponentError = calc->MxExponent(p0, p8, U);
      if (exponentError) {
        calc->MxMemoryCleanUp();
        delete calc;
        delete[] U;
        delete[] p8;
        delete[] p0;
        free(Data);

        exit(EXIT_FAILURE);
      }
    } else {
      if (opt->rrecover) {
        /* read pre-calculated eigenv{alues,ectors} */
        calc->MxRecover(&S, p8);
      } else {
        clck1 = clock();
        calc->MxDiagonalize(U, &S, p8);
        if (!opt->quiet)
          fprintf(stderr,
                  "Time to diagonalize: %.2f secs.\n",
                  (clock() - clck1) / (double)CLOCKS_PER_SEC);

        clck1 = clock();
      }

      calc->MxIterate(p0, p8, S);
      if (!opt->quiet)
        fprintf(stderr, "Time to iterate: %.2f secs.\n",
                (clock() - clck1) / (double)CLOCKS_PER_SEC);
    }
  }

  if (opt->pini != NULL)
    free(opt->pini);

  /* memory cleanup */
  calc->MxMemoryCleanUp();
  delete calc;
  delete[] U;
  delete[] p8;
  delete[] p0;
  free(Data);
  return 0;
}


int
main(int  argc,
     char **argv)
{
  Globals *globalParameters = Globals::initGlobals();

  globalParameters->parse_commandline(argc, argv);

  int     error;

  switch (globalParameters->opt.mlapackMethod) {
#ifdef WITH_MLAPACK_DOUBLE
    case MLAPACK_DOUBLE:
      error = treekin_main_precision<double>(globalParameters);
      break;
#endif
#ifdef WITH_MLAPACK_LD
    case MLAPACK_LD:
      error = treekin_main_precision<long double>(globalParameters);
      break;
#endif
#ifdef WITH_MLAPACK___FLOAT128
    case MLAPACK_FLOAT128:
      error = treekin_main_precision<__float128>(globalParameters);
      break;
#endif
#ifdef WITH_MLAPACK_MPFR
    case MLAPACK_MPFR:
      mpfr::mpreal::set_default_prec(globalParameters->opt.mlapackMethod_Bits);
      error = treekin_main_precision<mpreal_castable>(globalParameters);
      break;
#endif
#ifdef WITH_MLAPACK_DD
    case MLAPACK_DD:
      error = treekin_main_precision<dd_real_castable>(globalParameters);
      break;
#endif
#ifdef WITH_MLAPACK_QD
    case MLAPACK_QD:
      error = treekin_main_precision<qd_real_castable>(globalParameters);
      break;
#endif
#ifdef WITH_MLAPACK_GMP
    case MLAPACK_GMP:
      mpf_set_default_prec(globalParameters->opt.mlapackMethod_Bits);
      error = treekin_main_precision<mpf_real_castable>(globalParameters);
      break;
#endif

    default:
      error = treekin_main_precision<double>(globalParameters);
      break;
  }

  globalParameters->destroy();

  return error;
}


/* End of file */
