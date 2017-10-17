/*=================================================================*/
/*=   main.c                                                      =*/
/*=   main file for treekin                                       =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2017-07-27 18:40:23 ivo>          =*/
/*=   $Id: main.c,v 1.24 2006/11/27 13:47:57 mtw Exp $            =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#include "mxccm.h"
#include "treekinCastableTypes.h"
#include "calc.h"

using namespace std;

using namespace treekinCastableTypes;

template<class T>
T *convertRates(double *rates, size_t dim){
  size_t length = dim*dim;
  T *data = new T[length];
  for(size_t i = 0; i < length; i++){
   data[i] = (T)(rates[i]);
  }
  free(rates);
  return data;
}

template<class T>
int treekin_main_precision(Globals *globalParameters){

  clock_t clck1 = clock();

  treekin_options* opt = &globalParameters->opt;
  int lmins = globalParameters->lmins;
  SubInfo *E; //will be filled in Barparser::ParseInfile

  BarData *Data=NULL;
  /*  U - matrix (Q+I)^T, where Q is infitisimal generator (^T - transposed)
      S - eigenvectors of U
      p0 - distribution in the beginning
      p8 - stable (end) distribution
      R - tmp rates matrix (transposed) */
  T *U=NULL, *S=NULL, *p0=NULL, *p8=NULL, *R=NULL;
  double *tmpReadRates=NULL;
  int dim;

  Barparser *bp = new Barparser(opt,dim,E);

  switch (opt->method) {
    case 'F': dim = bp->ParseInfile(opt->INFILE, opt->RATFILE, &tmpReadRates); break;
    case 'I':
      if (opt->binrates){
        dim = bp->MxReadBinRates(opt->RATFILE, &tmpReadRates, opt->n, opt->max_decrease);
      }
      else {
        dim = bp->ParseRatesFile(opt->RATFILE, &tmpReadRates);
        // shorten the matrix from dimension my_dim to nstates:
        dim = bp->MxShorten(&tmpReadRates, opt->n, dim,  opt->max_decrease);
      }

      if (opt->INFILE) {
        bp->ParseBarfile(opt->INFILE, &Data);
      }
      if (dim == 0) {
        fprintf(stderr, "ERROR: Rate file empty!\n");
        globalParameters->free_gengetopt();
        exit(EXIT_FAILURE);
      }

      R = convertRates<T>(tmpReadRates,(size_t)dim);
      // visualize the graph:
      if (opt->vis_file) {
        bp->VisualizeRates(opt->vis_file, R, Data, dim);
      }
      break;
    case 'A': dim = bp->ParseBarfile (opt->INFILE, &Data); break;
  }

  delete bp;

  Calc<T> *calc = new Calc<T>(globalParameters,E,dim);

  // here we create the "almighty" rate matrix U which is actually only matrix needed for whole program
  U  = (T*)calc->MxBar2Matrix(Data, R);

  Calccpp *calccpp = new Calccpp(opt, E);

  // rescale if minimal_rate has been set:
  if (opt->hard_rescale != 1.0) calccpp->MxRescaleH(U, dim, opt->hard_rescale);
  else if (opt->minimal_rate > 0.0) calccpp->MxRescale(U, dim, opt->minimal_rate);

  // multiply if times was set:
  if (opt->times != 1.0) calccpp->MxTimes(U, dim, opt->times);

  delete calccpp;

  // create initial probability vector
  calc->MxStartVec(&p0);

  // check for ergodicity + adjust to that
  dim = calc->MxEgro(&U, &p0, dim);
  if (opt->want_verbose) calc->MxPrint(U, "Ergodic U", 'm');

  // allocate space for other matrices
  calc->MxGetSpace(&p8);

  if (!opt->quiet) fprintf(stderr, "Time to initialize: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
  clck1 = clock();

  // calculate equilibrium distribution
  if(opt->method == 'F')
    calc->MxEqDistrFULL (E, p8);
  else {
    if (opt->method == 'A')
      calc->MxEqDistr (Data, &p8);
    else 
      calc->MxEqDistrFromDetailedBalance(U, &p8);
    //MxEqDistrFromLinSys(U, &p8);
  }

  // write the equilibrium if we should
  if (opt->equil_file) {
    FILE *equil = fopen(opt->equil_file, "w");
    if (equil) {
      calc->MxFPrint(p8, "Equilibrium distribution", 'v', equil, 1);
      fclose(equil);
    }
  }

  // first passage times computation
  if(opt->fpt) {
    // output file?
    FILE *fpt = stderr;
    if (opt->fpt_file!=NULL) {
      fpt = fopen(opt->fpt_file, "w");
      if (fpt==NULL) {
        if (!opt->quiet) fprintf(stderr, "Cannot open file \"%s\" for FPT output, using \"stderr\" instead\n", opt->fpt_file);
        fpt = stderr;
      }
    }
    // all or only one state?
    if ((opt->fpt_num == -1) || opt->absrb) { // all ftp's
      calc->MxFPT(U, p8, fpt);
    } else {                                // only to state opt.ftp_num
      T *res = calc->MxFPTOneState(U, opt->fpt_num-1);  // starting from 0
      if (res != NULL) {
        if (fpt==stderr) fprintf(fpt, "First passage times to state number %d:\n", opt->fpt_num);
        int k;
        for (k = 0; k < dim; k++) fprintf(fpt, "%15.7f ", (double)res[k]);
        if (fpt==stderr) fprintf(fpt,"\n---\n");
        else fprintf(fpt, "\n");
        delete[] res; //free(res);
      }
    }
    fclose(fpt);
    if (!opt->quiet) fprintf(stderr, "Time to compute fpt: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  int exponentError = 0;
  if (!opt->just_sh) {
    // diagonalization + iteration
    if(opt->matexp) {
      exponentError = calc->MxExponent(p0,p8,U);
      if(exponentError){

        calc->MxMemoryCleanUp();
        delete calc;
        delete[] U;
        delete[] p8;
        delete[] p0;
        free(Data);

        exit(EXIT_FAILURE);
      }
    }
    else {
      if(opt->rrecover) /* read pre-calculated eigenv{alues,ectors} */
        calc->MxRecover(&S, p8);
      else {
        clck1 = clock();
        calc->MxDiagonalize (U, &S, p8);
        if (!opt->quiet) fprintf(stderr, "Time to diagonalize: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
        clck1 = clock();
      }
      calc->MxIterate (p0, p8, S);
      if (!opt->quiet) fprintf(stderr, "Time to iterate: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    }
  }


  if (opt->pini != NULL) free(opt->pini);
  // clean up the memory
  calc->MxMemoryCleanUp();
  delete calc;
  delete[] U;
  delete[] p8;
  delete[] p0;
  free(Data);
  return 0;
}


int
main (int argc, char **argv)
{
  Globals *globalParameters = Globals::initGlobals();
  globalParameters->parse_commandline(argc, argv);

  int error;

  switch(globalParameters->opt.mpackMethod){
    #ifdef WITH_MPACK_DOUBLE
    case MPACK_DOUBLE:
      error = treekin_main_precision<double>(globalParameters);
      break;
    #endif
    #ifdef WITH_MPACK_LD
    case MPACK_LD:
      error = treekin_main_precision<long double>(globalParameters);
      break;
    #endif
    #ifdef WITH_MPACK___FLOAT128
    case MPACK_FLOAT128:
      error = treekin_main_precision<__float128>(globalParameters);
      break;
    #endif
    #ifdef WITH_MPACK_MPFR
    case MPACK_MPFR:
      mpfr::mpreal::set_default_prec(globalParameters->opt.mpackMethod_Bits);
      error = treekin_main_precision<mpreal_castable>(globalParameters);
      break;
    #endif
    #ifdef WITH_MPACK_DD
    case MPACK_DD:
      error = treekin_main_precision<dd_real_castable>(globalParameters);
      break;
    #endif
    #ifdef WITH_MPACK_QD
    case MPACK_QD:
      error = treekin_main_precision<qd_real_castable>(globalParameters);
      break;
    #endif
    #ifdef WITH_MPACK_GMP
    case MPACK_GMP:
      mpf_set_default_prec(globalParameters->opt.mpackMethod_Bits);
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



