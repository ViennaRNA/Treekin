/*=================================================================*/
/*=   barparser.h                                                 =*/
/*=   header file for parsing bar-files etc.                      =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 11:43:07 mtw>          =*/
/*=   $Id: barparser.h,v 1.8 2006/03/15 11:08:15 mtw Exp $    =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#ifndef _BARPARSER_H_
#define _BARPARSER_H_

#define SADDLE_LIST 400

/* structures */ 
typedef struct _BarData { /* structure for bar-file */ 
  float energy;
  float ediff;
  int father;
  int number;
  float bsize;          /* # of structures in this lmin  == basin size */
  float eff_bsize;      /* bsize without children */
  float fathers_bsize;  /* bsize of father at merging-time */
  float F;              /* free energy of lmin i and its children */ 
  float F_eff;          /* free engery of lmin (alone) */
  float Z;              /* partition function of lmin */
  float Z_eff;          /* Z of lmin (alone) */
  float Gr_bsize;       /* gradient basin size */
  float FGr;            /* F of gradient basin */ 
} BarData;

typedef struct _inData { /* structure for data.out file for full process */ 
  int i;
  int j;
  double rate;
} InData;

typedef struct _TypeDegSaddle{ /* structure for processing degeneracy */
  float energy;   /* energy of saddle */
  int cc;        /* # of structures in that particular connected component */
  int list[100];     /* string containing lmins which are connected by the saddle */
} TypeDegSaddle;

typedef struct _SubInfo{ /* struct containing info on all subopt structures */
  double energy;         /* just needed in FULL process */
  int    ag;             /* associated gradient basin */
} SubInfo;

SubInfo *E; 
 
/* functions */ 
int  ParseBarfile (FILE *fp, BarData **lmin);
int  ParseInfile(FILE *fp, InData **InD);
int  ParseSaddleFile(TypeDegSaddle **my_saddle);
void ParseRatesFile(double **Raten, int dim);
#endif

/* End of file */
