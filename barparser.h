/* barparser.h */
/* Last changed Time-stamp: <2003-09-12 12:36:38 mtw> */
#ifndef _BARPARSER_H_
#define _BARPARSER_H_
#define SADDLE_LIST 400

/*  static char rcsid[] = "$Id: barparser.h,v 1.5 2003/09/25 13:51:01 mtw Exp $"; */

/* structures */ 
typedef struct _TypeBarData { /* structure for bar-file */ 
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
} TypeBarData;

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
int  ParseBarfile (FILE *fp, TypeBarData **lmin);
int  ParseInfile(FILE *fp, InData **InD, int *lmins);
int  ParseSaddleFile(TypeDegSaddle **my_saddle);
void ParseRatesFile(double **Raten, int dim);
#endif

/* End of file */
