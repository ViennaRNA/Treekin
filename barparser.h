/* barparser.h */
/* Last changed Time-stamp: <2001-08-28 17:11:02 mtw> */
#ifndef _BARPARSER_H_
#define _BARPARSER_H_
#define SADDLE_LIST 400

/*  static char rcsid[] = "$Id: barparser.h,v 1.1 2003/07/14 07:42:20 mtw Exp $"; */

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

/* functions */ 
int ParseBarfile (FILE *fp, TypeBarData **lmin);
int ParseInfile(FILE *fp, InData **InD, double **En, int **lmin_nr_so, int **assoc_gradbas);
int ParseSaddleFile(TypeDegSaddle **my_saddle);
#endif

/* End of file */
