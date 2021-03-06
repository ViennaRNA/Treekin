#ifndef _BARDATA_H_
#define _BARDATA_H_

/* structures */
typedef struct _BarData {
  /* structure for bar-file */
  double  energy;
  double  ediff;
  int     father;
  int     number;
  double  bsize;          /* # of structures in this lmin  == basin size */
  double  eff_bsize;      /* bsize without children */
  double  fathers_bsize;  /* bsize of father at merging-time */
  double  F;              /* free energy of lmin i and its children */
  double  F_eff;          /* free engery of lmin (alone) */
  double  Z;              /* partition function of lmin */
  double  Z_eff;          /* Z of lmin (alone) */
  double  Gr_bsize;       /* gradient basin size */
  double  FGr;            /* F of gradient basin */
} BarData;

typedef struct _TypeDegSaddle {
  /* structure for processing degeneracy */
  double  energy;     /* energy of saddle */
  int     cc;         /* # of structures in that particular connected component */
  int     list[100];  /* string containing lmins which are connected by the saddle */
} TypeDegSaddle;

typedef struct _SubInfo {
  /* struct containing info on all subopt structures */
  double  energy;         /* just needed in FULL process */
  int     ag;             /* associated gradient basin */
} SubInfo;

#endif

/* End of file */
