/*=================================================================*/
/*=   barparser.c                                                 =*/
/*=   routines for reading bar-files and other input for treekin  =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-05-25 18:24:46 mtw>          =*/
/*=   $Id: barparser.c,v 1.20 2006/06/09 15:49:35 mtw Exp $       =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "barparser.h"
#include "globals.h"
#include "limits.h"
#include "ctype.h"

#define LMINBASE 100
#define ARRAYSIZE 10000

static char *getline(FILE *fp);

/*==*/
int
ParseInfile(FILE *infile_fp, double **microrates)
{
  char *line = NULL, *mrfile = "microrates.out";
  int a,dim,indx=0, l=0, as1 = ARRAYSIZE, as2 = ARRAYSIZE;
  double *tmp_mr=NULL;    /* temp matrix containing microrates */
  SubInfo *tmp_subI=NULL; /* temp array for info on energies of all states */
  FILE *mr_FP=NULL;       /* file pointers 4 "microrates.out" */
  typedef struct {
    int i,j;
    double rate;
  } mr_t;
  mr_t *mr=NULL;

  /* FIRST: read (subopt)-infile with energies & gradient basins */
  tmp_subI = (SubInfo *)calloc(as1, sizeof(SubInfo));
  assert(tmp_subI != NULL);
  line = getline(infile_fp);
  opt.sequence = (char *)calloc(2048,sizeof(char));
  assert(opt.sequence != NULL);
  sscanf(line, "%s %*f", opt.sequence);
  free(line);
  while((line=getline(infile_fp)) != NULL){
    if(indx+1 >= as1){
      as1 *= 2;
      tmp_subI = (SubInfo *)realloc(tmp_subI,as1*sizeof(SubInfo));
    }
    sscanf(line, "%*s %f %d %*d", &tmp_subI[indx].energy, &tmp_subI[indx].ag);
    if (tmp_subI[indx].ag > l) l = tmp_subI[indx].ag;
    free(line);
    indx++;
  }
  dim = indx;
  indx = 0;
  
  /* SECOND: read microrates.out */
  mr_FP = fopen(mrfile, "r+");
  mr = (mr_t *)calloc(as2, sizeof(mr_t));
  assert(mr != NULL);
  while((line=getline(mr_FP)) != NULL){
    if(indx+1 >= as2){
      as2 *= 2;
      mr = (mr_t *)realloc(mr, as2*sizeof(mr_t));
    }
    sscanf(line, "%d %d %lf %*d", &mr[indx].i, &mr[indx].j, &mr[indx].rate);
    free(line);
    indx++;
  }
  fclose(mr_FP);

  /* THIRD: fill up transition matrix */
  tmp_mr  = (double *)calloc(dim*dim, sizeof(double));
  assert(tmp_mr != NULL);
  for(a=0;a<indx;a++){
    int i,j;
    i = mr[a].i;
    j = mr[a].j;
    tmp_mr[dim*(i-1)+(j-1)]=mr[a].rate;
    tmp_mr[dim*(j-1)+(i-1)]=1;
  }
  *microrates = tmp_mr;
 
  lmins = l;
  E = tmp_subI;
  fprintf(stderr, "dimension = %d, lmins = %d \n", dim, lmins);
  
  return dim;
}

/*==*/
void
ParseRatesFile(double **Raten, int dim)
{
  int i = 0, j = 0, read =0;
  char *cp, *raten_line = NULL, *rate_file = "rates.out";
  double *rate, *tmp_rates;

  tmp_rates = (double *)calloc(dim*dim,sizeof(double));
  if (tmp_rates == NULL){
    fprintf(stderr, "could not allocate tmp_rates in ParseRatesFile\n");
    exit(EXIT_FAILURE);
  }
  rate = (double *)calloc(dim,sizeof(double));
  if (rate == NULL){
    fprintf(stderr, "could not allocate rate in ParseRatesFile\n");
    exit(EXIT_FAILURE);
  }
  
  opt.RATENFILE = fopen(rate_file, "r+");
  while((raten_line = getline(opt.RATENFILE)) != NULL){ 
    cp = raten_line;
    while(sscanf(cp,"%lf%n", rate+j,&read) == 1){
      tmp_rates[dim*j+i] = *(rate+j);
      cp+=read;
      j++;
    }
    j=0;
    i++;
    memset(rate, 0, (size_t)(dim*sizeof(double)));
    free(raten_line);
  }

  *Raten = tmp_rates;
  fclose(opt.RATENFILE);
  free(rate);
  return;
}

/*==*/
int
ParseBarfile( FILE *fp, BarData **lmin)
{
  char *line = NULL, *tmpseq = NULL, *p, sep[] = " ";
  int count = 0, v = 0;
  int size = LMINBASE;
  BarData *tmp;

  tmp = (BarData *) calloc (LMINBASE, sizeof(BarData));
  tmpseq = (char *) calloc (500, sizeof(char));
  *lmin = NULL;
  line = getline(fp); /* get sequence from first line */
  if(sscanf(line, "%s", tmpseq) == 0){
    fprintf(stderr, "could not get sequence from first line of barfile\n");
    exit(EXIT_FAILURE);
  }
  opt.sequence = tmpseq;
  free(line);
  for (count = 0, line = getline(fp); line != NULL; count++, line = getline(fp)) {
    if (count >= size) {
      tmp = (BarData *) realloc (tmp, (size+LMINBASE)*sizeof(BarData));
      size += LMINBASE;
    }
    p = strtok(line, sep);
    v = 1;
    sscanf(p, "%d",  &tmp[count].number);
    while(p != NULL){
      p = strtok(NULL, sep);
      if (p == NULL) break; 
      if (isgraph(*p) && *p != '-' && !isdigit(*p)){ /* secondary structures */
	switch ((int)*p){
	case 40:      /* '(' */
	case 46:      /* '.' */
	case 70:      /* 'F' */
	case 126:     /* '~' */
	  continue;
	default:
	  fprintf(stderr, "%c does not seem to be a RNA secondary structure nor a SAW\n", *p);
	  continue;
	}
      }
   
      v++;
      switch (v){
      case 2:
	sscanf(p, "%f", &tmp[count].energy);
	break;
      case 3:
	sscanf(p, "%d", &tmp[count].father);
	break;
      case 4:
	sscanf(p, "%f", &tmp[count].ediff);
	break;
      case 5:
	sscanf(p, "%e", &tmp[count].bsize);
	break;
      case 6:
	sscanf(p, "%e", &tmp[count].fathers_bsize);
	break;
      case 7:
	sscanf(p, "%e", &tmp[count].F);
	break;
      case 8:
	sscanf(p, "%e", &tmp[count].Gr_bsize);
	break;
      case 9:
	sscanf(p, "%e",  &tmp[count].FGr);
	break;
      default:
	break;
      }
    }
    if (line != NULL) free(line);
  }
   
  if (line != NULL) free(line);

  tmp = (BarData *) realloc (tmp, count*sizeof(BarData));
  *lmin = tmp;
  return count;
}

/*==*/
/* parses file saddles.txt from barriers */
int
ParseSaddleFile(TypeDegSaddle **my_saddle)
{
  int count = 0, size = 100, newsize , p = 0;
  char *line = NULL, *remember_line, *filename = "saddles.txt";
  FILE *my_file;
  TypeDegSaddle *temp;
  temp = (TypeDegSaddle *) calloc (size, sizeof(TypeDegSaddle));
   
  my_file = fopen(filename, "r");
  while ((line = getline(my_file)) != NULL) {
    remember_line = line;
    if((count+1) == size){
      newsize = 2 * size;
      temp = (TypeDegSaddle *) realloc (temp, newsize*sizeof(TypeDegSaddle));
      if(!temp) {
	fprintf(stderr, "Could not realloc tmp structure array in ParseSaddleFile\n");
	exit(EXIT_FAILURE);
      }
      size = newsize;
    }
    sscanf(line, "%f %4d %*s %n", &temp[count].energy, &temp[count].cc, &p);
    line += p;
    temp[count].list[0] = 0; /* list[0] contains # of lmins connected by that saddle */
    while( *line != '\0') {
      if(temp[count].list[0] == 99) {
	fprintf (stderr, " too many lmins in saddles.txt\n");
	exit(EXIT_FAILURE);
      }
      p = 0;
      sscanf(line, "%d %n", &temp[count].list[temp[count].list[0]+1], &p);
      line += p;
      temp[count].list[0]++;
    }
    count++;
    free(remember_line); 
  }
  *my_saddle = temp;
  fclose(my_file);
  return count; /*return # of saddles read */
}

/*==*/
char *
getline(FILE *fp)
{ 
  char s[512], *line, *cp; 
  line = NULL;  
  do { 
    if(fgets(s, 512, fp) == NULL) break;
    cp = strchr(s, '\n');
    if(cp != NULL) *cp = '\0'; 
    if(line == NULL) line = (char *) calloc(strlen(s) + 1, sizeof(char)); 
    else line = (char *) realloc(line, strlen(s) + strlen(line) + 1); 
    strcat (line, s);
  } while (cp == NULL);
  return (line);
}
/* End of file */




