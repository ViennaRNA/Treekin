/*=================================================================*/
/*=   barparser.c                                                 =*/
/*=   routines for reading bar-files and other input for treekin  =*/
/*=   ---------------------------------------------------------   =*/
/*=   Last changed Time-stamp: <2006-03-15 10:34:59 mtw>          =*/
/*=   $Id: barparser.c,v 1.15 2006/03/15 11:08:15 mtw Exp $    =*/
/*=   ---------------------------------------------------------   =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "barparser.h"
#include "globals.h"
#include "limits.h"
#include "ctype.h"

#define LMINBASE 100

static char *getline(FILE *fp);

/*==*/
int
ParseInfile(FILE *fp, InData **transition)
{
  char *line = NULL, *line_tr = NULL, *grad_bas = "assoc_gradbas.out";
  int dimensione, c, i, l, newsize = 5000, limit;
  InData *tmp;        /* tmp array 4 rates between two states */ 
  SubInfo *tmp_subI;  /* tmp array 4 info on energies of all subopts */
  FILE *gb_FP;        /* file pointer 4 associated gradient basins */
  
  line = getline(fp);  /* read first line */
  sscanf(line, "%d %*s", &dimensione);
  /* HIER AUCH NOCH opt.sequence herausscannen */
  if(line != NULL) free(line);
  
  tmp =  (InData *) calloc (dimensione, sizeof(InData));
  if(tmp == NULL){
    fprintf(stderr, "tmp could not be allocated\n");
    exit(888);
  }
  tmp_subI = (SubInfo *) calloc (dimensione, sizeof(SubInfo));
  if(tmp_subI == NULL){
    fprintf(stderr, "struct tmp_subI could not be allocated in barparser.c\n");
    exit(888);
  }
  
  /* read energies of the different states into array tmp_subI from stdin (data.out) */ 
  for(i = 0; i < dimensione; i++){
    int o, p;
    line = getline(stdin);
    sscanf(line, "%d %d %lf", &o, &p, &tmp_subI[i].energy);
    if(o != p) {
      fprintf(stderr, "error while reading energies from data.out: %d != %d\n", o, p);
      exit(555);
    }
    if(line != NULL) free(line);
  }

  i = 0;
  limit = dimensione;
  /* read states and transition between them */ 
  while(line != NULL){
    line = getline(stdin);
    if (line == NULL) break;
    if(i+1 >= limit){ /* realloc tmp array */ 
      newsize *= 3;
      fprintf(stderr,"realloc data array: new size is %d\n", newsize);
      tmp =  (InData *) realloc (tmp, newsize*sizeof(InData));
      if(tmp == NULL){
	fprintf(stderr, "realloc of data array in barparser failed\n");
	exit(888);
      }
      limit = newsize;
    }
    sscanf(line, "%d %d %lf", &tmp[i].j, &tmp[i].i, &tmp[i].rate);
    i++;
    if(line != NULL) free(line);
  }

  tmp =  (InData *) realloc (tmp, i*sizeof(InData));
  in_nr = i-1;
  if(line) free(line);
  
  /* ==================================== */
  /* >>> begin read assoc_gradbas.out <<< */
  c = 0; l = 0;
  gb_FP = fopen(grad_bas, "r+");  /* file pointer for assoc_gradbas.out */
  line_tr = getline(gb_FP);       /* read first line containing info stuff */
  if(line_tr != NULL) free(line_tr);
  /* read the gradient basin of each entry from subopt */
  while((line_tr = getline(gb_FP)) != NULL){
    sscanf(line_tr, "%*d %5d", &tmp_subI[c].ag);
    if (tmp_subI[c].ag > l) l = tmp_subI[c].ag;
    c++;
    if(line_tr != NULL) free(line_tr);
  }
  if(c != dimensione){
    fprintf(stderr, " read more lines from gradient basin file than our dim is!\n");
    exit(777);
  }
  if(gb_FP) fclose(gb_FP);
  /* >>> end read assoc_gradbas.out <<< */
  /* ================================== */ 

  *transition = tmp;
  lmins = l;
  E = tmp_subI;
  fprintf(stderr, "read %d items, dimension = %d, lmins = %d \n", i, dimensione, lmins);
  
  if(line_tr != NULL) free(line_tr); 
  
  return dimensione;
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
    exit(888);
  }
  rate = (double *)calloc(dim,sizeof(double));
  if (rate == NULL){
    fprintf(stderr, "could not allocate rate in ParseRatesFile\n");
    exit(888);
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
    exit(888);
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
	exit(111);
      }
      size = newsize;
    }
    sscanf(line, "%f %4d %*s %n", &temp[count].energy, &temp[count].cc, &p);
    line += p;
    temp[count].list[0] = 0; /* list[0] contains # of lmins connected by that saddle */
    while( *line != '\0') {
      if(temp[count].list[0] == 99) {
	fprintf (stderr, " too many lmins in saddles.txt\n");
	exit(111);
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




