/* barparser.c */
/* Last changed Time-stamp: <2003-09-16 18:00:49 mtw> */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <barparser.h>
#include "globals.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define LMINBASE 100

/*  static char rcsid[] = "$Id: barparser.c,v 1.9 2003/09/25 13:51:01 mtw Exp $"; */

static char *getline(FILE *fp);

/*==*/
int ParseInfile(FILE *fp, InData **transition, int *lmins){

  char *line = NULL, *line_tr = NULL, *grad_bas = "assoc_gradbas.out";
  int dimensione, a, c, l,x, newsize = 5000, mem_inc = 10000, uhu;
  InData *tmp;        /* tmp array 4 rates between two states */ 
  SubInfo *tmp_subI;  /* tmp array 4 info on energies of all subopts */
  extern int in_nr;
  FILE *gb_FP;        /* file pointer 4 associated gradient basins */
  
  /* read first line */
  line = getline(fp);
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
  for(x = 0; x < dimensione; x++){
    int o, p;
    line = getline(stdin);
    sscanf(line, "%d %d %lf", &o, &p, &tmp_subI[x].energy);
    if(o != p) {
      fprintf(stderr, "error while reading energies from data.out: %d != %d\n", o, p);
      exit(555);
    }
    if(line != NULL) free(line);
  }

  a = 0;
  uhu = dimensione;
  /* read states and transition between them */ 
  while(line != NULL){
    line = getline(stdin);
    if (line == NULL) break;
    if(a+1 >= uhu){ /* realloc tmp array */ 
      newsize += mem_inc;
      fprintf(stderr,"realloc: new size: %d\n", newsize);
      tmp =  (InData *) realloc (tmp, newsize*sizeof(InData));
      if(tmp == NULL){
	fprintf(stderr, "realloc of tmp failed\n");
	exit(888);
      }
      uhu = newsize;
    }
    sscanf(line, "%d %d %lf", &tmp[a].j, &tmp[a].i, &tmp[a].rate);
  /*   fprintf(stderr, "a: %d\n", a); */
    a++;
    if(line != NULL) free(line);
  }

  tmp =  (InData *) realloc (tmp, a*sizeof(InData));
  in_nr = a-1;
  if(line) free(line);
  
  /* ==================================== */
  /* >>> begin read assoc_gradbas.out <<< */
  c = 0; l = 0;
  gb_FP = fopen(grad_bas, "r+");  /* file pointer for assoc_gradbas.out */
  line_tr = getline(gb_FP);       /* read first line containing info stuff */
  if(line_tr != NULL) free(line_tr);
  while((line_tr = getline(gb_FP)) != NULL){ /* read the gradient basin of each entry from subopt */  
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
  *lmins = l;
  E = tmp_subI;
  fprintf(stderr, "read %d items, dimension = %d, lmins = %d \n", a, dimensione, *lmins);
  
  if(line_tr != NULL) free(line_tr); 
  
  return dimensione;
}

/*==*/
void ParseRatesFile(double **Raten, int dim){
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
int ParseBarfile( FILE *fp, TypeBarData **lmin){

  char *line = NULL, *tmpseq = NULL, *p, sep[] = " ";
  int count = 0, v = 0;
  int size = LMINBASE;
  TypeBarData *tmp;

  tmp = (TypeBarData *) calloc (LMINBASE, sizeof(TypeBarData));
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
      tmp = (TypeBarData *) realloc (tmp, size + LMINBASE);
      size += LMINBASE;
    }
    p = strtok(line, sep);
    v = 1;
    sscanf(p, "%d",  &tmp[count].number);
    while(p != NULL){
      p = strtok(NULL, sep);
      if (p == NULL) break; 
      if (*p == '.' || *p == '(' ) continue; /* skip secondary structures */
      else if (*p == '~') continue;
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

  tmp = (TypeBarData *) realloc (tmp, count*sizeof(TypeBarData));
  *lmin = tmp;
  return count;
}

/*==*/
/* parses file saddles.txt from barriers */
int ParseSaddleFile(TypeDegSaddle **my_saddle) {

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
char *getline(FILE *fp) { 
  
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




