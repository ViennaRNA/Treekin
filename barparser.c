/* barparser.c */
/* Last changed Time-stamp: <2003-07-09 00:12:14 mtw> */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <barparser.h>
#include "globals.h"

#define LMINBASE 100

/*  static char rcsid[] = "$Id: barparser.c,v 1.1 2003/07/14 07:42:20 mtw Exp $"; */

static char *getline(FILE *fp);

/*==*/
int ParseInfile(FILE *fp, InData **transition, double **En, int **lmin_nr_so, int **assoc_gradbas){

  char *line = NULL, *line_tr = NULL, *f_trace = "tracelmins.out", *grad_bas = "assoc_gradbas.out";
  int dimensione, a, b, c, x, newsize = 500, mem_inc = 1000, uhu;
  int *tmp_trace;     /* tmp array 4 tracelmins */
  int *tmp_gradbas;   /* tmp array 4 associated gradient basins */
  InData *tmp;        /* tmp array 4 rates between two states */ 
  double *Etmp;       /* tmp array 4 Energies */ 
  extern int in_nr;
  FILE *lmins;        /* file pointer 4 tracelmins.out */
  FILE *gradient_basins;  /* file pointer 4 associated gradient basins */
  
  /* read first line */
  line = getline(fp);
  sscanf(line, "%d %*s", &dimensione);
  /* HIER AUCH NOCH opt.sequence herausscannen */
  
  tmp =  (InData *) calloc (dimensione, sizeof(InData));
  if(tmp == NULL){
    fprintf(stderr, "tmp could not be allocated\n");
    exit(888);
  }
  Etmp = (double *) calloc (dimensione, sizeof(double));
  if(Etmp == NULL){
    fprintf(stderr, "Etmp could not be allocated\n");
    exit(888);
  }
  
  /* read energies of the different states into array Etmp */ 
  for(x = 0; x < dimensione; x++){
    int o, p;
    line = getline(stdin);
    sscanf(line, "%d %d %lf", &o, &p, &Etmp[x]);
    if(o != p) {
      fprintf(stderr, "error while reading energies! %d != %d\n", o, p);
      exit(555);
    }
  }

  a = 0;
  uhu = dimensione;
  /* read states and transition between them */ 
  while(line != NULL){
    line = getline(stdin);
    if (line == NULL) break;
    if(a+1 >= uhu){ /* realloc tmp array */ 
      newsize += mem_inc;
      tmp =  (InData *) realloc (tmp, newsize*sizeof(InData));
      if(tmp == NULL){
	fprintf(stderr, "realloc of tmp failed\n");
	exit(888);
      }
      uhu = newsize;
    }
    sscanf(line, "%d %d %lf", &tmp[a].j, &tmp[a].i, &tmp[a].rate);
    a++;
  }

  tmp =  (InData *) realloc (tmp, a*sizeof(InData));
  in_nr = a-1;
  free(line);
  
  /* ================================= */
  /* >>> begin read tracelmins.out <<< */
  tmp_trace = (int*) calloc (LMINBASE, sizeof(int));
  if(tmp_trace == NULL){
    fprintf(stderr, "tmp_trace could not be allocated\n");
    exit(888);
  }

  b = 1;
  lmins = fopen(f_trace, "r+");    /* file pointer for tracelmins.out */
  line_tr = getline(lmins);        /* read first line containing info stuff */
  while((line_tr = getline(lmins)) != NULL){    /* read which lmin is which # in subopt */ 
    sscanf(line_tr, "%*d %6d", &tmp_trace[b]);
    b++;
    if(b >= LMINBASE){
      fprintf(stderr, " read more than 100 lmins..too many!\n");
      exit(777);
    }
  }
  tmp_trace[0] = b-1;
  tmp_trace = (int*) realloc (tmp_trace, b*sizeof(int));
  fclose(lmins);
  /* tmp_trace[0] = # of lmins which are now associated with a # from subopt */ 
  /* >>> end  read tracelmins.out <<< */
  /* ================================ */
  
  /* ==================================== */
  /* >>> begin read assoc_gradbas.out <<< */
  tmp_gradbas = (int *) calloc (dimensione, sizeof(int));
  if(tmp_gradbas == NULL){
    fprintf(stderr, "tmp_gradbas could not be allocated\n");
    exit(888);
  }
  
  c = 0;
  gradient_basins = fopen(grad_bas, "r+");  /* file pointer for assoc_gradbas.out */
  line_tr = getline(gradient_basins);       /* read first line containing info stuff */
  while((line_tr = getline(gradient_basins)) != NULL){ /* read the gradient basin of each entry from subopt */  
    sscanf(line_tr, "%*d %5d", &tmp_gradbas[c]);
    c++;
  }
  if(c != dimensione){
    fprintf(stderr, " read more than 100 lmins..too many!\n");
    exit(777);
  }
  fclose(gradient_basins);
  /* >>> end read assoc_gradbas.out <<< */
  /* ================================== */ 

  fprintf(stderr, "read %d items, dimension = %d\n", a, dimensione);

  *transition = tmp;
  *En = Etmp;
  *lmin_nr_so = tmp_trace;
  *assoc_gradbas = tmp_gradbas;
  
  /* fclose(lmins); */
/*   fclose(gradient_basins); */
/*   free(line_tr); */
  
  return dimensione;
}

/*==*/
int ParseBarfile( FILE *fp, TypeBarData **lmin){

  char *line = NULL, *p, sep[] = " ";
  int count = 0, v = 0, i;
  int size = LMINBASE;
 /*   float mfe; */
  TypeBarData *tmp;

  tmp = (TypeBarData *) calloc (LMINBASE, sizeof(TypeBarData));
  *lmin = NULL;
  
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
      if( *p == '.' || *p == '(' ) continue; /* skip secondary structures */
      else if ( *p == 'A' || *p == 'U' || *p == 'G' || *p == 'C' || *p == '~'){
	int count = 1;   /* read sequence from bar-file (saddle of 1st lmin) */
	while ( *(p+count) != '\0')
	  count++;
	opt.sequence = (char *) calloc(count, sizeof(char));
	for (i = 0; i < count; i++)
	  *(opt.sequence + i) = *(p + i);  
	continue;
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
  /* if we have an absorbing state lower the state's energy by 15kcal*/
  /*   if(opt.absrb > 0) */
  /*      tmp[opt.absrb-1].FGr -= 15; */
  
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




