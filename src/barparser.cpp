/*=================================================================*/
/*=   barparser.c                                                 =*/
/*=   routines for reading bar-files and other input for treekin  =*/
/*=                                                               =*/
/*=                 (c) Michael Thomas Wolfinger                  =*/
/*=                      mtw@tbi.univie.ac.at                     =*/
/*=                             treekin                           =*/
/*=================================================================*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>

#include "barparser.h"

const int Barparser::LMINBASE   = 100;
const int Barparser::ARRAYSIZE  = 10000;

/*==*/
char *
Barparser::my_getline(FILE *fp)
{
  char s[512], *line, *cp;

  line = NULL;
  do {
    if (fgets(s, 512, fp) == NULL)
      break;

    cp = strchr(s, '\n');
    if (cp != NULL)
      *cp = '\0';

    if (line == NULL)
      line = (char *)calloc(strlen(s) + 1, sizeof(char));
    else
      line = (char *)realloc(line, strlen(s) + strlen(line) + 1);

    strcat(line, s);
  } while (cp == NULL);
  return line;
}


/*==*/
/* parses file saddles.txt from barriers */
int
Barparser::ParseSaddleFile(TypeDegSaddle **my_saddle)
{
  int           count = 0, size = 100, newsize, p = 0;
  char          *line = NULL, *remember_line;
  const char    *filename = "saddles.txt";
  FILE          *my_file;
  TypeDegSaddle *temp;

  temp = (TypeDegSaddle *)calloc(size, sizeof(TypeDegSaddle));

  my_file = fopen(filename, "r");
  while ((line = my_getline(my_file)) != NULL) {
    remember_line = line;
    if ((count + 1) == size) {
      newsize = 2 * size;
      temp    = (TypeDegSaddle *)realloc(temp, newsize * sizeof(TypeDegSaddle));
      if (!temp) {
        fprintf(stderr, "Could not realloc tmp structure array in ParseSaddleFile\n");
        exit(EXIT_FAILURE);
      }

      size = newsize;
    }

    sscanf(line, "%lg %4d %*s %n", &temp[count].energy, &temp[count].cc, &p);
    line                += p;
    temp[count].list[0] = 0; /* list[0] contains # of lmins connected by that saddle */
    while (*line != '\0') {
      if (temp[count].list[0] == 99) {
        fprintf(stderr, " too many lmins in saddles.txt\n");
        exit(EXIT_FAILURE);
      }

      p = 0;
      sscanf(line, "%d %n", &temp[count].list[temp[count].list[0] + 1], &p);
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
int
Barparser::ParseBarfile(FILE    *fp,
                        BarData **lmin)
{
  char    *line = NULL, *tmpseq = NULL, *p, sep[] = " ";
  int     count = 0, v = 0;
  int     size = LMINBASE;
  BarData *tmp;

  if ((fp == NULL) || (fcntl(fileno(fp), F_GETFD) == -1)) {
    fprintf(stderr, "ERROR: Barrier file handle is either NULL or already closed!");
    exit(EXIT_FAILURE);
  }

  tmp     = (BarData *)calloc(LMINBASE, sizeof(BarData));
  tmpseq  = (char *)calloc(500, sizeof(char));
  *lmin   = NULL;
  line    = my_getline(fp); /* get sequence from first line */
  if (sscanf(line, "%s", tmpseq) == 0) {
    fprintf(stderr, "could not get sequence from first line of barfile\n");
    exit(EXIT_FAILURE);
  }

  if (_opt->sequence)
    free(_opt->sequence);

  _opt->sequence = tmpseq;
  free(line);
  count = 0;
  for (line = my_getline(fp); line != NULL; line = my_getline(fp)) {
    if (count >= _opt->n)
      break;

    if (count >= size) {
      tmp   = (BarData *)realloc(tmp, (size + LMINBASE) * sizeof(BarData));
      size  += LMINBASE;
    }

    p = strtok(line, sep);
    v = 1;
    sscanf(p, "%d", &tmp[count].number);
    while (p != NULL) {
      p = strtok(NULL, sep);
      if (p == NULL)
        break;

      if (isgraph(*p) && *p != '-' && !isdigit(*p)) {
        /* secondary structures */
        switch ((int)*p) {
          case 40:    /* '(' */
          case 46:    /* '.' */
          case 70:    /* 'F' */
          case 126:   /* '~' */
            continue;
          default:
            fprintf(stderr, "%c does not seem to be a RNA secondary structure nor a SAW\n", *p);
            continue;
        }
      }

      v++;
      switch (v) {
        case 2:
          sscanf(p, "%lg", &tmp[count].energy);
          break;
        case 3:
          sscanf(p, "%d", &tmp[count].father);
          break;
        case 4:
          sscanf(p, "%lg", &tmp[count].ediff);
          break;
        case 5:
          sscanf(p, "%lg", &tmp[count].bsize);
          break;
        case 6:
          sscanf(p, "%lg", &tmp[count].fathers_bsize);
          break;
        case 7:
          sscanf(p, "%lg", &tmp[count].F);
          break;
        case 8:
          sscanf(p, "%lg", &tmp[count].Gr_bsize);
          break;
        case 9:
          sscanf(p, "%lg", &tmp[count].FGr);
          break;
        default:
          break;
      }
    }
    if (line != NULL)
      free(line);

    count++;
  }

  if (line != NULL)
    free(line);

  tmp   = (BarData *)realloc(tmp, count * sizeof(BarData));
  *lmin = tmp;
  fclose(fp);
  return count;
}


int
Barparser::MxReadBinRates(FILE    *rate_file,
                          double  **rate_mx,
                          int     nstates,
                          int     max)
{
  int dimension = 0;
  int ref;

  if ((rate_file == NULL) || (fcntl(fileno(rate_file), F_GETFD) == -1)) {
    fprintf(stderr, "ERROR: Rates matrix (binary) file handle is either NULL or already closed!");
    exit(EXIT_FAILURE);
  }

  /* read dimension from file */
  ref = fread(&dimension, sizeof(int), 1, rate_file);
  double *data = (double *)calloc(dimension * dimension, sizeof(double));
  ref = fread((void *)data, sizeof(double), dimension * dimension, rate_file);

  *rate_mx = data;
  fclose(rate_file);

  /* decrease dimension */
  dimension = MxShorten(rate_mx, nstates, dimension, max);

  return dimension;
}


/*==*/
int
Barparser::ParseInfile(FILE   *infile_fp,
                       FILE   *mr_FP,
                       double **microrates)
{
  char    *line = NULL;
  int     a, dim, indx = 0, l = 0, as1 = ARRAYSIZE, as2 = ARRAYSIZE;
  double  *tmp_mr   = NULL; /* temp matrix containing microrates */
  SubInfo *tmp_subI = NULL; /* temp array for info on energies of all states */

  typedef struct {
    int     i, j;
    double  rate;
  } mr_t;
  mr_t *mr = NULL;

  /* FIRST: read (subopt)-infile with energies & gradient basins */
  tmp_subI = (SubInfo *)calloc(as1, sizeof(SubInfo));
  assert(tmp_subI != NULL);
  line            = my_getline(infile_fp);
  _opt->sequence  = (char *)calloc(2048, sizeof(char));
  assert(_opt->sequence != NULL);
  sscanf(line, "%s %*f", _opt->sequence);
  free(line);
  while ((line = my_getline(infile_fp)) != NULL && (indx < _opt->n)) {
    if (indx + 1 >= as1) {
      as1       *= 2;
      tmp_subI  = (SubInfo *)realloc(tmp_subI, as1 * sizeof(SubInfo));
    }

    int res = sscanf(line, "%*s %lg %d %*d", &tmp_subI[indx].energy, &tmp_subI[indx].ag);
    if (res == 0)
      res = sscanf(line, "%*d %*s %lg %d %*d", &tmp_subI[indx].energy, &tmp_subI[indx].ag);

    /*fprintf(stderr, "%s %lg %d %d\n", line, tmp_subI[indx].energy, tmp_subI[indx].ag, res); */
    if (tmp_subI[indx].ag > l)
      l = tmp_subI[indx].ag;

    /* fix allegiance status: */
    if (tmp_subI[indx].ag > 0)
      tmp_subI[indx].ag--;

    free(line);
    indx++;
  }
  dim   = indx;
  indx  = 0;
  fclose(infile_fp);

  /* SECOND: read microrates file */
  mr = (mr_t *)calloc(as2, sizeof(mr_t));
  assert(mr != NULL);
  while ((line = my_getline(mr_FP)) != NULL) {
    if (indx + 1 >= as2) {
      as2 *= 2;
      mr  = (mr_t *)realloc(mr, as2 * sizeof(mr_t));
    }

    if (sscanf(line, "%d %d %lf %*d", &mr[indx].i, &mr[indx].j, &mr[indx].rate) == 3)
      /*fprintf(stderr,"indx %i\n", indx); */
      indx++;

    free(line);
  }
  fclose(mr_FP);

  /* THIRD: fill up transition matrix */
  tmp_mr = (double *)calloc(dim * dim, sizeof(double));
  assert(tmp_mr != NULL);
  for (a = 0; a < indx; a++) {
    int i, j;
    i                               = mr[a].i;
    j                               = mr[a].j;
    tmp_mr[dim * (i - 1) + (j - 1)] = mr[a].rate;
    tmp_mr[dim * (j - 1) + (i - 1)] = 1.;
  }
  *microrates = tmp_mr;
  free(mr);
  _lmins  = l;
  _E      = tmp_subI;
  fprintf(stderr, "dimension = %d, lmins = %d \n", dim, _lmins);
  return dim;
}


/*==*/
int
Barparser::ParseRatesFile(FILE    *rates_FP,
                          double  **Raten)
{
  char    *cp = NULL, *raten_line = NULL;
  /*char *suffix = "rates.out"; */
  double  rate, *tmp_rates = NULL;

  /* valid rate file? */
  if ((rates_FP == NULL) || (fcntl(fileno(rates_FP), F_GETFD) == -1)) {
    fprintf(stderr, "ERROR: Rates matrix file handle is either NULL or already closed!");
    exit(EXIT_FAILURE);
  }

  /* actual reading: */
  int   my_dim = 0;
  /* get length of line */
  char  *tmp_line = my_getline(rates_FP);
  raten_line = (char *)calloc(strlen(tmp_line) + 1, sizeof(char));
  strcpy(raten_line, tmp_line);

  char  *p = strtok(tmp_line, " \t\n");
  while (p && sscanf(p, "%lf", &rate) == 1) {
    /*fprintf(stderr, "%d-%lf-\n", my_dim, rate); */
    my_dim++;
    p = strtok(NULL, " \t\n");
  }
  free(tmp_line);

  /* allocate space */
  tmp_rates = (double *)calloc(my_dim * my_dim, sizeof(double));
  assert(tmp_rates != NULL);

  /* read! */
  int i = 0, j = 0, read = 0;
  while (raten_line != NULL && i < my_dim) {
    cp = raten_line;
    while (cp != NULL && sscanf(cp, "%lf%n", &rate, &read) == 1 && j < my_dim) {
      tmp_rates[my_dim * j + i] = rate;
      cp                        += read;
      j++;
    }
    j = 0;
    i++;
    free(raten_line);
    raten_line = my_getline(rates_FP);
  }

  /* check dimensions: */
  if (j == 0) {
    j = my_dim - 1;
    i--;
  }                              /* if last line empty */

  if (i != my_dim - 1 || j != my_dim - 1) {
    if (!_opt->quiet) {
      fprintf(stderr,
              "WARNING: dimensions are corrupted lines: %d(%d); rows: %d(%d)\n",
              j,
              my_dim - 1,
              i,
              my_dim - 1);
    }
  }

  /*output rates */
  *Raten = tmp_rates;

  fclose(rates_FP);
  return my_dim;
}


/* End of file */
