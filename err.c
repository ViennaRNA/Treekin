
/**************************************************************************
**
** Copyright (C) 1993 David E. Stewart & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/*
  File with basic error-handling operations
  Based on previous version on Zilog
  System 8000 setret() etc.
  Ported to Pyramid 9810 late 1987
  */

#include	<stdio.h>
#include	<setjmp.h>
#include	<ctype.h>
#include        "err.h"
#include	<signal.h>



#define		FALSE	0
#define		TRUE	1

#define	EF_EXIT		0
#define	EF_ABORT	1
#define	EF_JUMP		2
#define	EF_SILENT	3

/* The only error caught in this file! */
#define	E_SIGNAL	16

static	char	*err_mesg[] = {
  "unknown error",			    /* 0 */
  "sizes of objects don't match",	    /* 1 */
  "index out of bounds",		    /* 2 */
  "can't allocate memory",		    /* 3 */
  "singular matrix",			    /* 4 */
  "matrix not positive definite",	    /* 5 */
  "incorrect format input",		    /* 6 */
  "bad input file/device",		    /* 7 */
  "NULL objects passed",		    /* 8 */
  "matrix not square",			    /* 9 */
  "object out of range",		    /* 10 */
  "can't do operation in situ for non-square matrix",   /* 11 */
  "can't do operation in situ",		    /* 12 */
  "excessive number of iterations",	    /* 13 */
  "convergence criterion failed",	    /* 14 */
  "bad starting value",			    /* 15 */
  "floating exception",			    /* 16 */
  "internal inconsistency (data structure)",/* 17 */
  "unexpected end-of-file",		    /* 18 */
  "shared vectors (cannot release them)",   /* 19 */  
  "negative argument",			    /* 20 */
  "cannot overwrite object",                /* 21 */
  "breakdown in iterative method"           /* 22 */
};

#define	MAXERR	(sizeof(err_mesg)/sizeof(char *))

static char *warn_mesg[] = {
   "unknown warning",				  /* 0 */
   "wrong type number (use macro TYPE_*)",	  /* 1 */
   "no corresponding mem_stat_mark",		  /* 2 */
   "computed norm of a residual is less than 0",  /* 3 */
   "resizing a shared vector"			  /* 4 */
};

#define MAXWARN  (sizeof(warn_mesg)/sizeof(char *))



#define	MAX_ERRS	100

jmp_buf	restart;


/* array of pointers to lists of errors */

typedef struct {
   char **listp;    /* pointer to a list of errors */
   unsigned len;    /* length of the list */
   unsigned warn;   /* =FALSE - errors, =TRUE - warnings */
}  Err_list;

static Err_list     err_list[ERR_LIST_MAX_LEN] = {
 {err_mesg,MAXERR,FALSE},	/* basic errors list */
 {warn_mesg,MAXWARN,TRUE}	/* basic warnings list */
};


static int err_list_end = 2;   /* number of elements in err_list */


/* other local variables */

static	int	err_flag = EF_EXIT, num_errs = 0, cnt_errs = 1;

/* set_err_flag -- sets err_flag -- returns old err_flag */
int	set_err_flag(flag)
int	flag;
{
   int	tmp;
   
   tmp = err_flag;
   err_flag = flag;
   return tmp;
}


/* ev_err -- reports error (err_num) in file "file" at line "line_num" and
   returns to user error handler;
   list_num is an error list number (0 is the basic list 
   pointed by err_mesg, 1 is the basic list of warnings)
 */
int ev_err(char *file, int err_num, int line_num, char *fn_name, int list_num) {
  int	num;
  
  if ( err_num < 0 ) err_num = 0;
  
  if (list_num < 0 || list_num >= err_list_end ||
      err_list[list_num].listp == (char **)NULL) {
    fprintf(stderr, "\n Not (properly) attached list of errors: list_num = %d\n", list_num);
    fprintf(stderr," Call \"err_list_attach\" in your program\n");
    if ( ! isatty(fileno(stdout)) ) {
      fprintf(stderr, "\n Not (properly) attached list of errors: list_num = %d\n",
	      list_num);
      fprintf(stderr," Call \"err_list_attach\" in your program\n");
    }
    printf("\nExiting program\n");
    exit(0);
  }
  
  num = err_num;
  if ( num >= err_list[list_num].len ) num = 0;
  
  if ( cnt_errs && ++num_errs >= MAX_ERRS ) {	/* too many errors */
    fprintf(stderr,"\n\"%s\", line %d: %s in function %s()\n",
	    file,line_num,err_list[list_num].listp[num],
	    isascii(*fn_name) ? fn_name : "???");
    if ( ! isatty(fileno(stdout)) )
      fprintf(stdout,"\n\"%s\", line %d: %s in function %s()\n",
	      file,line_num,err_list[list_num].listp[num],
	      isascii(*fn_name) ? fn_name : "???");
    printf("Sorry, too many errors: %d\n",num_errs);
    printf("Exiting program\n");
    exit(0);
  }
  if ( err_list[list_num].warn )
    switch ( err_flag )   {
    case EF_SILENT: break;
    default:
      fprintf(stderr,"\n\"%s\", line %d: %s in function %s()\n\n",
	      file,line_num,err_list[list_num].listp[num],
	      isascii(*fn_name) ? fn_name : "???");
      if ( ! isatty(fileno(stdout)) )
	fprintf(stdout,"\n\"%s\", line %d: %s in function %s()\n\n",
		file,line_num,err_list[list_num].listp[num],
		isascii(*fn_name) ? fn_name : "???");
      break;
    }
  else
    switch ( err_flag )   {
    case EF_SILENT:
      longjmp(restart,(err_num==0)? -1 : err_num);
      break;
    case EF_ABORT:
      fprintf(stderr,"\n\"%s\", line %d: %s in function %s()\n",
	      file,line_num,err_list[list_num].listp[num],
	      isascii(*fn_name) ? fn_name : "???");
      if ( ! isatty(fileno(stdout)) )
	fprintf(stdout,"\n\"%s\", line %d: %s in function %s()\n",
		file,line_num,err_list[list_num].listp[num],
		isascii(*fn_name) ? fn_name : "???");
      abort();
      break;
    case EF_JUMP:
      fprintf(stderr,"\n\"%s\", line %d: %s in function %s()\n",
	      file,line_num,err_list[list_num].listp[num],
	      isascii(*fn_name) ? fn_name : "???");
      if ( ! isatty(fileno(stdout)) )
	fprintf(stdout,"\n\"%s\", line %d: %s in function %s()\n",
		file,line_num,err_list[list_num].listp[num],
		isascii(*fn_name) ? fn_name : "???");
      longjmp(restart,(err_num==0)? -1 : err_num);
      break;
    default:
      fprintf(stderr,"\n\"%s\", line %d: %s in function %s()\n\n",
	      file,line_num,err_list[list_num].listp[num],
	      isascii(*fn_name) ? fn_name : "???");
      if ( ! isatty(fileno(stdout)) )
	fprintf(stdout,"\n\"%s\", line %d: %s in function %s()\n\n",
		file,line_num,err_list[list_num].listp[num],
		isascii(*fn_name) ? fn_name : "???");
      
      break;
    }
  
  /* ensure exit if fall through */
  if ( ! err_list[list_num].warn )  exit(0);
  
  return 0;
}
