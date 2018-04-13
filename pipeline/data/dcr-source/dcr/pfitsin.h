/*========================================================*/
/*                                                        */
/*  pfitsin.h       version 3.7.0         2015.05.11      */
/*                                                        */
/*  Copyright (C) 2010-2015 by Wojtek Pych, CAMK PAN      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS I/O library.                             */
/*                                                        */
/*========================================================*/

#include "pfitshead.h"
#include "swap.h"

#define EPSILON 1.0e-15

extern int   read_FITS_image(FILE *, char **, const int, void **);
extern float *read_FITS_1Dfile(char *, int *, char ***, int *);
extern float **read_FITS_2Dfile(char *, int *, char ***, int *, int *);
extern float *read_FITS_2D1file(char *, int *, char ***, int *, int *);

/*** END ***/
