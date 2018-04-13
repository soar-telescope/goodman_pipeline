/*========================================================*/
/*                                                        */
/*  pfitsio.h       version 3.3.0         2013.03.25      */
/*                                                        */
/*  Copyright (C) 2007-2013 by Wojtek Pych, CAMK PAN      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS I/O library.                             */
/*                                                        */
/*========================================================*/

#include "pfitsin.h"

#define EPSILON 1.0e-15

extern double scale(const int, float *, short *, double *);
extern int write_FITS_1Dimage(FILE *, const int, const int, void *);
extern int write_FITS_2Dimage(FILE *, const int, const int, const int,
  void **);
extern int write_FITS_1Dfile(char *, const int, char **, const int, const int,
  void *);
extern int write_FITS_2Dfile(char *, const int, char **, const int, const int,
  const int, void **);

/*** END ***/
