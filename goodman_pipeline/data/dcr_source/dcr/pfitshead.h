/*========================================================*/
/*                                                        */
/*  pfitshead.h     version 3.7.0         2015.05.11      */
/*                                                        */
/*  Copyright (C) 2008-2015 by Wojtek Pych, CAMK PAN      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS I/O library.                             */
/*                                                        */
/*  Definition of the FITS:                               */
/*  Hanisch, R. J.; et al. 2001, A&A, 376, 359-380        */
/*                                                        */
/*========================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Some FITS magic numbers */
#define RECORD_SIZE   2880
#define RECORD_CARDS  36
#define CARD_SIZE     80
#define KEYWORD_SIZE  8
#define VALUE_SIZE    70

extern int get_FITS_key(int, char **, char *, char*);
extern int FITS_image_size(const int, char **);
extern int read_FITS_header(FILE *, char ***);
extern int read_FITS_headers(FILE *, int **, char ****, int **);
extern int del_header_card(int, char ***, char *);
extern int add_header_card(char *, int *, char ***);
extern int write_FITS_header(FILE *, const int, char **);

/*** END ***/
