/*========================================================*/
/*                                                        */
/*  dcr.h           2015.08.21      version 2.2.1         */
/*                                                        */
/*  Copyright (C) 2006-2015 by Wojtek Pych, CAMK PAN      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/*  Detect cosmic rays on a CCD frame.                    */
/*                                                        */
/*--------------------------------------------------------*/
/*                                                        */
/*========================================================*/

#define LLEN 200

typedef struct
{
  char  *iname,
        *orname,
        *ocname;
  short vopt;
  int   xr, yr,
        rlr, rur,
        gr,
        npass,
        dispaxis;
  float ths;
} PARAMS;

void  usage(void);
int   readpar(char *, PARAMS *);
void  calc_mean(float **, int, int, double *, double *);
int   calc_submean(PARAMS, float **, int, int, int, int, double *, double *);
float max(float **data, int, int, int *, int *);
void  minmax(float **data, int, int, int, int, float *, float *);
int   make_hist(int, int, float **, int, int, float, int, float, int *);
int   detect(PARAMS, float **, int, int, int, int, int, int, char **);
int   make_map(PARAMS, float **, int, int, char **);
int   clean(PARAMS, float **, int, int, double, char **, float **);
int   modify_header(int, char ***, PARAMS);

/*** END ***/
