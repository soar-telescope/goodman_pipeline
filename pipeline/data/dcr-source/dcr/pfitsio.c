/*========================================================*/
/*                                                        */
/*  pfitsio.c       version 3.4.0         2015.05.11      */
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
/*                                                        */
/* Important assumptions:                                 */
/* short:	2 bytes = 16 bits                         */
/* int:		4 bytes = 32 bits                         */
/* long:	8 bytes = 64 bits                         */
/*                                                        */
/*========================================================*/

#include <math.h>

#include "pfitsio.h"

/*--------------------------------------------------------*/
double scale(const int size, float *data, short *idata, double *fitsbzero)
{
        int     i;
        float   min,
                max;
        double  FITS_bscale,
                FITS_bzero;

/** scalling factors **/
  min = max = data[0];
  for (i=1; i<size; i++)
  {
    if (data[i] > max) max=data[i];
    if (data[i] < min) min=data[i];
  }
  FITS_bzero=(max-min)/2.0;
  FITS_bscale=(max-min)/65535.0;

  for (i=0; i<size; i++)
    idata[i]=(short)rint((data[i]-FITS_bzero)/FITS_bscale);

  *fitsbzero=FITS_bzero;

  return(FITS_bscale);
}
/*--------------------------------------------------------*/
int write_FITS_1Dimage(FILE *outf, const int nx, const int bytepix,
  void *data)
{
        char  *tmp;
        int   size,   /* image size in bytes  */
              ls;

  if (needswap() != 0) swap_bytes(data, nx, bytepix);

  if (fwrite(data, bytepix, nx, outf) != nx)
  {
    perror("\n\tERROR! write_FITS_1Dimage(): fwrite(data)");
    return(EXIT_FAILURE);
  }

  size=nx*bytepix;
  if ((ls=size%RECORD_SIZE) != 0) ls=RECORD_SIZE-ls;
  if ((tmp=(char *)malloc(ls)) == NULL)
  {
    perror("\n\tERROR! write_FITS_1Dimage(): malloc(tmp)");
    return(EXIT_FAILURE);
  }

  memset(tmp, 0, ls);
  if (fwrite(tmp, 1, ls, outf) != ls)
  {
    printf("\n\tERROR! write_FITS_1Dimage(): writing zero padding\n");
    return(EXIT_FAILURE);
  }

  free(tmp);

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int write_FITS_2Dimage(FILE *outf, const int nx, const int ny,
  const int bytepix, void **data)
{
        char  *tmp;
        int   i,
              size,   /* image size in bytes  */
              ls;

  if (needswap() != 0)
  {
    for (i=0; i<ny; i++)
    {
      swap_bytes(data[i], nx, bytepix);
      if (fwrite(data[i], bytepix, nx, outf) != nx)
      {
        perror("\n\tERROR! write_FITS_2Dimage(): fwrite(data[i])");
        return(EXIT_FAILURE);
      }
    }
  }
  else
  {
    for (i=0; i<ny; i++)
    {
      if (fwrite(data[i], bytepix, nx, outf) != nx)
      {
        perror("\n\tERROR! write_FITS_2Dimage(): fwrite(data[i])");
        return(EXIT_FAILURE);
      }
    }
  }

  size=nx*ny*bytepix;
  if ((ls=size%RECORD_SIZE) != 0) ls=RECORD_SIZE-ls;
  if ((tmp=(char *)malloc(ls)) == NULL)
  {
    perror("\n\tERROR! write_FITS_2Dimage(): malloc(tmp)");
    return(EXIT_FAILURE);
  }

  memset(tmp, 0, ls);
  if (fwrite(tmp, 1, ls, outf) != ls)
  {
    printf("\n\tERROR! write_FITS_2Dimage(): writing zero padding\n");
    return(EXIT_FAILURE);
  }

  free(tmp);

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int write_FITS_1Dfile(char *oname, const int ncards, char **header,
  const int nx, const int bytepix, void *data)
{
        FILE  *outf;

  if ((outf=fopen(oname, "w")) == NULL)
  { perror(oname); return(EXIT_FAILURE); }

  if (write_FITS_header(outf, ncards, header) != ncards)
    return(EXIT_FAILURE);
  if (write_FITS_1Dimage(outf, nx, bytepix, data) != EXIT_SUCCESS)
    return(EXIT_FAILURE);

  if (fclose(outf) != 0)
    printf("\n\tERROR! write_FITS_1Dfile(): cannot close %s\n", oname);

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int write_FITS_2Dfile(char *oname, const int ncards, char **header,
  const int nx, const int ny, const int bytepix, void **data)
{
        FILE  *outf;

  if ((outf=fopen(oname, "w")) == NULL)
  { perror(oname); return(EXIT_FAILURE); }

  if (write_FITS_header(outf, ncards, header) != ncards)
    return(EXIT_FAILURE);
  if (write_FITS_2Dimage(outf, nx, ny, bytepix, data) != EXIT_SUCCESS)
    return(EXIT_FAILURE);

  if (fclose(outf) != 0)
    printf("\n\tERROR! write_FITS_2Dfile(): cannot close %s\n", oname);

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
/*** END ***/
