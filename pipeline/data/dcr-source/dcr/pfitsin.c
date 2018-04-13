/*========================================================*/
/*                                                        */
/*  pfitsin.c       version 3.7.0         2015.05.11      */
/*                                                        */
/*  Copyright (C) 2010-2015 by Wojtek Pych, CAMK PAN      */
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pfitsin.h"

/*--------------------------------------------------------*/
int read_FITS_image(FILE *inf, char **header, const int ncards,
  void **buf)
{
        void  *ptr;                     /* data table pointer   */
        char  keyword[KEYWORD_SIZE+1],  /* header entry keyword */
              value[VALUE_SIZE];
        int   i,              /* loop numerator               */
              naxes,          /* image dimension              */
              *naxis,         /* table of sizes of axes       */
              bitpix,         /* bits per pixel               */
              bytepix,        /* bytes per pixel              */
              pixnum;         /* number of image pixels       */

  memset(value, ' ', VALUE_SIZE);

  if (get_FITS_key(ncards, header, "NAXIS", value) == -1) return(0);
  if (sscanf(value, "%d", &naxes) < 1)
  {
    printf("\n\tERROR! read_FITS_image(): reading NAXIS\n");
    return(0);
  }

  if ((naxis=(int *)calloc(naxes, sizeof(int))) == NULL)
  { perror("\n\tERROR! read_FITS_image(): calloc(naxis)"); return(0); }

  for (i=0; i<naxes; i++)
  {
    if (snprintf(keyword, KEYWORD_SIZE, "NAXIS%d", i+1) < 6)
    {
      printf("\n\tERROR! read_FITS_image(): setting NAXIS_ keyword\n");
      free(naxis);
      return(0);
    }
    if (get_FITS_key(ncards, header, keyword, value) == -1)
    {
      free(naxis);
      return(0);
    }
    if (sscanf(value, "%d", &naxis[i]) < 1)
    {
      printf("\n\tERROR! read_FITS_image(): reading %s\n", keyword);
      free(naxis);
      return(0);
    }
  }

  pixnum=1;
  for (i=0; i<naxes; i++) pixnum*=naxis[i];

  free(naxis);

  if (get_FITS_key(ncards, header, "BITPIX", value) == -1) return(0);
  if (sscanf(value, "%d", &bitpix) < 1)
  {
    printf("\n\tERROR! read_FITS_image(): reading BITPIX from header\n");
    return(0);
  }
  bytepix=abs(bitpix)/8;

  if ((ptr=(void *)malloc(pixnum*bytepix)) == NULL)
  { perror("\n\tERROR! read_FITS_image(): malloc(ptr)"); return(0); }

  if (fread(ptr, bytepix, pixnum, inf) != pixnum)
  {
    printf("\n\tERROR! read_FITS_image(): reading the image\n");
    return(0);
  }

  if (needswap() != 0) swap_bytes(ptr, pixnum, bytepix);

  *buf=ptr;

  return(pixnum);
}
/*--------------------------------------------------------*/
float *read_FITS_1Dfile(char *iname, int *ncards, char ***header, int *n1)
{
        void    *buf;         /* FITS image buffer        */
        char    value[VALUE_SIZE],  /* card value buffer  */
		**lhead;            /* pointer to header  */
        int     i,            /* loop numerator           */
		hs,           /* number of header cards   */
                bitpix,       /* FITS BITPIX value        */
                naxes,        /* number of image axes     */
                nx,
                nbytes;       /* length of the image buffer */
        float   *ldata;       /* output data pointer      */
        double  FITS_bzero,
                FITS_bscale;
        FILE    *inf;

  memset(value, ' ', VALUE_SIZE);

  if ((inf=fopen(iname, "r")) == NULL)
  { perror(iname); return(NULL); }

  if ((hs=read_FITS_header(inf, &lhead)) < RECORD_CARDS)
  {
    printf("\n\tERROR! read_FITS_1Dfile(): incorrect header length\n");
    return(NULL);
  }

  if ((nbytes=read_FITS_image(inf, lhead, hs, &buf)) < 1)
  {
    printf("\n\tERROR! read_FITS_1Dfile(): incorrect image length\n");
    return(NULL);
  }

  if (fclose(inf) != 0)
    printf("\n\tWARNING! read_FITS_1Dfile: cannot close %s\n", iname);

  if (get_FITS_key(hs, lhead, "NAXIS", value) == -1)
  {
    printf("\n\tERROR! read_FITS_1Dfile(%s): NAXIS not found in the header\n",
            iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &naxes) < 1)
  {
    printf("\n\tERROR! read_FITS_1Dfile(): reading NAXIS\n");
    return(NULL);
  }

  if (naxes != 1)
  {
    printf("\n\tERROR! read_FITS_1Dfile(%s): NAXIS= %d (expected 1)\n",
            iname, naxes);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "NAXIS1", value) == -1)
  {
    printf("\n\tERROR! read_FITS_1Dfile(%s): NAXIS1 not found in the header\n",
            iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &nx) < 1)
  {
    printf("\n\tERROR! read_FITS_1Dfile(): reading NAXIS1\n");
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "BITPIX", value) == -1)
  {
    printf("\n\tERROR! read_FITS_1Dfile(%s): BITPIX not found in the header\n",
            iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &bitpix) < 1)
  { printf("\n\tERROR! read_FITS_1Dfile(): reading BITPIX\n"); return(NULL); }

  if ((ldata=(float *)calloc(nx, sizeof(float))) == NULL)
  { perror("\n\tERROR! read_FITS_1Dfile(): calloc(ldata)"); return(NULL); }

  if (get_FITS_key(hs, lhead, "BZERO", value) == -1)  FITS_bzero=0.0;
  else
  {
    if (sscanf(value, "%lg", &FITS_bzero) < 1)
    {
      printf("\n\tERROR! read_FITS_1Dfile(): reading BZERO\n");
      return(NULL);
    }
  }

  if (get_FITS_key(hs, lhead, "BSCALE", value) == -1) FITS_bscale=1.0;
  else
  {
    if (sscanf(value, "%lg", &FITS_bscale) < 1)
    {
      printf("\n\tERROR! read_FITS_1Dfile(): reading BSCALE\n");
      return(NULL);
    }
  }

  switch(bitpix)
  {
    case   8: for (i=0; i<nx; i++)  ldata[i]=(float)((unsigned char *)buf)[i];
              break;
    case  16: for (i=0; i<nx; i++)  ldata[i]=(float)((short *)buf)[i];
              break;
    case  32: for (i=0; i<nx; i++)  ldata[i]=(float)((int *)buf)[i];
              break;
    case -32: for (i=0; i<nx; i++)  ldata[i]=((float *)buf)[i];
              break;
    case -64: for (i=0; i<nx; i++)  ldata[i]=(float)((double *)buf)[i];
              break;
    default:  printf("\n\tERROR! read_FITS_1Dfile(%s): bitpix= %d",
                      iname, bitpix);
              printf(" not supported\n");
              return(NULL);
  }

  free(buf);

  if (fabs(FITS_bscale-1.0) > EPSILON)
    for (i=0; i<nx; i++) ldata[i]*=FITS_bscale;
  if (fabs(FITS_bzero) > EPSILON)
    for (i=0; i<nx; i++) ldata[i]+=FITS_bzero;

  *n1=nx;
  *ncards=hs;
  *header=lhead;

  return(ldata);
}
/*--------------------------------------------------------*/
float **read_FITS_2Dfile(char *iname, int *ncards, char ***header,
  int *n1, int *n2)
{
        void    *buf;
        char    value[VALUE_SIZE],
		**lhead;
        int     i, j,
                hs,
                bitpix,
                naxes,
                nx, ny,
                npix;
        float   **ldata;
        double  FITS_bzero,
                FITS_bscale;
        FILE    *inf;

  memset(value, ' ', VALUE_SIZE);

  if ((inf=fopen(iname, "r")) == NULL) { perror(iname); return(NULL); }

  if ((hs=read_FITS_header(inf, &lhead)) < RECORD_CARDS)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): incorrect header length\n",
      iname);
    return(NULL);
  }

  if ((npix=read_FITS_image(inf, lhead, hs, &buf)) < 1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): incorrect image length\n", iname);
    return(NULL);
  }

  if (fclose(inf) != 0)
    printf("\n\tWARNING! read_FITS_2Dfile(%s): cannot close file\n", iname);

  if (get_FITS_key(hs, lhead, "NAXIS", value) == -1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): NAXIS not found in the header\n",
      iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &naxes) < 1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): reading NAXIS\n", iname);
    return(NULL);
  }

  if (naxes != 2)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): NAXIS= %d (expected 2)\n",
      iname, naxes);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "NAXIS1", value) == -1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): NAXIS1 not found in header\n",
      iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &nx) < 1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): reading NAXIS1\n", iname);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "NAXIS2", value) == -1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): NAXIS2 not found in header\n",
      iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &ny) < 1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): reading NAXIS2\n", iname);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "BITPIX", value) == -1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): BITPIX not found in header\n",
      iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &bitpix) < 1)
  {
    printf("\n\tERROR! read_FITS_2Dfile(%s): reading BITPIX\n", iname);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "BZERO", value) == -1)  FITS_bzero=0.0;
  else
  {
    if (sscanf(value, "%lg", &FITS_bzero) < 1)
    {
      printf("\n\tERROR! read_FITS_2Dfile(%s): reading BZERO\n", iname);
      return(NULL);
    }
  }

  if (get_FITS_key(hs, lhead, "BSCALE", value) == -1) FITS_bscale=1.0;
  else
  {
    if (sscanf(value, "%lg", &FITS_bscale) < 1)
    {
      printf("\n\tERROR! read_FITS_2Dfile(%s): reading BSCALE\n", iname);
      return(NULL);
    }
  }

  if ((ldata=(float **)calloc(ny, sizeof(float *))) == NULL)
  { perror("\n\tERROR! read_FITS_2Dfile(): calloc(ldata)"); return(NULL); }

  for (i=0; i<ny; i++)
  {
    if ((ldata[i]=(float *)calloc(nx, sizeof(float))) == NULL)
    {
      perror("\n\tERROR! read_FITS_2Dfile(): calloc(ldata[i])");
      return(NULL);
    }
  }

  switch(bitpix)
  {
    case   8: for (i=0; i<ny; i++)
                for (j=0; j<nx; j++)
                  ldata[i][j]=(float)((unsigned char *)buf)[i*nx+j];
              break;
    case  16: for (i=0; i<ny; i++)
                for (j=0; j<nx; j++)
                  ldata[i][j]=(float)((short *)buf)[i*nx+j];
              break;
    case  32: for (i=0; i<ny; i++)
                for (j=0; j<nx; j++)
                  ldata[i][j]=(float)((int *)buf)[i*nx+j];
              break;
    case -32: for (i=0; i<ny; i++)
                for (j=0; j<nx; j++)
                  ldata[i][j]=((float *)buf)[i*nx+j];
              break;
    case -64: for (i=0; i<ny; i++)
                for (j=0; j<nx; j++)
                  ldata[i][j]=(float)((double *)buf)[i*nx+j];
              break;
    default:  printf("\n\tERROR! read_FITS_2Dfile(%s): bitpix= %d\n",
                iname, bitpix);
              return(NULL);
  }

  free(buf);

  if (fabs(FITS_bscale-1.0) > EPSILON)
    for (i=0; i<ny; i++)
      for (j=0; j<nx; j++)
        ldata[i][j]*=FITS_bscale;

  if (fabs(FITS_bzero) > EPSILON)
    for (i=0; i<ny; i++)
      for (j=0; j<nx; j++)
        ldata[i][j]+=FITS_bzero;

  *n1=nx;
  *n2=ny;
  *header=lhead;
  *ncards=hs;

  return(ldata);
}
/*--------------------------------------------------------*/
float *read_FITS_2D1file(char *iname, int *ncards, char ***header,
  int *n1, int *n2)
{
        void    *buf;
        char    value[VALUE_SIZE],
		**lhead;
        int     i,
                hs,
                bitpix,
                naxes,
                nx,
                ny,
                npix;
        float   *ldata;
        double  FITS_bzero,
                FITS_bscale;
        FILE    *inf;

  memset(value, ' ', VALUE_SIZE);

  if ((inf=fopen(iname, "r")) == NULL) { perror(iname); return(NULL); }

  if ((hs=read_FITS_header(inf, &lhead)) < RECORD_CARDS)
  {
    printf("\n\tERROR! read_FITS_2D1file(): incorrect header length\n");
    return(NULL);
  }

  if ((npix=read_FITS_image(inf, lhead, hs, &buf)) < 1)
  {
    printf("\n\tERROR! read_FITS_2D1file(): incorrect image length\n");
    return(NULL);
  }

  if (fclose(inf) != 0)
    printf("\n\tWARNING! read_FITS_2D1file(%s): cannot close file\n", iname);

  if (get_FITS_key(hs, lhead, "NAXIS", value) == -1)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): NAXIS not found in the header\n",
      iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &naxes) < 1)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): reading NAXIS\n", iname);
    return(NULL);
  }

  if (naxes != 2)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): NAXIS= %d (expected 2)\n",
      iname, naxes);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "NAXIS1", value) == -1)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): NAXIS1 not found in header\n",
      iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &nx) < 1)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): reading NAXIS1\n", iname);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "NAXIS2", value) == -1)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): NAXIS2 not found in header\n",
      iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &ny) < 1)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): reading NAXIS2\n", iname);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "BITPIX", value) == -1)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): BITPIX not found in header\n",
      iname);
    return(NULL);
  }

  if (sscanf(value, "%d", &bitpix) < 1)
  {
    printf("\n\tERROR! read_FITS_2D1file(%s): reading BITPIX\n", iname);
    return(NULL);
  }

  if (get_FITS_key(hs, lhead, "BZERO", value) == -1)  FITS_bzero=0.0;
  else
  {
    if (sscanf(value, "%lg", &FITS_bzero) < 1)
    {
      printf("\n\tERROR! read_FITS_2D1file(%s): reading BZERO\n", iname);
      return(NULL);
    }
  }

  if (get_FITS_key(hs, lhead, "BSCALE", value) == -1) FITS_bscale=1.0;
  else
  {
    if (sscanf(value, "%lg", &FITS_bscale) < 1)
    {
      printf("\n\tERROR! read_FITS_2D1file(%s): reading BSCALE\n", iname);
      return(NULL);
    }
  }

  if ((ldata=(float *)calloc(npix, sizeof(float))) == NULL)
  {
    perror("\n\tERROR! read_FITS_2D1file(): calloc(ldata)");
    return(NULL);
  }

  switch(bitpix)
  {
    case   8: for (i=0; i<npix; i++)  ldata[i]=(float)((unsigned char *)buf)[i];
              break;
    case  16: for (i=0; i<npix; i++)  ldata[i]=(float)((short *)buf)[i];
              break;
    case  32: for (i=0; i<npix; i++)  ldata[i]=(float)((int *)buf)[i];
              break;
    case -32: for (i=0; i<npix; i++)  ldata[i]=((float *)buf)[i];
              break;
    case -64: for (i=0; i<npix; i++)  ldata[i]=(float)((double *)buf)[i];
              break;
    default:  printf("\n\tERROR! read_FITS_2D1file(%s): bitpix= %d\n",
                iname, bitpix);
              return(NULL);
  }

  free(buf);

  if (fabs(FITS_bscale-1.0) > EPSILON)
    for (i=0; i<npix; i++)  ldata[i]*=FITS_bscale;

  if (fabs(FITS_bzero) > EPSILON)
    for (i=0; i<npix; i++)  ldata[i]+=FITS_bzero;

  *n1=nx;
  *n2=ny;
  *header=lhead;
  *ncards=hs;

  return(ldata);
}
/*--------------------------------------------------------*/
/*** END ***/
