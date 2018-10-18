/*========================================================*/
/*                                                        */
/*  dcr.c           2015.08.21      version 2.2.1         */
/*                                                        */
/*  Copyright (C) 2008-2015 by Wojtek Pych, CAMK PAN      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/*  Detect cosmic rays on a single CCD frame.             */
/*                                                        */
/*--------------------------------------------------------*/
/*                                                        */
/*  Only standard FITS files supported:                   */
/*  - no unsigned 16-bit and 32-bit images                */
/*                                                        */
/*========================================================*/

/***************************************************************************/
/*                                                                         */
/*   This program is free software; you can redistribute it and/or modify  */
/*   it under the terms of the GNU General Public License as published by  */
/*   the Free Software Foundation; either version 2 of the License, or     */
/*   (at your option) any later version.                                   */
/*                                                                         */
/*   This program is distributed in the hope that it will be useful,       */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*   GNU General Public License for more details.                          */
/*                                                                         */
/*   You should have received a copy of the GNU General Public License     */
/*   along with this program; if not, write to the                         */
/*   Free Software Foundation, Inc.,                                       */
/*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             */
/*                                                                         */
/***************************************************************************/

#include <math.h>
#include <strings.h>

#include "pfitsio.h"
#include "dcr.h"

/*--------------------------------------------------------*/
void usage()
{
  printf("\n\tThis is a modified version of DCR! for the Goodman Spectroscopic Pipeline");
  printf("\n\tPlease visit the author's site to get the original version:");
  printf("\n\tModification Version: 0.0.1\n");
  printf("\n\thttp://users.camk.edu.pl/pych/DCR/\n\n");
  printf("\n\tUSAGE:  dcr  input_file  cleaned_file  cosmicrays_file\n\n");
  printf("File 'dcr.par' must be present in the working directory.\n");
  printf("      ~~~~~~\n");
  return;
}
/*--------------------------------------------------------*/
int readpar(char *parfname, PARAMS *par)
{
        char  line[LLEN],
              key[LLEN],
              val[LLEN];
        int   npar;
        FILE  *inpf;

  npar=0;
  if ((inpf=fopen(parfname, "r")) == NULL)
  { perror(parfname); return(-1); }

  while (feof(inpf) == 0)
  {
    if (fgets(line, LLEN, inpf) == NULL)
    {
      printf("\n\tERROR! reading %s\n", parfname);
      return(-1);
    }
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (strcasecmp(key, "XRAD") == 0)
    {
      sscanf(val, "%d", &par->xr);
      npar++;
    }
    else if (strcasecmp(key, "YRAD") == 0)
    {
      sscanf(val, "%d", &par->yr);
      npar++;
    }
    else if (strcasecmp(key, "NPASS") == 0)
    {
      sscanf(val, "%d", &par->npass);
      npar++;
    }
    else if (strcasecmp(key, "THRESH") == 0)
    {
      sscanf(val, "%f", &par->ths);
      npar++;
    }
    else if (strcasecmp(key, "DIAXIS") == 0)
    {
      sscanf(val, "%d", &par->dispaxis);
      npar++;
    }
    else if (strcasecmp(key, "LRAD") == 0)
    {
      sscanf(val, "%d", &par->rlr);
      npar++;
    }
    else if (strcasecmp(key, "URAD") == 0)
    {
      sscanf(val, "%d", &par->rur);
      npar++;
    }
    else if (strcasecmp(key, "GRAD") == 0)
    {
      sscanf(val, "%d", &par->gr);
      npar++;
    }
    else if (strcasecmp(key, "VERBOSE") == 0)
    {
      sscanf(val, "%hd", &par->vopt);
      npar++;
    }
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->vopt)
  {
    printf("\tReading dcr.par file:\n");
    printf("\t------------------------\n");
    printf("Cleaning box x-radius:    %d\n", par->xr);
    printf("Cleaning box y-radius:    %d\n", par->yr);
    printf("Number of clean passes:   %d\n", par->npass);
    printf("Clean threshold:          %f\n", par->ths);
    printf("Dispersion axis:          %d\n", par->dispaxis);
    printf("Replacement lower radius: %d\n", par->rlr);
    printf("Replacement upper radius: %d\n", par->rur);
    printf("Growing radius:           %d\n", par->gr);
    printf("\t------------------------\n");
    printf("%d parameters read\n", npar);
  }

  return(npar);
}
/*--------------------------------------------------------*/
void calc_mean(float **dat, int nx, int ny, double *mean, double *sdev)
{
        int     i,            /* loop numerator           */
                j;            /* loop numerator           */
        double  s,
                ss;

  s = ss = 0.0;
  for (i=0; i<ny; i++)
  {
    for (j=0; j<nx; j++)
    {
      s+=dat[i][j];
      ss+=dat[i][j]*dat[i][j];
    }
  }
  *mean=s/nx/ny;
  *sdev=sqrt((ss-s*s/nx/ny)/nx/ny);

  return;
}
/*--------------------------------------------------------*/
int calc_submean(PARAMS par, float **data, int xs, int ys, int dx, int dy,
                 double *mean, double *sdev)
{
        int     i,
                j,
                k;
        double  s,
                ss,
                thu,
                thl;

  s = ss = 0.0;
  for (i=ys; i<ys+dy; i++)
  {
    for (j=xs; j<xs+dx; j++)
    {
      s+=data[i][j];
      ss+=data[i][j]*data[i][j];
    }
  }

  *mean=s/dx/dy;
  *sdev=sqrt((ss-s*s/dx/dy)/dx/dy);

  if (*sdev == 0.0) return(0);

  thu=*mean+par.ths*(*sdev);
  thl=*mean-par.ths*(*sdev);

  k=0;
  s = ss = 0.0;
  for (i=ys; i<ys+dy; i++)
  {
    for (j=xs; j<xs+dx; j++)
    {
      if ((data[i][j] < thu) && (data[i][j] > thl))
      {
        k++;
        s+=data[i][j];
        ss+=data[i][j]*data[i][j];
      }
    }
  }
  *mean=s/k;
  *sdev=sqrt((ss-s*s/k)/k);

  return(k);
}
/*--------------------------------------------------------*/
float max(float **data, int xs, int ys, int *xmax, int *ymax)
{
        int   i,              /* rows loop numerator      */
              j;              /* columns loop numerator   */
        float maxc;           /* maximum count            */

  maxc=data[0][0];
  *xmax=xs;
  *ymax=ys;
  for (i=0; i<ys; i++)
  {
    for (j=0; j<xs; j++)
      if (data[i][j] > maxc) { maxc=data[i][j]; *xmax=j; *ymax=i; }
  }

  return(maxc);
}
/*--------------------------------------------------------*/
void minmax(float **data, int xs, int ys, int dx, int dy,
            float *minc, float *maxc)
{
        int i,      /* rows loop numerator                */
            j;      /* columns loop numerator             */

  *minc = *maxc = data[ys][xs];
  for (i=ys; i<ys+dy; i++)
  {
    for (j=xs; j<xs+dx; j++)
    {
      if (data[i][j] < *minc) *minc=data[i][j];
      if (data[i][j] > *maxc) *maxc=data[i][j];
    }
  }

  return;
}
/*--------------------------------------------------------*/
int make_hist(int xs, int ys, float **data, int dx, int dy, float min,
  int hs, float bin_width, int *hbuf)
{
        int     i,            /* loop numerator           */
                j,            /* loop numerator           */
                hi;           /* histogram buffer index   */

  memset((void *)hbuf, 0, hs*sizeof(int));
  for (i=ys; i<ys+dy; i++)
  {
    for (j=xs; j<xs+dx; j++)
    {
      hi=(int)((data[i][j]-min)/bin_width);
      if (hi < 0)
      {
        printf("\n\tERROR! make_hist(): (hi= %ld) < 0\n", (long)hi);
        return(EXIT_FAILURE);
      }
      if (hi > hs-1)
      {
        printf("\n\tERROR! make_hist(): (hi= %d) > (hs= %d)\n", hi, hs-1);
        return(EXIT_FAILURE);
      }
      hbuf[hi]++;
    }
  }

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int detect(PARAMS par, float **data, int naxis1, int naxis2, int xs, int ys,
  int dx, int dy, char **map)
{
        int     ix, iy,       /* loop numerator           */
                kx, ky,       /* loop numerator           */
                x, y,         /* whole frame coordinates  */
                k,            /* num. pixels in submean   */
                *hbuf,        /* histogram buffer         */
                hmax,         /* hbuf maximum value       */
                n,            /* number of pixels cleaned */
                nmax,         /* limit of pixels cleaned  */
                hbs,          /* size of histogram buffer */
                mode,         /* histogram  mode          */
                i,            /* loop numerator           */
                j;            /* loop numerator           */
        float   minc,         /* minimum data             */
                maxc,         /* maximum data             */
                hw,           /* width of histogram bin   */
                th;           /* threshold                */
        double  mean,         /* data mean                */
                sdev;         /* std. deviation           */

  minmax(data, xs, ys, dx, dy, &minc, &maxc);
  if (minc == maxc) return(0);

  if ((k=calc_submean(par, data, xs, ys, dx, dy, &mean, &sdev)) <= 0)
  { printf("\n\tERROR! calc_submean() failed\n"); return(-1); }

  hw=1.0;	/* this is valid for the counts in the image */
  hbs=(int)(maxc-minc)/hw+1;

  if ((hbuf=(int *)calloc(hbs, sizeof(int))) == NULL)
  {
    printf("\n\tERROR! detect(): calloc(hbuf)\n");
    return(-1);
  }

  if (make_hist(xs, ys, data, dx, dy, minc, hbs, hw, hbuf) != EXIT_SUCCESS)
    return(-1);

/** find mode = maximum peak of histogram **/
  mode=0;
  hmax=hbuf[0];
  for (i=0; i<hbs; i++)
  {
    if (hbuf[i] > hmax)
    {
      hmax=hbuf[i];
      mode=i;
    }
  }

/** determine clean threshold **/
  if (mode == 0) mode=(int)((mean-minc)/hw);

  j=0;
  for (i=mode; i<hbs; i++)
  {
    if ((hbuf[i]) == 0) j++;
    else                j=0;
    if (j > (int)(par.ths*sdev/hw)) break;
  }
  th=minc+(float)i*hw;

/** count number of pixels to be cleaned **/
  n=0;
  for (j=i; j<hbs; j++) n+=hbuf[j];

  free(hbuf);

  nmax=(int)sqrt((double)dx*dy);
  if (n > nmax)
  {
    if (par.vopt)
      printf("\n\tWARNING: [%d:%d,%d:%d] number of pixels to be cleaned: %d > %d\n",
                        xs+1, xs+dx, ys+1, ys+dy, n, nmax);
    return(0);
  }

/** detect **/
  n=0;
  for (iy=0; iy<dy; iy++)
  {
    y=ys+iy;

    for (ix=0; ix<dx; ix++)
    {
      x=xs+ix;

      if (data[y][x] > th)
      {
  	n++;
        for (ky=-par.gr; ky<=par.gr; ky++)
        {
          if ((y+ky >= 0) && (y+ky < naxis2))
          {
            for (kx=-par.gr; kx<=par.gr; kx++)
            {
              if ((x+kx >= 0) && (x+kx < naxis1))
              {
                map[y+ky][x+kx]=1;
              }
            }
          }
        }
      }
    }
  }

  if ((n != 0) && (par.vopt > 1))
  {
    printf("  min= %.1f max= %.1f mean= %.1f+-%.1f (npix= %d/%d) mode= %.1f\n",
      minc, maxc, mean, sdev, k, dx*dy, minc+(float)mode*hw);
    printf("  threshold= %.1f -> %d pixels to be cleaned\n", th, n);
  }

  return(n);
}
/*--------------------------------------------------------*/
int make_map(PARAMS par, float **data, int naxis1, int naxis2, char **map)
{
        int     i,
                j,
                xs,
                ys,
                hx,
                hy,
                dx,
                dy,
                imax,
                jmax,
                nc,           /* number of pixels cleaned */
                n;            /* number of pixels cleaned */

  hx=par.xr;
  dx=2*hx;
  if (dx > naxis1)
  {
    printf("\n\tERROR! make_map(): x-radius of the box too large\n");
    return(-1);
  }
  imax=naxis1/hx-1;

  hy=par.yr;
  dy=2*hy;
  if (dy > naxis2)
  {
    printf("\n\tERROR! make_map(): y-radius of the box too large\n");
    return(-1);
  }
  jmax=naxis2/hy-1;

/** clean most of the frame **/
  n = nc = 0;
  for (j=0; j<jmax; j++)
  {
    ys=j*hy;
    for (i=0; i<imax; i++)
    {
      xs=i*hx;

      nc=detect(par, data, naxis1, naxis2, xs, ys, dx, dy, map);
      if (nc < 0) return(-1);
      if (nc > 0)
      {
        n+=nc;
        if (par.vopt > 1) printf("[%d:%d,%d:%d]\n", xs+1, xs+dx, ys+1, ys+dy);
      }
    }
  }

/** clean margins of the frame **/
/** left margin **/
  xs=naxis1-dx;
  for (j=0; j<jmax; j++)
  {
    ys=j*hy;

    nc=detect(par, data, naxis1, naxis2, xs, ys, dx, dy, map);
    if (nc < 0) return(-1);
    if (nc > 0)
    {
      n+=nc;
      if (par.vopt > 1) printf("[%d:%d,%d:%d]\n", xs+1, xs+dx, ys+1, ys+dy);
    }
  }

/** top margin **/
  ys=naxis2-dy;
  for (i=0; i<imax; i++)
  {
    xs=i*hx;

    nc=detect(par, data, naxis1, naxis2, xs, ys, dx, dy, map);
    if (nc < 0) return(-1);
    if (nc > 0)
    {
      n+=nc;
      if (par.vopt > 1) printf("[%d:%d,%d:%d]\n", xs+1, xs+dx, ys+1, ys+dy);
    }
  }

/** left-top margin **/
  xs=naxis1-dx;
  ys=naxis2-dy;

  nc=detect(par, data, naxis1, naxis2, xs, ys, dx, dy, map);
  if (nc < 0) return(-1);
  if (nc > 0)
  {
    n+=nc;
    if (par.vopt > 1) printf("[%d:%d,%d:%d]\n", xs+1, xs+dx, ys+1, ys+dy);
  }

  return(n);
}
/*--------------------------------------------------------*/
void clean_xdisp(PARAMS par, float **data, int i, int j, int naxis1, int naxis2,
                 double mean, char **map, float **cf)
{
        int     mx, ns;
        double  s;

  ns=0;
  s=0.0;

  for (mx=-par.rur; mx<=-par.rlr; mx++)
  {
    if (j+mx < 0) continue;
    if (map[i][j+mx]) continue;
    s+=data[i][j+mx];
    ns++;
  }

  for (mx=par.rlr; mx<=par.rur; mx++)
  {
    if (j+mx >= naxis1) break;
    if (map[i][j+mx]) continue;
    s+=data[i][j+mx];
    ns++;
  }

  if (ns) s/=ns;
  else    s=mean;

  cf[i][j]+=data[i][j]-s;
  data[i][j]=s;

  return;
}
/*--------------------------------------------------------*/
void clean_ydisp(PARAMS par, float **data, int i, int j, int nasix1, int naxis2,
                 double mean, char **map, float **cf)
{
        int     my, ns;
        double  s;

  ns=0;
  s=0.0;

  for (my=-par.rur; my<=-par.rlr; my++)
  {
    if (i+my < 0) continue;
    if (map[i+my][j]) continue;
    s+=data[i+my][j];
    ns++;
  }

  for (my=par.rlr; my<=par.rur; my++)
  {
    if (i+my >= naxis2) break;
    if (map[i+my][j]) continue;
    s+=data[i+my][j];
    ns++;
  }

  if (ns) s/=ns;
  else    s=mean;

  cf[i][j]+=data[i][j]-s;
  data[i][j]=s;

  return;
}
/*--------------------------------------------------------*/
void clean_nodisp(PARAMS par, float **data, int i, int j,
                  int naxis1, int naxis2,
                  double mean, char **map, float **cf)
{
        int     mx, my, ns;
        float   d;
        double  s;

  ns=0;
  s=0.0;

  for (my=par.rlr; my<=par.rur; my++)
  {
    for (mx=par.rlr; mx<=par.rur; mx++)
    {
      d=sqrt((float)mx*mx+(float)my*my);
      if ((d > (float)par.rur) || (d < (float)par.rlr)) continue;

      if ((j+mx < naxis1) && (i+my < naxis2))
      {
        if (!map[i+my][j+mx])
        {
          s+=data[i+my][j+mx];
          ns++;
        }
      }

      if ((j-mx >= 0) && (i+my < naxis2))
      {
        if (!map[i+my][j-mx])
        {
          s+=data[i+my][j-mx];
          ns++;
        }
      }

      if ((j+mx < naxis1) && (i-my >= 0))
      {
        if (!map[i-my][j+mx])
        {
          s+=data[i-my][j+mx];
          ns++;
        }
      }

      if ((j-mx >= 0) && (i-my >= 0))
      {
        if (!map[i-my][j-mx])
        {
          s+=data[i-my][j-mx];
          ns++;
        }
      }
    }
  }

  if (ns) s/=ns;
  else    s=mean;

  cf[i][j]+=data[i][j]-s;
  data[i][j]=s;

  return;
}
/*--------------------------------------------------------*/
int clean(PARAMS par, float **data, int naxis1, int naxis2, double mean,
          char **map, float **cf)
{
        int i, j;	/* loop numerator	*/

  for (i=0; i<naxis2; i++)
  {
    for (j=0; j<naxis1; j++)
    {
      if (map[i][j])
      {
        switch (par.dispaxis)
        {
          case 1:  clean_xdisp(par, data, i, j, naxis1, naxis2, mean, map, cf);
                   break;
          case 2:  clean_ydisp(par, data, i, j, naxis1, naxis2, mean, map, cf);
                   break;
          default: clean_nodisp(par, data, i, j, naxis1, naxis2, mean, map, cf);
                   break;
        }
      }
    }
  }

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int modify_header(int ncards, char ***header, PARAMS par)
{
        char  **nheader,
              **lheader,
              ncard[CARD_SIZE],
              ecard[CARD_SIZE],
              val[VALUE_SIZE];
        char gsp_fnam[CARD_SIZE];
        char gsp_dcr_reference[CARD_SIZE];
        size_t fnam_size,
               dcrr_size;

        int   i,
              ni,
              nncards,
              p;

  lheader=*header;

  for (i=0; i<CARD_SIZE; i++) ecard[i]=' ';

  if ((nheader=(char **)calloc(ncards, sizeof(char *))) == NULL)
  { printf("\n\tERROR! modify_header(): calloc(nheader)\n"); return(-1); }

  for (i=0, ni=0; i<ncards; i++)
  {
    if (strncmp(lheader[i], "END     ", KEYWORD_SIZE) == 0) break;

    if (strncmp(lheader[i], ecard, CARD_SIZE) != 0)
    {
      if ((nheader[ni]=(char *)calloc(CARD_SIZE, sizeof(char))) == NULL)
      {
	printf("\n\tERROR! modify_header(): calloc(nheader[ni])"); return(-1);
      }

      memcpy(nheader[ni], lheader[i], CARD_SIZE);
      ni++;
    }
  }
  nncards=ni;

  for (i=0; i<ncards; i++) free(lheader[i]);
  free(lheader);

  if ((p=get_FITS_key(nncards, nheader, "BITPIX", val)) == -1)
  {
    printf("\n\tERROR! modify_header(): BITPIX not found in the header\n");
    return(-1);
  }
  if (p != 1)
    printf("WARNING! header does not conform to FITS standard (BITPIX)\n");
  memcpy(ncard, ecard, CARD_SIZE);
  memcpy(ncard, "BITPIX  =                  -32  / Bits per pixel", 48);
  memcpy(nheader[p], ncard, CARD_SIZE);

  if ((p=get_FITS_key(nncards, nheader, "BZERO", val)) != -1)
  {
    memcpy(ncard, ecard, CARD_SIZE);
    memcpy(ncard, "BZERO   =                  0.0  / Bits per pixel", 48);
    memcpy(nheader[p], ncard, CARD_SIZE);
  }

  if ((p=get_FITS_key(nncards, nheader, "BSCALE", val)) != -1)
  {
    memcpy(ncard, ecard, CARD_SIZE);
    memcpy(ncard, "BSCALE  =                  1.0  /", 33);
    memcpy(nheader[p], ncard, CARD_SIZE);
  }

  nncards+=2;
  if ((nheader=(char **)realloc(nheader, nncards*sizeof(char *))) == NULL)
  {
    printf("\n\tERROR! modify_header(): realloc(nheader)");
    return(-1);
  }
  if ((nheader[nncards-2]=(char *)calloc(CARD_SIZE, sizeof(char))) == NULL)
  {
    printf("\n\tERROR! modify_header(): calloc(nheader[nncards-2])");
    return(-1);
  }
  if ((nheader[nncards-1]=(char *)calloc(CARD_SIZE, sizeof(char))) == NULL)
  {
    printf("\n\tERROR! modify_header(): calloc(nheader[nncards-1])");
    return(-1);
  }

  /* Update GSP_FNAM keyword which stores the name used to save the file */
  if ((p=get_FITS_key(nncards, nheader, "GSP_FNAM", val)) != -1)
  {
    memcpy(ncard, ecard, CARD_SIZE);

    if (strrchr(par.ocname, '/') == NULL)
    {
        /* if the name was parsed as file name only (not path)*/
        fnam_size = snprintf(NULL, 0, "GSP_FNAM= '%s' / Current file name", par.ocname) + 1;

        snprintf(gsp_fnam, fnam_size, "GSP_FNAM= '%s' / Current file name", par.ocname);
    }
    else
    {
        /* if the name was parsed as a full path*/
        fnam_size = snprintf(NULL, 0, "GSP_FNAM= '%s' / Current file name", &par.ocname[(strrchr(par.ocname, '/') -  par.ocname) + 1]) + 1;

        snprintf(gsp_fnam, fnam_size, "GSP_FNAM= '%s' / Current file name", &par.ocname[(strrchr(par.ocname, '/') -  par.ocname) + 1]);
    }

    memcpy(ncard, gsp_fnam , strlen(gsp_fnam));

    memcpy(nheader[p], ncard, CARD_SIZE);
  }


  /* Add DCR reference */
  memcpy(ncard, ecard, CARD_SIZE);

  dcrr_size = snprintf(NULL, 0, "GSP_DCRR= 'see Pych, W., 2004, PASP, 116, 148' / DCR Reference") + 1;

  snprintf(gsp_dcr_reference, dcrr_size, "GSP_DCRR= 'see Pych, W., 2004, PASP, 116, 148' / DCR Reference");

  memcpy(ncard, gsp_dcr_reference, strlen(gsp_dcr_reference));
  memcpy(nheader[nncards-2], ncard, CARD_SIZE);

  /* END card */
  memcpy(ncard, ecard, CARD_SIZE);
  memcpy(ncard, "END", 3);
  memcpy(nheader[nncards-1], ncard, CARD_SIZE);

  *header=nheader;

  return(nncards);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char    **header,     /* FITS header                */
                **map;        /* map of pixels to clean     */
        int     ncards,        /* number of header lines     */
                i, j, ipas,   /* loop numerator             */
                np,           /* number of pixels cleaned in one pass */
                nc,           /* number of pixels cleaned   */
                naxis1,       /* number of image columns    */
                naxis2,       /* number of image rows       */
                xmax,         /* x coordinate of max. peak  */
                ymax;         /* y coordinate of max. peak  */
        float   **data,       /* image data                 */
                **cf,         /* cosmic rays image data     */
                maxc;         /* maximum count              */
        double  mean,         /* mean value of pixel counts */
                sdev;         /* STDDEV of pixel counts     */
        PARAMS  par;

  if (argc != 4) { usage(); return(EXIT_FAILURE); }

  par.iname=argv[1];
  par.ocname=argv[2];
  par.orname=argv[3];

  if (readpar("dcr.par", &par) <= 0)
  { printf("\n\tERROR! readpar() failed\n"); return(EXIT_FAILURE); }

  if (par.vopt > 0)
  {
    printf("Input frame file:         %s\n", par.iname);
    printf("Cosmic rays file:         %s\n", par.orname);
    printf("Cleaned frame file:       %s\n", par.ocname);
  }

  if ((data=read_FITS_2Dfile(par.iname, &ncards, &header, &naxis1, &naxis2)) == NULL)
    return(EXIT_FAILURE);

  if (par.vopt)
  {
    printf("Whole frame before cleaning:\n");
    calc_mean(data, naxis1, naxis2, &mean, &sdev);
    printf("mean= %f +- %f\n", mean, sdev);
    maxc=max(data, naxis1, naxis2, &xmax, &ymax);
    printf("max count= %f (%d,%d)\n", maxc, xmax+1, ymax+1);
  }

  if ((map=(char **)calloc(naxis2, sizeof(char *))) == NULL)
  { printf("\n\tERROR! calloc(map)"); return(EXIT_FAILURE); }
  if ((cf=(float **)calloc(naxis2, sizeof(float *))) == NULL)
  { printf("\n\tERROR! calloc(cf)"); return(EXIT_FAILURE); }
  for (i=0; i<naxis2; i++)
  {
    if ((map[i]=(char *)calloc(naxis1, sizeof(char))) == NULL)
    { printf("\n\tERROR! calloc(map[i])"); return(EXIT_FAILURE); }
    if ((cf[i]=(float *)calloc(naxis1, sizeof(float))) == NULL)
    { printf("\n\tERROR! calloc(cf[i])"); return(EXIT_FAILURE); }
  }

  nc=0;
  for (ipas=1; ipas<=par.npass; ipas++)
  {
    for (i=0; i<naxis2; i++)
      for (j=0; j<naxis1; j++)
        map[i][j]=0;

    if ((np=make_map(par, data, naxis1, naxis2, map)) < 0)
      return(EXIT_FAILURE);
    clean(par, data, naxis1, naxis2, mean, map, cf);
    nc+=np;
    if (par.vopt) printf("Pass %d: %d pixels cleaned\n", ipas, np);
    if (par.vopt > 1) printf("--------------------------\n");
    if (np == 0) break;
  }

  for (i=0; i<naxis2; i++) free(map[i]);
  free(map);

  if (par.vopt)
  {
    printf("Total number of pixels/cleaned= %d/%d (%.1f%%)\n",
            naxis1*naxis2, nc, nc*100.0/naxis1/naxis2);
    printf("Whole frame after cleaning:\n");
    calc_mean(data, naxis1, naxis2, &mean, &sdev);
    printf("mean= %f +- %f\n", mean, sdev);
    maxc=max(data, naxis1, naxis2, &xmax, &ymax);
    printf("max count= %f (%d,%d)\n", maxc, xmax+1, ymax+1);
  }

  if ((ncards=modify_header(ncards, &header, par)) < 1)
  { printf("\n\tERROR! modify_header() failed\n"); return(EXIT_FAILURE); }

  if (par.vopt) printf("Cleaned frame file: %s\n", par.ocname);
  write_FITS_2Dfile(par.ocname, ncards, header, naxis1, naxis2, 4, (void **)data);

  for (i=0; i<naxis2; i++) free(data[i]);
  free(data);

  if (par.vopt) printf("Cosmic rays file:   %s\n", par.orname);
  write_FITS_2Dfile(par.orname, ncards, header, naxis1, naxis2, 4, (void **)cf);

  for (i=0; i<naxis2; i++) free(cf[i]);
  free(cf);
  for (i=0; i<ncards; i++) free(header[i]);
  free(header);

  if (par.vopt) printf("==============================\n");

  return(EXIT_SUCCESS);
}
/*** END ***/
