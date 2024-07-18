/*========================================================*/
/*                                                        */
/*  pfitshead.c     version 3.7.0         2015.05.11      */
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

#include "pfitshead.h"

/*--------------------------------------------------------*/
/*  Read header card value with comment                   */
/*  Return card number (starting from 0)                  */
/*  or -1 on error.                                       */
/*--------------------------------------------------------*/
int get_FITS_key(int ncards, char **header, char *keyword, char *value)
{
        char  lkeyword[KEYWORD_SIZE];
        int   i;              /*  loop numerator          */

  if (strlen(keyword) > KEYWORD_SIZE)
  {
    printf("\n\tWARNING! get_FITS_key(%s): invalid keyword\n", keyword);
    return(-1);
  }

  memcpy(lkeyword, keyword, strlen(keyword));
  for (i=(int)strlen(keyword); i<KEYWORD_SIZE; i++)  lkeyword[i]=' ';

  for (i=0; i<ncards; i++)
  {
    if (strncmp(header[i], lkeyword, KEYWORD_SIZE) == 0)
    {
/* value may start in column 11 */
      memcpy(value, header[i]+10, VALUE_SIZE);
      break;
    }
  }

  if (i == ncards) return(-1);

  return(i);
}
/*--------------------------------------------------------*/
/*  Return: image size in bytes (without the filling)     */
/*          0 no image                                    */
/*         -1 error                                       */
/*--------------------------------------------------------*/
int FITS_image_size(const int ncards, char **header)
{
        char  val[VALUE_SIZE],
              keyword[KEYWORD_SIZE];
        int   i,
              p,
              bitpix,
              bytepix,
              naxis,
              naxisn,
	      size;

  memset(val, ' ', VALUE_SIZE);

  if ((p=get_FITS_key(ncards, header, "BITPIX", val)) == -1)
  {
    printf("\n\tERROR! FITS_image_size(): BITPIX not found in the header\n");
    return(0);
  }
  if (p != 1)
    printf("\n\tWARNING! header does not conform to FITS standard (BITPIX)\n");
  if (sscanf(val, "%d", &bitpix) < 1)
  {
    printf("\n\tERROR! FITS_image_size(): reading BITPIX\n");
    return(0);
  }
  bytepix=abs(bitpix)/8;

  strcpy(keyword, "NAXIS");
  if ((p=get_FITS_key(ncards, header, keyword, val)) == -1)
  {
    printf("\n\tERROR! FITS_image_size(): NAXIS not found in the header\n");
    return(0);
  }
  if (p != 2)
    printf("\n\tWARNING! header does not conform to FITS standard (NAXIS)\n");

  if (sscanf(val, "%d", &naxis) < 1)
  {
    printf("\n\tERROR! FITS_image_size(): reading NAXIS\n");
    return(0);
  }

  if (naxis == 0) return(0);

  size=bytepix;
  for (i=1; i<=naxis; i++)
  {
    if (snprintf((keyword+5), KEYWORD_SIZE-5, "%d", i) < 1)
    {
      printf("\n\tERROR! FITS_image_size(): %s\n", keyword);
      return(0);
    }
    if ((p=get_FITS_key(ncards, header, keyword, val)) == -1)
    {
      printf("\n\tERROR! FITS_image_size(): %s not found in the header\n",
              keyword);
      return(0);
    }
    if (p != i+2)
      printf("\n\tWARNING! header does not conform to FITS standard (%s)\n",
        keyword);

    if (sscanf(val, "%d", &naxisn) < 1)
    {
      printf("\n\tERROR! FITS_image_size(): reading %s\n", keyword);
      return(0);
    }
    size*=naxisn;
  }

  return(size);
}
/*--------------------------------------------------------*/
/*  Read FITS header (primary or extension)               */
/*  Return number of card images.                         */
/*--------------------------------------------------------*/
int read_FITS_header(FILE *inf, char ***header)
{
        char  **lheader;
        short eoh;            /* end of header  */
        int   i,
              hc,             /* number of header cards   */
              hr;             /* number of header records */

  eoh=0;

  if ((lheader=(char **)calloc(RECORD_CARDS, sizeof(char *))) == NULL)
  {
    perror("\n\tERROR! read_FITS_header(): calloc(lheader)");
    return(0);
  }

  hc=0;
  for (hr=1; eoh==0; hr++)
  {
    if (!(lheader=(char **)realloc(lheader, hr*RECORD_CARDS*sizeof(char *))))
    {
      perror("\n\tERROR! read_FITS_header(): realloc(lheader)");
      return(0);
    }

    for (i=0; i<RECORD_CARDS; i++)
    {
      if ((lheader[hc+i]=(char *)calloc(CARD_SIZE, sizeof(char))) == NULL)
      {
        perror("\n\tERROR! read_FITS_header(): calloc(lheader[hc+i])");
        return(0);
      }

      if (fread(lheader[hc+i], sizeof(char), CARD_SIZE, inf) != CARD_SIZE)
      {
        printf("\n\tERROR! read_FITS_header(): header corrupted\n");
        return(0);
      }
    }

    for (i=0; i<RECORD_CARDS; i++)
      if (strncmp(lheader[hc+i], "END     ", KEYWORD_SIZE) == 0) eoh=1;

    hc+=RECORD_CARDS;
  }

  *header=lheader;

  return(hc);
}
/*--------------------------------------------------------*/
/*  Read primary and extensions headers                   */
/*  Return number of extensions (0=primary header only)   */
/*  Byte offsets for extension headers recorded.          */
/*--------------------------------------------------------*/
int read_FITS_headers(FILE *inf, int **ncards, char ****header, int **offsets)
{
        char    ***lheader,         /* local pointer to the header  */
                value[VALUE_SIZE],  /* card value buffer            */
                val[VALUE_SIZE],    /* card value buffer            */
                tmp;
        short   eoF;
        int     p,
                hen,        /* number of header extensions (0=primary header) */
                *lncards,   /* size of the header extension                 */
                offset,     /* FITS file offset                             */
                *loffsets,  /* table of FITS file offsets to the extensions */
                imsize;     /* image size                                   */

  memset(value, ' ', VALUE_SIZE);

/* read primary header */
  if ((lheader=(char ***)calloc(1, sizeof(char **))) == NULL)
  {
    perror("\n\tERROR! read_FITS_headers(): calloc(lheader)");
    return(-1);
  }

  if ((lncards=(int *)calloc(1, sizeof(int))) == NULL)
  {
    perror("\n\tERROR! read_FITS_headers(): calloc(lncards)");
    return(-1);
  }

  if ((loffsets=(int *)calloc(1, sizeof(int))) == NULL)
  {
    perror("\n\tERROR! read_FITS_headers(): calloc(loffsets)");
    return(-1);
  }

  loffsets[0]=0;
  if ((lncards[0]=read_FITS_header(inf, &lheader[0])) <= 0)
  {
    free(lheader);
    free(lncards);
    free(loffsets);
    return(-1);
  }

/* check conformance with the FITS standard */
  if ((p=get_FITS_key(lncards[0], lheader[0], "SIMPLE", value)) == -1)
  {
    printf("\n\tERROR! read_FITS_headers(): SIMPLE not found in the header\n");
    free(lheader);
    free(lncards);
    free(loffsets);
    return(-1);
  }
  if (p != 0)
    printf("\n\tWARNING! Header does not conform to FITS standard (SIMPLE)\n");
  if (sscanf(value, "%s", val) < 1)
  {
    printf("\n\tERROR! read_FITS_headers(): reading SIMPLE\n");
    free(lheader);
    free(lncards);
    free(loffsets);
    return(-1);
  }
  if ((strcmp(val, "T") != 0) && (strncmp(val, "T/", 2) != 0))
    printf("\n\tWARNING! File does not conform to FITS standard (SIMPLE)\n");

/* check for possible extensions */
  if ((p=get_FITS_key(lncards[0], lheader[0], "EXTEND", value)) == -1)
  {
    *ncards=lncards;
    *header=lheader;
    *offsets=loffsets;

    return(0);
  }

  if (p < 3)
    printf("\n\tWARNING! Header does not conform to FITS standard (EXTEND)\n");
  if (sscanf(value, "%s", val) < 1)
  {
    printf("\n\tERROR! read_FITS_headers(): reading SIMPLE\n");
    free(lheader);
    free(lncards);
    free(loffsets);
    return(-1);
  }
  if ((strcmp(val, "T") != 0) && (strncmp(val, "T/", 2) != 0))
  {
    *ncards=lncards;
    *header=lheader;
    *offsets=loffsets;

    return(0);
  }

/* read extension headers */
  hen=0;
  eoF=0;
  while (eoF == 0)
  {
/* skip the image */
    imsize=FITS_image_size(lncards[hen], lheader[hen]);
    offset=imsize;
    if (imsize%RECORD_SIZE != 0) offset+=(RECORD_SIZE-imsize%RECORD_SIZE);
    if (fseek(inf, offset, SEEK_CUR) != 0)
    {
      printf("\n\tERROR! read_FITS_headers(): fseek(%d)\n", offset);
      free(lheader);
      free(lncards);
      free(loffsets);
      return(-1);
    }

/* check for the existence of another extension */
    imsize=fread(&tmp, 1, 1, inf);
    if (feof(inf) != 0) break;
    if (fseek(inf, -1, SEEK_CUR) != 0)
    {
      printf("\n\tERROR! read_FITS_headers(): fseek(CUR-1)\n");
      free(lheader);
      free(lncards);
      free(loffsets);
      return(-1);
    }

/* read another header */
    hen++;

    if ((lheader=(char ***)realloc(lheader, (hen+1)*sizeof(char **))) == NULL)
    {
      perror("\n\tERROR! read_FITS_headers(): realloc(lheader)");
      return(-1);
    }

    if ((lncards=(int *)realloc(lncards, (hen+1)*sizeof(int))) == NULL)
    {
      perror("\n\tERROR! read_FITS_headers(): realloc(lncards)");
      return(-1);
    }

    if ((loffsets=(int *)realloc(loffsets, (hen+1)*sizeof(int))) == NULL)
    {
      perror("\n\tERROR! read_FITS_headers(): realloc(loffsets)");
      return(-1);
    }

    loffsets[hen]=ftell(inf);
    if ((lncards[hen]=read_FITS_header(inf, &lheader[hen])) <= 0)
    {
      free(lheader);
      free(lncards);
      free(loffsets);
      return(-1);
    }

/* check the XTENSION value - only IMAGE supported */
    if ((p=get_FITS_key(lncards[hen], lheader[hen], "XTENSION", value)) == -1)
    {
      printf("\n\tERROR! read_FITS_headers(): XTENSION not found\n");
      free(lheader);
      free(lncards);
      free(loffsets);
      return(-1);
    }
    if (p != 0)
      printf("\n\tWARNING! Header does not conform to FITS standard (XTENSION)\n");

    if (strncmp(value, "'IMAGE   '", 10) != 0)
      printf("\n\tWARNING! Unsupported extension type: %s\n", value);

    if (lncards[hen] == 0) eoF=1;
  }

  *ncards=lncards;
  *header=lheader;
  *offsets=loffsets;

  return(hen);
}
/*--------------------------------------------------------*/
int del_header_card(int ncards, char ***header, char *keyword)
{
        char  **lheader,
              val[VALUE_SIZE];
        int   p,
              i;

  lheader=*header;
  memset(val, ' ', VALUE_SIZE);

  if ((p=get_FITS_key(ncards, lheader, keyword, val)) == -1) return(ncards);

  for (i=p; i<ncards-1; i++)
    memcpy(lheader[i], lheader[i+1], CARD_SIZE);

  if ((p=get_FITS_key(ncards, lheader, "END", val)) == -1)
  {
    printf("\n\tERROR! del_header_card(): header corrupted (END)\n");
    return(-1);
  }

  if ((p+1)%RECORD_CARDS == 0)
  {
    ncards-=RECORD_CARDS;

    for (i=0; i<RECORD_CARDS; i++) free(lheader[ncards+i]);
    if ((lheader=(char **)realloc(lheader, ncards*sizeof(char *))) == NULL)
    {
      perror("\n\tERROR! del_header_card(): realloc(lheader)");
      return(-1);
    }

    *header=lheader;
  }

  return(ncards);
}
/*------------------------------------------------------*/
int add_header_card(char *new_card, int *ncards, char ***header)
{
        char  **lheader,        /* local pointer to the header  */
              card[CARD_SIZE];  /* header card                  */
        int   p,                /* pointer to the card          */
              i,                /* loop numerator               */
              lncards;          /* local number of header cards */

  lncards=*ncards;
  lheader=*header;

  if ((p=get_FITS_key(lncards, lheader, "END", card)) < 1)
  {
    printf("\n\tERROR! add_header_card(): END not found in the FITS header\n");
    return(-1);
  }

  memset(lheader[p], ' ', CARD_SIZE);

  memset(card, ' ', CARD_SIZE);
  for (i=0; i<p; i++)
  {
    if (strncmp(card, lheader[i], CARD_SIZE) == 0)
    {
      p=i;
      break;
    }
  }

  memcpy(lheader[p], new_card, strlen(new_card));

  if (p >= lncards-1)
  {
    if (lncards%RECORD_CARDS == 0)
    {
      lncards+=RECORD_CARDS;
      if ((lheader=(char **)realloc(lheader, lncards*sizeof(char *))) == NULL)
      {
        perror("\n\tERROR! add_header_card(): realloc(lheader)");
        return(-1);
      }

      for (i=1; i<=RECORD_CARDS; i++)
      {
        if ((lheader[lncards-i]=(char *)calloc(CARD_SIZE, sizeof(char))) == NULL)
        {
          perror("\n\tERROR! add_header_card(): malloc(lheader[lncards-i])");
          return(-1);
        }
        memset(lheader[lncards-i], ' ', CARD_SIZE);
      }
    }
    else
    {
      lncards++;
      if ((lheader=(char **)realloc(lheader, lncards*sizeof(char *))) == NULL)
      {
        perror("\n\tERROR! add_header_card(): realloc(lheader)");
        return(-1);
      }
      memset(lheader[lncards-1], ' ', CARD_SIZE);
    }
  }

  memcpy(lheader[p+1], "END", 3);

  *ncards=lncards;
  *header=lheader;

  return(p);
}
/*--------------------------------------------------------*/
int write_FITS_header(FILE *outf, const int ncards, char **header)
{
        char  tmp[CARD_SIZE];
        int   i;

  for (i=0; i<ncards; i++)
  {
    if (fwrite(header[i], 1, CARD_SIZE, outf) != CARD_SIZE)
    {
      printf("\n\tERROR! write_FITS_header(): writing FITS header\n");
      return(-1);
    }
  }

  if (ncards%RECORD_CARDS != 0)
  {
    memset(tmp, ' ', CARD_SIZE);
    for (i=0; i<RECORD_CARDS-ncards%RECORD_CARDS; i++)
    {
      if (fwrite(tmp, 1, CARD_SIZE, outf) != CARD_SIZE)
      {
        printf("\n\tERROR! write_FITS_header(): padding FITS header\n");
        return(-1);
      }
    }
  }

  return(ncards);
}
/*** END ***/
