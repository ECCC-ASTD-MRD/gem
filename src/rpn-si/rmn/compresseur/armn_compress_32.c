#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <rpnmacros.h>
#include <unistd.h>
#include "zfstlib.h"
#include <string.h>

int  c_armn_compress32(unsigned char *, float *, int, int, int, int);
int  c_armn_uncompress32(float *fld, unsigned char *zstream, int ni, int nj, int nk, int nchiffres_sign);
int  c_fstzip32(unsigned int *zfld, unsigned int *fld, int ni, int nj, int nk, int step, int nbits, int remaining_space);
int  c_fstzip_parallelogram32(unsigned int *zfld, int *zlng, unsigned int *fld, int ni, int nj, int step, int nbits, word *header);
int  f77name(armn_compress32)(unsigned char *, float *, int *, int *, int *, int *);
int  f77name(armn_uncompress32)(float *fld, unsigned char *zstream, int *ni, int *nj, int *nk, int *nchiffres_sign);
void c_fstunzip(unsigned int *fld, unsigned int *zfld, int ni, int nj, int nbits);    
void pack1bitRLE(unsigned int z[], unsigned int *zlng, unsigned char ufld[], int npts);
void packTokensParallelogram(unsigned int z[], unsigned int *zlng, unsigned short ufld[], int ni, int nj, int nbits, int istep, word *header);
void packTokensParallelogram32(unsigned int z[], int *zlng, unsigned int ufld[], int ni, int nj, int nbits, int istep, int remaining_space);
void packTokensParallelogram_8(unsigned int z[], unsigned int *zlng, unsigned char ufld[], int ni, int nj, int nbits, int istep);
void pack_stream_nbits_32(unsigned int z[], unsigned int *zlng, unsigned int ufld[], int npts, int nbits);
void unpackTokensParallelogram(unsigned short ufld[], unsigned int z[], int ni, int nj, int nbits, int istep, word *header);
void unpackTokensParallelogram32(unsigned int ufld[], unsigned int z[], int ni, int nj, int nbits, int istep);
void unpackTokensParallelogram_8(unsigned char ufld[], unsigned int z[], int ni, int nj, int nbits, int istep);


int compact_mask_char(unsigned int *dest, unsigned char *src, int npts);
int uncompact_mask_char(int *dest, unsigned int *src, int npts);
void pack1bitRLE(unsigned int z[], unsigned int *zlng, unsigned char ufld[], int npts);
void pack_stream_nbits_16(unsigned int z[], unsigned int *zlng, unsigned short ufld[], int npts, int nbits);
void pack_stream_nbits_32(unsigned int z[], unsigned int *zlng, unsigned int ufld[], int npts, int nbits);
void pack_stream_nbits_8(unsigned int z[], unsigned int *zlng, unsigned char ufld[], int npts, int nbits);
void unpack1bitRLE(unsigned char ufld[], unsigned int z[], unsigned int *zlng,  int npts);
void unpack_stream_nbits_16(unsigned short ufld[], unsigned int z[], int npts, int nbits);
void unpack_stream_nbits_32(unsigned int ufld[], unsigned int z[], int npts, int nbits);
void unpack_stream_nbits_8(unsigned char ufld[], unsigned int z[], int npts, int nbits);


#define MEME_SIGNE_POSITIF   0x00
#define MEME_SIGNE_NEGATIF   0x10
#define DIFF_SIGNE_PACKED    0x20
#define DIFF_SIGNE_STREAM    0x30
#define MEME_EXPOSANT        0x00
#define DIFF_EXPOSANT_PACKED 0x08
#define DIFF_EXPOSANT_STREAM 0x0C
#define MANTISSE_PACKED      0x00
#define MANTISSE_STREAM      0x01  

#define SEQUENCE 0
#define COUNT    1
#define ZERO     0
#define UN       1
  

static unsigned char fastlog[256];
static int once = 0;
extern int zfst_msglevel;


int f77name(armn_compress32)(unsigned char *zstream, float *fld, int *ni, int *nj, int *nk, int *nbits)
  {
  return c_armn_compress32(zstream, fld, *ni, *nj, *nk, *nbits);
  }
  
int c_armn_compress32(unsigned char *zstream, float *fld, int ni, int nj, int nk, int znbits)
  {
  float *p_fld;
  int i, meme_signe, remaining_space;
  int nbits;
  unsigned char *exposant, *exposant2;
  unsigned char *le_pointeur, *pos_lng_signe, *pos_lng_exposant, *pos_lng_mantisse;
  unsigned char *signe, *zsigne, code_signe, code_exposant, code_mantisse;
  unsigned char codes;
  unsigned int *mantisse, *mantisse_stream, *la_mantisse,*zmantisse, exp_base;
  unsigned int *temp;
  unsigned int exp_min, exp_max, *zexposant;
  unsigned int le_signe_or, le_signe_and;
  unsigned int npts, zlng, lng_signe, lng_exposant,lng_mantisse,nbits_needed;
  unsigned int zieee_info;
 

  _fstzip zfstzip; 
  _floatint r_exp_max;
  
  if (ni < 16 || nj < 16)
    {
    zlng = -1;
    fprintf(stderr, "*** <armn_compress32> : The dimensions of NI and NJ have to be > 16\n");
    return zlng;
    }
  
  nbits = znbits - 9;
  memset(&zfstzip, 0, sizeof(_fstzip));
  
  zfstzip.predictor_type = PARALLELOGRAM32;
  zfstzip.step           = 3;
  zfstzip.degree         = 1;
  zfstzip.nbits          = nbits;
  zfstzip.levels         = 1;
  zfstzip.version        = 2;

  mantisse_stream = NULL;
  pos_lng_signe = NULL;
  
  npts = ni * nj;
  signe    = (unsigned char *)   malloc(2*npts * sizeof(unsigned char));
  zsigne    = (unsigned char *)  malloc(2*npts * sizeof(unsigned char));
  exposant = (unsigned char *)  malloc(2*npts * sizeof(unsigned char));
  exposant2 = (unsigned char *) malloc(2*npts * sizeof(unsigned char));
  zexposant = (unsigned int *)  malloc(npts*sizeof(unsigned int));
  mantisse = (unsigned int *)    malloc(2*npts * sizeof(unsigned int));
  zmantisse = (unsigned int *)    malloc(2*npts * sizeof(unsigned int));
  le_signe_or = 0;
  le_signe_and = 0xFFFFFFFF;
  
  p_fld = fld;
  for (i=0; i < npts; i++)
    {
    temp = (unsigned int *) p_fld;
    signe[i]    = *temp >> 31;
    le_signe_or |= *temp;
    le_signe_and &= *temp;
    exposant[i] = (*temp >> 23) & 0xff;
    exposant2[i] = exposant[i];
    mantisse[i] = (*temp & 0x07fffff);
    p_fld++;
    }    

  if (nbits < 23)
    {
    for (i=0; i < npts; i++)
      {
      mantisse[i] = (mantisse[i] >> (23-nbits));
      }
    }
  meme_signe=0;
  if ((le_signe_or>>31) == (le_signe_and>>31))
    {
    meme_signe=1;
    }
  
  exp_min = exposant2[0];
  exp_max = exp_min;
  for (i=0; i < npts; i++)
    {
    if (exposant2[i] < exp_min) exp_min = exposant2[i];
    if (exposant2[i] > exp_max) exp_max = exposant2[i];
    }      
    
 
  exp_base = exp_min;
  
  exp_max = exp_max - exp_min;
  exp_min = 0;
  for (i=0; i < npts; i++)
    {
    exposant2[i] = exposant2[i]- exp_base;
    }        
    
  if (exp_max == exp_min)
    {
    nbits_needed =0;
    }
  else
    {
    nbits_needed =(int)(1+log(exp_max+0.5)/log(2.0));
    r_exp_max.f = (float)(exp_max);
    nbits_needed = (r_exp_max.i >> 23) - 126;
    }
  

/* ------------------- Encodage du stream de signe ------------------------- */  
  le_pointeur = &(zstream[8]);
  if (meme_signe == 1)
    {
    lng_signe = 0;

    if ((le_signe_or >> 31) == 1)
      {
      code_signe = MEME_SIGNE_NEGATIF;
      }
    else
      {
      code_signe = MEME_SIGNE_POSITIF;
      }
    }
  else
    {
    pos_lng_signe = le_pointeur;
    le_pointeur += sizeof(unsigned int);
    pack1bitRLE((unsigned int *)le_pointeur, &lng_signe, signe, npts);
    code_signe = DIFF_SIGNE_PACKED;
    }
  
  if (lng_signe > (npts/4))
    {
    compact_mask_char((unsigned int *) le_pointeur, signe, npts);
    lng_signe = 1 + (npts >> 5);
    code_signe = DIFF_SIGNE_STREAM;

    }
    
    if (0 != (lng_signe%4))
      {
      lng_signe += (4 - (lng_signe%4));
      }
    
    le_pointeur+= lng_signe;
    
  
   if (code_signe == DIFF_SIGNE_PACKED || code_signe == DIFF_SIGNE_STREAM )
    {
    memcpy(pos_lng_signe, &lng_signe, sizeof(unsigned int));
    }

    
/* ------------------- Encodage du stream d'exposants  ------------------------- */  
  pos_lng_exposant = le_pointeur;
  if (nbits_needed == 0)
    {
    lng_exposant = 0;
    code_exposant = MEME_EXPOSANT;
    }
  else
    {
    code_exposant = DIFF_EXPOSANT_PACKED;
    le_pointeur += sizeof(unsigned int);
    packTokensParallelogram_8((unsigned int *)le_pointeur, &lng_exposant, exposant2, ni, nj, nbits_needed, 3);
    if (lng_exposant > ni*nj)
      {
         zlng = -1;
         fprintf(stderr, "*** <armn_compress32> : Exponent range too large\n");
         fprintf(stderr, "*** <armn_compress32> : Original field left uncompressed\n");
         return zlng;
         /*pack_stream_nbits_8((unsigned int *)le_pointeur, &lng_exposant, exposant2, npts, nbits_needed);
      code_exposant = DIFF_EXPOSANT_STREAM;*/
      }
    
     if (0 != (lng_exposant%4))
      {
      lng_exposant += 4 - (lng_exposant%4);
      }
    memcpy(pos_lng_exposant, &lng_exposant, sizeof(unsigned int));
    le_pointeur += lng_exposant;
    }
/* ------------------- Encodage du stream de mantisse ------------------------- */  
  pos_lng_mantisse = le_pointeur;
  le_pointeur += sizeof(unsigned int);
  code_mantisse = MANTISSE_PACKED;
  remaining_space = (ni*nj*znbits)/32;
  remaining_space = remaining_space - (le_pointeur - zstream)/sizeof(unsigned int);
  lng_mantisse   = c_fstzip32((unsigned int *)le_pointeur, mantisse, ni, nj, nk, zfstzip.step, nbits, remaining_space);
  la_mantisse = zmantisse;

  if (lng_mantisse == 0)
    {
    free(signe);
    free(zsigne);
    free(exposant);
    free(exposant2);
    free(zexposant);
    if (la_mantisse != mantisse)
      {
      free(mantisse_stream);
      }
    free(mantisse);
    free(zmantisse);
    return -1;
    }

  if (0 != (lng_mantisse%4))
    {
    lng_mantisse += 4 - (lng_mantisse%4);
    }

     memcpy(pos_lng_mantisse, &lng_signe, sizeof(unsigned int));
 
/* ------------------- ------------------------------- ------------------------- */  
  

/* ------------------- ----------- Assemblage final -------------------          */

  
/* ------------------- Entete et codes ------------------- */  

  zstream[0] = (unsigned char)'\0';
  
  i = 0;
  memcpy(&zstream[0], &zfstzip, sizeof(unsigned int));
  zieee_info = 0;
  i+= sizeof(unsigned int);
  codes = (unsigned char)(code_signe | code_exposant | code_mantisse);
  zieee_info = ((unsigned char) exp_base) << 16;
  zieee_info = zieee_info | ((unsigned char)nbits_needed) << 8;
  zieee_info = zieee_info  | (unsigned char) codes;
  memcpy(&zstream[i], &zieee_info, sizeof(unsigned int));
  i+=sizeof(unsigned int);
  
  le_pointeur += lng_mantisse;
  zlng = le_pointeur - zstream;
 
/* ----------------  Menage avant de s'en aller  ------------------- */
   
  free(signe);
  free(zsigne);
  free(exposant);
  free(exposant2);
  free(zexposant);
  if (la_mantisse != mantisse)
    {
    free(mantisse_stream);
    }
  free(mantisse);
  free(zmantisse);
  return zlng;
  }
 

  
int f77name(armn_uncompress32)(float *fld, unsigned char *zstream, int *ni, int *nj, int *nk, int *nbits)
  {
  return c_armn_uncompress32(fld, zstream, *ni, *nj, *nk, *nbits);
  }
  
int c_armn_uncompress32(float *fld, unsigned char *zstream, int ni, int nj, int nk, int znbits)
  {
  _fstzip zfstzip;
  int bitPackInWord, saut, zlng;
  int i, nbits_mantisse, nbits;
  unsigned char *exposant, *exposant2;
  unsigned char *signe, *zsigne, code_signe, code_exposant, code_mantisse;
  unsigned char codes;
  unsigned int *cur, curword, nbits_needed_exposant;
  unsigned int *mantisse;
  unsigned int *temp;
  unsigned int exp_min;
  unsigned int npts, lng_signe, lng_exposant,lng_mantisse, zieee_info;
      
  npts = ni * nj;
  signe    = (unsigned char *) malloc(2*npts * sizeof(unsigned char));
  zsigne    = (unsigned char *) malloc(2*npts * sizeof(unsigned char));
  exposant = (unsigned char *) malloc(2*npts * sizeof(unsigned char));
  exposant2 = (unsigned char *) malloc(2*npts * sizeof(unsigned char));
  mantisse = (unsigned int *)   malloc(2*npts * sizeof(unsigned int));
  
  bitPackInWord = 32;
  
  cur = (unsigned int *)zstream;
  curword = *cur;
  memcpy(&zfstzip, cur, sizeof(unsigned int));
  cur++;
  memcpy(&zieee_info, cur, sizeof(unsigned int));
  cur++;
  
  
  exp_min = zieee_info >> 16;
  nbits_needed_exposant = (zieee_info >> 8) & 0xFF;
  nbits_mantisse = zfstzip.nbits;
  nbits = nbits_mantisse;
  codes = zieee_info& 0xFF;
  
  code_signe =    (codes & 0x30);
  code_exposant = (codes & 0xC);
  code_mantisse = (codes & 0x3);
  
 
  if (code_signe == DIFF_SIGNE_PACKED || code_signe == DIFF_SIGNE_STREAM)
    {
    memcpy(&lng_signe, cur, sizeof(unsigned int));
    cur++;
    unpack1bitRLE(signe, cur, (unsigned int *)&zlng,npts);
    saut = (lng_signe >> 2);
    cur += saut;
    }
  else
    {
    if (code_signe == MEME_SIGNE_POSITIF)
      {
      for (i=0; i < npts; i++)
        {
        signe[i] = 0;
        }
      }
    else
      {
      for (i=0; i < npts; i++)
        {
        signe[i] = 1;
        }
      
      }
    }

  if (code_exposant == DIFF_EXPOSANT_PACKED || code_exposant == DIFF_EXPOSANT_STREAM)
    {
    memcpy(&lng_exposant, cur, sizeof(unsigned int));
    cur++;
    unpackTokensParallelogram_8(exposant2, cur, ni, nj, nbits_needed_exposant, 3);
/*    lng_exposant = armn_compress(exposant2, *ni, *nj, *nk, nbits_needed_exposant, 2);*/
    saut = (lng_exposant >> 2);
    cur += saut;
    for (i=0; i < npts; i++)
      {
      exposant2[i] = exposant2[i] + exp_min;
      }    
    }
  else
    {
    for (i=0; i < npts; i++)
      {
      exposant2[i] = exp_min;
      }    

    }
  
  
  memcpy(&lng_mantisse, cur, sizeof(unsigned int));
  cur++;
  if (code_mantisse == MANTISSE_PACKED)
    {
    unpackTokensParallelogram32(mantisse, cur, ni, nj, nbits_mantisse, 3);
    }
  else
    {
    unpack_stream_nbits_32(mantisse, cur, npts, nbits_mantisse);
    }
  
  switch(code_signe)
    {
    case MEME_SIGNE_POSITIF:
    for (i=0; i < npts; i++)
      {
      temp = (unsigned int *)&fld[i];
      *temp = 0;
      *temp |= (exposant2[i] << 23);
      *temp |=  (mantisse[i] << (23-nbits));
      }
    break;
    
    case MEME_SIGNE_NEGATIF:
    for (i=0; i < npts; i++)
      {
      temp = (unsigned int *)&fld[i];
       *temp = 0x80000000;
      *temp = *temp | (exposant2[i] << 23);
      *temp = *temp | (mantisse[i] << (23-nbits));
      }
    break;
    
    default:
    for (i=0; i < npts; i++)
      {
      temp = (unsigned int *)&fld[i];
       *temp = signe[i] << 31;
      *temp |= (exposant2[i] << 23);
      *temp |= (mantisse[i] << (23-nbits));
      }
    break;
    }  
     
  free(signe);
  free(zsigne);
  free(exposant);
  free(exposant2);
  free(mantisse);

  return ni*nj;
  }
  
int c_fstzip32(unsigned int *zfld, unsigned int *fld, int ni, int nj, int nk, int step, int nbits, int remaining_space)
  {
  int zlng, lng_origin;
  
  lng_origin = (1+(ni*nj*nk*1.0*nbits)/8);

  if (ni == 1 || nj == 1)
    {
    return lng_origin;
    }
      
  packTokensParallelogram32(zfld, &zlng, fld, ni, nj, step, nbits, remaining_space);
  
 
  if (zlng == 0 && zfst_msglevel <= 2)
    {
    fprintf(stdout, "IEEE compressed field is larger than original... Returning original\n\n");
    return 0;
    }
  
  return zlng;
}  


  
void packTokensParallelogram32(unsigned int z[], int *zlng, unsigned int ufld[], int ni, int nj, int istep, int nbits, int remaining_space)
{

  _floatint r_lmax;
  int *ufld_dst;
  int k22, nbits2;
  int lcl_m, lcl_n;
  int local_max;
  unsigned int *cur;
  unsigned int i, j, k, m, n;
  unsigned int lastWordShifted, spaceInLastWord, lastSlot;
  unsigned int lcl, nbits_needed;
  unsigned int nbits_req_container, token;
  
  lastSlot = 0;
  cur = z;
  
  ufld_dst=(int *) malloc(ni*nj*sizeof(int));
  
  for (j=1; j <= nj; j++)
   {
   k = FTN2C(1,j,ni);
   ufld_dst[k] = 0;
   } 
   
  for (i=1; i <= ni; i++)
   {
   k = FTN2C(i,1,ni);
   ufld_dst[k] = 0;
   }
   
  for (j=2; j<=nj; j++)
   {
   for (i=2; i <=ni; i++)
      {
      k22 = FTN2C(i,  j,  ni);
      ufld_dst[k22] = ufld[k22] - (ufld[k22-ni]+ufld[k22-1]-ufld[k22-1-ni]);
      }
   }  
  

  nbits_req_container = 5;

  lastWordShifted = 0;
  spaceInLastWord = 32;
  *cur = 0;
  
  stuff(nbits_req_container, cur, 32, istep, lastWordShifted, spaceInLastWord);
  
  for (i=1; i <= ni; i++)
   {
   k = FTN2C(i,1,ni);
   stuff(ufld[k], cur, 32, nbits, lastWordShifted, spaceInLastWord);
   }
   
  for (j=2; j <= nj; j++)
   {
   k = FTN2C(1,j,ni);
   stuff(ufld[k], cur, 32, nbits, lastWordShifted, spaceInLastWord);
   }
   
  for (j=2; j <= nj; j+=istep)
    {
    lcl_n = ((j + istep - 1) >= nj ? nj - j : istep - 1);
    for (i=2; i <= ni; i+=istep)
      {
      k = FTN2C(i,j,ni);
      local_max = ufld_dst[k];
      lcl_m = ((i + istep - 1) >= ni ? ni - i : istep - 1);
      for (n=0; n <= lcl_n; n++)
        {
        for (m=0; m <= lcl_m; m++)
          {
          k = FTN2C(i+m,j+n,ni);
          if (local_max < abs(ufld_dst[k])) local_max = abs(ufld_dst[k]);
          }
        }
      if (local_max == 0)
        {
        nbits_needed = 0;
        }
      else
        {
        r_lmax.f = (float)local_max;
        nbits_needed = (r_lmax.i >> 23) - 126;
        }
      stuff(nbits_needed, cur, 32, nbits_req_container, lastWordShifted, spaceInLastWord);
      switch (nbits_needed)
        {
        case 0:
        break;
              
        default:
        nbits2 = nbits_needed + 1;
        for (n=0; n <= lcl_n; n++)
          {
          for (m=0; m <= lcl_m; m++)
            {
            k = FTN2C(i+m,j+n,ni);
            token = (unsigned int) (ufld_dst[k] & ~((-1)<<nbits2));
            stuff(token, cur, 32, nbits2, lastWordShifted, spaceInLastWord);
            }
          }
        if (remaining_space < ((cur - z)+(1+((nbits_needed+9*nbits)>>5))))
          {
          *zlng = 0;
          return;
          }
        break;
        } 
       }
    }

  
    lcl = 0;
    stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord);
    stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord); 

   *zlng = 1 + (int) (cur-z) * 4;
    free(ufld_dst);

  }


void unpackTokensParallelogram32(unsigned int ufld[], unsigned int z[], int ni, int nj, int nbits, int istep)
{

  int *ufld_tmp;
  int bitPackInWord;
  int i, j, k, m, n;
  int k11, k12, k21, k22;
  int lcl_m, lcl_n;
  unsigned int  nbits_needed, curword;
  unsigned int *cur;
  unsigned int nbits_req_container, token, nbits2;

  bitPackInWord = 32;
  
  cur = z;
  curword = *cur;
  ufld_tmp = (int *) malloc(ni*nj*sizeof(int));
  
  extract(nbits_req_container, cur, 32, istep, curword, bitPackInWord); 
  
  for (i=1; i <= ni; i++)
   {
   k = FTN2C(i,1,ni);
   extract(token, cur, 32, nbits, curword, bitPackInWord); 
   ufld[k] = token;
   }
   
  for (j=2; j <= nj; j++)
   {
   k = FTN2C(1,j,ni);
   extract(token, cur, 32, nbits, curword, bitPackInWord); 
   ufld[k] = token;
   }
  
  
  for (j=2; j <= nj; j+=istep)
    {
    lcl_n = ((j + istep - 1) >= nj ? nj - j : istep - 1);
    for (i=2; i <= ni; i+=istep)
      {
      lcl_m = ((i + istep - 1) >= ni ? ni - i : istep - 1);
      extract(nbits_needed, cur, 32, nbits_req_container, curword, bitPackInWord); 
      switch (nbits_needed)
        {
        case 0:
        for (n=0; n <= lcl_n; n++)
          {
          for (m=0; m <= lcl_m; m++)
            {
            k = FTN2C(i+m,j+n,ni);
            ufld_tmp[k] = 0;
            }
          }
        break;
        
        default:
        nbits2 = nbits_needed + 1;
        for (n=0; n <= lcl_n; n++)
          {
          for (m=0; m <= lcl_m; m++)
            {
            k = FTN2C(i+m,j+n,ni);
            extract(token, cur, 32, nbits2, curword, bitPackInWord);
            ufld_tmp[k] = token; 
            ufld_tmp[k] = (ufld_tmp[k] << (32-nbits2)) >> (32-nbits2);
            }  
          }
         } 
                
        }
      }
      
  for (j=2; j<=nj; j++)
   {
   for (i=2; i <=ni; i++)
      {
      k11 = FTN2C(i-1,j-1,ni);
      k12 = FTN2C(i-1,j  ,ni);
      k21 = FTN2C(i,  j-1,ni);
      k22 = FTN2C(i,  j,  ni);
      ufld[k22] =  ufld_tmp[k22] + (ufld[k21]+ufld[k12]-ufld[k11]);
      }
   }  
  
      free(ufld_tmp);
    }   

void packTokensParallelogram_8(unsigned int z[], unsigned int *zlng, unsigned char ufld[], int ni, int nj, int nbits, int istep)
{

  float rlog2;
  int *ufld_dst;
  int k22, nbits2;
  int lcl_m, lcl_n;
  int local_max;
  unsigned int *cur;
  unsigned int i, j, k, m, n;
  unsigned int lastWordShifted, spaceInLastWord, lastSlot;
  unsigned int lcl, nbits_needed;
  unsigned int nbits_req_container, token;
  
  lastSlot = 0;
  cur = z;
  
if (once == 0)
   {
   rlog2 = 1.0/log(2.0);
   for (i=0; i < 256; i++)
      {
      fastlog[i] = (int)(1+log(i+0.5)*rlog2);
      }
   once = 1;
   }
  
  ufld_dst=(int *) malloc(ni*nj*sizeof(int));
  
  for (j=1; j <= nj; j++)
   {
   k = FTN2C(1,j,ni);
   ufld_dst[k] = 0;
   } 
   
  for (i=1; i <= ni; i++)
   {
   k = FTN2C(i,1,ni);
   ufld_dst[k] = 0;
   }
   
  for (j=2; j<=nj; j++)
   {
   for (i=2; i <=ni; i++)
      {
      k22 = FTN2C(i,  j,  ni);
      ufld_dst[k22] = ufld[k22] - (ufld[k22-ni]+ufld[k22-1]-ufld[k22-1-ni]);
      }
   }  
   
  nbits_req_container = 4;


  lastWordShifted = 0;
  spaceInLastWord = 32;
  *cur = 0;
  
  stuff(nbits_req_container, cur, 32, istep, lastWordShifted, spaceInLastWord);
  
  for (i=1; i <= ni; i++)
   {
   k = FTN2C(i,1,ni);
   stuff(ufld[k], cur, 32, nbits, lastWordShifted, spaceInLastWord);
   }
   
  for (j=2; j <= nj; j++)
   {
   k = FTN2C(1,j,ni);
   stuff(ufld[k], cur, 32, nbits, lastWordShifted, spaceInLastWord);
   }
   
  for (j=2; j <= nj; j+=istep)
    {
    lcl_n = ((j + istep - 1) >= nj ? nj - j : istep - 1);
    for (i=2; i <= ni; i+=istep)
      {
      k = FTN2C(i,j,ni);
      local_max = ufld_dst[k];
      lcl_m = ((i + istep - 1) >= ni ? ni - i : istep - 1);
      for (n=0; n <= lcl_n; n++)
        {
        for (m=0; m <= lcl_m; m++)
          {
          k = FTN2C(i+m,j+n,ni);
          if (local_max < abs(ufld_dst[k])) local_max = abs(ufld_dst[k]);
          }
        }
      if (local_max == 0)
        {
        nbits_needed = 0;
        }
      else
        {
        if (local_max < 256)
         {
         nbits_needed = fastlog[local_max];
         }
        else
         {
         nbits_needed = 8 + fastlog[local_max>>8];
         }
        }
      stuff(nbits_needed, cur, 32, nbits_req_container, lastWordShifted, spaceInLastWord);
      switch (nbits_needed)
        {
        case 0:
        break;
        
    
        default:
        nbits2 = nbits_needed + 1;
        for (n=0; n <= lcl_n; n++)
          {
          for (m=0; m <= lcl_m; m++)
            {
            k = FTN2C(i+m,j+n,ni);
            token = (unsigned int) (ufld_dst[k] & ~((-1)<<nbits2));
            stuff(token, cur, 32, nbits2, lastWordShifted, spaceInLastWord);
            }
          }
        break;
        } 
       
       }
    }

  
    lcl = 0;
    stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord);
    stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord); 

   *zlng = 1 + (int) (cur-z) * 4;
    free(ufld_dst);

}


void unpackTokensParallelogram_8(unsigned char ufld[], unsigned int z[], int ni, int nj, int nbits, int istep)
{

  int *ufld_tmp;
  int bitPackInWord;
  int i, j, k, m, n;
  int k11, k12, k21, k22;
  int lcl_m, lcl_n;
  unsigned int  nbits_needed, curword;
  unsigned int *cur;
  unsigned int nbits_req_container, token, nbits2;

  bitPackInWord = 32;
  
  cur = z;
  curword = *cur;
  ufld_tmp = (int *) malloc(ni*nj*sizeof(int));
  
  extract(nbits_req_container, cur, 32, istep, curword, bitPackInWord); 
  
  for (i=1; i <= ni; i++)
   {
   k = FTN2C(i,1,ni);
   extract(token, cur, 32, nbits, curword, bitPackInWord); 
   ufld[k] = token;
   }
   
  for (j=2; j <= nj; j++)
   {
   k = FTN2C(1,j,ni);
   extract(token, cur, 32, nbits, curword, bitPackInWord); 
   ufld[k] = token;
   }
  
  
  for (j=2; j <= nj; j+=istep)
    {
    lcl_n = ((j + istep - 1) >= nj ? nj - j : istep - 1);
    for (i=2; i <= ni; i+=istep)
      {
      lcl_m = ((i + istep - 1) >= ni ? ni - i : istep - 1);
      extract(nbits_needed, cur, 32, nbits_req_container, curword, bitPackInWord); 
      switch (nbits_needed)
        {
        case 0:
        for (n=0; n <= lcl_n; n++)
          {
          for (m=0; m <= lcl_m; m++)
            {
            k = FTN2C(i+m,j+n,ni);
            ufld_tmp[k] = 0;
            }
          }
        break;
        
        default:
        nbits2 = nbits_needed + 1;
        for (n=0; n <= lcl_n; n++)
          {
          for (m=0; m <= lcl_m; m++)
            {
            k = FTN2C(i+m,j+n,ni);
            extract(token, cur, 32, nbits2, curword, bitPackInWord);
            ufld_tmp[k] = token; 
            ufld_tmp[k] = (ufld_tmp[k] << (32-nbits2)) >> (32-nbits2);
            }  
          }
         } 
                
        }
      }
      
  for (j=2; j<=nj; j++)
   {
   for (i=2; i <=ni; i++)
      {
      k11 = FTN2C(i-1,j-1,ni);
      k12 = FTN2C(i-1,j  ,ni);
      k21 = FTN2C(i,  j-1,ni);
      k22 = FTN2C(i,  j,  ni);
      ufld[k22] =  ufld_tmp[k22] + (ufld[k21]+ufld[k12]-ufld[k11]);
      }
   }  
  
      free(ufld_tmp);
    }   
        
void pack1bitRLE(unsigned int z[], unsigned int *zlng, unsigned char ufld[], int npts)
{
  unsigned int i, j;
  unsigned int lastWordShifted, spaceInLastWord, lastSlot;

  unsigned int *cur;
  unsigned int lcl, indx, last_indx;
  unsigned char lcl_count;
  int count, limite, repeat;

  lastSlot = 0;
  cur = z;
   
  lastWordShifted = 0;
  spaceInLastWord = 32;
  *cur = 0;
  last_indx=0;
  indx=1;
  while (indx <= npts)
    {
    while (ufld[indx] == ufld[last_indx] && indx < npts)
      {
      indx++;
      }
    count = indx - last_indx; 
    if (count < 8)
      {
      stuff(SEQUENCE, cur, 32, 1, lastWordShifted, spaceInLastWord);
      limite = (last_indx + 7) > npts ? (npts - last_indx) : 7;
      for (i=0; i < limite; i++)
        {
        stuff(ufld[last_indx+i], cur, 32, 1, lastWordShifted, spaceInLastWord);
        }
      last_indx +=7;
      indx = last_indx + 1;
      }
    else
      {
      i = 0;
      repeat = 0;
      while (i < count)
        {
        lcl_count = (count-i) >= 63 ? 62 : (count - i);
        if (lcl_count < 8)
          {
          stuff(SEQUENCE, cur, 32, 1, lastWordShifted, spaceInLastWord);
          limite = (last_indx + 7) > npts ? (npts - last_indx) : 7;
          for (j=0; j < limite; j++)
            {
            stuff(ufld[last_indx+j], cur, 32, 1, lastWordShifted, spaceInLastWord);
            }
          last_indx +=7;
          indx = last_indx + 1;
          }
        else
          {
          if (lcl_count == 62)
            {
            if ((count - i) > 256 && (repeat == 1))
              {
              lcl_count = 0xFF;
              stuff(lcl_count, cur, 32, 8,lastWordShifted, spaceInLastWord);
              last_indx+=lcl_count;
              indx = last_indx + 1;
              }
            else
              {
              stuff(COUNT, cur, 32, 1, lastWordShifted, spaceInLastWord);
              stuff(ufld[last_indx], cur, 32, 1,lastWordShifted, spaceInLastWord);
              stuff(lcl_count, cur, 32, 6,lastWordShifted, spaceInLastWord);
              last_indx+=lcl_count;
              indx = last_indx + 1;
              repeat = 1;
              }
            }
          else
            {
            stuff(COUNT, cur, 32, 1, lastWordShifted, spaceInLastWord);
            stuff(ufld[last_indx], cur, 32, 1,lastWordShifted, spaceInLastWord);
            stuff(lcl_count, cur, 32, 6,lastWordShifted, spaceInLastWord);
            last_indx+=lcl_count;
            indx = last_indx + 1;
            }
          }
        i += lcl_count;

        }
      }  
    }
    

  lcl = 0;
  stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord);
  stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord); 

   *zlng = 1 + (int) (cur-z) * 4;
}

void unpack1bitRLE(unsigned char ufld[], unsigned int z[], unsigned int *zlng,  int npts)
{
  unsigned int i, j;

  unsigned int *cur, seq_type, val, last_val;
  int count, limite;

  int bitPackInWord;

  unsigned int curword;
  unsigned int token;

  bitPackInWord = 32;
  
  cur = z;
  curword = *cur;
  last_val = 0xFFFFFFFF;
  
  i = 0;
  while (i < npts)
    {
    extract(seq_type, cur, 32, 1, curword, bitPackInWord); 
    switch(seq_type)
      {
      case SEQUENCE:
      limite = (i + 7) > npts ? (npts - i) : 7;
      for (j=0; j < limite; j++)
        {
        extract(token, cur, 32, 1, curword, bitPackInWord); 
        ufld[i+j] = (unsigned char) token;
        }
      i+=limite;
      break;
      
      case COUNT:
      extract(val, cur, 32, 1, curword, bitPackInWord); 
      extract(count, cur, 32, 6, curword, bitPackInWord); 
      switch (count)
        {
        case 63:
        for (j=0; j < 255; j++)
          {
          ufld[i+j] = (unsigned char) last_val;
          }
        i+=255;
        break;
        
        default:
        for (j=0; j < count; j++)
          {
          ufld[i+j] = (unsigned char) val;
          }
        i+=count;
        last_val = val;
        break;
        }
      break;
      }
    }
}


int compact_mask_char(unsigned int *dest, unsigned char *src, int npts)
  {
  int i,entier, fraction,npts32;

  npts32 = 1 + (npts >> 5);

  for (i=0; i < npts32; i++)
    {
    dest[i] = 0;
    }

  for (i=0; i < npts; i++)
    {
    entier = i >> 5;
    fraction = i - (entier << 5);
    dest[entier] |= (src[i] << fraction);
    }
  return 0;
  }

int uncompact_mask_char(int *dest, unsigned int *src, int npts)
  {
  int i,entier, fraction;

  for (i=0; i < npts; i++)
    {
    entier = i >> 5;
    fraction = i - (entier << 5);
    dest[i] = (src[entier]  & (1 << fraction)) >> fraction;
    }
  return 0;
  }

void pack_stream_nbits_32(unsigned int z[], unsigned int *zlng, unsigned int ufld[], int npts, int nbits)
{
  unsigned int i;
  unsigned int lastWordShifted, spaceInLastWord, lastSlot;
  unsigned int *cur;

  unsigned int lcl;
  lastSlot = 0;
  cur = z;
   
  lastWordShifted = 0;
  spaceInLastWord = 32;
  *cur = 0;
  for (i=0; i < npts; i++)
    {
    stuff(ufld[i], cur, 32, nbits, lastWordShifted, spaceInLastWord);
    }
  lcl = 0;
  stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord);
  stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord); 

  *zlng = 1 + (int) (cur-z) * 4;
}

void unpack_stream_nbits_32(unsigned int ufld[], unsigned int z[], int npts, int nbits)
{
  unsigned int i;
  unsigned int lastSlot;

  unsigned int *cur;

  int bitPackInWord;
  unsigned int curword;

  lastSlot = 0;
  cur = z;
   
  bitPackInWord = 32;
  
  curword = *cur;
  for (i=0; i < npts; i++)
    {
    extract(ufld[i], cur, 32, nbits, curword, bitPackInWord); 
    }
}

void pack_stream_nbits_16(unsigned int z[], unsigned int *zlng, unsigned short ufld[], int npts, int nbits)
{
  unsigned int i;
  unsigned int lastWordShifted, spaceInLastWord, lastSlot;
  unsigned int *cur;

  unsigned int lcl;
  lastSlot = 0;
  cur = z;
   
  lastWordShifted = 0;
  spaceInLastWord = 32;
  *cur = 0;
  for (i=0; i < npts; i++)
    {
    stuff(ufld[i], cur, 32, nbits, lastWordShifted, spaceInLastWord);
    }
  lcl = 0;
  stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord);
  stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord); 

  *zlng = 1 + (int) (cur-z) * 4;
}

void unpack_stream_nbits_16(unsigned short ufld[], unsigned int z[], int npts, int nbits)
{
  unsigned int i;
  unsigned int lastSlot;

  unsigned int *cur;

  int bitPackInWord;
  unsigned int curword;

  lastSlot = 0;
  cur = z;
   
  bitPackInWord = 32;
  
  curword = *cur;
  for (i=0; i < npts; i++)
    {
    extract(ufld[i], cur, 32, nbits, curword, bitPackInWord); 
    }
}
void pack_stream_nbits_8(unsigned int z[], unsigned int *zlng, unsigned char ufld[], int npts, int nbits)
{
  unsigned int i;
  unsigned int lastWordShifted, spaceInLastWord, lastSlot;
  unsigned int *cur;

  unsigned int lcl;
  lastSlot = 0;
  cur = z;
   
  lastWordShifted = 0;
  spaceInLastWord = 32;
  *cur = 0;
  for (i=0; i < npts; i++)
    {
    stuff(ufld[i], cur, 32, nbits, lastWordShifted, spaceInLastWord);
    }
  lcl = 0;
  stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord);
  stuff(lcl, cur, 32, 16, lastWordShifted, spaceInLastWord); 

  *zlng = 1 + (int) (cur-z) * 4;
}

void unpack_stream_nbits_8(unsigned char ufld[], unsigned int z[], int npts, int nbits)
{
  unsigned int i;
  unsigned int lastSlot;

  unsigned int *cur;

  int bitPackInWord;
  unsigned int curword;

  lastSlot = 0;
  cur = z;
   
  bitPackInWord = 32;
  
  curword = *cur;
  for (i=0; i < npts; i++)
    {
    extract(ufld[i], cur, 32, nbits, curword, bitPackInWord); 
    }
}
