#include <stdlib.h>
#include <stdint.h>

#include "crc16_.h"
#include "crc24_.h"
#include "crc32_.h"


static int one = 1;
static char *little_endian = (char *) &one;

uint32_t update_crc_ne(uint32_t seed, int crclen, void *data, int datasiz, int datalen, int mode);

/* datasiz = 1/2/4/8 (in bytes) */
/* datalen = number of data elements (datalen*datasiz bytes) */
/* crclen  = 16/24/32 crc type desired */
/* data    = pointer to data area for which crc is desired */
/* seed    = old crc (must be 0 at the very beginning) */
/* mode    = 0/1/2  0=init, 1=update, 2=finalize */
uint32_t update_crc_ne(uint32_t old_crc, int crclen, void *data, int datasiz, int datalen, int mode)
{
  uint32_t seed=old_crc;
  uint32_t new_crc;
  crc16_t crc16;
  crc24_t crc24;
  crc32_t crc32;
  unsigned int mask=0 ;
  int data_len=datasiz*datalen ;
  
  new_crc = seed ;
    
  if( (*little_endian) & (datasiz>1)){  /* multibyte tokens on little endian machine */
    mask = datasiz-1 ; /* 0/1/3/7 magic mask to get stitching order count on little endian machines */
  }
  if (crclen==16) {
    crc16 = seed & 0xFFFF ;
    if(mode==0) crc16 = crc16_init() ;
    crc16 = crc16_update_le(crc16, (const unsigned char *) data, data_len, mask) ;
    if(mode==2) crc16 = crc16_finalize(crc16) ;
    new_crc = crc16 ;
  } else if (crclen==24) {
    crc24 = seed & 0xFFFFFF ;
    if(mode==0) crc24 = crc24_init() ;
    crc24 = crc24_update_le(crc24, (const unsigned char *) data, data_len, mask) ;
    if(mode==2) crc24 = crc24_finalize(crc24) ;
    new_crc = crc24 ;
  } else if (crclen==32) {
    crc32 = seed & 0xFFFFFFFF ;
    if(mode==0) crc32 = crc32_init() ;
    crc32 = crc32_update_le(crc32, (const unsigned char *) data, data_len, mask) ;
    if(mode==2) crc32 = crc32_finalize(crc32) ;
    new_crc = crc32 ;
  } else {
    new_crc = seed ;
  }
  return new_crc;
}

#pragma weak f_update_crc_ne_  = f_update_crc_ne
#pragma weak f_update_crc_ne__ = f_update_crc_ne
uint32_t f_update_crc_ne_(uint32_t *old_crc, int *crclen, void *data, int *datasiz, int *datalen, int *mode);
uint32_t f_update_crc_ne__(uint32_t *old_crc, int *crclen, void *data, int *datasiz, int *datalen, int *mode);
uint32_t f_update_crc_ne(uint32_t *old_crc, int *crclen, void *data, int *datasiz, int *datalen, int *mode)
{
  return update_crc_ne(*old_crc, *crclen, data, *datasiz, *datalen, *mode);
}

#ifdef SELFTEST
#include <stdio.h>
int main(int argc, char** argv)
{
  unsigned int i,mask ;
  unsigned long long lmask ;
  unsigned int iarray[32];
  unsigned long long larray [64];
  crc16_t crc16=0;
  crc24_t crc24=0;
  crc32_t crc32=0;
  unsigned int p16, p24, p32;

  mask = 1;
  lmask = 1;
  for (i=0;i<32;i++) {iarray[i] = mask  ; mask <<= 1 ; }
  for (i=0;i<64;i++) {larray[i] = lmask ; lmask <<= 1 ; }
  p16=crc16 ; p24=crc24 ; p32=crc32 ;
  printf("crc16=%8.8x, crc24=%8.8x, crc32=%8.8x\n",p16,p24,p32) ;
  crc16=update_crc_ne(crc16, 16, iarray, 4, 32, 0);
  crc24=update_crc_ne(crc24, 24, iarray, 4, 32, 0);
  crc32=update_crc_ne(crc32, 32, iarray, 4, 32, 0);
  p16=crc16 ; p24=crc24 ; p32=crc32 ;
  printf("iarray , %8.8x, %8.8x\n",iarray[0],iarray[31]);
  printf("crc16=0000f39a, crc24=006b7cac, crc32=e3b736e4 is the correct answer\n") ;
  printf("crc16=%8.8x, crc24=%8.8x, crc32=%8.8x\n",p16,p24,p32) ;
  crc16=update_crc_ne(crc16, 16, larray, 8, 64, 2);
  crc24=update_crc_ne(crc24, 24, larray, 8, 64, 2);
  crc32=update_crc_ne(crc32, 32, larray, 8, 64, 2);
  p16=crc16 ; p24=crc24 ; p32=crc32 ;
  printf("larray , %16.16Lx, %16.16Lx\n",larray[0],larray[63]);
  printf("crc16=00000fc6, crc24=0061f67c, crc32=5560ac91 is the correct answer\n") ;
  printf("crc16=%8.8x, crc24=%8.8x, crc32=%8.8x\n",p16,p24,p32) ;
  return(0);
}
#endif
