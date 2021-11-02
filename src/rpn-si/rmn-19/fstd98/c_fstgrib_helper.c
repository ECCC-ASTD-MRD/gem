#include <stdio.h>
#include <string.h>
#include <rpnmacros.h>

/* 

Routines commodes pour ecrire des records GRIB encodes dans des fichiers standars
Yves Chartier - Juin 2006

*/

extern ftnword calc_crc(unsigned char *p, ftnword *flen, ftnword *fseed, ftnword stride);
extern long long f77name(f_gettimeofday_micro)();
extern void f77name(f_bits_put)(ftnword *bit_array, ftnword *bits_per_slice, ftnword *slices, ftnword *nslices);
extern void f77name(f_bits_get)(ftnword *bit_array, ftnword *bits_per_slice, ftnword *slices, ftnword *nslices);
void grb_84bits_to_ig1234(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *the_84bit_token);
void grb_84bits_to_ip123(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *the_84bit_token);
void c_84bits_ig_get(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *the_84bit_token);
void c_84bits_ig_put(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *the_84bit_token);
void c_84bits_ip_get(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *the_84bit_token);
void c_84bits_ip_put(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *the_84bit_token);
void c_igaip84(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4);
void c_ipaig84(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *ip1, ftnword *ip2, ftnword *ip3);
void f77name(f_igaip84)(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4);
void f77name(f_ipaig84)(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *ip1, ftnword *ip2, ftnword *ip3);
void c_84bits_token(ftnword *the_84bit_token, unsigned char *grib_header, ftnword length_grib_header);
void f77name(f_def_84bitkey)(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *ig1, ftnword *ig2, 
  ftnword *ig3, ftnword *ig4, unsigned char *grib_header, int *len_grib_header);
void c_def_84bitkey(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *ig1, ftnword *ig2, 
  ftnword *ig3, ftnword *ig4, unsigned char *grib_header, int len_grib_header);


void f77name(f_def_84bitkey)(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *ig1, ftnword *ig2, 
  ftnword *ig3, ftnword *ig4, unsigned char *grib_header, int *len_grib_header)
  {
  c_def_84bitkey(ip1, ip2, ip3, ip1, ig2, ig3, ig4, grib_header, *len_grib_header);
  }

void c_def_84bitkey(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *ig1, ftnword *ig2, 
  ftnword *ig3, ftnword *ig4, unsigned char *grib_header, int len_grib_header)
  {
  ftnword the_84bit_token[3];
  c_84bits_token(the_84bit_token, grib_header, len_grib_header);
  c_84bits_ip_get(ip1, ip2, ip3, the_84bit_token);
  c_84bits_ig_get(ig1, ig2, ig3, ig4, the_84bit_token);
  }


void c_84bits_token(ftnword *the_84bit_token, unsigned char *grib_header, ftnword length_grib_header)
  {
  long long time_of_day_micro;
  unsigned ftnword time_of_day, micro_secs;
  ftnword slices[3];
  ftnword nslices = 3;
  ftnword fseed = 0;
  ftnword bits_per_slice[] = {32, 32, 20};
  ftnword header_crc;
  
  time_of_day_micro = f77name(f_gettimeofday_micro)();
  
  header_crc = calc_crc(grib_header, &length_grib_header, &fseed, 1);
/*  printf("%d crc: %d timeofday: %lld\n", strlen(grib_header), header_crc, time_of_day_micro);*/
  
  micro_secs = time_of_day_micro % 1000000;
  time_of_day = (time_of_day_micro - micro_secs) / 1000000;
/*   printf("%d %d %lld %lld\n", time_of_day, micro_secs, time_of_day_micro, (long long)(time_of_day * 1000000 + micro_secs)); */

  slices[0] = header_crc;
  slices[1] = time_of_day;
  slices[2] = micro_secs;
  f77name(f_bits_put)((ftnword *)the_84bit_token, (ftnword *)bits_per_slice, (ftnword *)slices, (ftnword *)&nslices);

  }

void grb_84bits_to_ip123(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *the_84bit_token)
  {
  ftnword ip_bits_per_slice[] = {28, 28, 28};
  ftnword ip_slices[3];
  ftnword nslices = 3;
   
  f77name(f_bits_get)((ftnword *)the_84bit_token, (ftnword *)ip_bits_per_slice, (ftnword *)ip_slices, (ftnword *)&nslices);
  *ip1 = ip_slices[0];
  *ip2 = ip_slices[1];
  *ip3 = ip_slices[2];
   
  }

void grb_84bits_to_ig1234(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *the_84bit_token)
  {
  ftnword ig_bits_per_slice[] = {21, 21, 21, 21};
  ftnword ig_slices[4];
  ftnword nslices = 4;
   
  f77name(f_bits_get)((ftnword *)the_84bit_token, (ftnword *)ig_bits_per_slice, (ftnword *)ig_slices, (ftnword *)&nslices);
  *ig1 = ig_slices[0];
  *ig2 = ig_slices[1];
  *ig3 = ig_slices[2];
  *ig4 = ig_slices[3];
   
  }
  
void c_84bits_ig_get(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *the_84bit_token)
  {
  ftnword ig_bits_per_slice[] = {21, 21, 21, 21};
  ftnword ig_slices[4];
  ftnword nslices = 4;
   
  f77name(f_bits_get)((ftnword *)the_84bit_token, (ftnword *)ig_bits_per_slice, (ftnword *)ig_slices, (ftnword *)&nslices);
  *ig1 = ig_slices[0];
  *ig2 = ig_slices[1];
  *ig3 = ig_slices[2];
  *ig4 = ig_slices[3];
  }

void c_84bits_ig_put(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *the_84bit_token)
  {
  ftnword ig_bits_per_slice[] = {21, 21, 21, 21};
  ftnword ig_slices[4];
  ftnword nslices = 4;
   
  memset(the_84bit_token, 0, 3*sizeof(ftnword));
  ig_slices[0] = *ig1;
  ig_slices[1] = *ig2;
  ig_slices[2] = *ig3;
  ig_slices[3] = *ig4;
  f77name(f_bits_put)((ftnword *)the_84bit_token, (ftnword *)ig_bits_per_slice, (ftnword *)ig_slices, (ftnword *)&nslices);
  }

void c_84bits_ip_get(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *the_84bit_token)
  {
  ftnword ip_bits_per_slice[] = {28, 28, 28};
  ftnword ip_slices[3];
  ftnword nslices = 3;
   
  f77name(f_bits_get)((ftnword *)the_84bit_token, (ftnword *)ip_bits_per_slice, (ftnword *)ip_slices, (ftnword *)&nslices);
  *ip1 = ip_slices[0];
  *ip2 = ip_slices[1];
  *ip3 = ip_slices[2];
  }

void c_84bits_ip_put(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *the_84bit_token)
  {
  ftnword ip_bits_per_slice[] = {28, 28, 28};
  ftnword ip_slices[3];
  ftnword nslices = 3;
   
  memset(the_84bit_token, 0, 3*sizeof(ftnword));
  ip_slices[0] = *ip1;
  ip_slices[1] = *ip2;
  ip_slices[2] = *ip3;
  f77name(f_bits_put)((ftnword *)the_84bit_token, (ftnword *)ip_bits_per_slice, (ftnword *)ip_slices, (ftnword *)&nslices);
  }

void c_igaip84(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4)
  {
  ftnword the_84bit_token[3];
  
  c_84bits_ig_put(ig1, ig2, ig3, ig4, the_84bit_token);
  c_84bits_ip_get(ip1, ip2, ip3, the_84bit_token);
  }
  
void c_ipaig84(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *ip1, ftnword *ip2, ftnword *ip3)
  {
  ftnword the_84bit_token[3];
  
  c_84bits_ip_put(ip1, ip2, ip3, the_84bit_token);
  c_84bits_ig_get(ig1, ig2, ig3, ig4, the_84bit_token);
  }

void f77name(f_igaip84)(ftnword *ip1, ftnword *ip2, ftnword *ip3, ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4)
  {
  c_igaip84(ip1,ip2, ip3, ig1, ig2, ig3, ig4);
  }
  
void f77name(f_ipaig84)(ftnword *ig1, ftnword *ig2, ftnword *ig3, ftnword *ig4, ftnword *ip1, ftnword *ip2, ftnword *ip3)
  {
  c_ipaig84(ig1, ig2, ig3, ig4, ip1, ip2, ip3);
  }


