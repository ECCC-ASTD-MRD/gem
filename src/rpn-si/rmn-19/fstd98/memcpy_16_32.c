#include <stdio.h>
void memcpy_16_32(int *p32, short *p16, int nbits, int n)
{
  int i;
  int mask;
  
  mask = ~ (-1 << nbits);
  for (i=0; i < n; i++)
    *p32++ = *p16++ & mask;
}

void memcpy_32_16(short *p16, int *p32, int nbits, int n)
{
  int i;
  int mask;
  
  mask = ~ (-1 << nbits);
  for (i=0; i < n; i++)
    *p16++ = *p32++ & mask;
}


