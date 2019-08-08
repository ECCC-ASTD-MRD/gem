void memcpy_8_16(short *p16, char *p8, int nb)
{
  int i;
  for (i=0; i < nb; i++)
    *p16++ = *p8++;
}

void memcpy_16_8(char *p8, short *p16, int nb)
{
  int i;
  for (i=0; i < nb; i++)
    *p8++ = *p16++;
}


