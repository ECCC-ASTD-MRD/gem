static int value = 1;
static unsigned char *endian = (unsigned char *)&value;

#pragma weak is_little_endian__=is_little_endian
#pragma weak is_little_endian_=is_little_endian
int is_little_endian__();
int is_little_endian_();
int is_little_endian() { return(*endian) ; }

#pragma weak is_big_endian__=is_big_endian
#pragma weak is_big_endian_=is_big_endian
int is_big_endian__();
int is_big_endian_();
int is_big_endian() { return(1 - *endian) ; }

#if defined(SELF_TEST)
#include <stdio.h>
int main()
{
  printf("this machine is %s little endian \n",is_little_endian()?"":"not");
  printf("this machine is %s big endian \n",is_big_endian()?"":"not");
  return(0);
}
#endif