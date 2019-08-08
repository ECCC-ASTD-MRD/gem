#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

#pragma weak generate_unique_name__=generate_unique_name
int32_t generate_unique_name__(char *str, int32_t *sz);
#pragma weak generate_unique_name_=generate_unique_name
int32_t generate_unique_name_(char *str, int32_t *sz);
int32_t generate_unique_name(char *str, int32_t *sz){
  int32_t hid = gethostid();
  int32_t pid = getpid();
  int32_t t1;
  int64_t tf;
  struct timeval t;
  size_t szt = *sz;

  t1 = gettimeofday(&t, NULL);
  tf = t.tv_sec;
  tf *= 1000000;
  tf += t.tv_usec;
  t1 = snprintf(str, szt, "%8.8x%8.8x%16.16x" , hid, pid, tf);
  return t1 > *sz ? *sz : t1;
}
