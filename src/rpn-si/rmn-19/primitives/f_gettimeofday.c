#ifndef WIN32	/*CHC/NRC*/
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#else
#include <time.h>
#include <windows.h>
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif
 
struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};
 
int gettimeofday(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;
 
  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);
 
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
 
    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    tmpres /= 10;  /*convert into microseconds*/
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }
 
  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
    tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }
 
  return 0;
}
#endif

#include <rpnmacros.h>

/* time of day in seconds.microseconds */
double f77name(f_gettimeofday)()
{
  struct timeval tvnow;
  struct timezone tz;
  int ier;
  double t1;

  ier = gettimeofday(&tvnow,&tz);
  if (ier != 0) printf("gettimeofday error: ier=%d\n",ier);
  t1 = tvnow.tv_usec * 1.0e-6 ;
  t1 = tvnow.tv_sec + t1;
  return(t1);
}
/* time of day in microseconds */
long long f77name(f_gettimeofday_micro)()
{
  struct timeval tvnow;
  struct timezone tz;
  int ier;
  long long t1;

  ier = gettimeofday(&tvnow,&tz);
  if (ier != 0) printf("gettimeofday error: ier=%d\n",ier);
  t1 = tvnow.tv_sec ;
  t1 = t1 * 1000000;
  t1 = t1 + tvnow.tv_usec ;
  return(t1);
}
