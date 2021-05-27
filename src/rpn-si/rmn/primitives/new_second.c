#include <sys/time.h>
#include <sys/resource.h>

#if !defined(WIN32)
float second()
{
struct rusage mydata;
float temp;

getrusage(RUSAGE_SELF,&mydata);
/*fprintf(stderr,"%d %d %d %d\n",mydata.ru_utime.tv_sec,mydata.ru_stime.tv_sec, mydata.ru_utime.tv_usec, mydata.ru_stime.tv_usec); */
temp = mydata.ru_utime.tv_sec + mydata.ru_stime.tv_sec + (.000001 * (mydata.ru_utime.tv_usec + mydata.ru_stime.tv_usec));
return (temp);
}

float second_()
{
return (second());
}

float second__()
{
return (second());
}
#endif