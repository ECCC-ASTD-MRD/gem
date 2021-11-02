#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <rpnmacros.h>

ftnword f77name(get_max_rss)()
{
struct rlimit limits;
struct rusage mydata;

getrusage(RUSAGE_SELF,&mydata);
/* printf("maxrss=%d\n",mydata.ru_maxrss);  */
return (mydata.ru_maxrss);
}
