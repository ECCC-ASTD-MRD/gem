#include <rpnmacros.h>
int c_logical_smt_cpus()
{
  char my_host_name[128];
  if (gethostname(my_host_name,128)) return (1);
  if((my_host_name[0]=='c') && (my_host_name[2]=='f') && (my_host_name[5]=='p') && (my_host_name[1]=='3')) return(2);
  if((my_host_name[0]=='c') && (my_host_name[2]=='f') && (my_host_name[5]=='p') && (my_host_name[1]=='4')) return(2);
  if((my_host_name[0]=='c') && (my_host_name[2]=='f') && (my_host_name[5]=='p') && (my_host_name[1]=='5')) return(2);
  return(1);
}
wordint f77_name(logical_smt_cpus)()
{
  return(c_logical_smt_cpus());
}
