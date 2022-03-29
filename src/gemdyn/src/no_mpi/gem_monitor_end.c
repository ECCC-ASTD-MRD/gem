#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <signal.h>
#include <fcntl.h>

void main_gem_monitor_end (int argc, char **argv) {
   int i, fd;
   int count=30;
   char buffer[32768];
   pid_t pp=getppid();

   if(argc-1 != 1) {
      printf("argument count must be 1 \n"); exit(1);
   }

   while(count-- >=0) {
      if(kill(pp,0)) exit(0);
      if( (fd=open(argv[1],O_RDONLY )) < 0 ) exit (0);
      close(fd);
      sleep(1);
   }
}
