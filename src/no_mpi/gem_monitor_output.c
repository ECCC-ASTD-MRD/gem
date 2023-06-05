#include <unistd.h> 
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

void main_gem_monitor_output(int argc, char **argv){
int i, fd;
char buffer[32768];
pid_t pp=getppid();
time_t now;

if(argc-1 & 1) { printf("argument count must be even \n"); exit(1); }
if(fork()) exit(0);
while(1){
  if(kill(pp,0)) exit(0);
  for(i=1 ; i<argc ; i+=2){
    if( (fd=open(argv[i],O_RDONLY )) >= 0 ) {
      close(fd);
      printf("file=%s,cmd=%s\n",argv[i],argv[i+1]);
      snprintf(buffer,sizeof(buffer)-1,"%s %s",argv[i+1],argv[i]);
      now = time(NULL);
      printf("Executing:%s - %s",buffer,asctime(gmtime(&now)));
      system(buffer);
      now = time(NULL);
      printf("Deleting:%s - %s",argv[i],asctime(gmtime(&now)));
      fflush(stdout); 
      unlink(argv[i]);
      }
    }
  sleep(2);
  }
}
