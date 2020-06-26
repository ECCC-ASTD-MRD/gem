/* mgi.h */

#define MAX_CHANNELS 24
#define MAX_NAME 125
#define MAX_STR 1024
#define BUFSIZE 40960

#ifndef FALSE
#define FALSE            0
#define TRUE        !FALSE
#endif

typedef struct 
{
  int fd_data;
  int fd_sig;
  int msgno_W;
  int msgno_R;
  int nblks;
  char name[MAX_NAME];
  char mode;
  int *buffer;
  int pos;
  int gchannel;
} channel;
