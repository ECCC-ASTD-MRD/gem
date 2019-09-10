/* mgi.h */

#define MAX_CHANNELS 10
#define MAX_NAME 41
#define MAX_STR 1024
#define BUFSIZE 40960

typedef struct {
        int fd_data;
        int fd_sig;
        int msgno_W;
        int msgno_R;
        int nblks;
        char name[MAX_NAME];
        char mode;
        int *buffer;
        int pos;
        } channel;

