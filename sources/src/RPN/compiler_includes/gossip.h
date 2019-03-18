#include <rpnmacros.h>

#define INT_16 short

#define swap_2(mot) { register unsigned INT_16 tmp =(unsigned INT_16)mot; \
   mot = (tmp>>8) | (tmp<<8) ; }

#define swap_4(mot) { register unsigned INT_32 tmp =(unsigned INT_32)mot; \
   mot = (tmp>>24) | (tmp<<24) | ((tmp>>8)&0xFF00) | ((tmp&0xFF00)<<8); }

#define swap_8(mot) { register unsigned INT_64 tmp1; register unsigned INT_64 tmp2; register unsigned INT_32 mot1;\
   tmp1 = ((unsigned INT_64) mot << 32) >> 32 ; \
   tmp2 = (unsigned INT_64) mot >> 32 ; \
   mot  = ( (tmp1>>24) | (tmp1<<24) | ((tmp1>>8)&0xFF00) | ((tmp1&0xFF00)<<8) ) << 32 ; \
   mot1 = ( (tmp2>>24) | (tmp2<<24) | ((tmp2>>8)&0xFF00) | ((tmp2&0xFF00)<<8) ); \
   mot |= mot1; }
   

#define swap_4_4(mot1,mot2) { register unsigned INT_32 tmp1 = (unsigned INT_32)mot1; \
                              register unsigned INT_32 tmp2 = (unsigned INT_32)mot2; \
     mot2 = (tmp1>>24) | (tmp1<<24) | ((tmp1>>8)&0xFF00) | ((tmp1&0xFF00)<<8); \
     mot1 = (tmp2>>24) | (tmp2<<24) | ((tmp2>>8)&0xFF00) | ((tmp2&0xFF00)<<8); }

/* macro to swap the bytes in a 32-bit variable */
#define swapbytes32(x) \
{ \
    unsigned int data = *(unsigned int*)&(x); \
    data = ((data & 0xff000000) >> 24) |    \
           ((data & 0x00ff0000) >>  8) |    \
           ((data & 0x0000ff00) <<  8) |    \
           ((data & 0x000000ff) << 24);     \
    *(unsigned int*)&(x) = data;            \
}

/* macro to swap the bytes in a 64-bit variable */
#define swapbytes64(x) \
{ \
    unsigned int *words = (unsigned int *)&(x); \
    unsigned int temp0;  \
    unsigned int temp1;  \
    temp0 = words[0];    \
    swapbytes32(temp0);  \
    temp1 = words[1];    \
    swapbytes32(temp1);  \
    words[1] = temp0;    \
    words[0] = temp1;    \
}

#define ONE_BYTE    1
#define TWO_BYTES   2
#define FOUR_BYTES  4
#define EIGHT_BYTES 8


#define MAX_CLIENTS 128
#define MAX_EXTENDED_CLIENTS 128

typedef struct {
  int uid;
  int pid;
  int socket;
  int client_id;
  char * command;
  } CLIENT_SLOT ;

typedef struct {
  int uid;
  int pid;
  int socket;
  int client_id;
  char * command;
  void *data;
  void (*user_function)(void *);
  } EXTENDED_CLIENT_SLOT ;

typedef struct {
  char *name;
  void (*function)(CLIENT_SLOT *);
  } TABLE_SLOT;

typedef struct {
  char *name;
  void (*function)(EXTENDED_CLIENT_SLOT *);
  } EXTENDED_TABLE_SLOT;

typedef struct {
  int fd;
  unsigned char *buffer;
  unsigned char *in;
  unsigned char *out;
  unsigned char *limit;
  unsigned int BufferSize;
  unsigned int RecLen;
  unsigned int Log2Siz;
  unsigned char flags[4];
  } gossip_stream;

int set_host_and_port(char *channel_file, char *host_and_port)  ;
char *get_host_and_port(char *channel_file)  ;
char *get_broker_Authorization()  ;
void set_broker_Authorization(int auth_token)  ;
int accept_from_sock(int fserver)  ;
int bind_sock_to_port(int s)  ;
int get_sock_net()  ;
int set_sock_opt(int s)  ;
int get_ip_address(char *hostname)  ;
int get_own_ip_address()  ;
int connect_to_hostport(char *target)  ;
int connect_to_localport(int port)  ;
int bind_to_localport(int *port, char *buf, int maxbuf)  ;
void send_ack_nack(int fclient,int status)  ;
int get_ack_nack(int fserver)  ;
int send_command_to_server(int fserver, char *buf)  ;
INT_32 get_int32_from_channel(int channel)  ;
void put_int32_to_channel(int channel, INT_32 to_send)  ;
int connect_to_channel_by_name(char *name)  ;
void gossip_fork_server(char *LOGFILE, char *channel, int (*user_client)(int,int,int,int,char *), int PING_INTERVAL, int from_inetd)   ;
int init_gossip_stream(gossip_stream *s,int fd,int bufsz)  ;
int fill_gossip_stream(gossip_stream *s)  ;
int flush_gossip_stream(gossip_stream *s)  ;
int gossip_record_head(gossip_stream *s)  ;
int gossip_record_read(gossip_stream *s, unsigned char *where, int ToRead)  ;
int gossip_record_skip(gossip_stream *s)  ;
int gossip_record_write(gossip_stream *s, unsigned char *where, int ToWrite, int Log2Size)  ;
void gossip_copy(unsigned char *from, unsigned char *to, int many, int log2siz)  ;
void set_exit_requested()   ;
void exit_from_client_thread(EXTENDED_CLIENT_SLOT *client)  ;
void start_client_module(void (*client_address)(CLIENT_SLOT *), int client_uid, int client_pid, int fclient,char *command)  ;
void start_client_thread(void (*client_address)(CLIENT_SLOT *), int client_uid, int client_pid, int fclient,char *command)   ;
void increment_client_count()  ;
void decrement_client_count()  ;
int get_client_count()  ;
void reset_timeout_counter()  ;
void increment_timeout_counter()  ;
void decrement_timeout_counter()  ;
int set_timeout_counter(int timeout_value)  ;
int get_timeout_counter()  ;
