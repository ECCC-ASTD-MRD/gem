/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2000  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h> 
#include <strings.h> 
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <time.h>
#include <rpnmacros.h>

#include <gossip.h>


static int exit_requested = 0;   /* flag to indicate pending SHUTDOWN request */

void set_exit_requested()   /*   %ENTRY%   */
{
  exit_requested = 1;
}
int is_exit_requested()     /*   %ENTRY%   */
{
  return(exit_requested);
}


static EXTENDED_CLIENT_SLOT clients[MAX_EXTENDED_CLIENTS];
static int client_table_initialized = 0 ;  /* initialization flag */

/* initialize client table by setting everything to 0*/
static void initialize_client_table() 
{
  if(client_table_initialized) return;
  memset(&clients[0], 0, sizeof(clients));
  client_table_initialized = 1;
}

/* find a free slot in the client table. If table full, return -1 */
static int find_client_slot() 
{
  int i;
  for(i = 0;  clients[i].client_id != 0 && i<MAX_EXTENDED_CLIENTS ; i++);
  return( i < MAX_EXTENDED_CLIENTS ? i : -1 );
}

/* =================================================================== */

#ifdef HP
static int mutex_initialized = 0;
static pthread_mutex_t mutex ;
#else
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
#endif
#define MAX_LOCKS 256

static pthread_mutex_t locks[MAX_LOCKS] ;
static int locks_initialized = 0;

static int client_no = 0;        /* number of active clients */

/* initialize locks table from already initialized mutex */
static void initialize_locks()
{
  int i;
#ifdef HP
 /* mutex must be initialized by a call to pthread_mutex_init under HP-UX
     under IRIX, LINUX, SUPER-UX, a static initializer can be used */

  if( mutex_initialized == 0 )
    {
      pthread_mutex_init(&mutex, pthread_mutexattr_default);
      mutex_initialized == 1;
    }
#endif
  if(locks_initialized) 
    return;
  for (i=1 ; i<MAX_LOCKS ; i++) 
    {
#ifdef HP
      pthread_mutex_init(&locks[i], pthread_mutexattr_default);
#else
      memcpy(&locks[i], &mutex, sizeof(mutex));
#endif
    }
  locks_initialized = 1;
}

/* exit from a client serving thread, reinitialize client structure */
void exit_from_client_thread(EXTENDED_CLIENT_SLOT *client)  /*   %ENTRY%   */
{
  if(client->command)
    free(client->command);
  memset(client, 0, sizeof(EXTENDED_CLIENT_SLOT));
  pthread_exit(0);

}

static int client_ord = 0;       /* total number of clients served */

/* start a client module, non threaded */
void start_client_module_2(void (*client_address)(), int client_uid, int client_pid, int fclient, char *command, void *data)  /*   %ENTRY%   */
{
  EXTENDED_CLIENT_SLOT client;
  client.uid = client_uid;
  client.pid = client_pid;
  client.socket = fclient;
  client.client_id = -1;
  client.command = command;
  client.data = NULL;
  client.user_function = NULL; 
  (*client_address)(&client);
  
}

void start_client_thread_2(void (*client_address)(EXTENDED_CLIENT_SLOT *), int client_uid, int client_pid, int fclient, char *command, void *data, void (*user_server)())   /*   %ENTRY%   */

{
  int slot;
  pthread_t client_thread[1];
  pthread_attr_t client_attr;

  initialize_locks();
  initialize_client_table();
  
  pthread_mutex_lock(&mutex);
  client_no++;
  client_ord++;
  pthread_mutex_unlock(&mutex);
  
  
  slot = find_client_slot();
  clients[slot].uid = client_uid;
  clients[slot].pid = client_pid;
  clients[slot].socket = fclient;
  
  clients[slot].client_id = client_ord;
  clients[slot].command = (char *)malloc(strlen(command) + 1);
  strncpy(clients[slot].command, command, strlen(command) + 1);
  
  clients[slot].data = data;
  
  clients[slot].user_function = user_server;
  
  /* get_client(); */

#ifdef HP
   
  pthread_create(&client_thread[0],
		 client_attr = pthread_attr_default,
		 (void *) client_address,
		 (void *) &clients[slot]);
  
  pthread_detach(client_thread);

#else
  pthread_attr_init(&client_attr);
  pthread_create(&client_thread[0],
		 &client_attr,
		 (void *(*)(void *))client_address,
		 (void *) &clients[slot]);
  
  pthread_detach(client_thread[0]);
#endif

}

/* =================================================================== */
/* manage number of active clients                       */

/* increment (under lock) the count of connected clients */
void increment_client_count()  /*   %ENTRY%   */
{
     pthread_mutex_lock(&mutex);
     client_no++ ;
     pthread_mutex_unlock(&mutex);
}

/* decrement (under lock) the count of connected clients */
void decrement_client_count()  /*   %ENTRY%   */
{
     pthread_mutex_lock(&mutex);
     client_no-- ;
     client_ord--;
     pthread_mutex_unlock(&mutex);
}

/* get the count of connected clients */
int get_client_count()  /*   %ENTRY%   */
{
  return(client_no > MAX_EXTENDED_CLIENTS ? -1 : client_no);
}

/* EXTENDED_CLIENT_SLOT get_client(int fclient) */
int exit_from_client(int fclient)
{
  int i;

  for (i = 0; i < get_client_count(); i++)
    {
      if(clients[i].socket == fclient)
	{
	  exit_from_client_thread((void *)&clients[i]);
	  return 0;
	}
    }
  
  
  return -1;
}
