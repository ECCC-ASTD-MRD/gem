/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* MGILIB.C
   
   Functions written for programs that run in "parallel" and need 
   to exchange information simultaneously. MGI stands for Model Gossip
   Interface. Each program using these functions will have to compile 
   with this library. There are 6 functions:

   MGI_INIT
   to initialize a channel using the name given and open the socket.
   It will also allocate a writing buffer dynamically and if all is
   successful, it will return a channel number (socket descriptor).

   MGI_OPEN
   to open the channel in a certain mode:
   'R' for read: returns 0 to signal that data been written.
                 or returns nblks to be read
   'W' for write: returns 1 if open is ok.
   'S' for storing a restart file:returns 1 if open is ok.

   MGI_READ
   to read from a channel that is open for READ mode.
   It accepts the following type of data:
   'C' for CHARACTER
   'I' for INTEGER
   'R' for REAL
   'D' for REAL*8
   It returns the number of blocks left to read from channel.

   MGI_WRITE
   to write to a channel that is open for WRITE mode.
   It accepts the same type of data as MGI_READ.
   It returns the number of blocks written to channel.

   MGI_CLOS
   to close the mode of a channel and check to make sure all is
   transmitted as requested. It returns the status of the data
   file after it is closed.

   MGI_TERM
   to delete the PID file that was created in the beginning and 
   to release all the memory allocated dynamically. It closes all
   the filepipes therefore, breaking all the pipe connections with
   the other programs.

   ***NOTE: These functions are written to keep enough bits for the
   equivalent of an integer/float in C or integer/real*4 in FORTRAN.
   In other words, you will lose some precision with the 64-bit
   compilation unless real*8 is used.

*/

#include <stdio.h>
#include <stdlib.h>
#include <rpnmacros.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h> 
#include <string.h>
#include <ctype.h>
#include "mgi.h"
#include <gossip.h>
//JMB
#include <sys/resource.h>

/* error codes header file */
#include "cgossip.h"

#define  CLOSE      -5
#define  TIMEOUT    -5

static channel chn[MAX_CHANNELS];

static int ichan = 0;
static int init = 0;
static int SIG_ACTIVE = 1;
static char PID_file[MAX_STR];
static char *mgidir;
static int *intBuffer; 

static void getmgidir ();
static int makepidfile ();
static void removepidfile ();
static void strcopy (char *s, char *t, int charlen);
static int validchan (int chan);
static int bwrite (int chan, void *buffer, int nelem, char *dtype);
ftnword f77name (mgi_init) (char *channel_name, F2Cl lname);
ftnword f77name (mgi_open) (ftnword *f_chan, char *mode, F2Cl lmode);
ftnword f77name (mgi_read) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
ftnword f77name (mgi_write) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
ftnword f77name (mgi_clos) (ftnword *f_chan);
ftnword f77name (mgi_term) ();
void f77name (mgi_set_timeout) (ftnword *chan, ftnword *timeout);

extern int connect_to_subchannel_by_name (char *channel, char *subchannel, char *mode);
extern int get_ack_nack (int socket);
extern int write_record (int fclient, void *buf, int longueur, int tokensize);
extern void *read_record (int fclient, void *buf, int *longueur, int maxlongueur, int tokensize);
extern char *get_gossip_dir (int display);

extern void init_client_table ();
extern void set_client_timeout (int fclient, int timeout);
extern int get_client_timeout (int fclient);
extern int close_channel (int fclient, char *channel);

ftnword f77name (mgi_read_oob) ();
ftnword f77name (mgi_write_oob) ();

/* --------------------------------------------------------------------------- */
/* #define DEBUG */

/*********************************************************************************************/

void f77name (mgi_nosig) ()
     /* to disable the signals between filepipes */
{
  /* SIG_ACTIVE = 0; */
  fprintf(stderr,"MGI_NOSIG: deprecated call\n");
}

/* to copy a string given by a fortran routine, 
   the space character is ignored, it is taken 
   here as the end of the string */
static void strcopy_( char *s1, char *s2, int lengths1, int lengths2 )    
{
  int i = 0;

  while ( (*s1++ = *s2++) != ' ' && i++ < lengths2 );
  
}


/* to compare two strings given by a fortran routine, considering 
   the space character as the end of the string*/
static int f_strcmp( unsigned char *s1, unsigned char *s2, int s1length, int s2length )
{
  int i = 0;	 
  int length;
  
  if(s1length <= s2length)
    length = s1length;
  else
    length = s2length;

  while( i < length && *s1 != ' ' && *s2 != ' ' && *s1 == *s2 && *s1 != '\0' && *s2 != '\0')
    {
      s1++;
      s2++;
      i++;
    }

  
  if(*s1 == ' ')
    {
      fprintf(stderr, "mgilib2::f_strcmp(), before return if(*s1 == ' '), s1 => %s\n ", s1);
      return check_ends(s1, s2, s1length, s2length, i);
     
    }
  else if(*s2 == ' ')
    {
      fprintf(stderr, "mgilib2::f_strcmp(), before return if(*s2 == ' '), s2 => %s\n ", s2);
      return 2 + check_ends(s2, s1, s2length, s1length, i);
    }

  return (*s1 - *s2);

}

int check_ends( char *s1, char *s2, int s1length, int s2length, int i )
{
  if(*s2 == ' ' )
    {
      return (*(s1 - 1) - *(s2 - 1));
    }
  else
    {
      if(i == s2length)
	    return (*(s1 - 1) - *(s2 - 1));
      else
	return -1;
    }
}
/***********************************************************************************************/


static void strcopy( char *s, char *t, int charlen )
     /* to copy a string given by a fortran routine and place the NULL 
	character at the end of the true (charlen) length of the string */
{
  int i;
  
  i = 0;
  while ( (*s++ = *t++) != ' ' && i++ < charlen);
  if (*s-- == ' ') *s = '\0';
  else *s++ = '\0';
}

static int validchan( int chan )
     /* to validate the channel number; it must be greater than
	zero and less than or equal to ICHAN*/
{
  if ( chn[chan].buffer == NULL )
    return (-1);
  return(0);
}

static void getmgidir()
     /* to get the value of the environment variable "MGI_DIR" */
{
  if ( (mgidir = getenv("MGI_DIR")) == NULL)
    {
      fprintf(stderr,"Environment variable \"MGI_DIR\" undefined --\n");
      /* exit(1); */
    }
}

static int makepidfile()
     /* to make the PID file */
{
  char stuff[MAX_STR];
  sprintf( PID_file, "%s/%d", mgidir, getpid() );
  sprintf( stuff, "%s/%s", mgidir, "PROCS" );
  fprintf(stderr, "linking :%s: to :%s:\n", PID_file, stuff );
  return ( link ( stuff, PID_file ) );
}

static void removepidfile()
     /* to remove the PID file */
{
  fprintf(stderr, "removing %s\n", PID_file );
  unlink( PID_file );
}

static int bwrite ( int chan, void *buffer, int nelem, char *dtype )
     /* To fill the write buffer of initialized channel "chan" on the server */
{
  int ier;

#ifdef TRACE_LEVEL1
  fprintf(stderr, "\n bwrite: COMM1-W send_command_to_server() channel: %d",chn[chan].gchannel); 
  fflush(stderr);
#endif

   ier = send_command_to_server(chn[chan].gchannel, "WRITE");

   if (ier < 0)
    { 
#ifdef DEBUG
        fprintf(stderr, "\n bwrite send_command_to_server returns ier= %d\n",ier);
	fprintf(stderr, "\n bwrite: COMM1-W send_command_to_server() error on channel: %d",chn[chan].gchannel); 
        fprintf(stderr,"\n bwrite, unable to send write command\n");
        fflush(stderr);
#endif
	return(-1);
    }


  if(*dtype == 'I' || *dtype == 'R')
    {
      ier = write_record(chn[chan].gchannel, (unsigned char *)buffer, nelem, sizeof(int));
    }

  else if(*dtype == 'D')
    {
      ier = write_record(chn[chan].gchannel, (char *)buffer, nelem, sizeof(double));
    }
  
  else if(*dtype == 'C')
    {
      ier = write_record(chn[chan].gchannel, buffer, nelem, 1);
    }

#ifdef TRACE_LEVEL1  
    fprintf(stderr, "\n bwrite: COMM8-R get_ack_nack() \n");
    fflush(stderr);
#endif

    if (get_ack_nack(chn[chan].gchannel)) 
      { 
#ifdef DEBUG
	fprintf(stderr, "\n bwrite: COMM8-R get_ack_nack() error on channel: %d",chn[chan].gchannel); 
        fflush(stderr);
#endif
	return(-1);
      }
    else 
      { 
	return(0);
      }

  return ier;
}

ftnword f77name (mgi_clos) (ftnword *f_chan)
     /* close a channel and signal that it can be opened in another mode */
{
  int ier = 0, chan;
  char buf[1024];
  chan = (int) *f_chan;

  if(chn[chan].gchannel != 0)
    {
      snprintf(buf, 1023, "%s %s", "END", chn[chan].name);
      ier = send_command(buf);
      fprintf(stderr,"MGI_CLOS: subchannel \"%s\" is closed \n", chn[chan].name);
    }
  
   if(chn[chan].buffer)
    {
      free(chn[chan].buffer);
      chn[chan].buffer = NULL;
    }  
  return ier;
  
}

ftnword f77name (mgi_term) ()
{
  /* close all channels */
  int chan, ier = -1;

  for (chan = 0; chan <= ichan; chan++)
    {
      if(chn[chan].name && strcmp((char *)chn[chan].name, "") && chn[chan].gchannel > 0)
	{
	  ier = send_command("END");
	  fprintf(stderr,"MGI_TERM: subchannel \"%s\" has been closed!\n", chn[chan].name);
	  
	  if(chn[chan].buffer)
	    {
	      free(chn[chan].buffer);
	      chn[chan].buffer = NULL;
	    }
	}
    }

  return ier;
}

ftnword f77name (mgi_init) (char *channel_name, F2Cl lname)
     /* To initialize a channel given a channel_name.
	It will return a number to represent this channel (1 to MAX_CHANNELS-1 */
{
  int chan;

  if (init == 0)
    {
      init = 1;
    }

#ifdef DEBUG
  fprintf(stderr,"MGI_INIT ** \n"); 
#endif
  ichan++;
  if (ichan >= MAX_CHANNELS)
    {
      fprintf(stderr,"MGI_INIT: ERROR, Too many channels assigned; MAX = %d\n", MAX_CHANNELS);
      /* return -1; */
      return INIT_ERROR;
    }
  else
    {
      chan = ichan;
      if (lname < MAX_NAME)
      {
	strcopy(chn[chan].name, channel_name, lname);
      }
      else 
      {
	fprintf(stderr,"MGI_INIT: ERROR, Length of channel name > %d chars.\n",
	       MAX_NAME-1);
	/* return -1; */
	return INIT_ERROR;
      }
      chn[chan].fd_data = -1;
      if (SIG_ACTIVE)
	{
	  fprintf(stderr,"MGI_INIT: Opening channel: \"%s\" \n", chn[chan].name);
	}
    
      /* initialize channel */
      chn[chan].msgno_W = 0;
      chn[chan].msgno_R = 0;
      chn[chan].nblks = 0;
      chn[chan].mode = ' ';
      chn[chan].pos = 0;
      chn[chan].gchannel = 0;




    if ((intBuffer = (int *) malloc(BUFSIZE * sizeof(int))) == NULL)
      {
	fprintf(stderr,"MGI_INIT: ERROR on channel %s: Cannot allocate memory for intBuffer\n",
	       chn[chan].name);
	/* return -1; */
	return INIT_ERROR;
      }

    chn[chan].buffer = intBuffer;
    }

  return(chan);
}

ftnword f77name (mgi_open) (ftnword *f_chan, char *mode, F2Cl lmode)
     /* to open a channel in mode "mode"; where mode can be:
	'R' for reading
	'W' for writing
	'S' for storing
     */
{
  int chan;
  chan = (int) *f_chan;

  if (*mode == 'W') 
    {
      chn[chan].gchannel = connect_to_subchannel_by_name( get_gossip_dir(0), chn[chan].name, "write" );

      if( chn[chan].gchannel < 0 )
	chn[chan].gchannel = retry_connect( chan );

    }
  else if (*mode == 'R') 
    {
      chn[chan].gchannel = connect_to_subchannel_by_name( get_gossip_dir(0), chn[chan].name, "read" );

      if( chn[chan].gchannel < 0 )
	chn[chan].gchannel = retry_connect( chan );
    }
  else if (*mode == 'S') 
    { /* store mode (for restart files)*/
      chn[chan].mode = 'S';
      chn[chan].nblks = 0;
      chn[chan].msgno_W++;
      chn[chan].pos = 0;
    }

  if(chn[chan].gchannel < 0)
    {
      fprintf(stderr, "MGI_OPEN, Connection Failed, the Server may be down !!\n" );
      /* exit(-1); */
      return CONNECTION_ERROR; 
    }

  /* initialize timeout table */
  init_client_table( chn[chan].gchannel );

  return chan;
}


/* if connection to server fails          */
/* default:  retry 10 times after a sleep */
/* else use user value                    */

int USER_TRY_CONNECT = 10;
void f77name (mgi_set_retry_connect) (ftnword *try_nbr)
{
  printf( "MGI_OPEN, setting try to connect USER_TRY_CONNECT: \"%d\" times\n", (int) *try_nbr );  
  if((int) *try_nbr > 0 && (int) *try_nbr < 10)
    USER_TRY_CONNECT = (int) *try_nbr;
}

int mgi_get_retry_connect(int chan)
{
  
  return USER_TRY_CONNECT;

}

/* if connection to server fails          */
/* default: retry 10 times after a sleep  */
/* interval of 10 secs                    */
int retry_connect( int chan )
{
  int PING_INTERVAL = 10;
  int ping_ord0 = mgi_get_retry_connect(chan);
  int ping_ord =  mgi_get_retry_connect(chan);

  while( chn[chan].gchannel < 0 && ping_ord > 0 )
    {
      sleep( PING_INTERVAL );
      fprintf(stderr, "MGI_OPEN, Connection to Server Failed,  retry to connect: \"%d/%d\" \n", ping_ord0 - ping_ord + 1, ping_ord0 );
	  chn[chan].gchannel = connect_to_subchannel_by_name( get_gossip_dir(0), chn[chan].name, "write" );
	  ping_ord--;
    }
  return chn[chan].gchannel; 
  
}
ftnword f77name (mgi_write) (ftnword *f_chan, void *buffer, ftnword *f_nelem, char *dtype, F2Cl ltype)
     /* to write elements from "buffer" into the specified channel
	opened for WRITEMODE. It actually writes
	
	The following data types (dtype) are accepted:
	'C': character
	'I': integer
	'R': real
	'D': real*8 ; note that only the precision of a real would be kept
     */
{
  int ier, chan, nelem;
  int lnblnk_();

  chan = (int) *f_chan;
  nelem = (int) *f_nelem;
  char *tmpstr;

#ifdef DEBUG
      fprintf(stderr,"\nMGI_WRITE JMB: data type = %c, elts Nbr = %d, subchannel = %s\n", dtype[0], nelem, chn[chan].name);
#endif
  
  if( nelem <= 0 )
    {
      fprintf(stderr,"\nMGI_WRITE, Error, cannot write data with length = %d\n", nelem);
      
      return WRITE_ERROR;
    }

  if( chn[chan].gchannel < 0 )
    {
      fprintf(stderr,"\nMGI_WRITE, Error, cannot connect to server using descriptor: \"%d\"!!!\n", chn[chan].gchannel);
      
      return WRITE_ERROR;
    }

  if ( *dtype == 'C' )
    {
      nelem = ( *f_nelem < ltype ) ? (int) *f_nelem:ltype;

      tmpstr = (char *)malloc(nelem + 1);

      strncpy( tmpstr, (char *)buffer, nelem);
      tmpstr[nelem] = '\0';

#ifdef DEBUG
      fprintf(stderr,"\nMGI_WRITE CHARACTER JMB: data type = %c, elts Nbr = %d, strlen = %d,  subchannel = %s\n", dtype[0], nelem, ltype, chn[chan].name);
#endif

      if ((ier = bwrite(chan, (unsigned char *)tmpstr, nelem, dtype)) < 0)
	{
	  fprintf(stderr,"\nMGI_WRITE (C): ERROR on %s\n", chn[chan].name);
	  free( tmpstr );
	  /* return number of bytes not sent */
	  return WRITE_ERROR;
	}
      free( tmpstr );

    }


  else if (*dtype == 'I' || *dtype == 'R' || *dtype == 'D' ) 
    {
      chn[chan].nblks++;

      if ((ier = bwrite(chan, (unsigned char *)buffer, nelem, dtype)) < 0)
	{
	  fprintf(stderr,"\nMGI_WRITE(I || R || D) : ERROR on %s\n", chn[chan].name);
	  /* return number of bytes not sent */
	  return WRITE_ERROR;
	}
     
    }

  else 
    {
      fprintf(stderr,"\nMGI_WRITE: ERROR on channel %s: Unknown data type: %c\n", chn[chan].name, *dtype);
      /* return -1; */
      return WRITE_TYPE_ERROR;
    }


  if(ier < 0)
    {
      if(get_timeout_signal(chn[chan].gchannel))
	{
	  if (*dtype == 'C')
	    fprintf(stderr, "\nMGI_WRITE: TIMEOUT for write \"%d of Character data\" \n", nelem);

	  else if(*dtype == 'I')
	    fprintf(stderr, "\nMGI_WRITE: TIMEOUT for write \"%d of Integer data\" \n", nelem);

	  else if(*dtype == 'R')
	    fprintf(stderr, "\nMGI_WRITE: TIMEOUT for write \"%d of Real data\" \n", nelem);

	  else if(*dtype == 'D')
	    fprintf(stderr, "\nMGI_WRITE: TIMEOUT for write \"%d of Double data\" \n", nelem);

	  return signal_timeout(chn[chan].gchannel);
	}
    }

  return ier;
}

void f77name (mgi_set_timeout) (ftnword *chan, ftnword *timeout)
{
  set_client_timeout(chn[(int) *chan].gchannel, (int) *timeout);

}

ftnword f77name (mgi_read) (ftnword *f_chan, void *buffer, ftnword *f_nelem, char *dtype, F2Cl ltype)

     /* to read elements directly from the data file related to the 
	specified channel into "buffer". The channel must be opened for 
	READMODE only.
	The following data types (dtype) are accepted:
	'C': character
	'I': integer (int)
	'R': real    (float)
	'D': real*8  (double)
     */
{
  int ier, chan, nelem;
  
  chan = (int) *f_chan;
  nelem = (int) *f_nelem;

#ifdef DEBUG
  fprintf(stderr,"MGI_READ JMB: data type = %c, elts Nbr = %d, strlen = %d,  subchannel = %s\n", dtype[0], nelem, ltype, chn[chan].name);
#endif

  if(nelem <= 0)
    {
#ifdef DEBUG
      fprintf(stderr,"MGI_READ, Error: cannot read data with length = %d\n", nelem);
#endif
      /* return -1; */
      return DATA_LENGTH_ERROR;
    }

  bzero(buffer, nelem);

  ier = send_command_to_server(chn[chan].gchannel, "READ");

  if(ier < 0)
    {
#ifdef DEBUG   
      fprintf(stderr,"MGI_READ, Error: unable to send write command for channel: \"%s\"\n", chn[chan].name);
#endif    
      return SEND_COMMAND_ERROR;
    }
  
  if (*dtype == 'I')
    { /* integer */
  
      buffer = (int *)read_record( chn[chan].gchannel, (int *)buffer, &nelem, nelem, sizeof(int) );

      if(buffer != NULL)
	{
	  get_ack_nack( chn[chan].gchannel );
	
	  return ier = nelem;
	}
      else
	{
	  if( get_timeout_signal(chn[chan].gchannel) )
	    {
	      fprintf(stderr, "MGI_READ: TIMEOUT for read \"Integer\" \n" );
	   
	      ier = READ_TIMEOUT;
	    }
	  else
	    {
	      fprintf( stderr, "MGI_READ: Problem read Integer\n" );
	   
	      return READ_ERROR;
	    }
	}

    }

  else if (*dtype == 'R')
    { /* float */

      //JMB      
      struct rusage mydata;

      //      getrusage(RUSAGE_SELF,&mydata);
      //      printf(" MGI_READ: before read_record 'R', maxrss=%d\n",mydata.ru_maxrss);

      buffer = (float *)read_record(chn[chan].gchannel, (float *)buffer, &nelem, nelem, sizeof(int));

      //      getrusage(RUSAGE_SELF,&mydata);
      //printf(" MGI_READ: after read_record 'R', maxrss=%d\n",mydata.ru_maxrss);

      
      if(buffer != NULL)
	{
	  get_ack_nack(chn[chan].gchannel);
	
	  return ier = nelem;

	}
      else
	{
	  
	  if( get_timeout_signal( chn[chan].gchannel ) )
	    {
	      fprintf(stderr, "MGI_READ:  TIMEOUT for read \"Real\" \n");

	      ier = READ_TIMEOUT;
	    }
	  else
	    {
	      fprintf( stderr, "MGI_READ: problem read Real data\n" );

	      return READ_ERROR;
	    }
	}
      
    }
  else if (*dtype == 'D')
    { /* double */


      buffer = (double *)read_record(chn[chan].gchannel, (double *)buffer, &nelem, nelem, sizeof(double));

      if(buffer != NULL)
	{
	  get_ack_nack(chn[chan].gchannel);
	
	  return ier = nelem;
	}
      else
	{
	  
	  if( get_timeout_signal( chn[chan].gchannel ) )
	    {
	      fprintf(stderr, "MGI_READ: TIMEOUT for read \"Double\"\n");

	      ier = READ_TIMEOUT;
	    }
	  else
	    {
	      fprintf( stderr, "MGI_READ: Problem read Double data\n" );

	      return READ_ERROR;
	    }
 	}
     }
  
  else if (*dtype == 'C')
    { /* character */
      int i;
      char *temp = (char *)buffer;

      for(i = 0; i < ltype ; i++ ) 
	{
	  temp[i] = ' ';
	}

      buffer = (char *)read_record(chn[chan].gchannel, (char *)buffer, &nelem, nelem, sizeof(char));

      for(i = nelem+1 ; i < ltype ; i++ ) 
	{
	  temp[i] = ' ';
	}
      
      if(buffer != NULL)
	{
	  get_ack_nack(chn[chan].gchannel);
	  
	  return ier = nelem;
	} 
      else
	{
	  
	  if( get_timeout_signal( chn[chan].gchannel ) )
	    {
	      fprintf(stderr, "MGI_READ: TIMEOUT for read \"Character\"\n");
	      
	      ier = READ_TIMEOUT;
	    }
	  else
	    {
	      fprintf( stderr, "MGI_READ: Problem read Character data\n" );
	      
	      return READ_ERROR;
	    }
	}

    }
  
  else
    {
      fprintf(stderr,"MGI_READ: ERROR on channel %s: Unknown data type: %c\n", chn[chan].name, *dtype);
      
      return READ_TYPE_ERROR;
    }
  
  if(ier == CLOSE)
    {
      close_channel(chn[chan].gchannel, chn[chan].name);
    }


  return ier;
}


