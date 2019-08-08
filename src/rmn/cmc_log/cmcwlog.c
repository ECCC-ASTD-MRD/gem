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

#include <rpnmacros.h>
#include "envrecv.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <string.h>

static int setfp();

/*
  c_cmcwlog3 is the lowest interface level
  filen is the log file name
  cl is the message class, 2 chars max
  msgno is the message number, 4 digits max
  id is the message id, MAXID-1 chars max
  txt is the message text, MAXTEXT-1 chars max
*/
int c_cmcwlog3(char *filen,char *cl,int msgno,char *id,char *txt)
{
FILE *fp;
int fd;
int i;
int cllen,idlen, txtlen;
char msgcl[3];
char msgid[MAXID];
char msgtxt[MAXTEXT];
char tstamp[25];

struct tm *ttm;
time_t cloc,time();
clock_t cpuclock;
clock_t clock();
mode_t mask=umask(0);

if(msgno>9999) msgno=9999;
if(msgno<0) msgno=0;


/*
  fprintf(stderr,"file=%s, id=%s, msgno=%d, cl=%s, txt=%s, umask=%d\n",filen,id,msgno,cl,txt,mask);
*/
/*
  Make sure that there are valid id and text arguments
*/

if ( id == NULL || txt == NULL ) {
	return(1);
}

cloc=time(0L);                       /* get cpu consumption */
ttm = localtime(&cloc);
cpuclock = clock();
#if defined(HP)
cpuclock = cpuclock / CLOCKS_PER_SEC;
#else
cpuclock = cpuclock / 1000;
#endif
if (cpuclock < 1) cpuclock = 1;

/*
  Open the log file and if it does not exist 
  then create it. Follow this by converting
  it to a stream
*/
if((fd = open(filen,O_CREAT | O_RDWR, S_IREAD | S_IWRITE | S_IRGRP | S_IROTH )) == -1){
  perror(filen);
  return(2);
}
fp=fdopen(fd,"r+");

/* 
 Lock the log file since more than one process may be writing to it
*/
#ifndef __CYGWIN__
lockf(fileno(fp),F_LOCK,0L); 
#endif

/*
  set up the write pointer and check that the log file is not full
  setfp will also update the write pointer
*/
if( setfp(fp) == -1){
  fprintf(stderr,"log file %s is full\n",filen);
#ifndef __CYGWIN__
  lockf(fileno(fp),F_ULOCK,0L);
#endif
  fclose(fp);
  return(3);
}

cllen = strlen (strncpy( msgcl,cl,2 ));
for ( i=cllen; i < 2; i++)
  msgid[i] = ' ';
msgcl[2] = '\0';

for ( i=0; i < 24; i++)
  tstamp[i] = '%';
tstamp[24] = '\0';

/*
  Determine the length of the message id which may be a maximum of MAXID 
  message id is blank filled except for element MAXID-1 which is \0
*/
idlen = strlen (strncpy( msgid,id,MAXID-1 ));
for ( i=idlen; i < MAXID -1; i++)
	msgid[i] = ' ';
msgid[MAXID -1] = '\0';

/*
  Determine the length of the message  which
  may be a maximum of MAXTEXT and will be blank filled
  except for element MAXTEXT-1 which is \0
*/
txtlen = strlen(strncpy( msgtxt,txt,MAXTEXT-1));
for ( i=txtlen; i < MAXTEXT -1; i++)
	msgtxt[i] = ' ';
msgtxt[MAXTEXT - 1] = '\0';

/*
  write the message id and text to the log
*/
sprintf(tstamp,"%4i%02i%02i-%02i%02i%02i-%08lu",ttm->tm_year+1900,
        ttm->tm_mon+1,ttm->tm_mday,ttm->tm_hour,
        ttm->tm_min,ttm->tm_sec,cpuclock);
fprintf(fp,"%s%04i%s%s%s\n",msgcl,msgno,msgid,tstamp,msgtxt);
/*
  clean up
*/
#ifndef __CYGWIN__
lockf(fileno(fp),F_ULOCK,0L);
#endif
fclose(fp);
return(0);  /* SUCCESS, return 0 */
}
/* ===========================================================================*/
int c_cmcwlog2B(char *filen,char *cl,int msgno,char *id,char *txt)
{
char *pivot=strchr(filen,'@');
char buffer[1024];
int fserver, errcode;

if(pivot!=NULL){  /* filename@@host:port or filename@channel_name, this is a network log call */
  *pivot='\0'; pivot++;

  /*   connect to server, specified by channel_name or @host:port */
  snprintf(buffer,(sizeof buffer) -1,"%s",pivot);
  fserver=connect_to_channel_by_name(pivot); 
  if(fserver < 0) {
    fprintf(stderr,"Error while connecting to %s\n",buffer);
    return(4);
    }
  snprintf(buffer,(sizeof buffer) -1,"%s;%s;%d;%s;%s;",filen,cl,msgno,id,txt);
  errcode=send_command_to_server(fserver,buffer);
/*
  fprintf(stderr,"Closing server connection and exiting\n");
*/
  close(fserver);
  return(errcode ? 6 : 0);
  }
else { /* simple case straight file name  contained in environment variable CMC_LOGFILE */
  return(c_cmcwlog3(filen,cl,msgno,id,txt));
  }
}

/* ===========================================================================*/
/*
    2nd level interface, c_cmcwlog2 will process the multiple targets
    supplied by filen containing a string of the following  type
    "target_number_1,target_number_2,target_number_3" (comma separated list )

    upon success, c_cmcwlog2 returns 0 to the caller, upon failure
    c_cmcwlog2 will try targets in succession. If all targets fail,
    an error code is returned
*/
int c_cmcwlog2(char *filen,char *cl,int msgno,char *id,char *txt)
{
char *next_target=strchr(filen,',');
int retcode;
char buffer[1024];

if(next_target!=NULL) *next_target='\0';

while(*filen){                 /* LOOP OVER POTENTIAL TARGETS UNTIL SUCCESFUL */
  snprintf(buffer,(sizeof buffer) -1,"%s",filen);
  retcode=c_cmcwlog2B(filen,cl,msgno,id,txt);

  if(retcode==0) return(0);   /* SUCCESS, return immediately with the good news */

  if(next_target!=NULL) {     /* OOPS, anything more to try */
    fprintf(stderr,"LOG TARGET %s failed, trying next one \n",buffer);
    }
  else {                      /* nope */
    fprintf(stderr,"LOG TARGET %s failed \n",filen);
    return(retcode);
    }

  filen=next_target+1;
  next_target=strchr(filen,',');
  if(next_target!=NULL) *next_target='\0';
  }
return(7); /* everything failed sigh !! */
}
/* ===========================================================================*/
/*
    c_cmcwlog(CLASS,msgno,ID,TEXT);

    1st level interface, c_cmcwlog will check for environment
    variable CMC_LOGFILE and complain if not found
    if found, the 2nd level interface will be called
*/
int c_cmcwlog(char *cl,int msgno,char *id,char *txt)
{
  char *filen;
  extern char *getenv();

/*
  Get environmental variable CMC_LOGFILE which is the pathname of the target log
  file and make sure it is declared.
*/
  if ( (filen=getenv("CMC_LOGFILE_PLUS")) == NULL) {
   if ( (filen=getenv("CMC_LOGFILE")) == NULL) {
    return(5);
   }
  }
  return(  c_cmcwlog2(filen,cl,msgno,id,txt)  ) ;
}

/*
  Find the position of the next write in the logfile either by initializing the
  write pointer in a new file or reading from an existing one
*/
/* ===========================================================================*/
static int setfp(fp)
FILE *fp;
{
long int frptr;
long int fwptr,fwptr2;
long int flptr;

/*
	  if the file is empty initialize the read
	  write and limit pointers,writing them to the
	  file.
*/
if(fscanf(fp,"%ld\n",&frptr) == EOF) {
  fwptr=DATA_POS;
  frptr=DATA_POS;
  flptr=LIMT_PTR;
  fseek(fp,READ_PTR,SEEK_SET);
  fprintf(fp,"%-19ld\n",frptr);
  fseek(fp,WRIT_PTR,SEEK_SET);
  fprintf(fp,"%-19ld\n",fwptr);
  fseek(fp,LIMT_PTR,SEEK_SET);
  fprintf(fp,"%-19ld\n",ENTRYSIZE);
}

/* 
	 file exists so read the pointers
*/
fseek(fp,READ_PTR,SEEK_SET);
fscanf(fp,"%ld\n",&frptr);

fseek(fp,WRIT_PTR,SEEK_SET);
fscanf(fp,"%ld\n",&fwptr);

fseek(fp,LIMT_PTR,SEEK_SET);
fscanf(fp,"%ld\n",&flptr);

/*
  update the value of the write pointer
*/
fwptr2 = fwptr + ENTRYSIZE;
fseek(fp,WRIT_PTR,SEEK_SET);
fprintf(fp,"%-19ld\n",fwptr2);
/*
  leave file positioned at writing position
*/
fseek(fp,fwptr,SEEK_SET);
return(0);
}
/* ===========================================================================*/
/*
   FORTRAN interface
   call cmcwlog(CLASS,msgno,ID,TEXT)
   character *(*) CLASS,ID,TEXT
   integer msgno
   ALL ARGUMENTS ARE INPUT ONLY
*/
wordint f77name(cmcwlog)(fclstr,msgno,fidstr,ftxtstr,lclstr,lidstr,ltxtstr)
char *fidstr, *ftxtstr, *fclstr;
ftnword *msgno;
F2Cl lidstr,ltxtstr,lclstr;
{
char msgcl[30];
char msgid[MAXID];
char msgtxt[MAXTEXT];

if (lclstr >= 2) {
   strncpy( msgcl,fclstr,2 );
   msgcl[2] = '\0';
   }
else {
   strncpy( msgcl,fclstr,lclstr );
   msgcl[lclstr] = '\0';
   }

if (lidstr > MAXID-1) {
   strncpy( msgid,fidstr,MAXID-1 );
   msgid[MAXID-1] = '\0';
   }
else {
   strncpy( msgid,fidstr,lidstr );
   msgid[lidstr] = '\0';
   }

if (ltxtstr > MAXTEXT-1) {
   strncpy( msgtxt,ftxtstr,MAXTEXT-1 );
   msgtxt[MAXTEXT-1] = '\0';
   }
else {
   strncpy( msgtxt,ftxtstr,ltxtstr );
   msgtxt[ltxtstr] = '\0';
   }

return(c_cmcwlog(fclstr,*msgno,msgid,msgtxt));
}
