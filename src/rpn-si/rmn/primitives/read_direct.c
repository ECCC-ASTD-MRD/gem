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

/*

author: Michel Valin    August 2001

revisions: 
      V.Lee  July 2002   -increased value of MAXARGS (256 to 51200)
                         -increased value of MAXARGLEN (8192 to 1638400)
      S. Chamberland - Dec 2011 - add option to suppress informative messages
*/

#include <stdio.h>
#include <stdlib.h>
#include <rpnmacros.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

/*   structure containing callback entries  */
#define MAX_CMD_LEN 17
typedef struct{
	char command_name[MAX_CMD_LEN];  /* command name */
	int (*command)();                /* pointer to callback function */
	void *private_data;              /* pointer to private data */
	void *private_data_2;            /* pointer to private data */
	int max_args;                    /* maximun number of arguments permitted */
	wordint *actual_args;            /* actual number of arguments */
	int is_ftn;                      /* 1 if FORTRAN routine, 0 if C routine */
	}callback_entry;

/*   buffer used to process input file 1024 + 128 */
#define MAXBUFFER 1152
static unsigned char buffer[MAXBUFFER];
static unsigned char *buffer_start=buffer;
static unsigned char *buffer_end=buffer+MAXBUFFER-1;
static unsigned char *buffer_in=buffer+1024;
static unsigned char *buffer_out=buffer+1024;

/*   table containing character syntactic flags */
static unsigned char char_type[256];

/*   misc local flags and variables */
static int INIT=1;
static FILE *streamd=NULL;     /* input file descriptor, STDIN by default */
static int cur_char=0;         /* current input character */
static int cur_char_typ=0;     /* type of current input character */

/*  table of callback entries, number of entries in table */
#define MAXKEYS 1024
static callback_entry callback_table[MAXKEYS];
static int callbacks=0;

/* current token (verb or argument) being collected */
#define MAXTOKEN 80
static char token[MAXTOKEN];

/* argument collection table and related poointers */
#define MAXARGS 51200
static char *ARGV[MAXARGS];
static int ARGC=0;
#define MAXARGLEN 1638400
static char ARGTXT[MAXARGLEN];
static int ARGLEN;

static char abort_on_error=0;

#define RPNCB_VERBOSE 1
#define RPNCB_QUIET   0
static int rpnCBverbose=RPNCB_VERBOSE;

/* macro used to skip blanks, tabs, newlines, control characters from input file */
/* will not skip beyond a NEWLINE */
#define Skip_Blanks(some_dummy) { while( (cur_char_typ = char_type[cur_char = *buffer_out]) == 0 ) buffer_out++; }

/* macro used to set current_char to next input character */
#define Next_Char(some_dummy) {\
if(buffer_in<=buffer_out)fill_buffer();  \
buffer_out++; \
}

/*  register a C callback function */
int rpn_c_callback(char *VERB, void *callback, char *OPTIONS,
		void *Private_data, void *Private_data_2) {

 strncpy(callback_table[callbacks].command_name,VERB,MAX_CMD_LEN-1);
 callback_table[callbacks].command_name[MAX_CMD_LEN-1]='\0';
 callback_table[callbacks].command=callback;
 callback_table[callbacks].is_ftn=0;
 callback_table[callbacks].private_data=Private_data;
 callback_table[callbacks].private_data_2=Private_data_2;

 if(callbacks<MAXKEYS-1)callbacks++;

 return (callbacks);
}

/*  set verbose status */
void rpn_c_callback_setverbose(int verbose) {
  if(verbose>0) 
	 rpnCBverbose=RPNCB_VERBOSE;
  else 
	 rpnCBverbose=RPNCB_QUIET;
}

void f77name(rpn_f_callback_setverbose)(wordint *verbose) {
  rpn_c_callback_setverbose((int) *verbose);
}

/*  register a FORTRAN callback function */
wordint f77name(rpn_fortran_callback)(char *VERB, void *callback, char *OPTIONS,
		void *Private_data,void *Private_data_2,  F2Cl l_VERB, F2Cl l_OPTIONS) {

 strncpy(callback_table[callbacks].command_name,VERB,l_VERB<MAX_CMD_LEN-1?l_VERB:MAX_CMD_LEN-1);
 callback_table[callbacks].command_name[MAX_CMD_LEN-1]='\0';
 callback_table[callbacks].command=callback;
 callback_table[callbacks].is_ftn=1;
 callback_table[callbacks].private_data=Private_data;
 callback_table[callbacks].private_data_2=Private_data_2;

 if(callbacks<MAXKEYS-1)callbacks++;

 return (callbacks);
}

/* find command in callback table, return it's index if found, -1 otherwise */
static int Find_Callback(char *name){
 int i=callbacks;

 while(i--) {
   if(strcmp(callback_table[i].command_name,name)==0) return(i);
 }
 return(-1);
}

/* input buffer is empty, read as many characters as possible */
static void fill_buffer(){
 int nc;
 char *temp;
redo:
 buffer_in=buffer_out=buffer+1024;   /* make sure that it is possible to push back some chars */
 *buffer_in=0xFF;

 if(NULL != fgets(buffer_in,buffer_end-buffer_in-1,streamd))
   nc=strlen(buffer_in);
 else
   nc=0;

 temp=buffer_in;
 while(*temp==' ' || *temp=='\t')temp++;
 if(*temp=='\n') goto redo;

 if(nc>0) {
   buffer_in+=nc;
   if(rpnCBverbose==RPNCB_VERBOSE) fprintf(stderr,"%s",buffer_out);
 }
 cur_char = *buffer_out ;
 cur_char_typ = char_type[cur_char] ;
}

/* get current input character, do not bump input pointer, if buffer empty, fill it */
static unsigned int Current_Char(){

 if(buffer_in<=buffer_out) fill_buffer();
 if(buffer_in<=buffer_out) return(0xFF);    /* no input available, return EOF */

 cur_char = *buffer_out ;
 cur_char_typ = char_type[cur_char] ;

 return (cur_char);
}

/* initialize table containing syntactic character types */
static init_char_table(){
int i;

char_type[0] = 0xFF;
for (i=1;i<31;i++) char_type[i]=0;        /* control characters */
for (i=32;i<126;i++) char_type[i]=128;    /* regular ASCII chars */
for (i=127;i<254;i++) char_type[i]=0;     /* accented and special characters, ignore */

for (i='a';i<='z';i++) char_type[i] = 4+1;  /* TOKEN + LETTER */
for (i='A';i<='Z';i++) char_type[i] = 4+1;  /* TOKEN + LETTER */
for (i='0';i<='9';i++) char_type[i] = 4+2;  /* TOKEN + DIGIT */

char_type[255] = 0xFF;
char_type[' '] = 0;                         /* SPACE */
char_type['('] = 32;                        /* START OF CMD */
char_type['='] = 32;                        /* START OF CMD */
char_type['.'] = 4+16;                      /* TOKEN + NUM */
char_type['-'] = 4+16;                      /* TOKEN + NUM */
char_type['+'] = 4+16;                      /* TOKEN + NUM */
char_type['_'] = 4+1;                       /* TOKEN + LETTER */

char_type['%'] = 8; char_type['@'] = 8;     /* OPER */
char_type['"'] = 4+64;                      /* TOKEN + SEPAR */
char_type['<'] = 4+64;                      /* TOKEN + SEPAR */
char_type['['] = 4+64;                      /* TOKEN + SEPAR */
char_type['\''] = 4+64;                     /* TOKEN + SEPAR */
}

/* collect an argument for a command ,return token and length */
static int Argument(char *token){
 int len=0; char tmp_char;

 /* Skip blanks, even through a newline */
 Skip_Blanks(0); Current_Char(); Skip_Blanks(0);

 token[0]=Current_Char();                          /* get first character */
 if( (char_type[token[0]] & 4) == 0) {             /* check that charset is OK */
   fprintf(stderr,"Bad char in token:%c:\n",token[0]);
   token[1]='\0';
   return (-1);
 }
 len = 1;

 Next_Char(0);

 if(token[0]=='[' || token[0]==']') {
  token[1]='\0';
  return(len);
 }

 if(token[0]!='"' && token[0]!='\'' && token[0]!='<') {    /* not a delimited string */
  while ( (len < MAXTOKEN) && (char_type[token[len]=Current_Char()] & 4) ) {
   len++;
   Next_Char(0);
  }
 }else{                                   /* delimited string, collect to matching delimiter */
  if(token[0]=='<') token[0]='>';
  while ( (len < MAXTOKEN-1) && ((token[len]=Current_Char()) != token[0]) ) {
   token[MAXTOKEN-1]='\0';                /* make sure that string terminator is always present */
   if(token[len]=='\n') return (-1);      /* end of line, OOPS */
   len++;
   Next_Char(0);
  }
  token[len]=Current_Char();
  token[MAXTOKEN-1]='\0';                 /* make sure that string terminator is always present */
  if(token[len] != token[0]) return(-1);  /* did not find correct string delimiter */
  len++;
  Next_Char(0);
 }
 token[MAXTOKEN-1]='\0';                  /* make sure that string terminator is always present */

 if(len >=MAXTOKEN) return (-1);   /* token is too large */
 token[len]='\0';
 return(len);
}

/* collect the verb for a command ,return token and length */
static int Verb(char *token){
 int len=0;
 token[0]=Current_Char();
 Skip_Blanks(0);
 token[0]=Current_Char();  /* get first character of verb */

 if((token[0]&0xFF)==0xFF){         /* EOF ? */
  if(rpnCBverbose==RPNCB_VERBOSE) fprintf(stderr,"END of directives\n");
  fclose(streamd);
  return (-1);
  }
 if( (char_type[token[0]] & 1) == 0) {  /* bad first character ? (must be letter) */
  fprintf(stderr,"Bad starting character for Verb:%c:\n",token[0]);
  return (-1);
  }

 len = 1;
 Next_Char(0);  /* colect rest of verb, letters + digits only */
 while ( (len < MAXTOKEN) && (char_type[token[len]=Current_Char()] & (1+2)) ) {
   len++;
   Next_Char(0);
 }
 token[MAXTOKEN-1]='\0'; 

 if(len >=MAXTOKEN) return (-1);   /* token is too large */
 token[len]='\0';
 return(len);
}

/* test function for C callbacks, prints arguments and private data */
static int Print_Args(int argc , char **argv,char cmd_strt,  char *Private_Data, char *Private_Data_2) {
int cmd=0xFF&cmd_strt;
int i;
printf("Private_Data=%s\n",Private_Data);
printf("Private_Data_2=%s\n",Private_Data_2);
printf("Command Type=%s\n",cmd=='('?"EXEC":"ASSIGN");
for (i=0;i<=argc;i++) printf("%s\n",argv[i]);
return(0);
}

/* process directives from file */
int process_c_callback(char *filename){
int oo, errors=0;
int i, nrange, start_of_list=0;
int status;
wordint wstatus;
int cmd_end = ')'; unsigned char cmd_strt;

if(INIT) init_char_table();
INIT = 0;

/* open input file if a name was specified, use stdin otherwise */
streamd=stdin;
if(filename != NULL) streamd=fopen(filename,"r");
if (streamd == NULL) {
  fprintf(stderr,"process_c_callback: Cannot Open File -%s-\n",filename);
  goto error;
}

while(Verb(token)>0){ /* collect "verb" (command name) */
  ARGLEN=0;
  ARGC=0;
  ARGV[0]=token;
  Skip_Blanks(0); 
  if(cur_char_typ != 32 ) {
    fprintf(stderr,"Badly formed command \n");
    goto error;
   }
  cmd_strt = cur_char;
  cmd_end = ( cur_char == '(' ) ? ')' : ';' ;   /* assign proper terminator to command */
  Next_Char(0);

  /* collect comma separated command arguments */ 
  while(1){

    collect_arg:

    ARGC++;
    ARGV[ARGC]=&ARGTXT[ARGLEN];

    if((oo=Argument(ARGTXT+ARGLEN)) <= 0 ) {
      fprintf(stderr,"TOKEN error in argument\n"); 
      goto error;
    }

    if((ARGC>=MAXARGS-1) || (ARGLEN+oo>=MAXARGLEN-1)) {
      fprintf(stderr,"Too many arguments or arguments too long \n");
      fprintf(stderr,"ARGC=%d, ARGLEN=%d, oo=%d \n",ARGC,ARGLEN,oo);
      goto error;
    }
    if(*ARGV[ARGC]=='[') {        /* Argument is '[   ', will become [nnn */
     /* printf("Start of list at ARGC=%d\n",ARGC); */
     ARGTXT[ARGLEN+1] = '0';
     ARGTXT[ARGLEN+2] = '0';
     ARGTXT[ARGLEN+3] = '0';
     ARGTXT[ARGLEN+4] = ']';
     ARGTXT[ARGLEN+5] = '\0';
     ARGLEN=ARGLEN+6;
     if(start_of_list!=0) {
      fprintf(stderr,"List already open \n");
      goto error;
     }
     start_of_list=ARGC;
     goto collect_arg;
    }
    ARGTXT[ARGLEN+oo]='\0';
    ARGLEN=ARGLEN+oo+1;

    if(*ARGV[ARGC]=='>'){          /* range detected */
     double dble0,dble1,dble2;
     char term;
     char pad[1024];
     int nargs;
     nargs=sscanf(ARGV[ARGC]," > %lf , %lf , %lf%[%> ]%c",&dble0,&dble1,&dble2,pad);
     /* printf("Range detected, nargs=%d,pad=:%s:,from %f to %f by %f\n",nargs,pad,dble0,dble1,dble2);  */
     if(nargs!=4 || pad[strlen(pad)-1]!='>') { fprintf(stderr,"bad range\n"); goto error;}

     ARGLEN=ARGLEN-oo-1;ARGC--;
     nrange=(dble1-dble0)/dble2;
     if(nrange<0) nrange=0;
     while(nrange-->=0) {       /* expand range into argument list */
       if((ARGC>=MAXARGS-1) || (ARGLEN+oo>=MAXARGLEN-30)) {
         fprintf(stderr,"Too many arguments or arguments too long \n");
         fprintf(stderr,"nrange:ARGC=%d, ARGLEN=%d, oo=%d \n",ARGC,ARGLEN,oo);
         goto error;
       }
       ARGC++;
       ARGV[ARGC]=ARGTXT+ARGLEN;
       oo=sprintf(ARGV[ARGC],"%f",dble0);
       ARGTXT[ARGLEN+oo]='\0';
       ARGLEN=ARGLEN+oo+1;
       dble0=dble0+dble2;
     }
    }

    again:

    Skip_Blanks(0); Current_Char();

    if(cur_char == ']') {
     if(start_of_list==0) {
      fprintf(stderr,"Cannot close non existent list \n");
      goto error;
     }
     /* printf("End of  %d element list at ARGC=%d\n",ARGC-start_of_list,ARGC); */
     sprintf(ARGV[start_of_list],"[%3d]",ARGC-start_of_list);
     start_of_list=0;
     Next_Char(0);
     goto again;
    }
    if(cur_char == cmd_end) break; /* command terminator found */
    if(cur_char != ',') {fprintf(stderr,"bad separator :%c:\n",cur_char); goto error;}

    Next_Char(0);
  }
  Next_Char(0); Skip_Blanks(0);

  if(start_of_list!=0) {
   fprintf(stderr,"Unmatched [ \n");
   goto error;
  }

  if( (oo=Find_Callback(token)) >= 0 ) {

   if(callback_table[oo].is_ftn){
 
    int max_arg_lng=0; int i; int arg_lng; char *F_ARGV; int j; char *temp;
 
    for (i=0 ; i<=ARGC ; i++) {  /* find the length of the longest argument */
     arg_lng=strlen(ARGV[i]); 
     max_arg_lng=arg_lng>max_arg_lng?arg_lng:max_arg_lng;
    }
    F_ARGV=(char *)malloc(max_arg_lng*(1+ARGC));  /* allocate FORTRAN string space */
 
    for (i=0 ; i<max_arg_lng*(1+ARGC) ; i++) F_ARGV[i]=' ';
    temp=F_ARGV;
    for (j=0;j<=ARGC;j++) {              /* copy strings from C strings to FORTRAN strings */
     arg_lng=strlen(ARGV[j]);
     for (i=0 ; i<max_arg_lng ; i++) {
      if(i<arg_lng)*temp=ARGV[j][i];
      temp++;
     }
    }
    wstatus=callback_table[oo].command(&ARGC,F_ARGV,&cmd_strt,
                              callback_table[oo].private_data,
                              callback_table[oo].private_data_2,max_arg_lng,1);
    status=wstatus;
    free(F_ARGV);
   }
   else {
 
    status=callback_table[oo].command(ARGC,ARGV,cmd_strt,
                              callback_table[oo].private_data,
                              callback_table[oo].private_data_2);
   }
   if(status!=0) goto error;
  }else{
   fprintf(stderr,"Command %s NOT FOUND\n",token);
  }
continue;
error: fprintf(stderr,"skipping rest of line\n");
  errors++;
  buffer_out=buffer_in;
  start_of_list=0;
  if(abort_on_error) return(errors);
}  /* while verb */
return (errors);
}

/* equivalent of process_c_callback but FORTRAN callable */
wordint f77name(process_f_callback)(char *filename, F2Cl l_filename){
char the_name[1024];
strncpy(the_name, filename, l_filename<1023 ? l_filename : 1023);
the_name[l_filename<1023 ? l_filename : 1023]='\0';
return(process_c_callback(the_name));
}

#ifdef TEST
f77name(c_callback_test)()
{
int errs;
rpn_c_callback("verb1",Print_Args,"","VERB1","verb1");
rpn_c_callback("verb2",Print_Args,"","VERB2","verb2");
errs=process_c_callback("test_directives_c");
printf("Number of errors=%d\n",errs);
}
#endif
