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

/*CoMpIlAtIoN_OpTiOnS ::SX4=-Onooverlap::SX5=-Onooverlap::*/
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64

#include <rpnmacros.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>

#ifdef WIN32    /*CHC/NRC*/
#include <stdio.h>
#include <sys/stat.h>
#include <Winsock2.h>
#define L_SET SEEK_SET
#define L_INCR SEEK_CUR
#define L_XTND SEEK_END
#define S_IRUSR _S_IREAD
#define S_IWUSR _S_IWRITE
#define S_IRGRP _S_IREAD
#define S_IWGRP _S_IWRITE
#define S_IROTH _S_IREAD
#define S_IWOTH _S_IWRITE
#define WIN32_O_BINARY O_BINARY
#else
#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/signal.h>
#define WIN32_O_BINARY 0
#endif

#include <fcntl.h>
#include <errno.h>

/* ETG Normally in <sys/file.h> but forgotten on VPP */
#if defined (__uxpv__) || defined(SUN)
#define L_INCR          1       /* relative to current offset */
#define L_XTND          2       /* relative to end of file */
#endif

#define FNOM_OWNER
#include "../INTRALIB_INCLUDES/fnom.h"
#include "wafile.h"

#if defined (NEC)
#define tell(fd) lseek(fd,0,1)
#endif

#if defined(__linux__) || defined(__AIX__)
#define tell64(fd) lseek64(fd,0,1)
#endif

#ifdef __CYGWIN__
#define lseek64 lseek
#define open64 open
#define off64_t off_t
#define tell(fd) lseek(fd,0,1)
#define tell64(fd) lseek(fd,0,1)
#endif

void static dump_file_entry(int i);
static void reset_file_entry(int i);
static int find_file_entry(char *caller, int iun);
static int c_qqqfscr(char *type);
static int qqcclos(int indf);
static void scrap_page(int ind0,int ind1);
static void process_decay();
static void get_new_page(int ind);
static void wa_pages_flush(int ind);
static long long filepos(int indf);
static int qqcopen(int indf);
static void wa_page_read(int fd,word *buf,unsigned int adr,int nmots,int indf);
static void wa_page_write(int fd,word *buf,unsigned int adr,int nmots,int indf);
static void qqcwawr(word *buf,unsigned int ladr,int lnmots,int indf);
static void qqcward(word *buf,unsigned int ladr,int  lnmots,int indf);
static INT_32 qqcnblk(int lfd,int indf);
static void MOVE (word *src, word *dest, int nwords);
static void ZERO ( word *dest, int nwords);
static int fnom_rem_connect(int ind, char* remote_host);


static ENTETE_CMCARC cmcarc;
static ENTETE_CMCARC_V5 cmcarc64;

static FILEINFO wafile[MAXWAFILES];
static word *free_list[MAXWAFILES*MAXPAGES];
static int dastat[MAXWAFILES] = {MAXWAFILES * 0};

static int BLKSIZE = 512;

static int WA_PAGE_SIZE = 0;
static int WA_PAGE_NB   = 1;
static int WA_PAGE_LIMIT = 0;

static int global_count = 0;
static int nfree = -1;
static int init = 0;
static int debug_mode = 0;
static int subfile_length=0;
static int fnom_initialized=0;
static int endian_int=1;
static int stdoutflag=0;
static int stdinflag=0;
static char *little_endian=(char *)&endian_int;

static char *AFSISIO=NULL;
static char *ARMNLIB=NULL;
static char *LOCALDIR="./";

/****************************************************************************
*                   C _ F R E T O U R ,   F R E T O U R                     *
*****************************************************************************
*
***function c_fretour, fretour
*
*OBJECT: Kept only for backward compatibility.
*
*ARGUMENTS: in iun   unit number, ignored
*
*RETURNS: zero.
*
*/
int c_fretour(int iun){ return(0) ; }
ftnword f77name(fretour)(ftnword *fiun){ return(0) ; }


/****************************************************************************
*                         D U M P _ F I L E _ E N T R Y                     *
*****************************************************************************
*
***function dump_file_entry
*
*OBJECT: Prints the characteristics and attributes of one file 
*        of the master file table(for debugging use).
*
*ARGUMENTS: in  i   index of the entry to dump
*
*/
void static dump_file_entry(int i)
{
      fprintf(stderr,"FGFDT[%d] ",i);
      fprintf(stderr,"file_name=%s subname=%s file_type=%s\n",
              FGFDT[i].file_name,FGFDT[i].subname,FGFDT[i].file_type);
      fprintf(stderr,"iun=%d,fd=%d,size=%d,esize=%d,lrec=%d,flags=%s%s%s%s%s%s%s%s%s%s%s%s\n",
              FGFDT[i].iun,
              FGFDT[i].fd,
              FGFDT[i].file_size,
              FGFDT[i].eff_file_size,
              FGFDT[i].lrec,
              FGFDT[i].open_flag?"OPEN":"",
              FGFDT[i].attr.stream?"+STREAM":"",
              FGFDT[i].attr.std?"+STD":"",
              FGFDT[i].attr.burp?"+BURP":"",
              FGFDT[i].attr.rnd?"+RND":"+SEQ",
              FGFDT[i].attr.wa?"+WA":"",
              FGFDT[i].attr.ftn?"+FTN":"",
              FGFDT[i].attr.unf?"+UNF":"+FMT",
              FGFDT[i].attr.read_only?"+R/O":"+R/W",
              FGFDT[i].attr.old?"+OLD":"",
              FGFDT[i].attr.notpaged?"+NOT PAGED":"",
              FGFDT[i].attr.scratch?"+SCRATCH":"");
      fprintf(stderr,"\n");
}

/****************************************************************************
*                                D _ F G F D T                              *
*****************************************************************************
*
***function d_fgfdt
*
*OBJECT: Prints file characteristics and attributes of in use files in 
*        the master file table(for debugging use).
*
*ARGUMENTS: none
*
*/

void f77name(d_fgfdt)()
{
   int i=0;
   fprintf(stderr,"********** DUMP of MASTER FILE TABLE **********\n");
   for (i=0 ; i<MAXFILES ; i++) {
      if(FGFDT[i].iun != 0) dump_file_entry(i);
   }
}

/****************************************************************************
*                     R E S E T _ F I L E _ E N T R Y                       *
*****************************************************************************
*
***function reset_file_entry
*
*OBJECT: Resets a file entry in the master file table to "not in use" values.
*
*ARGUMENTS: in  i  index of the file to reset
*
*/
static void reset_file_entry(int i){
   if (FGFDT[i].file_name) free(FGFDT[i].file_name);
   if (FGFDT[i].subname)   free(FGFDT[i].subname);
   if (FGFDT[i].file_type) free(FGFDT[i].file_type);

   FGFDT[i].file_name      = (char *) NULL;
   FGFDT[i].subname        = (char *) NULL;
   FGFDT[i].file_type      = (char *) NULL;
   FGFDT[i].iun            = 0;
   FGFDT[i].fd             = -1;
   FGFDT[i].file_size      = 0;
   FGFDT[i].eff_file_size  = 0;
   FGFDT[i].lrec           = 0;
   FGFDT[i].open_flag      = 0;
   FGFDT[i].attr.stream    = 0;
   FGFDT[i].attr.std       = 0;
   FGFDT[i].attr.burp      = 0;
   FGFDT[i].attr.rnd       = 0;
   FGFDT[i].attr.wa        = 0;
   FGFDT[i].attr.ftn       = 0;
   FGFDT[i].attr.unf       = 0;
   FGFDT[i].attr.read_only = 0;
   FGFDT[i].attr.old       = 0;
   FGFDT[i].attr.scratch   = 0;
   FGFDT[i].attr.notpaged  = 0;
   FGFDT[i].attr.write_mode= 0;
   FGFDT[i].attr.remote    = 0;    /* remote file, socket wa file */
}
/****************************************************************************
*                      F I N D _ F I L E _ E N T R Y                        *
*****************************************************************************
*
***function find_file_entry
*
*OBJECT: Returns the index in the master file table of the 
*        file having iun as unit number.
*
*ARGUMENTS: in  caller   name of the function calling find_file_entry
*           in  iun      unit number
*
*RETURNS: the index in the master file table of the 
*         file having iun as unit number.
*/
static int find_file_entry(char *caller, int iun)
{
   int i;
   for (i=0; i<MAXFILES; i++)
      if (FGFDT[i].iun == iun)
         return(i);
/*   if (i == MAXFILES) { */
   fprintf(stderr,"%s error: unit %d is not associated with any file\n",caller,iun);
   return(-1);
/*   } */
}

/****************************************************************************
*                            C _ F N O M , F N O M                          *
*****************************************************************************
*
***function c_fnom, fnom
*
*OBJECT: Open a file and make the connection with a unit number.
*        Process record file attributes.
*
*ARGUMENTS:in  iun    unit number
*          in  nom    string containing the name of the file
*          in  type   string that contains the desired file attributes
*          in  lrec   length of record(must be 0 except if type contains D77)
*
*RETURNS: zero if connection is successful, non-zero otherwise.
*
*NOTES: If name is all in upper case it will be converted to lower case.
*       c_fnom is intended to be called from C.
*       fnom is intended to be called from FORTRAN.
*
*AUTHOR: Mario Lepine - RPN - nov 1995
* 
*Revision - mars 1999 - Bug fix allocation pour filename
*           avril 2008 - Correction pour reconnaissance de iun=6,output iun=5,input
*           sept 2008 - Correction du nom de fichier passe pour fichier cmcarc remote 
*
*/
int c_fnom(int *iun,char *nom,char *type,int lrec)
{
  int liun,ier = 0, minus = 0, majus = 0, lng, i, j, pid, rndflag, unfflag, lngt, junk, mode;
  char *c, *c2, *tmpdir, *cmcarc, *pos2p;
  char remote_mach[256];
  char nom2[1024];
  unsigned INT_32 hid;
  PTR_AS_INT ptr_as_int;

  if(fnom_initialized == 0) {
    /* Make sure that file descriptor 0 (stdin) will not be returned by open for use with a regular file */
    /* This is a workaround for a particular case on Linux in batch mode with PBS */
    mode = O_RDONLY;
#ifndef WIN32
    junk = open("/dev/null", mode);
    if (junk != 0)
      close(junk);
    /*    else
          printf("Debug junk associe a /dev/null\n"); */
#endif
     ARMNLIB=getenv("ARMNLIB");
     if( ARMNLIB == NULL ) ARMNLIB = LOCALDIR;
     AFSISIO=getenv("AFSISIO");
     if( AFSISIO == NULL ) AFSISIO = LOCALDIR;
     for (i=0; i<MAXFILES; i++) reset_file_entry(i);
     hid = f77name(check_host_id)();
     fnom_initialized=1;
     }

  if (((PTR_AS_INT) iun > 0) && ((PTR_AS_INT) iun < 1000)) {
     /* an integer value has been passed to c_fnom as iun */
     ptr_as_int = (PTR_AS_INT) iun ;
     liun = ptr_as_int;
     }
  else {
     /* a pointer has been passed to c_fnom as iun */
     if (*iun == 0)
        *iun = c_qqqfscr(type);
        if (*iun == -1) {
	  fprintf(stderr,"C_FNOM ERROR: no more units available\n");
	  return(-1);
	}
     liun = *iun;
     }

  if ((liun ==6) && ((strcmp(nom,"$OUT") == 0) || (strcmp(nom,"$OUTPUT") == 0) || (strcmp(nom,"OUTPUT") == 0) ||                     (strcmp(nom,"output") == 0))) {
    stdoutflag=1;
/*    fprintf(stderr,"C_FNOM DEBUG already connected: iun=%d filename=%s\n",liun,nom);   */
    return(0);
  }
  if ((liun ==5) && ((strcmp(nom,"$IN") == 0) || (strcmp(nom,"$INPUT") == 0) || (strcmp(nom,"INPUT") == 0) ||                     (strcmp(nom,"input") == 0))) {
    stdinflag=1;
/*    fprintf(stderr,"C_FNOM DEBUG already connected: iun=%d filename=%s\n",liun,nom);   */
    return(0);
  }

  if ((liun == 6) || (liun == -2)) {
    c = nom;
    c2 = nom2;
    nom2[strlen(nom)] = '\0';
    for (j = 0; (j < strlen(nom) && j < 1024); j++, c++, c2++) {
       if (islower(*c)) {
         minus=1;
         *c2 = *c;
       }
       else if (isupper(*c)) {
         majus=1;
         *c2 = tolower(*c);
       }
       else
         *c2 = *c;
    }
  }
  if (liun == 6) {
    fclose(stdout);
    if (minus && majus)
      freopen(nom,"a",stdout);
    else
      freopen(nom2,"a",stdout);
/*    fprintf(stderr,"C_FNOM DEBUG: freopen %s pour stdout\n",nom2) ; */
    stdoutflag=1;
    return(0);
  }
  else if (liun == -2) {
    fclose(stderr);
    if (minus && majus)
      freopen(nom,"a",stderr);
    else
      freopen(nom2,"a",stderr);
    return(0);
  }
  for (i=0; i<MAXFILES; i++)
     if (FGFDT[i].iun == liun) {
        fprintf(stderr,"c_fnom error: unit %d is already in use\n", liun);
        return(-1);
        }
  for (i=0; i<MAXFILES; i++)
     if (FGFDT[i].iun == 0) {
        FGFDT[i].iun = liun;
        break;
        }
  if (i == MAXFILES) {
     fprintf(stderr,"c_fnom error: too many files, file table is full\n");
     return(-1);
     }
/*
 *   record file attributes
 */ 
  lngt = strlen(type) + 1;
  FGFDT[i].file_type = malloc(lngt+1);
  strncpy(FGFDT[i].file_type,type,lngt);
  FGFDT[i].attr.stream=0;
  FGFDT[i].attr.std=0;
  FGFDT[i].attr.burp=0;
  FGFDT[i].attr.rnd=0;
  FGFDT[i].attr.wa=0;
  FGFDT[i].attr.ftn=0;
  FGFDT[i].attr.unf=0;
  FGFDT[i].attr.read_only=0;
  FGFDT[i].attr.old=0;
  FGFDT[i].attr.notpaged=0;
  FGFDT[i].attr.scratch = 0;
  FGFDT[i].attr.pipe = 0;
  FGFDT[i].attr.remote=0;

  if (strstr(type,"STREAM") || strstr(type,"stream")){ FGFDT[i].attr.stream=1;
                                                       FGFDT[i].attr.rnd=1; }
  if (strstr(type,"STD")    || strstr(type,"std"))   { FGFDT[i].attr.std=1;
                                                       FGFDT[i].attr.rnd=1; }
  if (strstr(type,"BURP")   || strstr(type,"burp"))  { FGFDT[i].attr.burp=1;
                                                       FGFDT[i].attr.rnd=1; }
  if (strstr(type,"RND")    || strstr(type,"rnd"))     FGFDT[i].attr.rnd=1;
  if (strstr(type,"WA")     || strstr(type,"wa"))      FGFDT[i].attr.rnd=1;   /* wa attribute will be set by waopen */
  if (strstr(type,"FTN")    || strstr(type,"ftn"))   { FGFDT[i].attr.ftn=1;
                                                       FGFDT[i].attr.rnd=0; }
  if (strstr(type,"UNF")    || strstr(type,"unf"))   { FGFDT[i].attr.unf=1;
                                                       FGFDT[i].attr.ftn=1;
                                                       FGFDT[i].attr.rnd=0; }
  if (strstr(type,"OLD")    || strstr(type,"old"))     FGFDT[i].attr.old=1;
  if (strstr(type,"R/O")    || strstr(type,"r/o"))   { FGFDT[i].attr.read_only=1;
                                                       FGFDT[i].attr.old=1; }
  if (strstr(type,"R/W")    || strstr(type,"r/w"))   { FGFDT[i].attr.read_only=0;
                                                       FGFDT[i].attr.write_mode=1; }
  if (strstr(type,"D77")    || strstr(type,"d77"))   { FGFDT[i].attr.ftn=1;
                                                       FGFDT[i].attr.rnd=1; }
  if (strstr(type,"SCRATCH") || strstr(type,"scratch")) FGFDT[i].attr.scratch=1;
  if (strstr(type,"REMOTE") || strstr(type,"remote")) { FGFDT[i].attr.remote=1; }
    
  if (!FGFDT[i].attr.std && !FGFDT[i].attr.burp && 
      !FGFDT[i].attr.wa && !FGFDT[i].attr.rnd  && !FGFDT[i].attr.stream)
     FGFDT[i].attr.ftn = 1;
  FGFDT[i].lrec = lrec;
  FGFDT[i].open_flag = 0;

/*
 *   if scratch file, add tmpdir directory and pid to the file name 
 */
  if (FGFDT[i].attr.scratch) {
     if(strstr(nom,"/")){
        fprintf(stderr,"c_fnom error: / is illegal in scratch file name\n");
        fprintf(stderr,"              specified name was %s\n",nom);
        return(-1);
     }
     pid = getpid();
     tmpdir = getenv("TMPDIR");
     if (tmpdir == NULL) {
       fprintf(stderr,
"c_fnom warning: TMPDIR environment variable is not defined, /tmp is used\n");
       tmpdir = "/tmp";
     }
     lng = strlen(nom) + strlen(tmpdir) + 10 + 3 + 128;    /* espace tampon supplementaire */
     if ((FGFDT[i].file_name = malloc(lng)) == NULL) {
       fprintf(stderr,"c_fnom error: can't allocate memory for file name\n");
       perror("c_fnom");
       exit(1);
     }
     sprintf(FGFDT[i].file_name, "%s/%d_%s", tmpdir, pid, nom);
  }
  else {
/*
 *   convert file name to lower case unless it contains a mix of 
 *   upper case / lower case
 */
     lng = strlen(nom);
     FGFDT[i].file_name = malloc(lng+1);
     if ((FGFDT[i].attr.remote) && (strchr(nom, ':'))) {
       if (FGFDT[i].attr.rnd) {
/*          remote_mach = strtok(nom,":"); */
          pos2p = strchr(nom, ':');
          if (pos2p != NULL) {
            strncpy(remote_mach, nom, pos2p-nom);
            remote_mach[pos2p-nom]='\0';
            nom = ++pos2p;
            printf("Debug+ remote_mach=%s file name=%s\n", remote_mach, nom);
            lng = strlen(nom);
          }
       }  
       else {
          /* code to remote copy the file on local machine and change file name */
          FGFDT[i].attr.remote=0;
          FGFDT[i].attr.read_only = 0;     /* file becomes read only */
       }
     }
     else 
        FGFDT[i].attr.remote=0;
     c = nom;
     if (nom[0] == '@') {  /* name is of the form @some_file_name */
       c++;              /* skip the @, scan later under        */
       lng--;            /* AFSISIO & ARMNLIB if not local file */
       }
     if (nom[0] == '%') {  /* name is of the form %[%@]some_pipe_ file */
       c++;              /* skip the %  */
       lng--;
       FGFDT[i].attr.pipe = 1;
       }
     if (nom[0] == '+') {  /* name is of the form +some_file_name */
       c++;              /* skip the +, do not convert to lowercase */
       lng--;            
       strncpy(FGFDT[i].file_name, nom+1, lng);
       c2 = FGFDT[i].file_name ;
       *(c2+lng)  = '\0';
       }
     else {
       c2 = FGFDT[i].file_name ;
       *(c2 + lng) = '\0';
       for (j = 0; j < lng; j++, c++, c2++) {
         if (islower(*c)) {
           minus = 1;
           *c2 = *c;
         }
         else if (isupper(*c)) {
           majus=1;
           *c2 = tolower(*c);
         }
         else
           *c2 = *c;
       }
       if (majus && minus)
         strncpy(FGFDT[i].file_name,nom,lng);
     }
  }
/*
 *   verify for a cmcarc type of file (filename@subfilename)
 */
  if ((cmcarc = strchr(FGFDT[i].file_name,'@')) && !(FGFDT[i].attr.remote)) {
     FGFDT[i].subname = malloc(lng+1);
     strcpy(FGFDT[i].subname,cmcarc+1);
     *cmcarc = '\0';
     FGFDT[i].attr.old = 1;  /* file must exist and be read/only */
     FGFDT[i].attr.read_only = 1;
     }
  else
     FGFDT[i].subname = NULL;
/*
 *   check for @file (@filename)
 */
  if (nom[0]=='@') {
     char filename[1024];
     struct stat etat,etat2,etat3;
     strcpy(filename,FGFDT[i].file_name);

     if (access(filename,F_OK) == -1) {         /* no local file */
       sprintf(filename,"%s/datafiles/constants/%s",AFSISIO,FGFDT[i].file_name);

       if (access(filename,F_OK)  == -1) {                /* not under AFSISIO */
          sprintf(filename,"%s/data/%s",ARMNLIB,FGFDT[i].file_name);

          if (access(filename,F_OK)  == -1) {             /* not under ARMNLIB either */
             return(-1);
             }
          }
       }
     free(FGFDT[i].file_name);
     FGFDT[i].file_name = malloc(strlen(filename)+10);
     strcpy(FGFDT[i].file_name,filename);
     lng=strlen(filename);
     }

  if ((FGFDT[i].attr.old || FGFDT[i].attr.read_only) && ! FGFDT[i].attr.remote)
     if (!f77name(existe)(FGFDT[i].file_name,(F2Cl) strlen(FGFDT[i].file_name))) {
        fprintf(stderr,"c_fnom error: file %s should exist and does not\n",FGFDT[i].file_name);
        junk=c_fclos(liun);
        return(-1);
        }
/*
 *   FORTRAN files must be opened by a FORTRAN module
 */
  ier = 0;
  if (FGFDT[i].attr.ftn) {
     ftnword iun77=liun;
     ftnword lrec77=lrec;
     ftnword rndflag77 = FGFDT[i].attr.rnd;
     ftnword unfflag77 = FGFDT[i].attr.unf;
     ftnword lmult = D77MULT;
     ier = open64(FGFDT[i].file_name,O_RDONLY | WIN32_O_BINARY);
     if (ier <=0) {
        FGFDT[i].file_size = -1;
        FGFDT[i].eff_file_size = -1;
        }
     else {
        LLSK dimm=0;
        dimm = LSEEK(ier,dimm,L_XTND);
        FGFDT[i].file_size = dimm / sizeof(word);
        FGFDT[i].eff_file_size = dimm / sizeof(word);
        close(ier);
        }
     ier = f77name(qqqf7op)(&iun77,FGFDT[i].file_name,&lrec77,&rndflag77,&unfflag77,&lmult,(F2Cl) lng);
  }
  else if (FGFDT[i].attr.stream || FGFDT[i].attr.std || FGFDT[i].attr.burp || FGFDT[i].attr.wa ||
          (FGFDT[i].attr.rnd && !FGFDT[i].attr.ftn) ) {
     ier = c_waopen2(liun);
     FGFDT[i].attr.wa = 0;   /* will be set by waopen */
  }
  
  if (FGFDT[i].attr.remote) ier = fnom_rem_connect(i,remote_mach);
    
  if (ier == 0) FGFDT[i].open_flag = 1;
  if (ier < 0) junk=c_fclos(liun);
  return(ier<0?-1:0);
}  

ftnword f77name(fnom)(ftnword *iun,char *nom,char *type,ftnword *flrec,F2Cl l1,F2Cl l2)
{
   int lrec,lng,tmp,liun=*iun;
   char filename[1025],filetype[257];

   lrec = *flrec;

   lng = (l1 <= 1024) ? l1 : 1024;
   strncpy(filename,nom,lng);        /*  copy filename into a C string  */
   filename[lng] = '\0';

   while ((filename[lng-1] == ' ') && (lng > 1)) {  /* strip trailing blanks */
      lng--;
      filename[lng] = '\0';
      }

   lng = (l2 <= 256) ? l2 : 256;
   strncpy(filetype,type,lng);   /*  copy file type into a C string  */
   filetype[lng] = '\0';

   while ((filetype[lng-1] == ' ') && (lng > 1)) { /* strip trailing blanks */
      lng--;
      filetype[lng] = '\0';
      }

   tmp=(c_fnom(&liun,filename,filetype,lrec));
   if(*iun==0) *iun = liun;
   return (tmp);
}


/****************************************************************************
*                            C _ F C L O S,   F C L O S                     *
*****************************************************************************
*
***functions c_fclos and fclos
*
*OBJECT: Close file associated with unit iun.
*        Returns zero 
*        and non-zero otherwise.
*
*ARGUMENTS: in iun   unit number
*
*RETURNS: zero if the connection is successful, non-zero 
*         otherwise
*
*/
int c_fclos(int iun)
{
   int i,ier;
   ftnword iun77;

/*
   for (i=0; i<MAXFILES; i++)
      if (FGFDT[i].iun == iun)
         break;
   if (i == MAXFILES) {
   fprintf(stderr,"c_fclos error: unit %d is not associated with any file\n",iun);
   return(-2);
   }
*/

   if ((iun == 6) && (stdoutflag)) return(0);
   if ((iun == 5) && (stdinflag)) return(0);

   if ((i=find_file_entry("c_fclos",iun)) < 0) return(i);
   iun77 = iun;
   ier=0;
   if (FGFDT[i].open_flag){
      if (FGFDT[i].attr.ftn)
         ier = f77name(ftnclos)(&iun77);
      else
         ier = close(FGFDT[i].fd);
      }

   reset_file_entry(i);
   return(ier);
}

ftnword f77name(fclos)(ftnword *fiun)
{
   int iun,ier;
   iun = *fiun;
   ier=c_fclos(iun);
   return(ier);
}


/****************************************************************************
*                          C _ Q Q Q F S C R                                *
*****************************************************************************
*
***function c_qqqfscr
*
*OBJECT: Generates a scratch unit number. 
*
*ARGUMENTS: in type  string that contains the file attributes (see FNOM)
*
*RETURNS: the scratch unit number
*
*/
static int c_qqqfscr(char *type)          
{
   int iun,i,j,inused,start;

   iun = -1;
   if (strstr(type,"FTN") || strstr(type,"ftn") || strstr(type,"D77") || strstr(type,"d77")) 
      start = 99;
   else
      start = 999;
   for (j=start; j>10; j--) {
      inused = 0;
      for (i=0; i<MAXFILES; i++)
         if (FGFDT[i].iun == j) {
            inused = 1;
            break;
            }
      if (! inused) {
         iun = j;
         break;
         }
      }
   return(iun);
}

/****************************************************************************
*                             Q Q Q F N O M                                 *
*****************************************************************************
*
***function qqqfnom
*
*OBJECT: Obtain some information about a unit: file name, file type,
*        lenght of record.
*
*ARGUMENTS:in  iun     unit number
*          out nom     file name(blank padded to have a total length of l1)
*          out type    type of file(blank padded to have a total length of l2)
*          out flrec   length of record
*          in  l1      length of filename
*          in  l2      length of type 
*
*/
ftnword f77name(qqqfnom)(ftnword *iun,char *nom,char *type,ftnword *flrec,F2Cl l1,F2Cl l2)
{
   int i,j;

   j=*iun;
   if ((i=find_file_entry("qqqfnom",j)) < 0) return(i);

   strncpy(nom,FGFDT[i].file_name,l1);
   for(j=strlen(FGFDT[i].file_name);j<l1;j++) nom[j]=' ';
   strncpy(type,FGFDT[i].file_type,l2);
   for(j=strlen(FGFDT[i].file_type);j<l2;j++) type[j]=' ';
   *flrec=FGFDT[i].lrec;
   return(0);
}
/****************************************************************************
*                              Q Q C C L O S                                *
*****************************************************************************
*
***function qqcclos
*
*OBJECT: Closes a file, given its index. 
*
*ARGUMENTS: in indf   index of the file in the master file table
*
*RETURNS: zero if file correctly closed, non-zero otherwise
*
*/
static int qqcclos(int indf)
{
  int lfd=FGFDT[indf].fd;
  int ind;

  ind = 0;
  while ((wafile[ind].file_desc != lfd) && (ind < MAXWAFILES))
     ind++;
  if (ind == MAXWAFILES) {
     fprintf(stderr,"qqcclos error: file is not open, fd=%d, name=%s\n",
                    lfd,FGFDT[indf].file_name);
     return(1);
     }
     
  if (FGFDT[indf].attr.remote) {
    int *s_ID, *addr, *nw, *RW_mode, *checksum;
    int sock_comm_ID=0xBABE;
    int demande[5];
    int nc;

    s_ID = &(demande[0]);
    addr = &(demande[1]);
    nw = &(demande[2]);
    RW_mode = &(demande[3]);
    checksum = &(demande[4]);
    *s_ID = sock_comm_ID;
    *addr = 0;
    *nw = 0;
    *RW_mode = 3;  /* close request */
    *checksum = *s_ID ^ *addr ^ *nw ^ *RW_mode;
    check_swap_records(demande,5,sizeof(int));
    nc=write_stream(FGFDT[indf].fd,demande,5*sizeof(int));
    if (nc == 0) {
      printf("socket qqcclos wrote to server OK\n");
      fflush(stdout);
      }
    else {
      fprintf(stderr,"socket qqcclos error: wrote only %d bytes to server\n",nc);
      fflush(stderr);
      }
  } /* end remote */
	else {
    if (WA_PAGE_SIZE != 0) {
      wa_pages_flush(ind);
      if (wafile[ind].nb_page_in_use != 0)
          fprintf(stderr,"qqcclos error: nb_page_in_use = %d\n",wafile[ind].nb_page_in_use);
      FGFDT[indf].file_size = 0;
      FGFDT[indf].eff_file_size = 0;
      wafile[ind].nb_page_in_use = 0;
      if (debug_mode > 1) {
        fprintf(stderr,"Debug fermeture du fichier ind=%d, fd=%d\n",ind,lfd);
        }
      }    
	}
  wafile[ind].file_desc = -1;
  FGFDT[indf].fd = -1;
  FGFDT[indf].open_flag = 0;
  close(lfd);
  return(0);
}

/****************************************************************************
*  C _ W A O P E N ,   C _ W A O P E N 2 ,   W A O P E N 2 ,   W A O P E N  *
*****************************************************************************
*
***function c_waopen
*
*OBJECT: Opens a word addressable file.
*
*ARGUMENTS: in iun   unit number
*
*RETURNS: (only for c_waopen2 and waopen2) zero if file correctly opened, non-zero otherwise
*
*/
void c_waopen(int iun) { int scrap=c_waopen2(iun); if(scrap<=0) exit(1); }
int c_waopen2(int iun)   /* open unit iun for WORD ADDRESSABLE access */
{
   int i,ier;

   for (i=0; i < MAXFILES; i++)
      if (FGFDT[i].iun == iun)
         break;
   if (i == MAXFILES) {
      for (i=0; i < MAXFILES; i++)
         if (FGFDT[i].iun == 0) {
            FGFDT[i].iun = iun;
            break;
            }
      if (i == MAXFILES) {
        fprintf(stderr,"c_waopen error: file table is full\n");
        return(-1);
        }
/*
 *  file is not associated with fnom, file name is set to Wafileiun
 */
      FGFDT[i].file_name = malloc(10);
      sprintf(FGFDT[i].file_name,"%s%d","Wafile",iun);
      FGFDT[i].attr.wa = 1;
      FGFDT[i].attr.rnd = 1;
      }
   else {
      if (FGFDT[i].attr.rnd == 0) {
         fprintf(stderr,"c_waopen error: waopen needs a file with the RND or WA type\n");
         return(-1);
      }
      if (FGFDT[i].open_flag) {
         if (FGFDT[i].attr.wa == 1){ /* fnom opens the file but does not set wa flag */
            fprintf(stderr,"c_waopen warning: unit %d already open as %s\n",iun,FGFDT[i].file_name); 
            }
         FGFDT[i].attr.wa = 1;
         return(FGFDT[i].fd);
         }
      }
   ier = qqcopen(i);
   if (ier >=0) { FGFDT[i].open_flag = 1; FGFDT[i].attr.wa = 1; FGFDT[i].attr.rnd = 1; }
   return ier;
}

ftnword f77name(waopen2)(ftnword *fiun)
{
   int iun;
   iun = *fiun;
   return(c_waopen2(iun));
}

void f77name(waopen)(ftnword *fiun)
{
   int iun;
   iun = *fiun;
   c_waopen(iun);
}
         
/****************************************************************************
* C _ W A C L O S ,   C _ W A C L O S 2 ,   W A C L O S ,   W A C L O S 2   *
*****************************************************************************
*
***function c_waclos
*
*OBJECT: Closes a word addressable file.
*
*ARGUMENTS: in iun   unit number
*
*RETURNS: (only for c_waclos2 and waclos2) zero if closure is successful,
*         non-zero otherwise.
*
*/
void c_waclos(int iun) { int scrap =  c_waclos2(iun) ; }
int c_waclos2(int iun)
{
   int i,ier;

   if ((i=find_file_entry("c_waclos",iun)) < 0) return(i);

   if (! FGFDT[i].open_flag) {
      fprintf(stderr,"c_waclos error: unit %d is not open\n",iun);
      return(-1);
      }

   ier = qqcclos(i);
   FGFDT[i].open_flag = 0;
   FGFDT[i].attr.wa = 0;
   return(ier);
}
ftnword f77name(waclos2)(ftnword *fiun)
{
   int iun;
   iun = *fiun;
   return(c_waclos2(iun));
}
void f77name(waclos)(ftnword *fiun)
{
   int iun;
   iun = *fiun;
   iun=c_waclos2(iun);
}


/****************************************************************************
*   C _ W A W R I T ,   C _ W A W R I T 2 ,   W A W R I T ,   W A W R I T 2 *
*****************************************************************************
*
***function c_wawrit, c_wawrit2, wawrit, wawrit2
*
*OBJECT: Writes in a word addressable. c_wawrit2 and wawrit2 are functions
*        returning an error code.
*
*ARGUMENTS: in iun    unit number
*           in buf    will contain data
*           in adr    file address
*           in nmots  number of words
*
*RETURNS: (only for c-wawrit2 and wawrit2) the number of words written.
*
*/

void c_wawrit(int iun,void *buf,unsigned int adr,int nmots)
{
  c_wawrit2(iun,buf,adr,nmots);
}
int c_wawrit2(int iun,void *buf,unsigned int adr,int nmots)
{
#define WA_HOLE 2048
   int i,ier;
   word scrap[WA_HOLE];
   word *bufswap = (word *) buf;

   if ((i=find_file_entry("c_wawrit",iun)) < 0) return(i);

   if (! FGFDT[i].open_flag) {
      fprintf(stderr,"c_wawrit error: unit %d is not open\n",iun);
      return(-1);
      }
   if ( FGFDT[i].attr.read_only != 0 ) {
      fprintf(stderr,"c_wawrit error: unit %d ,file= %s is READ ONLY\n",
                     iun,FGFDT[i].file_name);
      return(-1);
      }
   if ( adr > FGFDT[i].file_size+WA_HOLE ) {
      fprintf(stderr,"c_wawrit error: attempt to write beyond EOF+%d\n",WA_HOLE);
      fprintf(stderr,"                unit = %d, adr=%u > file_size=%d\n",
                     iun,adr,FGFDT[i].file_size);
      fprintf(stderr,"                filename=%s\n",FGFDT[i].file_name);
      exit(1);
      }
   if ( adr > FGFDT[i].file_size+1 ){
      qqcwawr(scrap,FGFDT[i].file_size+1,adr-FGFDT[i].file_size,i);
      }
   if (*little_endian) swap_buffer_endianness(bufswap,nmots)
   qqcwawr((word *)buf,adr,nmots,i);
   if (*little_endian) swap_buffer_endianness(bufswap,nmots)
   return( nmots>0 ? nmots : 0);
}
void f77name(wawrit)(ftnword *fiun,void *buf,unsigned ftnword *fadr,ftnword *fnmots){
     f77name(wawrit2)(fiun,buf,fadr,fnmots);
     }
ftnword f77name(wawrit2)(ftnword *fiun,void *buf,unsigned ftnword *fadr,ftnword *fnmots)
{
   int iun,adr,nmots;
   iun = *fiun; adr = *fadr; nmots = *fnmots;
#if defined (ALL64)
   if ( adr > 0 )
      return(c_wawrit2(iun,buf,(2*adr)-1,nmots*2));
   else
      return(c_wawrit2(iun,buf,adr,nmots));
#else
   return(c_wawrit2(iun,buf,adr,nmots));
#endif
}

/****************************************************************************
* C _ W A R E A D ,   C _ W A R E A D 2 ,   W A R E A D ,   W A R E A D 2   *
*****************************************************************************
*
***function c_waread, c_waread2, waread, waread2
*
*OBJECT: Read from a word addressable file.
*
*ARGUMENTS: in  iun     unit number
*           out buf     will contain the data read
*           in  adr     adress to start from
*           in  nmots   number of words to read
*
*RETURNS: (only for c-waread2 and waread2) the number of words read.
*
*/
void c_waread(int iun,void *buf,unsigned int adr,int nmots)
{
  int ier, i;
  ier = c_waread2(iun,buf,adr,nmots);
  if (ier == -2) {
    i = find_file_entry("c_waread",iun);
    fprintf(stderr,
            "c_waread error: attempt to read beyond EOF, of file %s\n",
            FGFDT[i].file_name);
    fprintf(stderr,"                addr = %u, EOF = %d\n",
            adr,FGFDT[i].eff_file_size);
  }

}
int c_waread2(int iun,void *buf,unsigned int adr,int nmots)
{
   int i,ier;
   word *bufswap = (word *) buf;

   if ((i=find_file_entry("c_waread",iun)) < 0) return(i);
   
   if (! FGFDT[i].open_flag) {
      fprintf(stderr,"c_waread error: unit %d is not open\n",iun);
      return(-1);
      }

   if ( adr > FGFDT[i].eff_file_size+2 ) return(-2);

   if ( FGFDT[i].eff_file_size == 0 ) return(0);

   if ( adr+nmots-1 > FGFDT[i].eff_file_size ) {
      nmots -= (adr+nmots-1-FGFDT[i].eff_file_size);
      }
   if ( nmots == 0 ) return(0);
   qqcward((word *)buf,adr,nmots,i);
   if (*little_endian) swap_buffer_endianness(bufswap,nmots)
   return(nmots);
}
void f77name(waread)(ftnword *fiun,void *buf,unsigned ftnword *fadr,
                     ftnword *fnmots)
{
  f77name(waread2)(fiun,buf,fadr,fnmots);
}
ftnword f77name(waread2)(ftnword *fiun,void *buf,unsigned ftnword *fadr,
                         ftnword *fnmots)
{
   int iun,adr,nmots;
   iun = *fiun; adr = *fadr; nmots = *fnmots;
#if defined (ALL64)
   if ( adr > 0 )
      return(c_waread2(iun,buf,(2*adr)-1,nmots*2));
   else
      return(c_waread2(iun,buf,adr,nmots));
#else
   return(c_waread2(iun,buf,adr,nmots));
#endif
}

/****************************************************************************
*                      C _ W A S I Z E ,   W A S I Z E                      *
*****************************************************************************
*
***function c_wasize, wasize
*
*OBJECT: Returns the size (in words) of a file, given its unit number.
*
*ARGUMENTS: in iun   unit number
*
*RETURNS: the size of a file in words.
*
*NOTE: Gives the size of any file that can be passed to fnom.
*
*/
INT_32 c_wasize(int iun)
{
   int i,ier;
   word n;

   if ((i=find_file_entry("c_wasize",iun)) < 0) return(i);

   if (! FGFDT[i].open_flag) {
      ier = qqcopen(i);
      n = FGFDT[i].eff_file_size;
      ier = qqcclos(i);
      }
   else
      n = FGFDT[i].eff_file_size;

   return(n);
}
ftnword f77name(wasize)(ftnword *fiun)  /* return file size in FORTRAN WORDS */
{
   int iun;
   iun = *fiun;
#if defined (ALL64)
   return(c_wasize(iun)/2);
#else
   return(c_wasize(iun));
#endif
}

/****************************************************************************
*                    C _ N U M B L K S ,   N U M B L K S                    *
*****************************************************************************
*
***function c_numblks, numblks
*
*OBJECT: Returns the size of a file in kilobytes.
*
*ARGUMENTS: in iun   unit number
*
*RETURNS: the size of a file in kilobytes.
*
*NOTE: It uses c_wasize.
*
*/
INT_32 c_numblks(int iun)
{
   int i,ier;
   int n;

   n = c_wasize(iun);
   if( n<0 ) return (n);
   i = 1024 / sizeof(word);
   return ( (n+i-1) / i );
}
ftnword f77name(numblks)(ftnword *fiun)     /* return file size in KiloBytes */
{
   int iun;
   iun = *fiun;
   return(c_numblks(iun));
}


/****************************************************************************
*                             E X I S T E                                   *
*****************************************************************************
*
***function existe
*
*OBJECT: Checks if file exists.
*
*ARGUMENTS: in  nom  name of the file
*           in  lng  length of name 
*
*RETURNS: one if the file exists,
*         zero if it doesn't.
*         
*
*/
ftnword f77name(existe)(char *nom,F2Cl llng) 
{
   struct stat etat;
   int l2, lng=llng;
   char filename[257];

   l2 = (lng <= 256) ? lng : 256;
   strncpy(filename,nom,l2);
   filename[lng] = '\0';

   while ((filename[lng-1] == ' ') && (lng > 1)) {
      lng--;
      filename[lng] = '\0';
      }
   
   if (access(filename, F_OK) == -1)
      return(0);                     /* file does not exist */
   else
      return(1);                     /* file exists */
}

/****************************************************************************
*                    C _ O P E N D A ,   O P E N D A                        *
*****************************************************************************
*
***function c_openda, openda
*
*OBJECT: Opens a direct access file.
*
*ARGUMENTS: in iun  unit number
*
*/
void c_openda(int iun)
{
   c_waopen(iun);
}
void f77name(openda)(ftnword *iun)
{
   int liun;
   liun = (int) *iun;
   c_waopen(liun);
}

/****************************************************************************
*                   C _ C L O S D A ,   C L O S D A                         *
*****************************************************************************
*
***function c_closda, closda
*
*OBJECT: Closes a direct access file.
*
*ARGUMENTS: in iun unit number
*
*/
void c_closda(int iun)
{
   c_waclos(iun);
}
void f77name(closda)(ftnword *iun)
{
   int liun;
   liun = (int) *iun;
   c_closda(liun);
}

/****************************************************************************
*                     C _ C H E C D A ,   C H E C D A                       *
*****************************************************************************
*
***function c_checda, checda
*
*OBJECT: Checks that I/O is done.
*
*ARGUMENTS: in iun   unit number
*
*/
void c_checda(int iun)
{
   int *pt;
   for (pt = dastat; pt < &dastat[MAXWAFILES]; pt++)
      if (*pt == iun) {
         *pt = 0;
         break ;
         }
}
void f77name(checda)(ftnword *iun)
{
   int liun;
   liun = (int) *iun;
   c_checda(liun);
}

/****************************************************************************
*                     C _ R E A D D A ,   R E A D D A                       *
*****************************************************************************
*
***function c_readda, readda
*
*OBJECT: Reads a direct access file.
*
*ARGUMENTS: in  iun      unit number
*           out bufptr   will contain the data read
*           in  ns       number of words to read
*           in  is       word to start from
*
*/
void c_readda(int iun,int *bufptr,int ns,int is)
{
   int *pt;
   for (pt = dastat; pt < &dastat[MAXWAFILES]; pt++)
     if( *pt == iun ) break ;

   if ( pt <= &dastat[MAXWAFILES] )
   if (*pt == iun)
      {
      fprintf(stderr,"readda error: consecutive calls to readda without call to checda, iun=%d \n",iun);
      return;
      }
   c_waread(iun,bufptr,(is-1)*BLKSIZE+1,ns*BLKSIZE);
   for (pt = dastat; pt < &dastat[MAXWAFILES]; pt++)
      if( *pt == 0 ) break ;
   if (  pt >= &dastat[MAXWAFILES] ) {
      fprintf(stderr," *** readda error: recompiler avec MAXWAFILES++ \n");
      return;
      }
   *pt = iun;
}
void f77name(readda)(ftnword *iun,ftnword *bufptr,ftnword *ns,ftnword *is)
{
   int liun,lns,lis,save=BLKSIZE;
   liun = (int) *iun;
   lns = (int) *ns;
   lis = (int) *is;
   BLKSIZE = BLKSIZE * (sizeof(ftnword)/sizeof(word));
   c_readda(liun,bufptr,lns,lis);
   BLKSIZE=save;
}
/****************************************************************************
*                      C _ W R I T D A ,   W R I T D A                      *
*****************************************************************************
*
***function c_writda, writda
*
*OBJECT: Writes to a direct access file.
*
*ARGUMENTS: in  iun     unit number
*           in  bufptr  will contain the data read
*           in  ns      number of words to write
*           in  is      word to start from 
*
*/
void c_writda(int iun,int *bufptr,int ns,int is)
{
   int *pt;
   for (pt = dastat; pt < &dastat[MAXWAFILES]; pt++)
      if( *pt == iun ) break ;
   if ( pt <= &dastat[MAXWAFILES] )
   if (*pt == iun)
      {
      fprintf(stderr,"writda error: consecutive calls to writda without call to checda, iun=%d \n",iun);
      return;
      }
   c_wawrit(iun,bufptr,(is-1)*BLKSIZE+1,ns*BLKSIZE);
   for (pt = dastat; pt < &dastat[MAXWAFILES]; pt++)
      if( *pt == 0 ) break ;
   if (  pt >= &dastat[MAXWAFILES] ) {
      fprintf(stderr," *** writda error: recompiler avec MAXWAFILES++ \n");
      return;
      }
   *pt = iun;
}

void f77name(writda)(ftnword *iun,ftnword *bufptr,ftnword *ns,ftnword *is)
{
   int liun,lns,lis,save=BLKSIZE;
   liun = (int) *iun;
   lns = (int) *ns;
   lis = (int) *is;
   BLKSIZE = BLKSIZE * (sizeof(ftnword)/sizeof(word));
   c_writda(liun,bufptr,lns,lis);
   BLKSIZE=save;
}

/***************************************************************************
*                     C _ G E T F D S C ,   G E T F D S C                  *
****************************************************************************
*
***function c_getfdsc,getfdsc
*
*OBJECT: Get file descriptor associated to unit iun. 
*
*ARGUMENTS: in iun   unit number
*
*RETURNS: a file descriptor, if there is no error,
*         or a negative number if there is an error.
*
*/
int c_getfdsc(int iun) { 
   int i;

   if ((i=find_file_entry("c_getfdsc",iun)) < 0) return(i);

   if (! FGFDT[i].attr.stream) {
      fprintf(stderr,"c_getfdsc error: unit %d does not have the STREAM attribute\n",iun);
      return(-1);
      }
   if (! FGFDT[i].open_flag) {
      fprintf(stderr,"c_getfdsc error: unit %d is not open\n",iun);
      return(-1);
      }

   return(FGFDT[i].fd) ;
   }
ftnword f77name(getfdsc)( ftnword *iun) { return(c_getfdsc((int) *iun)) ;}

/***************************************************************************
*                     C _ S Q O P E N ,   S Q O P E N                      *
****************************************************************************
*
***function c_sqopen, sqopen
*
*OBJECT: Opens a stream.
*
*ARGUMENTS: in iun   unit number
*
*/

void c_socket_open()
{
}
void c_sqopen(int iun)
{
  int i;

  if ((i=find_file_entry("c_sqopen",iun)) < 0) return;

  if (FGFDT[i].attr.pipe) {
    if ((FGFDT[i].file_name[0] == '@') || (FGFDT[i].file_name[0] == '%'))
      c_socket_open();
    else 
      if ((FGFDT[i].file_name[1] == '0') || (FGFDT[i].file_name[1] == '1'))
        i = sscanf(&(FGFDT[i].file_name[1]),"%d",&(FGFDT[i].fd));
      else {
        FGFDT[i].file_name++;
        c_waopen(iun);
      }
  }
  else
    c_waopen(iun) ;
}
void f77name(sqopen)(ftnword *iun) { c_sqopen((int) *iun) ; }

/***************************************************************************
*                     C _ S Q C L O S ,   S Q C L O S                      *
****************************************************************************
*
***function c_sqclos, sqclos
*
*OBJECT: Closes a stream.
*
*ARGUMENTS: in iun  unit number
*
*/
void c_sqclos(int iun) 
{
   int i;

   if ((i=find_file_entry("c_sqclos",iun)) < 0) return;
   if (FGFDT[i].attr.wa == 1) c_waclos(iun) ;
}
void f77name(sqclos)(ftnword *iun) { c_sqclos((int) *iun) ; }

/***************************************************************************
*                     C _ S Q R E W ,   S Q R E W                          *
****************************************************************************
*
***function c_sqrew, sqrew
*
*OBJECT: Rewinds a stream
*
*ARGUMENTS: in iun  unit number
*
*/
void c_sqrew(int iun) {
   int fd ;
   int i;

   if ((i=find_file_entry("c_sqrew",iun)) < 0) return;
   
   if (FGFDT[i].attr.pipe) return;

   fd = c_getfdsc(iun);
   if (fd <= 0) return;
   lseek(fd,0,L_SET);
}
void f77name(sqrew)(ftnword *iun) { c_sqrew((int) *iun) ; }

/***************************************************************************
*                     C _ S Q E O I ,   S Q E O I                          *
****************************************************************************
*
***function c_sqeoi, sqeoi
*
*OBJECT: Goes to the end of stream(information).
*
*ARGUMENTS: in iun  unit number
*
*/
void c_sqeoi(int iun) {
   int fd ;
   int i;

   if ((i=find_file_entry("c_sqeoi",iun)) < 0) return;
   
   if (FGFDT[i].attr.pipe) return;

   fd = c_getfdsc(iun);
   if (fd <= 0) return;
   lseek(fd,0,L_XTND);
}
void f77name(sqeoi)(ftnword *iun) { c_sqeoi((int) *iun) ; }

/**************************************************************************
*                     C _ S Q G E T W ,   S Q G E T W                     * 
***************************************************************************
*
***function c_sqgetw, sqgetw
*
*OBJECT: Gets a certain number of words from a stream.
*
*ARGUMENTS: in  iun      unit number
*           out bufptr   will contain the data read
*           in  nmots    number of words to read
*
*RETURNS: the number of words read if the read is successful,
*         or a negative or null number otherwise.
*
*/
int c_sqgetw(int iun, word *bufptr, int nmots) {
   int nlu, alu, alire;
   int fd ;

   alu=0; nlu=1;
   alire = nmots * sizeof(word);

   fd = c_getfdsc(iun);
   if (fd <= 0) return(fd);
   
   while (alire && (nlu > 0)) {
     nlu = read(fd,bufptr,alire);
     alire -= nlu;
     alu += nlu;
     bufptr += (nlu / sizeof(word));
   }
   return( (alire == 0) ? alu/sizeof(word) : -1);
}
ftnword f77name(sqgetw)(ftnword *iun, ftnword *bufptr, ftnword *nmots) {
   int mult = sizeof(ftnword) / sizeof(word);
   return(c_sqgetw((int) *iun, (word *) bufptr, (int) (*nmots * mult)));
}

/***************************************************************************
*                     C _ S Q P U T W ,   S Q P U T W                      *
****************************************************************************
*
***function c_sqputw, sqputw
*
*OBJECT: Writes in a stream.
*
*ARGUMENTS: in iun      unit number
*           in bufptr   what to write
*           in nmots    number of words to write
*
*RETURNS: the number of words written, if the write is successful,
*         or a negative or null number otherwise.
*
*/
int c_sqputw(int iun, word *bufptr, int nmots) {
   int necrit, aecrit, aecrire ;
   int fd ;

   aecrit = 0;
   aecrire = sizeof(word) * nmots;
   necrit = 1;
   fd = c_getfdsc(iun);
   if (fd <= 0) return(fd);

   while (aecrire && (necrit > 0)) {
     necrit = write(fd,bufptr,aecrire);
     aecrire -= necrit;
     aecrit += necrit;
     bufptr += (necrit / sizeof(word));
   }
   return( (aecrire == 0) ? necrit/sizeof(word) : -1);
}
ftnword f77name(sqputw)(ftnword *iun, ftnword *bufptr, ftnword *nmots) {
   int mult = sizeof(ftnword) / sizeof(word);
   return(c_sqputw((int) *iun, (word *) bufptr, (int) (*nmots * mult)));
}
/***************************************************************************
*                     C _ S Q G E T S ,   S Q G E T S                      *
****************************************************************************
*
***function c_sqgets, sqgets
*
*OBJECT: Gets a certain number of character from a stream.
*
*ARGUMENTS: in  iun     unit number
*           out bufptr  will contain the data read
*           in  nchar   number of char to get
*
*RETURNS: the number of characters read, if the read is successful,
*         or a negative or null number otherwise.
*
*/
int c_sqgets(int iun, char *bufptr, int nchar) {
   int nlu ;
   int fd ;
   fd = c_getfdsc(iun);
   if (fd <= 0) return(fd);
   nlu = read(fd,bufptr,nchar);
   return( (nlu > 0) ? nlu : -1);
}
ftnword f77name(sqgets)(ftnword *iun, char  *bufptr, ftnword *nchar, F2Cl llbuf) {
   int lbuf=llbuf;
   if (lbuf >= *nchar)
      return( c_sqgets(*iun, bufptr , *nchar));
   else
      return( c_sqgets(*iun, bufptr , lbuf));
}
/***************************************************************************
*                     C _ S Q P U T S ,   S Q P U T S                      *
****************************************************************************
*
***function c_sqputs, sqputs
*
*OBJECT: Writes a certain number of character to a stream.
*
*ARGUMENTS: in iun     unit number 
*           in bufptr  contains data to write
*           in nchar   number of char to write
*
*RETURNS: the number of characters written in the stream, if the read is successful
*         or a negative or null number otherwise.
*
*/
int c_sqputs(int iun, char *bufptr, int nchar) {
   int nlu ;
   int fd ;
   fd = c_getfdsc(iun);
   if (fd <= 0) return(fd);
   nlu = write(fd,bufptr,nchar);
   return( (nlu > 0) ? nlu : -1);
}
ftnword f77name(sqputs)(ftnword *iun, char  *bufptr, ftnword *nchar, F2Cl llbuf) {
   int lbuf=llbuf;
   if (lbuf >= *nchar)
      return( c_sqputs(*iun, bufptr , *nchar));
   else
      return( c_sqputs(*iun, bufptr , lbuf));
}


/****************************************************************************
*                           S C R A P _ P A G E                             *
*****************************************************************************
*
***function scrap_page
*
*OBJECT: Gets rid of the least useful page used by wafiles ind0 
*        up to (and including) ind1.
*
*ARGUMENTS: in ind0  index of first file
*           in ind1  index of last file
*
*/
static void scrap_page(int ind0,int ind1)
{
   int fl0, pg0, age0, i, j, nm, ier, found_a_page=0;;

   fl0 = ind0;
   pg0 = 0;
   age0 = 0x7fffffff;
   if (debug_mode > 4) {
      fprintf(stderr,"Debug SCRAP_PAGE ind0=%d, ind1=%d \n",ind0,ind1);
      }
/*
 * TROUVER LA PAGE LA MOINS "UTILE"
 */
   for (j = ind0; j <= ind1; j++){
      for (i = 0; i < wafile[j].nb_page_in_use; i++) {
         if (debug_mode > 5) {
            fprintf(stderr,"Debug SCRAP_PAGE j=%d, i=%d age0=%d\n",
                           j,i,wafile[j].page[i].access_count);
            }
         if (wafile[j].page[i].access_count < age0) {
            age0 = wafile[j].page[i].access_count;
            found_a_page=1;
            pg0 = i;
            fl0 = j;
            }
         }
      }
   if(found_a_page==0){
      fprintf(stderr,"scrap_page error: cannot find a page to scrap\n");
      exit(1);
      }
   if (debug_mode > 4) {
      fprintf(stderr,"Debug SCRAP_PAGE fl0=%d, pg0=%d, age0=%d\n",fl0,pg0,age0);
      }
/*
  REECRIRE LA PAGE SI CE N'EST PAS UNE PAGE READ/ONLY
*/
   if (wafile[fl0].page[pg0].touch_flag) {
      nm = wafile[fl0].page[pg0].walast - wafile[fl0].page[pg0].wa0 + 1;
      WSEEK(wafile[fl0].file_desc,wafile[fl0].page[pg0].wa0-1,L_SET);
      ier = write(wafile[fl0].file_desc,wafile[fl0].page[pg0].page_adr,sizeof(word)*nm);
      if (ier != sizeof(word)*nm) {
        fprintf(stderr,"scrap_page error: cannot write page, fd=%d\n",wafile[fl0].file_desc);
        fprintf(stderr,"scrap_page error: trying to write %d words buffer=%x, fileadr=%d\n",nm,
                             wafile[fl0].page[pg0].page_adr,wafile[fl0].page[pg0].wa0-1);
        fprintf(stderr,"scrap_page error: ier=%d,fl0=%d,ind0=%d,ind1=%d\n",ier,fl0,ind0,ind1);
        perror("FATAL WA ERROR");
        exit(1);
        }
      if (debug_mode > 4) {
         fprintf(stderr,"Debug SCRAP_PAGE ecriture disque adr=%d, nmots=%d \n",
                        wafile[fl0].page[pg0].wa0,nm);
         }
      }
   wafile[fl0].nb_page_in_use--;
   free_list[++nfree] = wafile[fl0].page[pg0].page_adr;
   if (pg0 != wafile[fl0].nb_page_in_use)
      memcpy(&wafile[fl0].page[pg0],&wafile[fl0].page[pg0+1],
             sizeof(PAGEINFO)*(wafile[fl0].nb_page_in_use - pg0));
   pg0 = wafile[fl0].nb_page_in_use;
   wafile[fl0].page[pg0].wa0 = 0;
   wafile[fl0].page[pg0].walast = 0;
   wafile[fl0].page[pg0].access_count = 0;
   wafile[fl0].page[pg0].last_access = 0;
   wafile[fl0].page[pg0].touch_flag = 0;
   }

/****************************************************************************
*                        P R O C E S S _ D E C A Y                          *
*****************************************************************************
*
***function process_decay
*
*OBJECT: Updates the age of the pages.  
*
*ARGUMENTS: none
*
*/
static void process_decay()
{
  int i,j;
  for (j = 0; j < MAXWAFILES; j++)
    for (i = 0; i < wafile[j].nb_page_in_use; i++) {
      wafile[j].page[i].access_count = decay(wafile[j].page[i].access_count);
      wafile[j].page[i].last_access++;
    }
}

/****************************************************************************
*                        G E T _ N E W _ P A G E                            *
*****************************************************************************
*
***function get_new_page
*
*OBJECT: Gets a new page for the wafile of index ind.
*
*ARGUMENTS: in ind   index of wafile
*
*/
static void get_new_page(int ind)
{
   word *new_adr; int pg0;

   if (wafile[ind].nb_page_in_use >= WA_PAGE_NB) {
      scrap_page(ind,ind);
      if (wafile[ind].nb_page_in_use >= WA_PAGE_NB) {
        fprintf(stderr,"get_new_page error: no page left !??");
        exit(1);
        }
      }

   if (nfree < 0) {
      if (global_count < WA_PAGE_LIMIT) {
         global_count++;
         free_list[++nfree] = (word *) malloc(WA_PAGE_SIZE * sizeof(word));
         if (free_list[nfree] == NULL) {
           fprintf(stderr,
           "WA get_new_page error: can't allocate (not enough memory)\n");
           exit(1);
           }
   if (debug_mode > 4) {
      fprintf(stderr,"Debug GET_NEW_PAGE nfree=%d malloc=%d \n",nfree,free_list[nfree]);
      }
         }
      }
   if (nfree < 0)
      scrap_page(0,MAXWAFILES-1);
   
   if (nfree < 0) {
      fprintf(stderr,"get_new_page error: no page left !??");
      exit(1);
      }

   pg0 = wafile[ind].nb_page_in_use++;
   
   wafile[ind].page[pg0].page_adr = free_list[nfree--];
   wafile[ind].page[pg0].wa0 = 0;
   wafile[ind].page[pg0].walast = 0;
   wafile[ind].page[pg0].access_count = 0;
   wafile[ind].page[pg0].last_access = 0;
   wafile[ind].page[pg0].touch_flag = 0;
   if (debug_mode > 4) {
      fprintf(stderr,"Debug GET_NEW_PAGE pg0=%d, page_adr=%d \n",
                     pg0,wafile[ind].page[pg0].page_adr);
      }
   }

/****************************************************************************
*                       W A _ P A G E S _ F L U S H                         *
*****************************************************************************
*
***function wa_pages_flush
*
*OBJECT: Flushes all pages of a word addressable file.
*
*ARGUMENTS: in ind   index of the wafile
*
*/
static void wa_pages_flush(int ind)

{
   while (wafile[ind].nb_page_in_use > 0) {
      scrap_page(ind,ind);
      }
}

/****************************************************************************
*                              F I L E P O S                                *
*****************************************************************************
*
***function filepos
*
*OBJECT: Returns the position of the start of the data of a subfile 
*        in a CMCARC file.
*
*ARGUMENTS: in indf  index of the subfile in the master file table
*
*RETURNS: the position of the start of the data of a subfile 
*         in a CMCARC file.
*
*/
static long long filepos(int indf)
{
  char sign[25];

  typedef struct {
    unsigned char ntotal[4];
    unsigned char ndata[4];
    char code;
    char header[MAX_NAME];
  } HEADER_CMCARC;
  
  typedef struct {
    unsigned char ntotal[8];
    unsigned char ndata[8];
    char code;
    char header[MAX_NAME];
  } HEADER_CMCARC_V5;

  HEADER_CMCARC *cmcarc_file;
  int nblu,lng,found=0,version=0,tail_offset;
  unsigned int nt,nd;
  INT_64 nt64, nd64, lng64, nblu64, pos64, retour;
  
  
  lseek(FGFDT[indf].fd,(off_t) 0,L_SET);
  nblu = read(FGFDT[indf].fd,sign,8);
  if (strncmp(sign,CMCARC_SIGN,8) != 0) {
    nblu = read(FGFDT[indf].fd,&sign[8],17);
    if (strncmp(&sign[9],CMCARC_SIGN,8) == 0) {                   /* skip to beginning of next file */
      version=4;
/*      printf("Debug+ signature version 4 trouvee\n"); */
      }
    else {
      if (strncmp(&sign[17],CMCARC_SIGN_V5,8) == 0) {
        version=5;
/*        printf("Debug+ signature version 5 trouvee\n"); */
        }
      else {
        fprintf(stderr,"%s is not a CMCARC type file\n",FGFDT[indf].file_name);
        return(-1);
        }
    }
    cmcarc_file = (HEADER_CMCARC *) &sign[0];
    nt = (cmcarc_file->ntotal[0] << 24) |
    (cmcarc_file->ntotal[1] << 16) |
    (cmcarc_file->ntotal[2] <<  8) |
    (cmcarc_file->ntotal[3]);
    
    nd = (cmcarc_file->ndata[0] << 24) |
    (cmcarc_file->ndata[1] << 16) |
    (cmcarc_file->ndata[2] <<  8) |
    (cmcarc_file->ndata[3]);
    
    if (version == 5)
      nt = nd;
    else
      if (nd != 0) {
        fprintf(stderr,
                "%s is a CMCARC file but nd=%d\n",FGFDT[indf].file_name,nd);
        return(-1);
      }
    lng = (nt *8) - 25;
    if (lseek(FGFDT[indf].fd,(off_t)lng,L_INCR) == (off_t)(-1)) {
      return (-1);
    }
  }
  subfile_length = 0;
  do {
    
    nblu = read(FGFDT[indf].fd,&cmcarc,8); /* lire nt et nd */
    if (nblu != 8) return -2;
    
    nt = (cmcarc.ntc[0] << 24) |
      (cmcarc.ntc[1] << 16) |
      (cmcarc.ntc[2] <<  8) |
      (cmcarc.ntc[3]);
    
    nd = (cmcarc.ndc[0] << 24) |
      (cmcarc.ndc[1] << 16) |
      (cmcarc.ndc[2] <<  8) |
      (cmcarc.ndc[3]);
    
    if (nt >= nd+4) {
      nt64 = nt;
      nd64 = nd;
      lng64 = (nt64 - nd64 - 2) * 8;
      tail_offset = 1;
      }
    else {
      tail_offset = 2;
      nt64 = nt;
      nt64 = (nt64 << 32) | nd;
      nblu = read(FGFDT[indf].fd,&cmcarc,8);
      nd64 = cmcarc.ntc[0];
      nd64 = (nd64 << 8) | cmcarc.ntc[1];
      nd64 = (nd64 << 8) | cmcarc.ntc[2];
      nd64 = (nd64 << 8) | cmcarc.ntc[3];
      nd64 = (nd64 << 8) | cmcarc.ndc[0];
      nd64 = (nd64 << 8) | cmcarc.ndc[1];
      nd64 = (nd64 << 8) | cmcarc.ndc[2];
      nd64 = (nd64 << 8) | cmcarc.ndc[3];
      lng64 = (nt64 - nd64 - 4) * 8;
      if (nt64 < nd64+6) {
        fprintf(stderr,
                "%s is a CMCARC file but nt=%d nd=%d\n",FGFDT[indf].file_name,nt64,nd64);
        return(-1);
        }
      }
/*    printf("Debug+ nt64=%Ld nd64=%Ld lng64=%Ld\n",nt64,nd64,lng64); */
    nblu64 = read(FGFDT[indf].fd,cmcarc.cmcarc_name,lng64);
/*    printf("Debug cmcarc.cmcarc_name=%s\n",&cmcarc.cmcarc_name[1]); */
    if (nblu64 != lng64) return -3;
    if (strcmp(FGFDT[indf].subname,&cmcarc.cmcarc_name[1]) == 0) {
      found = 1;
    }
    else {              /* sauter les donnees */
      lng64 = (nd64+tail_offset) * 8;
      if (lseek64(FGFDT[indf].fd,(off_t)lng64,L_INCR) == (off_t)(-1)) {
        return (-1);
      }
    }
  } while(!found);
  subfile_length = (nd*8)/sizeof(word);
  pos64=tell64(FGFDT[indf].fd);
  retour = pos64/sizeof(word);
  return(retour);
/*  return((tell(FGFDT[indf].fd))/sizeof(word)); */
}


/****************************************************************************
*                              Q Q C O P E N                                *
*****************************************************************************
*
***function qqcopen
*
*OBJECT: Opens a non-fortran file.
*        Active part of c_waopen2.
*
*ARGUMENTS: in indf  index of the file in the master file table 
*
*NOTE: Processes files contained in a CMCARC file.
*
*RETURNS: the file descriptor, if the open is successful,
*         or a negative number otherwise.
*
*/
static int qqcopen(int indf)
{
  struct stat etat;
  int fd, mode;
  char name[MAX_NAME], *pt;
  int ind, i, nset, n1, n2, n3, n4;
  char *ptsz, *ptnb, *ptlm, *ptenv, *getenv();
  int pagesz, pagenb, atoi();
  LLSK dim;
  int iii;
  char *errmsg="";
  
  /*    beginning of INITIALIZATION section    */
if (! init) {
  ptenv = getenv("WA_CONFIG");
  nset = 0;
  if (ptenv != NULL)
    nset = sscanf(ptenv,"%d %d %d %d",&n1,&n2,&n3,&n4);
  
  switch (nset)
    {
    case 4:
      debug_mode = n4;
      
    case 3:
      WA_PAGE_LIMIT = n3;
      
    case 2:
      WA_PAGE_NB = n2;
      
    case 1:
      WA_PAGE_SIZE = n1 * 1024 * (sizeof(ftnword) / sizeof(word));
      break;
      
    default:
      WA_PAGE_SIZE = 0;
      break;
      
    }
  
  
  WA_PAGE_NB   = (WA_PAGE_NB < MAXPAGES) ? WA_PAGE_NB : MAXPAGES;
  
  if (WA_PAGE_LIMIT == 0)
    WA_PAGE_LIMIT = WA_PAGE_NB * MAXWAFILES;
  if (WA_PAGE_SIZE > 0) {
    fprintf(stderr,"WA_PAGE_SZ = %d Bytes ",WA_PAGE_SIZE*sizeof(word));
    fprintf(stderr,"WA_PAGE_NB = %d ",WA_PAGE_NB);
    fprintf(stderr,"WA_PAGE_LIMIT = %d\n",WA_PAGE_LIMIT);
  }
  for (ind = 0; ind < MAXWAFILES; ind++) {
    wafile[ind].file_desc = -1;
    wafile[ind].nb_page_in_use = 0;
    wafile[ind].offset = 0;
  }
  init = 1;
}
/*    end of INITIALIZATION section    */
if (FGFDT[indf].attr.remote) return(0);  /* file will be open by fnom_rem_connect */

FGFDT[indf].fd = -1;
ind = 0;
while ((wafile[ind].file_desc != -1) && (ind < MAXWAFILES))
ind++;
if (ind == MAXWAFILES) {
  fprintf(stderr,"qqcopen error: too many open files\n");
  return(-1);
}

if (FGFDT[indf].subname) {    /* fichier de type cmcarc */
  if (debug_mode > 4) {
    fprintf(stderr,"Debug opening subfile %s from file %s\n",
            FGFDT[indf].subname,FGFDT[indf].file_name);
  }
  FGFDT[indf].attr.read_only = 1;
  mode = O_RDONLY;
  if ((fd = open64(FGFDT[indf].file_name,mode | WIN32_O_BINARY)) == -1) {
    fprintf(stderr,"qqcopen error: cannot open file %s\n",FGFDT[indf].file_name);
    return(-1);
  }
  wafile[ind].file_desc = fd;
  FGFDT[indf].fd = fd;
  if ((wafile[ind].offset = filepos(indf)) <= 0) {
    fprintf(stderr,"qqcopen error: subfile %s not found in %s\n",
            FGFDT[indf].subname,FGFDT[indf].file_name);
    return(-1);
  }
  FGFDT[indf].open_flag = 1;
  if (debug_mode > 4) {
    fprintf(stderr,"Debug subfile found at position %u\n",wafile[ind].offset);
  }
}

else {  /* not a CMCARC type file */
  if (access(FGFDT[indf].file_name, F_OK) == -1)
    {
      if (errno == ENOENT)     /* nouveau fichier, creation */
        {
          fd = open64(FGFDT[indf].file_name, O_RDWR | O_CREAT | WIN32_O_BINARY,
                    S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
          FGFDT[indf].attr.read_only = 0;
          errmsg="cannot create file";
        }
    }
  else
    if (! FGFDT[indf].attr.read_only)      /* tentative d'ouverture en mode R/W */
      {
        mode = O_RDWR;
        fd = open64(FGFDT[indf].file_name, mode | WIN32_O_BINARY);
        if (fd == -1) {
          if (!FGFDT[indf].attr.write_mode)
            {
              mode = O_RDONLY;
              FGFDT[indf].attr.read_only = 1;
              fd = open64(FGFDT[indf].file_name, mode | WIN32_O_BINARY);
              errmsg="cannot open file";
            }
          else                  /* ouverture demande en mode R/W */
            errmsg="cannot open in write mode";
        }
      }
    else if (FGFDT[indf].attr.read_only)  /* ouverture en mode R/O */
      {
        mode = O_RDONLY;
        fd = open64(FGFDT[indf].file_name, mode | WIN32_O_BINARY);
        errmsg="cannot open file";
      }
  if (fd == -1)
    {
      fprintf(stderr, "qqcopen error: %s filename=(%s) !\n",errmsg,FGFDT[indf].file_name);
      return(-1);
    }
  wafile[ind].file_desc = fd;
  FGFDT[indf].fd = fd;
  FGFDT[indf].open_flag = 1;
}

dim = 0;
dim = LSEEK(fd, dim, L_XTND);
FGFDT[indf].file_size = dim / sizeof(word);
FGFDT[indf].eff_file_size = dim / sizeof(word);
dim = 0;
dim = LSEEK(fd, dim, L_SET);
if (subfile_length > 0) 
FGFDT[indf].eff_file_size = subfile_length;
subfile_length = 0;

if (WA_PAGE_SIZE != 0) {
  for (i = 0; i < WA_PAGE_NB; i++) {
    wafile[ind].page[i].page_adr = NULL;
    wafile[ind].page[i].wa0 = 0;
    wafile[ind].page[i].walast = 0;
    wafile[ind].page[i].access_count = 0;
    wafile[ind].page[i].last_access = 0;
    wafile[ind].page[i].touch_flag = 0;
  }
  wafile[ind].nb_page_in_use = 0;
  if (debug_mode > 1) {
    fprintf(stderr,"Debug ouverture du fichier %s ind=%d, fd=%d\n",
            FGFDT[indf].file_name,ind,fd);
#if defined (NEC)
    fprintf(stderr,"Debug longueur du fichier =%lld Bytes\n",dim);
#else
    fprintf(stderr,"Debug longueur du fichier =%d Bytes\n",dim);
#endif
  }      
    }
return(fd);
}


/****************************************************************************
*                               D _ W A F D T                               *
*****************************************************************************
*
***function d_wafdt
*
*OBJECT: Prints the wa control table.
*
*ARGUMENTS: none.
*
*/


void f77name(d_wafdt)()
{
int i;
fprintf(stderr,"\n   DUMP OF WA CONTROL TABLE \n");
for (i=0;i<MAXWAFILES;i++){
   if(wafile[i].file_desc != -1)
   fprintf(stderr,"waindex=%d, fd=%d, npages=%d, offset=%Ld\n",
    i,wafile[i].file_desc,wafile[i].nb_page_in_use,wafile[i].offset);
   }
}

/****************************************************************************
*                         W A _ P A G E _ R E A D                           *
*****************************************************************************
*
***function wa_page_read
*
*OBJECT: Reads a word addressable file page.
*
*ARGUMENTS: in  fd     file descriptor
*           out buf    will contain the data read
*           in  adr    word to start at
*           in  nmots  number of words to read
*           in  indf   index of wafile in the master file table
*
*/
static void wa_page_read(int fd,word *buf,unsigned int adr,int nmots,int indf)

{
   int ind, j, wa0, offset, i=0, found=0, nbytes;
   word lnmots, zero=0;
   word readbytes;

   process_decay();
/*
  TROUVER ind DU FICHIER
*/
   ind = 0;
   while ((wafile[ind].file_desc != fd) && (ind < MAXWAFILES))
      ind++;
   if (debug_mode > 4) {
      fprintf(stderr,"Debug WA_PAGE_READ requete adr=%u, nmots=%d ind=%d\n",adr,nmots,ind); 
      }
   if (ind == MAXWAFILES) {
      fprintf(stderr,"wa_page_read error: file is not open\n");
      exit(1);
      }

/*
  TROUVER LA PAGE i QUI CONTIENT [adr,adr+nmots]
*/
   found = 0;
   while((! found) && (i < wafile[ind].nb_page_in_use))
      if((adr >= wafile[ind].page[i].wa0) && (adr+nmots <= wafile[ind].page[i].wa0+WA_PAGE_SIZE))
         found = 1;
      else
         i++;

    if (! found) {
/*
  OBTENIR UNE NOUVELLE PAGE
*/
       get_new_page(ind);
       i = wafile[ind].nb_page_in_use -1;
       wa0 = adr - (adr % WA_PAGE_SIZE) + 1;
       wafile[ind].page[i].wa0 = (wa0 > 1) ? wa0 : 1;
       wafile[ind].page[i].access_count++;
       if (debug_mode > 4) {
          fprintf(stderr,"Debug WA_PAGE_READ obtention d'une page %d\n",i);
          }
       WSEEK(fd,wafile[ind].page[i].wa0-1,L_SET);
       if (WA_PAGE_SIZE+wafile[ind].page[i].wa0 > FGFDT[indf].file_size)
          readbytes =  sizeof(word)*(FGFDT[indf].file_size+1-wafile[ind].page[i].wa0);
       else
          readbytes = sizeof(word)*WA_PAGE_SIZE;
       nbytes = read(fd,wafile[ind].page[i].page_adr,sizeof(word)*WA_PAGE_SIZE);
       if ( nbytes < readbytes ) {
          fprintf(stderr,"wa_page_read error: cannot read page from file %d,fd=%d\n",ind,fd);
          fprintf(stderr,"  tried to get %d bytes, got %d\n",sizeof(word)*WA_PAGE_SIZE,nbytes);
          perror("WA_PAGE_READ");
          exit(1);
          }
       if (nbytes < sizeof(word)*WA_PAGE_SIZE) {
         lnmots = WA_PAGE_SIZE - (nbytes/sizeof(word));
/*
         f77name(afix)(wafile[ind].page[i].page_adr+(nbytes/sizeof(word)),&zero,&lnmots);
         { word *move_dest=wafile[ind].page[i].page_adr+(nbytes/sizeof(word)),
           nwords=lnmots;
           while(lnmots--) *move_dest++ = 0;
         }
*/
         ZERO(wafile[ind].page[i].page_adr+(nbytes/sizeof(word)),lnmots);
         }
       wafile[ind].page[i].walast = wafile[ind].page[i].wa0 + nbytes / sizeof(word) -1;
       if (debug_mode > 4) {
          fprintf(stderr,"Debug WA_PAGE_READ lecture disque adr=%d\n",wafile[ind].page[i].wa0);
          }
       }

/*
  COPIER DANS BUF LA SECTION DESIREE
*/
   offset = adr - wafile[ind].page[i].wa0;
   lnmots = nmots;
/*
   { word *move_src=wafile[ind].page[i].page_adr+offset,
          *move_dest=buf , nwords=lnmots;
     while(lnmots--) *move_dest++ = *move_src++;
   }
   f77name(movlev)(wafile[ind].page[i].page_adr+offset,buf,&lnmots); 
*/
   MOVE(wafile[ind].page[i].page_adr+offset,buf,lnmots); 
   wafile[ind].page[i].last_access = 0;
   wafile[ind].page[i].access_count = new_age_rd(wafile[ind].page[i].access_count);

/*
  VERIFICATION DES CROISEMENTS POSSIBLES DES PAGES
*/
   for (j=0; j< wafile[ind].nb_page_in_use; j++)
      for (i=0; i < wafile[ind].nb_page_in_use; i++)
         if (j != i)
            if ((wafile[ind].page[j].wa0 >= wafile[ind].page[i].wa0) &&
                (wafile[ind].page[j].wa0 <= wafile[ind].page[i].wa0+WA_PAGE_SIZE-1)) {
               fprintf(stderr,"WA_PAGE_READ error: overlapping pages i=%d\n",i);
               fprintf(stderr,"page[j].wa0 =%d, page[i].wa0 =%d, page[i].wa0+WA_PAGE_SIZE =%d\n",
                              wafile[ind].page[j].wa0,wafile[ind].page[i].wa0,
                              wafile[ind].page[i].wa0+WA_PAGE_SIZE-1);
               f77name(tracebck)();
               exit(1);
               }
   }

/*****************************************************************************
*                      H R J U S T                                           *
******************************************************************************
*
***function hrjust
*
*OBJECT: Right justifies a group of ncar characters(8 bits each).
*
*ARGUMENTS: in moth   word to justify
*           in ncar number of characters
*
*RETURNS: the group of characters right justified.
*
*/
unsigned ftnword f77name(hrjust) (unsigned ftnword *moth, ftnword *ncar)
{
   int sc;
   sc = 8 * ( sizeof(ftnword) - *ncar );
   return (sc<=0 ? *moth : (*moth) >> sc);
}

/*****************************************************************************
*                      H L J U S T                                           *
******************************************************************************
*
***function hljust
*
*OBJECT: Left justifies a group of ncar characters(8 bits each).
*
*ARGUMENTS: in moth   word to justify 
*           in ncar   number of characters
*
*RETURNS: the group of characters left justified.
*
*/
unsigned ftnword f77name(hljust) (unsigned ftnword *moth, ftnword *ncar)
{
   int sc;
   sc = 8 * ( sizeof(ftnword) - *ncar );
   return (sc<=0 ? *moth : (*moth) << sc);
}

/****************************************************************************
*                       W A _ P A G E _ W R I T E                           *
*****************************************************************************
*
***function wa_page_write
*
*OBJECT: Writes on a wa page.
*
*ARGUMENTS: in fd    file descriptor
*           in buf   contains data to write
*           in adr   where to start writing
*           in nmots number of words to write
*           in indf  index of the wafile in the master file table
*
*/
static void wa_page_write(int fd,word *buf,unsigned int adr,int nmots,int indf)

{
   int ind, j, wa0, offset, i=0, found=0, nbytes, readbytes;
   word lnmots, zero=0;

   process_decay();
/*
  TROUVER ind DU FICHIER
*/
   ind = 0;
   while ((wafile[ind].file_desc != fd) && (ind < MAXWAFILES))
      ind++;
   if (debug_mode > 4) {
      fprintf(stderr,"Debug WA_PAGE_WRITE requete adr=%u, nmots=%d ind=%d\n",adr,nmots,ind); 
      }
   if (ind == MAXWAFILES) {
      fprintf(stderr,"wa_page_write error: file is not open\n");
      exit(1);
      }
/*
  TROUVER LA PAGE i QUI CONTIENT [adr,adr+nmots]
*/
   found = 0;
   while((! found) && (i < wafile[ind].nb_page_in_use))
      if((adr >= wafile[ind].page[i].wa0) && (adr+nmots <= wafile[ind].page[i].wa0+WA_PAGE_SIZE))
         found = 1;
      else
         i++;

   if (! found) {
/*
  OBTENIR UNE NOUVELLE PAGE
*/
      get_new_page(ind);
      i = wafile[ind].nb_page_in_use -1;
      if (debug_mode > 4) {
         fprintf(stderr,"Debug WA_PAGE_WRITE obtention d'une page %d\n",i); 
         }
      wa0 = adr - (adr % WA_PAGE_SIZE) + 1;
      wafile[ind].page[i].wa0 = (wa0 > 1) ? wa0 : 1;
      wafile[ind].page[i].access_count++;
/*
  SI LA REQUETE D'ECRITURE NE COUVRE PAS LA PAGE EN ENTIER, OU
  SI LA REQUETE NE PART PAS DU DEBUT DE LA PAGE POUR SE TERMINER
  A LA FIN DE FICHIER OU PLUS, ALORS IL FAUT RELIRE LA PAGE EN QUESTION
*/
     if ((adr > wafile[ind].page[i].wa0) || 
          ((adr+nmots != wafile[ind].page[i].wa0+WA_PAGE_SIZE) && 
           (adr+nmots < FGFDT[indf].file_size))) {
         WSEEK(fd,wafile[ind].page[i].wa0-1,L_SET);
         if (WA_PAGE_SIZE+wafile[ind].page[i].wa0 > FGFDT[indf].file_size)
            readbytes =  sizeof(word)*(FGFDT[indf].file_size+1-wafile[ind].page[i].wa0);
         else
            readbytes = sizeof(word)*WA_PAGE_SIZE;
         nbytes = read(fd,wafile[ind].page[i].page_adr,readbytes);
         if ( nbytes < readbytes ) {
            fprintf(stderr,
                    "wa_page_write error: cannot read page on file %s\n",
                    FGFDT[indf].file_name);
            fprintf(stderr," tried to get %d bytes, got %d\n",
                    readbytes,nbytes);
            fprintf(stderr,"WA_PAGE_SIZE=%d wa0=%d file_size=%d\n",
                    WA_PAGE_SIZE,wafile[ind].page[i].wa0 ,
                    FGFDT[indf].file_size);
            perror("WA_PAGE_WRITE");
            exit(1);
            }
         if (nbytes < sizeof(word)*WA_PAGE_SIZE) {
           lnmots = WA_PAGE_SIZE - (nbytes/sizeof(word));
/*
           { word *move_dest=wafile[ind].page[i].page_adr+(nbytes/sizeof(word)),
             nwords=lnmots;
             while(lnmots--) *move_dest++ = 0;
           }
*/
           ZERO(wafile[ind].page[i].page_adr+(nbytes/sizeof(word)),lnmots);
           }
         wafile[ind].page[i].walast = wafile[ind].page[i].wa0 + nbytes / sizeof(word) -1;
         if (debug_mode > 4) {
            fprintf(stderr,"Debug WA_PAGE_WRITE relecture disque de la page %d a l'adresse %d\n",
                           i,wafile[ind].page[i].wa0); 
            }
         }
      }
/*
  COPIER DANS LA PAGE LA SECTION buf
*/
   offset = adr - wafile[ind].page[i].wa0;
   lnmots = nmots;
/*
   { word *move_dest=wafile[ind].page[i].page_adr+offset,
          *move_src=buf , nwords=lnmots;
     while(lnmots--) *move_dest++ = *move_src++;
   }
   f77name(movlev)(buf,wafile[ind].page[i].page_adr+offset,&lnmots);
*/
   MOVE(buf,wafile[ind].page[i].page_adr+offset,lnmots);
   wafile[ind].page[i].last_access = 0;
   wafile[ind].page[i].access_count = new_age_wr(wafile[ind].page[i].access_count);
   wafile[ind].page[i].touch_flag = 1;
   if (adr > FGFDT[indf].file_size+1) {
      fprintf(stderr,"WA_PAGE_WRITE error: ind = %d, adr=%u > file_size=%d\n",
                     ind,adr,FGFDT[indf].file_size);
      fprintf(stderr,"                     filename=%s\n",FGFDT[indf].file_name);
      exit(1);
      }

   if (adr+nmots-1 > FGFDT[indf].file_size) {
      FGFDT[indf].file_size = adr+nmots-1;
      FGFDT[indf].eff_file_size = adr+nmots-1;
      }
   if (wafile[ind].page[i].walast < adr+nmots-1)
      wafile[ind].page[i].walast = adr+nmots-1;

   for (j=0; j< wafile[ind].nb_page_in_use; j++)
      for (i=0; i < wafile[ind].nb_page_in_use; i++)
         if (j != i)
            if ((wafile[ind].page[j].wa0 >= wafile[ind].page[i].wa0) && 
                (wafile[ind].page[j].wa0 <= wafile[ind].page[i].wa0+WA_PAGE_SIZE-1)) {
               fprintf(stderr,"WA_PAGE_WRITE error: overlapping pages i=%d\n",i);
               fprintf(stderr,"page[j].wa0 =%d, page[i].wa0 =%d, page[i].wa0+WA_PAGE_SIZE =%d\n",
                              wafile[ind].page[j].wa0,wafile[ind].page[i].wa0,
                              wafile[ind].page[i].wa0+WA_PAGE_SIZE-1);
               f77name(tracebck)();
               exit(1);
               }
      }


/****************************************************************************
*                              Q Q C W A W R                                *
*****************************************************************************
*
***function qqcwawr
*
*OBJECT: Writes in a word adressable file. 
*        Active part of c_wawrit2.
*
*ARGUMENTS: in lfd    file descriptor
*           in buf    contains data to write
*           in wadr   file address in words 
*           in nmots  number of words to write
*           in indf   index in the master file table
*
*/

static void qqcwawr(word *buf,unsigned int wadr,int lnmots,int indf)
{

int offset,i,adr0,nwritten,togo;
int lng, l, lastadr, ind, statut;
int lfd=FGFDT[indf].fd;
long long ladr=wadr;
char *cbuf;

ind = 0;
while ((wafile[ind].file_desc != lfd) && (ind < MAXWAFILES))
   ind++;
if (ind == MAXWAFILES) {
  fprintf(stderr,"qqcwawr error: filename=%s , fd=%d not found in table\n",
                 FGFDT[indf].file_name,lfd);
  exit(1);
  }

if (ladr != 0) 
  ladr += wafile[ind].offset;

if (FGFDT[indf].attr.read_only) {
  fprintf(stderr,"qqcwawr error: no write permission for file %s\n",FGFDT[indf].file_name);
  exit(1);
  }
  
if (FGFDT[indf].attr.remote) {
  int *s_ID, *addr, *nw, *RW_mode, *checksum;
  int sock_comm_ID=0xBABE;
  int demande[5];
  int nc, nelm;

  s_ID = &(demande[0]);
  addr = &(demande[1]);
  nw = &(demande[2]);
  RW_mode = &(demande[3]);
  checksum = &(demande[4]);
  *s_ID = sock_comm_ID;
  *addr = ladr;
  *nw = lnmots;
  *RW_mode = 2;  /* write request */
  *checksum = *s_ID ^ *addr ^ *nw ^ *RW_mode;
  check_swap_records(demande,5,sizeof(int));
  nc=write_stream(FGFDT[indf].fd,demande,5*sizeof(int));
  if (nc != 0) {
    fprintf(stderr,"socket qqcwawr error: wrote only %i bytes to server\n",nc);
    fflush(stderr);
    }
  nelm=write_stream(FGFDT[indf].fd,buf,lnmots*sizeof(int));
#if defined (DEBUG)  
  if (nelm == 0) printf("socket qqcwawr wrote %d bytes\n",lnmots*sizeof(int));
#endif
  if (ladr+lnmots-1 > FGFDT[indf].file_size) {
    FGFDT[indf].file_size = ladr+lnmots-1;
    FGFDT[indf].eff_file_size = ladr+lnmots-1;
    }
} /* end remote */
else {
  
  if ((WA_PAGE_SIZE == 0) || (ladr == 0)) {
    if(ladr!=0) WSEEK(lfd,ladr - 1, L_SET);
    if ((nwritten=write(lfd, buf, sizeof(word) * lnmots)) != sizeof(word) * lnmots)
        {
          if (errno == 14)
            {
              fprintf(stderr, "qqcwawr error: write error for file %s\n",FGFDT[indf].file_name);
              fprintf(stderr,"qqcwawr: filename=%s, buf=%0x adr=%u, nmots=%d, nwritten=%d, errno=%d\n",
                      FGFDT[indf].file_name,buf,ladr,lnmots,nwritten,errno);
              fprintf(stderr, "*** Contactez un membre de la section informatique de RPN ***\n");
              fprintf(stderr, "*** Seek support from RPN informatic section ***\n");
              /*            memorymap(1); */
              perror("qqcwawr");
              exit(1);
            }
          if (nwritten >= 0) {
            cbuf = (char *) buf;
            cbuf += nwritten;
            togo = (lnmots * sizeof(word)) - nwritten;
            nwritten = write(lfd,buf,togo);
            fprintf(stderr,"qqcwawr WARNING: multiple write attempt of file %s last write=%d bytes, total needed=%d bytes\n",
                    FGFDT[indf].file_name,togo,lnmots*sizeof(word));
            if (nwritten != togo) {
              fprintf(stderr, "qqcwawr error: write error for file %s\n",FGFDT[indf].file_name);
              fprintf(stderr,"qqcwawr: filename=%s, buf=%0x adr=%u, nmots=%d, nwritten=%d, errno=%d\n",
                      FGFDT[indf].file_name,buf,ladr,lnmots,nwritten,errno);
              perror("qqcwawr");
              exit(1);
            }
          }
          else {
            fprintf(stderr, "qqcwawr error: write error or file not open for write!\n");
            fprintf(stderr,"qqcwawr: filename=%s, buf=%0x adr=%u, nmots=%d, nwritten=%d, errno=%d\n",
                    FGFDT[indf].file_name,buf,ladr,lnmots,nwritten,errno);
            perror("qqcwawr");
            exit(1);
          }
        }
    if (ladr+lnmots-1 > FGFDT[indf].file_size) {
        FGFDT[indf].file_size = ladr+lnmots-1;
        FGFDT[indf].eff_file_size = ladr+lnmots-1;
        }
    }
  else {
    lng = lnmots;
    adr0 = ladr;
    offset = 0;
    lastadr = (adr0 + WA_PAGE_SIZE -1)/WA_PAGE_SIZE * WA_PAGE_SIZE;
    while (lng > 0) {
        if (lng > lastadr-adr0) {
          l = lastadr - adr0 +1;
          wa_page_write(lfd,buf+offset,adr0,l,indf);
          offset = offset + l;
          adr0 = adr0 +l;
          lng = lng -l;
          lastadr = (adr0 + WA_PAGE_SIZE -1)/WA_PAGE_SIZE * WA_PAGE_SIZE;
          }
        else {
          wa_page_write(lfd,buf+offset,adr0,lng,indf);
          lng = 0;
          }
        }
    }
} /* end else remote */
}

/****************************************************************************
*                              Q Q C W A R D                                *
*****************************************************************************
*
***function qqcward
*
*OBJECT: Reads a word addressable file.
*        Active part of c_waread2.
*
*ARGUMENTS: in  lfd     file descriptor
*           out buf     will contain data read
*           in  wadr    file address in words
*           in  lnmots  number of words to read
*           in  indf    index of file in the master file table
*
*/
static void qqcward(word *buf,unsigned int wadr,int  lnmots,int indf)
{
int offset,i,wa0,adr0,lng,l,lastadr;
int npages,reste,ind;
int lfd=FGFDT[indf].fd;
long long ladr=wadr;

ind = 0;
while ((wafile[ind].file_desc != lfd) && (ind < MAXWAFILES))
   ind++;
if (ind == MAXWAFILES) {
  fprintf(stderr,"qqcward error: fd=%d not found in table\n",lfd);
  exit(1);
  }
if (FGFDT[indf].attr.remote) {
  int *s_ID, *addr, *nw, *RW_mode, *checksum;
  int sock_comm_ID=0xBABE;
  int demande[5];
  int nc, nelm;

  s_ID = &(demande[0]);
  addr = &(demande[1]);
  nw = &(demande[2]);
  RW_mode = &(demande[3]);
  checksum = &(demande[4]);
  *s_ID = sock_comm_ID;
  *addr = ladr;
  *nw = lnmots;
  *RW_mode = 1;  /* read request */
  *checksum = *s_ID ^ *addr ^ *nw ^ *RW_mode;
  check_swap_records(demande,5,sizeof(int));
  nc=write_stream(FGFDT[indf].fd,demande,5*sizeof(int));
  if (nc != 0) {
    fprintf(stderr,"socket qqcward error: wrote only %d bytes to server\n",nc);
    fflush(stderr);
    }
  nelm=read_stream(FGFDT[indf].fd,buf,lnmots*sizeof(int));
#if defined (DEBUG)
  printf("qqcward read %d bytes\n",nelm);
#endif
} /* end remote */
else {
  if (ladr != 0) 
    ladr += wafile[ind].offset;

  if ((WA_PAGE_SIZE == 0) || (ladr == 0)) {
    if(ladr!=0) WSEEK(lfd, ladr - 1, L_SET);
    reste=read(lfd, buf, sizeof(word) * lnmots);
    if(reste != sizeof(word)*lnmots) {
        fprintf(stderr,"qqcward error: tried to read %d words, only read %d\n",
                      sizeof(word)*lnmots,reste);
        fprintf(stderr,"qqcward: wafile[ind].offset=%d ladr=%Ld\n",wafile[ind].offset,ladr);
        f77name(tracebck)();
        exit(1);
    }
    }
  else {
    lng = lnmots;
    adr0 = ladr;
    offset = 0;
    lastadr = (adr0 + WA_PAGE_SIZE -1)/WA_PAGE_SIZE * WA_PAGE_SIZE;
    while (lng > 0) {
        if (lng > lastadr-adr0) {
          l = lastadr - adr0 +1;
          wa_page_read(lfd,buf+offset,adr0,l,indf);
          offset = offset + l;
          adr0 = adr0 +l;
          lng = lng -l;
          lastadr = (adr0 + WA_PAGE_SIZE -1)/WA_PAGE_SIZE * WA_PAGE_SIZE;
          }
        else {
          wa_page_read(lfd,buf+offset,adr0,lng,indf);
          lng = 0;
          }
    } /* end while */
  } /* end else */
} /* end else */
}
/****************************************************************************
*                              fnom_rem_connect                             *
*****************************************************************************
*
***function fnom_rem_connect
*
*OBJECT: Establish socket connection with host server
*
*ARGUMENTS: in remote_host remote host server name
*           in ind         file index in fnom general file table index
*           
*Revisions
*           nov 2009 - Utilisation de ssh dans tous les cas pour demarrer wa_server
*           fev 2011 - Utilisation de system(r.remote_wa_server) au lieu de popen
*
*/
int fnom_rem_connect(int ind, char* remote_host)
{
  char buf[1024];
  char cbuf[1024];
  char pbuf[1024];
  char remote_command[1024];
  int wa_buf_lu[10];
  int fserver;
  int fclient=-1;
  int server_port = -1;
  FILE *comm = NULL, *file_ptr=NULL;
  fd_set rfds;
  fd_set wfds;
  fd_set efds;
  struct timeval tv;
  int nc, i, indx, isel, ier;
  int *s_ID, *addr, *nw, *RW_mode, *checksum, new_checksum;
  int sock_comm_ID=0xBABE;
  int demande[5];
  
  fserver = bind_to_localport(&server_port, cbuf, sizeof(cbuf)-1);
  listen(fserver,5);
  printf("bound to #%s#\n", cbuf);
  fflush(stdout);

/*
  file_ptr = fopen("ECssm/all/bin/remote_exec.sh","r");

  if (file_ptr != NULL) {
    printf("Debug+ passe par remote_exec\n");
    snprintf(pbuf,sizeof(pbuf)-1,"ssh %s -n %s %s %s @%s ",remote_host,"ECssm/all/bin/remote_exec.sh wa_server",
           FGFDT[ind].file_name,(FGFDT[ind].attr.read_only == 1) ? "R/O" : "R/W",cbuf);
    printf("Debug+ commande passee =\n%s\n",pbuf);
    }
  else {
    snprintf(pbuf,sizeof(pbuf)-1,"ssh %s -n %s/bin/%s %s %s @%s ",remote_host,armnlibpath,"wa_server",
           FGFDT[ind].file_name,(FGFDT[ind].attr.read_only == 1) ? "R/O" : "R/W",cbuf);
    printf("Debug+ ancienne commande passee =\n%s\n",pbuf);
    }
*/
/*
  snprintf(pbuf,sizeof(pbuf)-1,"rsh %s -n %s %s %s @%s ",remote_host,"/users/dor/armn/mlp/tests/SOCKETS/wa_server",
           FGFDT[ind].file_name,(FGFDT[ind].attr.read_only == 1) ? "R/O" : "R/W",cbuf);
*/
                            /*  rsh host -n commande fichier mode socket */
  snprintf(pbuf,sizeof(pbuf)-1,"echo wa_server %s %s @%s | ssh %s 'bash --login 1>/dev/null 2>/dev/null'",
           FGFDT[ind].file_name,(FGFDT[ind].attr.read_only == 1) ? "R/O" : "R/W",cbuf,remote_host);

/*  comm = popen(pbuf,"r"); */
/*  ier = system(pbuf); */
  
  snprintf(remote_command,sizeof(remote_command)-1,"r.remote_wa_server %s %s %s %s",
	  FGFDT[ind].file_name,(FGFDT[ind].attr.read_only == 1) ? "R/O" : "R/W",cbuf,remote_host);
  printf("Debug+ commande passee =\n%s\n",remote_command);
  ier = system(remote_command);
/*  printf("Debug+ ier de remote_command =%d\n",ier); */
/*  
  if(comm == NULL) {
    printf("fnom_rem_connect error: popen error !!\n");
    return(-1);
  }else{
    int nelm;
    buf[0] = '\0';
    nelm = fread(buf, 1, sizeof(buf)-1, comm);
    if(nelm > 0) buf[nelm] ='\0';
    pclose(comm);
    printf("read %d bytes, ::%s::\n", nelm, buf);
  }
*/

  fflush(stdout);
  FD_ZERO(&rfds);
  FD_SET(fserver, &rfds);
  tv.tv_sec = 5;
  tv.tv_usec = 0;
  isel = select(fserver+1, &rfds, NULL, NULL, &tv);
  /* if (select(fserver+1, &rfds, NULL, NULL, &tv)) */
/*  printf("Debug+ retour de select=%d\n",isel); */
  if (isel)
  {
    fclient = accept_from_sock(fserver);
    printf("connected to server\n");
    fflush(stdout);
    FGFDT[ind].fd = -1;
    indx = 0;
    while ((wafile[indx].file_desc != -1) && (indx < MAXWAFILES))
      indx++;
    if (indx == MAXWAFILES) {
      fprintf(stderr, "fnom_rem_connect error: too many open files\n");
      return(-1);
    }
   
    s_ID = &(demande[0]);
    addr = &(demande[1]);
    nw = &(demande[2]);
    RW_mode = &(demande[3]);
    checksum = &(demande[4]);
    *s_ID = sock_comm_ID;
    *addr = 0;
    *nw = 0;
    *RW_mode = 4;  /* wasize request */
    *checksum = *s_ID ^ *addr ^ *RW_mode;
    check_swap_records(demande, 5, sizeof(int));
    nc = write_stream(fclient, demande, 5*sizeof(int));

    if (nc == 0) {
      printf("fnom_rem_connect wrote to server OK\n");
      fflush(stdout);
      }
    else {
      fprintf(stderr, "fnom_rem_connect error: wrote only %d bytes to server\n", nc);
      fflush(stderr);
      close(fclient);
      return(-1);
      }

    demande[0] = 0; demande[1] = 0; demande[2] = 0; demande[3] = 0; demande[4] = 0;
    nc = read_stream(fclient, demande, 5*sizeof(int));
    if (nc !=  5*sizeof(int)) {
      fprintf(stderr, "fnom_rem_connect error: read only %d bytes from server\n", nc);
      fflush(stderr);
      close(fclient);
      return(-1);
      }

    check_swap_records(demande, 5, sizeof(int));
    new_checksum = *s_ID ^ *addr ^ *RW_mode;
    if (new_checksum != *checksum) {
      fprintf(stderr,"fnom_rem_connect error: invalid checksum=%X not %X\n",new_checksum,checksum);
      fflush(stderr);
      close(fclient);
      return(-1);
      }
      printf("Debug+ fnom_rem_connect wasize=%d\n",*nw);
      fflush(stdout);
    FGFDT[ind].file_size = *nw;
    FGFDT[ind].eff_file_size = *nw;
  }
  else {
    fprintf(stderr,"fnom_rem_connect error: cannot connect to server\n");
    return(-1);
  }

  wafile[indx].file_desc = fclient;
  FGFDT[ind].fd = fclient;
  FGFDT[ind].open_flag = 1;
  return(0);
}

/****************************************************************************
*                                  MOVE                                     *
*****************************************************************************
*
***function move
*
*OBJECT: Moves nwords words from src to dest.
*
*ARGUMENTS: in src     source of information
*           in dest    destination of move
*           in nwords  number of words to copy
*
*/
static void MOVE (word *src, word *dest, int nwords)
{
int i;
for (i=0 ; i<nwords ; i++) {dest[i]=src[i];};
}

/****************************************************************************
*                                   ZERO                                    *
*****************************************************************************
*
***function zero
*
*OBJECT: Puts zeros in nwords words at dest.
*
*ARGUMENTS: in dest    destination in memory
*           in nwords  number of words to put to zero
*
*/
static void ZERO ( word *dest, int nwords)
{
int i;
for (i=0 ; i<nwords ; i++) {dest[i]=0;};
}
/****************************************************************************
*                              check_host_id                                *
****************************************************************************/
/*
 check that RMNLIB license file (node locked) is valid
 check_host_id is FORTRAN callable
 check_host_id returns the HOST id as obtained by gethostid
*/

unsigned INT_32 f77name(check_host_id)()
{
#if defined NEC || !defined CHECK_RMNLIB_LIC
return(0);
#else
FILE *id_file;
unsigned INT_32 sysid, key, domain_ok , junk;
char ypdomain[200];
char *ARMNLIB;

/* find YP(NIS) domain name */
junk=getdomainname(ypdomain,19);
/* find HOST id */
sysid=gethostid();
/* check that ARMNLIB is an environment variable */
ARMNLIB=getenv("ARMNLIB");
if (ARMNLIB==NULL){
  printf("ERROR: ARMNLIB environment variable not defined\n");
  exit(1);
}

/* if NIS domain name is cmcnet, no further check */
domain_ok= (ypdomain[0]=='c') && (ypdomain[1]=='m') &&
           (ypdomain[2]=='c') && (ypdomain[3]=='n') &&
           (ypdomain[4]=='e') && (ypdomain[5]=='t');

/* in test mode, ignore the NIS domain name */
#if defined(TEST)
domain_ok=0;
#endif

if(domain_ok)return(sysid);

/* license file name is $ARMNLIB/data/.LIC */
sprintf(ypdomain,"%s/data/.LIC",ARMNLIB);
id_file=fopen(ypdomain,"r");
if (id_file == NULL) {
  printf(" ERROR: RMNLIB LICENSE FILE IS NOT VALID\n");
  exit(1);
}

/* check all numeric tokens found in license file */
while( EOF != fscanf(id_file,"%u",&key)){
  domain_ok = domain_ok || (sysid ^ 0xCAFEFADE)==key ;
}

fclose(id_file);

if ( domain_ok) {
  /*  printf(" LICENSE is VALID\n"); */
  return(sysid);
}else{
  printf(" ERROR: RMNLIB LICENSE FILE IS NOT VALID\n");
  exit(1);
}
#endif
}
