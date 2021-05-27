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
/* write_slab.c
   original author: Karim Tegguiche
   Revision V2.0 - V.Lee July 1997
   Revision V3.0 - V.Lee March 1998 (to catch more bugs)
   Revision V4.0 - V.Lee November 1998 (slabini will not allow append)
   Revision V5.0 - V.Lee January 1999  (slabopt added for parallel code)
   Revision V6.0 - V.Lee Dec 2000(convip1 determines IP1, remove IP1 check)
   Revision V7.0 - M.Valin Jan 2001 added # type grids and little_endian patch
                                     bugfix: added nbytes to SLB9 marker
                                     added: capability to write to pipes/sockets
   Revision V8.0 - V.Lee May 2001 (introduced slabig34 to determine IG3 and IG4 
                            for # type grids, and open file even if it exists)
   Revision V9.0 - V.Lee Jan 2003(28bit check for IP1,IP2,IP3;24bit for NIO,NJO)
*/

#include<rpnmacros.h>
#include<fcntl.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <unistd.h>
#include <errno.h>
#include "slab.h"


#ifdef DEBUG
#define BUFSIZE 64
#else
#define BUFSIZE 1048576
#endif

#define GRTYPES "ABGLNSXYZ#EH"

/* macro to insert 32 bit integers into buffer. if buffer is full
   it gets written to appropriate file / pipe / socket */

#define put_in_buffer(file_desc,from,Ibuffer,pos,nbre) \
                 {unsigned INT_32 *ptemp=(unsigned INT_32 *)from;\
                  unsigned INT_32 *buffer=(unsigned INT_32 *)Ibuffer;\
                  int nbmots=nbre;\
                  while(nbmots--){ \
		    if(pos==BUFSIZE) {\
                      write_buf(file_desc,buffer,BUFSIZE);\
                      pos=0;\
                    }\
	            buffer[pos++] = *ptemp++;\
                  }\
                 }

/* character pointer aliased to integer. if value returned is 1
   (low byte) machine is little endian */

static int entier_quelconque=1;
static char *little_endian=(char *)&entier_quelconque;

static int ERROR_LEVEL=-2;   /* tolerance level for abort */

static file_table_desc file_table[MAX_SLAB_FILES];
static int f_index[MAX_SLAB_FILES];

static Id_Block_file  id_block;         /* SLB0 block */
static Slab_Descrt_file slab_descrt;    /* SLB1 block */
static Data_Block_file data_block;      /* SLB2 block */
static Slab_End slab_end;               /* SLB9 block */

static int *intBuffer, pos;
static int *iVal;
static ftnword *pVal;

static float *fBuffer;
static ftnfloat *fVal;

static int init=0;     /* global initialization flag */
static int proc0=1;    /* am I processor 0 */
static int numproc=1;  /* number of processors in run */

/* prototypes declaration */

static void init_index();
static int slab_exit(int level);
static int get_free_index(int fd);
static int get_file_index(int fd);
static int write_buf(int fd, void *buffer, int nitems);

ftnword f77name(slabini)(char *f_name, ftnword dateo[2], ftnword *f_npas,
		     ftnword *f_deet, char *f_etiket, F2Cl l1, F2Cl l2);
ftnword f77name(slabopt)(ftnword *f_proc, ftnword *f_numproc );
ftnword f77name(slabig34)(unsigned ftnword *f_ig3, unsigned ftnword *f_ig4,
                          ftnword *f_xmin, ftnword *f_xmax, 
                          ftnword *f_ymin, ftnword *f_ymax);
ftnword f77name(slabdsc)(ftnword *f_hand, ftnword *f_snum,char *f_gxtyp,
		     ftnword *f_ixyg1,ftnword *f_ixyg2,
		     ftnword *f_ixyg3, ftnword *f_ixyg4,ftnword *f_nio,
		     ftnword *f_njo,ftnword *f_nxgrid, 
		     ftnword *f_nygrid, ftnfloat *f_xgrid,ftnfloat *f_ygrid,
		     char *f_grtyp,ftnword *f_ig1,
		     ftnword *f_ig2,ftnword *f_ig3,ftnword *f_ig4,
		     ftnword *f_mtout, ftnword *f_np,
		     char *f_typvar,char *f_nomvar,ftnword *f_ip1,
		     ftnword *f_ip2, ftnword *f_ip3,
		     ftnword *f_datyp,ftnword *f_nbits,ftnword *iflt,
		     ftnfloat *f_xp, F2Cl l1, F2Cl l2, F2Cl l3, F2Cl l4);
ftnword f77name(slabxtr)(ftnword *f_hand, ftnword *f_snum, ftnword *f_nx,
		      ftnword *f_xnio,ftnword *f_mt,ftnword *f_mtas,
		      ftnfloat *f_mtadd, ftnfloat *f_mtmult, ftnfloat *f_mtval);
ftnword f77name(slabend)(ftnword *f_hand, char *f_sf_hand, F2Cl l1);


/*****************************************************************************
 *                            W R I T E _ B U F                              *
 *                                                                           *
 *Object                                                                     *
 * write integer buffer to file / socket / pipe                              *
 *                                                                           * 
 *****************************************************************************/
static int write_buf(int fd, void *Ibuffer, int nitems)
{
  unsigned INT_32 *buffer=(unsigned INT_32 *)Ibuffer;
  int n=nitems;
  int nwritten;
  char *cbuf=(char *)buffer;

  if(*little_endian){  /* slab files are BIG ENDIAN */
    unsigned INT_32 *tmpbuf=buffer;
    while(n--){
     unsigned INT_32 temp=*tmpbuf;
     *tmpbuf=SWAP32(temp);
     tmpbuf++;
    }
  }
  n=nitems*sizeof(INT_32);
  while(n>0){
    nwritten=write(fd,cbuf,n);
    /*
      if (errno != 0) 
      {
      int junk;
      char hname[64];
      perror("error type");
      junk=gethostname(hname,sizeof(hname));
      fprintf(stderr,"\n***ERROR in write_buf, hostname=%s errno=%d\n",hname,errno);
      }
    */
    if(nwritten<=0) return(nwritten);
    n-=nwritten;
    cbuf+=nwritten;
  }
  return(nitems*sizeof(INT_32));
}
/***************************************************************************** 
 *                            S L A B _ E X I T                              *
 *                                                                           * 
 *Object                                                                     *
 *  To check level of error. If less or equal to ERROR_LEVEL, EXIT program   *
 *                                                                           * 
 *****************************************************************************/

static int slab_exit(int level)
  {
  
   if (level <= ERROR_LEVEL)
       exit(level);
   return(level);
 }


/***************************************************************************** 
 *                            I N I T _ I N D E X                            *
 *                                                                           * 
 *Object                                                                     *
 *  Initialize f_index values to all -1. Should only be used once.           *
 *                                                                           * 
 *****************************************************************************/

static void init_index()
  {
   int i, j;
   char *slab_config;
  
   if ( (slab_config=getenv("SLAB_CONFIG")) != NULL){
    ERROR_LEVEL=atoi(slab_config);
   printf("NOTE: ERROR_LEVEL set to %d\n",ERROR_LEVEL);
    }
   for(i=0;i<MAX_SLAB_FILES;i++) {
       f_index[i] = -1;
       file_table[i].file_name[0]='\0';
       for(j=0;j<MAX_SLAB_TYPES;j++) {
         file_table[i].nrows[j] = 0;
         file_table[i].count[j] = 0;
         file_table[i].nio[j] = 0;
         file_table[i].ni[j] = 0;
         file_table[i].i1[j] = 0;
         file_table[i].njo[j] = 0;
         file_table[i].nj[j] = 0;
         file_table[i].j1[j] = 0;
       }
       file_table[i].buffer = NULL;
       file_table[i].pos = 0;
   }
   init = 1;
 }/* end init_index */

/***************************************************************************** 
 *                          G E T _ F R E E _ I N D E X                      *
 *                                                                           *
 *Object                                                                     *
 *  Return a free index entry in file table and record file information.     *
 *****************************************************************************/

static int get_free_index(int fd)
  {
   int i;
  
   for(i=0;i<MAX_SLAB_FILES;i++)
       if(f_index[i] == -1){
          f_index[i] = fd;
	  return(i);
          }

   fprintf(stderr,"\n***ERROR in GET_FREE_INDEX: slab file table is full\n");
   fprintf(stderr,"   MAX_SLAB_FILES = %d\n",MAX_SLAB_FILES);
   return(slab_exit(ERR_TAB_FULL));
   }

/***************************************************************************** 
 *                          G E T _ F I L E _ I N D E X                      *
 *                                                                           *
 *Object                                                                     *
 *  Return a file index entry in file table given the file handler           *
 *****************************************************************************/

static int get_file_index(int fd)
  {
   int i;
  
   for(i=0;i<MAX_SLAB_FILES;i++)
       if(f_index[i] == fd){
	  return(i);
          }

   fprintf(stderr,"\n***ERROR in GET_FILE_INDEX: slab file not initialized\n");
   return(slab_exit(ERR_NO_FILE));
   }

/***************************************************************************** 
 *                          S L A B O P T                                    *
 *                                                                           *
 *Object                                                                     *
 *  Get number of processors in run and my processor number                  *
 *                                                                           * 
 *Arguments                                                                  *
 *  IN   f_proc     my processor number                                      *
 *  IN   f_numproc  number of processors in run                              *
 *****************************************************************************/
ftnword f77name(slabopt)(ftnword *f_proc, ftnword *f_numproc)
{
     numproc = *f_numproc; 
     if ( *f_proc != 0)  proc0=0;
     else proc0 = 1;
     return(proc0);
}
/***************************************************************************** 
 *                              S L A B I N I                                * 
 *                                                                           *
 *Object                                                                     *
 *  Initialize slab file. Upon completion the function returns a file        * 
 *  descriptor related to the file                                           *
 *                                                                           * 
 *Arguments                                                                  *
 *  IN   f_name   file name                                                  * 
 *  IN   dateo    origin date, dimension 2.                                  *
 *                dateo[0] AAAAMMJJ                                          * 
 *                dateo[1] HHMMSS00                                          *
 *  IN   npas     time step number                                           * 
 *  IN   deet     length of time step                                        * 
 *  IN   etiket   record identificator                                       * 
 *****************************************************************************/

ftnword f77name(slabini)(char *f_name, ftnword dateo[2], ftnword *f_npas,
			ftnword *f_deet, char *f_etiket, F2Cl l1, F2Cl l2)
{                                               
 int fd, ix, i, j, taille, npas;
 char name[MAX_LEN], etiket[MAX_ETIKET];


 if (init == 0) init_index();

 l1 = (l1 < MAX_LEN) ? l1 : MAX_LEN-1;   /* get desired file name */
 strncpy(name,f_name,l1);
 name[l1] = '\0';

 while ((name[l1-1] == ' ') && (l1 > 1)) {
   l1--;
   name[l1] = '\0';
 }

 for(i=0;i<MAX_ETIKET;i++) etiket[i]='\0'; /* initialize etiket to nulls */

 l2 = (l2 < 12) ? l2 : MAX_ETIKET-1; 
 strncpy(etiket,f_etiket,l2);              /* copy etiket */
 etiket[l2] = '\0';
   
 while ((etiket[l2-1] == ' ') && (l2 > 1)) {
   l2--;
   etiket[l2] = '\0';
 }

 if((fd = open(name, O_RDWR | O_CREAT ,0744)) == ERR_NO_FILE)
    {
     int junk;
     char hname[64];
     junk=gethostname(hname,sizeof(hname));
     fprintf(stderr,"\n***ERROR in SLABINI: error opening file %s errno=%d hostname=%s\n",name,errno,hname);
     slab_exit(-3);
     }

 ix = get_free_index(fd);
 if (ix == ERR_TAB_FULL) {
     fprintf(stderr,"\n***ERROR in SLABINI(%s): slab file table is full\n",name);
     return(slab_exit(-2));
     }

 strcpy(file_table[ix].file_name,name);  /* initialize tables associated to file */
 for (j=0;j < MAX_SLAB_TYPES; j++){
   file_table[ix].count[j] = 0;
   file_table[ix].nrows[j] = 0;
   file_table[ix].nio[j] = 0;
   file_table[ix].i1[j] = 0;
   file_table[ix].ni[j] = 0;
   file_table[ix].njo[j] = 0;
   file_table[ix].j1[j] = 0;
   file_table[ix].nj[j] = 0;
 }

 if (file_table[ix].buffer) {
     fprintf(stderr,"\n***ERROR in SLABINI(%s): memory for buffer already allocated\n",name);
     return(slab_exit(-3));
 }

 if ((intBuffer = (int *) malloc(BUFSIZE * sizeof(int)))==NULL)
    {
     fprintf(stderr,"\n***ERROR in SLABINI(%s): Cannot allocate memory for buffer\n",name);
     return(slab_exit(-3));
     }
 
 file_table[ix].buffer = (unsigned INT_32*)intBuffer;
   
 id_block.slb0    = 'SLB0';
 id_block.nBytes  = 32;
 id_block.deet    = (int ) *f_deet;
 id_block.npas    = (int ) *f_npas;
                     /* stuff etiket into id_block.Ietiket in endian proof fashion */
 id_block.Ietiket[0] = etiket[0]<<24 | etiket[1]<<16 | etiket[2]<<8 | etiket[3];
 id_block.Ietiket[1] = etiket[4]<<24 | etiket[5]<<16 | etiket[6]<<8 | etiket[7];
 id_block.Ietiket[2] = etiket[8]<<24 | etiket[9]<<16 | etiket[10]<<8 | etiket[11];
 id_block.dateo1  = (int ) dateo[0];
 id_block.dateo2  = (int ) dateo[1];
 id_block.val15   = 1.5;

 pos=0;

 iVal = (int *) &id_block;
 taille = (sizeof(id_block) / sizeof(int));

 put_in_buffer(fd,iVal,intBuffer,pos,taille);  /* write SLB0 block */

 file_table[ix].pos = pos;
 return (ftnword) fd;  /* return file descriptor to caller */
}

/*****************************************************************************
 *                              SLABIG34                                     *
 *Object                                                                     *
 *  Determine the values of IG3 and IG4 for '#' grids.                       *
 *                                                                           * 
 *Arguments                                                                  *
 *  OUT   f_ig3    :   grid descriptor 3 for '#' grids                       *
 *                       istart + nipoints<<20                               *
 *  OUT   f_ig4    :   grid descriptor 4 for '#' grids                       *
 *                       jstart + njpoints<<20                               *
 *  IN   f_xmin    :   istart                                                * 
 *  IN   f_xmax    :   iend                                                  * 
 *  IN   f_ymin    :   jstart                                                * 
 *  IN   f_ymax    :   jend                                                  * 
 *                                                                           *
 ****************************************************************************/
ftnword f77name(slabig34)(unsigned ftnword *f_ig3, unsigned ftnword *f_ig4,
                          ftnword *f_xmin, ftnword *f_xmax, 
                          ftnword *f_ymin, ftnword *f_ymax)
{
 unsigned ftnword niout, njout;
 
 if (*f_xmin > 0xFFFFF) {
    fprintf(stderr,"\n***ERROR in SLABIG34: (XMIN=%d) > 2**19\n",*f_xmin);
    return(slab_exit(-2));
    }
 if (*f_xmin < 0) {
    fprintf(stderr,"\n***ERROR in SLABIG34: (XMIN=%d) < 0\n",*f_xmin);
    return(slab_exit(-2));
    }
 if (*f_ymin > 0xFFFFF) {
    fprintf(stderr,"\n***ERROR in SLABIG34: (YMIN=%d) > 2**19\n",*f_ymin);
    return(slab_exit(-2));
    }
 if (*f_ymin < 0) {
    fprintf(stderr,"\n***ERROR in SLABIG34: (YMIN=%d) < 0\n",*f_ymin);
    return(slab_exit(-2));
    }

 niout = *f_xmax - *f_xmin + 1;
 njout = *f_ymax - *f_ymin + 1;
  
 if (niout > 0xFFF) {
    fprintf(stderr,"\n***ERROR in SLABIG34: (XMAX-XMIN=%d) > 2**11\n",niout);
    return(slab_exit(-2));
    }
 if (njout > 0xFFF) {
    fprintf(stderr,"\n***ERROR in SLABIG34: (YMAX-YMIN=%d) > 2**11\n",njout);
    return(slab_exit(-2));
    }

 *f_ig3 = (niout<<20) | *f_xmin;
 *f_ig4 = (njout<<20) | *f_ymin;
  
  return(0);
}


/*****************************************************************************
 *                              S L A B D S C                                *
 *Object :                                                                   *
 *  Description slab file.                                                   *
 *Arguments :                                                                *
 *     IN   f_hand   :   file handler of slab file                           *
 *     IN   f_snum   :   slab number                                         *
 *     IN   f_gxtyp  :   grid type for >> and ^^                             *
 *     IN   f_ixyg1  :   grid descriptor 1 for >> and ^^                     *
 *     IN   f_ixyg2  :   grid descriptor 2 for >> and ^^                     *
 *     IN   f_ixyg3  :   grid descriptor 3 for >> and ^^                     *
 *     IN   f_ixyg4  :   grid descriptor 4 for >> and ^^                     *
 *     IN   f_nio    :   dimension x of grid                                 *
 *     IN   f_njo    :   dimension y of grid                                 *
 *     IN   f_nxgrid :   dimension of f_xgrid                                *
 *     IN   f_nygrid :   dimension of f_ygrid                                *
 *     IN   f_xgrid  :   >> fields (dim : f_nxgrid)                          *
 *     IN   f_ygrid  :   ^^ fields (dim : f_nygrid)                          *
 *     IN   f_grtyp  :   grid type (A,B,G,L,N,S,X,Y,Z,E,H)                   *
 *     IN   f_ig1    :   grid descriptor 1                                   *
 *     IN   f_ig2    :   grid descriptor 2                                   *
 *     IN   f_ig3    :   grid descriptor 3 (except for '#' grids)            *
 *                       istart + nipoints<<20  (for '#' grids)              *
 *     IN   f_ig4    :   grid descriptor 4 (except for '#' grids)            *
 *                       jstart + njpoints<<20  (for '#' grids)              *
 *                                                                           *
 *     IN   f_mtout  :   dimension of ip(1..3),nomvar,datyp,nbits and first  *
 *                       dimension of f_xp                                   *
 *     IN   f_np     :   second dimension of f_xp                            *
 *     IN   f_typvar :   type field                                          *
 *     IN   f_nomvar :   variable name                                       *
 *     IN   f_ip1    :   level descriptor                                    *
 *     IN   f_ip2    :   time descriptor                                     *
 *     IN   f_ip3    :   descriptor no 3                                     *
 *     IN   f_datyp  :   data type                                           *
 *     IN   f_nbits  :   bits number for each data field                     *
 *     IN   f_iflt   :   number of filter passes for each data field         *
 *     IN   f_xp     :   optional variables (dim : mtout,np)                 * 
 *                                                                           *
 *****************************************************************************/

ftnword f77name(slabdsc)(ftnword *f_hand, ftnword *f_snum,char *f_gxtyp,
		     ftnword *f_ixyg1,ftnword *f_ixyg2,
                     ftnword *f_ixyg3, ftnword *f_ixyg4,ftnword *f_nio,
		     ftnword *f_njo, ftnword *f_nxgrid,
                     ftnword *f_nygrid, ftnfloat *f_xgrid,ftnfloat *f_ygrid,
		     char *f_grtyp,ftnword *f_ig1,
		     ftnword *f_ig2,ftnword *f_ig3,ftnword *f_ig4,
		     ftnword *f_mtout, ftnword *f_np, 
                     char *f_typvar, char *f_nomvar,ftnword *f_ip1,
		     ftnword *f_ip2, ftnword *f_ip3,
                     ftnword *f_datyp,ftnword *f_nbits,ftnword *f_iflt,
		     ftnfloat *f_xp, F2Cl l1, F2Cl l2, F2Cl l3, F2Cl l4)
{
 int nrows, nio, njo, nxtra, 
     taille, i, x, MAX_GRTYP;
 int lng, n, nn, ix, snum;
 char grtyp[MAX_LEN], grtyp_[MAX_LEN];
 char *p_nomvar, *p_typvar;
 int *typvars_0, *nomvars_0, *nomvars, *typvars;
 ftnword *p_ip1,*p_ip2,*p_ip3,*p_datyp,*p_nbits,*p_iflt;
 int fd;

/* for debug only 
  int ip1[20],ip2[20],ip3[20],datyp[20],nbits[20],iflt[20];
  ftnword *p_ip1,*p_ip2,*p_ip3,*p_datyp,*p_nbits,*p_iflt; 
*/

 MAX_GRTYP = strlen(GRTYPES);
 nio = (int ) *f_nio;
 njo = (int ) *f_njo;
 nrows = (int ) *f_mtout;
 nxtra = (int ) *f_np;

 if ( (ix = get_file_index( (int ) *f_hand)) < 0 ) return(slab_exit(-3));
 fd=f_index[ix];

 intBuffer = (int *)file_table[ix].buffer;
 pos = file_table[ix].pos;

 snum = *f_snum;
 
 if (snum < MAX_SLAB_TYPES && snum >=0){  /* valid slab type number ? */
       if (file_table[ix].nrows[snum] == 0){  /* slab type already defined ? */
           file_table[ix].nrows[snum] = nrows;
           file_table[ix].count[snum] = nio*njo;
           file_table[ix].nio[snum] = nio;
           file_table[ix].ni[snum] = nio;
           file_table[ix].i1[snum] = 1;
           file_table[ix].njo[snum] = njo;
           file_table[ix].nj[snum] = njo;
           file_table[ix].j1[snum] = 1;
           if(*f_grtyp == '#') {
               file_table[ix].ni[snum] = (*f_ig3 >> 20) & 0xFFF;  /* number of points along X axis */
               file_table[ix].i1[snum] = *f_ig3 & 0xFFFFF;        /* X offset of first point */
               file_table[ix].nj[snum] = (*f_ig4 >> 20) & 0xFFF;  /* number of points along Y axis */
               file_table[ix].j1[snum] = *f_ig4 & 0xFFFFF;        /* Y offset of first point */
               file_table[ix].count[snum] = file_table[ix].ni[snum] * file_table[ix].nj[snum];
               }
           }
       else{  /* already defined, OOPS */
         fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d is already defined\n",
                       file_table[ix].file_name,snum);
         fprintf(stderr,"   mtout=%d, nio*njo=%d\n",
                       file_table[ix].nrows[snum],file_table[ix].count[snum]);
         fprintf(stderr,"   set to: mtout=%d, nio=%d njo=%d ??\n",
                       nrows,nio,njo);
         return(slab_exit(-2));
       }
       }
 else{  /* slab type number out of valid range, OOPS */
       fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d is out of range\n",file_table[ix].file_name,snum);
       fprintf(stderr,"  slabid MUST be from 0 to %d\n",MAX_SLAB_TYPES-1);
       return(slab_exit(-2));
     }
 if (nio < 1 || nio > 16777000){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: invalid NIO=%d\n",file_table[ix].file_name,snum,nio); 
      return(slab_exit(-2));
      }
 if (njo < 1 || njo > 16777000){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: invalid NJO=%d\n",file_table[ix].file_name,snum,njo); 
      return(slab_exit(-2));
      }
 if ((nomvars_0 = (int *)malloc(4 * nrows))==NULL)
    {
     fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: Cannot allocate memory for nomvars\n",file_table[ix].file_name,snum);
     return(slab_exit(-3));
     }

 if ((typvars_0 = (int *)malloc(4 * nrows))==NULL)
    {
     fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: Can not allocate memory for typvars\n",file_table[ix].file_name,snum);
     return(slab_exit(-3));
     }

 nomvars = nomvars_0;
 typvars = typvars_0;
 p_nomvar = f_nomvar;
 p_typvar = f_typvar;

 lng = (l4 < 4) ? l4 : 4;        /* copy nomvar into an integer using an endian proof method */
 for (n=0; n < nrows; n++) {
   int temp=0;
   for (i=0; i < lng; i++) {
     temp = (temp << 8) | (*p_nomvar & 0xFF);
     p_nomvar++;
   }
   for (i=lng; i < 4; i++) {
     temp = (temp <<= 8) | ' ';
   }
   *nomvars = temp;
   nomvars++;
 }

 lng = (l3 < 4) ? l3 : 4;        /* copy typvar into an integer using an endian proof method */
 for (n=0; n < nrows; n++) {
   int temp=0;
   for (i=0; i < lng; i++) {
     temp = (temp << 8) | (*p_typvar & 0xFF);
     p_typvar++;
   }
   for (i=lng; i < 4; i++) {
     temp = (temp <<= 8) | ' ';
   }
   *typvars = temp;
   typvars++;
 }

/* check ip1s */
     p_ip1 = f_ip1;
     for (n=0; n< nrows; n++){
      if ( (*p_ip1 & 0xF0000000) != 0){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: ip1[%d]=%d\n",file_table[ix].file_name,snum,n+1,(int) *p_ip1);
      return(slab_exit(-2));
      }
      p_ip1++;
      }

/* check ip2s */
     p_ip2 = f_ip2;
     for (n=0; n< nrows; n++){
      if ( (*p_ip2 & 0xF0000000) != 0){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: ip2[%d]=%d\n",file_table[ix].file_name,snum,n+1,(int) *p_ip2);
      return(slab_exit(-2));
      }
      p_ip2++;
      }

/* check ip3s */
     p_ip3 = f_ip3;
     for (n=0; n< nrows; n++){
      if ( (*p_ip3 & 0xF0000000) != 0){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: ip3[%d]=%d\n",file_table[ix].file_name,snum,n+1,(int) *p_ip3);
      return(slab_exit(-2));
      }
      p_ip3++;
      }

/* check datyps */
     p_datyp = f_datyp;
     for (n=0; n< nrows; n++){
      if ( (int) *p_datyp < 0  || (int) *p_datyp > 5){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: datyp[%d]=%d\n",file_table[ix].file_name,snum,n+1,(int) *p_datyp);
      return(slab_exit(-2));
      }
      p_datyp++;
      }

/* check nbits */
     p_nbits = f_nbits;
     for (n=0; n< nrows; n++){
      if ( (int) *p_nbits < 0 || (int) *p_nbits >  32767){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: nbits[%d]=%d\n",file_table[ix].file_name,snum,n+1,(int) *p_nbits);
      return(slab_exit(-2));
      }
      p_nbits++;
      }

/* check iflts */
     p_iflt = f_iflt;
     for (n=0; n< nrows; n++){
      if ( (int) *p_iflt < 0 || (int) *p_iflt >  32767){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: iflt[%d]=%d\n",file_table[ix].file_name,snum,n+1,(int) *p_iflt);
      return(slab_exit(-2));
      }
      p_iflt++;
      }

/* for debug only 
     printf("INSIDE SLABDSC ix=%d\n",ix);
     printf("l1=%d,l2=%d,l3=%d,l4=%d\n",l1,l2,l3,l4);
     printf("nio=%d,njo=%d,snum=%d,nrows=%d, nxtra=%d\n",
                           nio,njo,snum,nrows,nxtra);
     printf("ip1= ");
     p_ip1 = f_ip1;
     for (n=0; n< nrows; n++){
      ip1[n] = (int )*p_ip1;
      p_ip1++;
      printf("%d,",ip1[n]);
      }
     printf("\n");
     printf("ip2= ");
     p_ip2 = f_ip2;
     for (n=0; n< nrows; n++){
      ip2[n] = (int )*p_ip2;
      p_ip2++;
      printf("%d,",ip2[n]);
      }
     printf("\n");
     printf("ip3= ");
     p_ip3 = f_ip3;
     for (n=0; n< nrows; n++){
      ip3[n] = (int )*p_ip3;
      p_ip3++;
      printf("%d,",ip3[n]);
      }
     printf("\n");
     printf("datyp= ");
     p_datyp = f_datyp;
     for (n=0; n< nrows; n++){
      datyp[n] = (int )*p_datyp;
      p_datyp++;
      printf("%d,",datyp[n]);
      }
     printf("\n");
     printf("nbits= ");
     p_nbits = f_nbits;
     for (n=0; n< nrows; n++){
      nbits[n] = (int )*p_nbits;
      p_nbits++;
      printf("%d,",nbits[n]);
      }
     printf("\n");
     printf("iflt= ");
     p_iflt = f_iflt;
     for (n=0; n< nrows; n++){
      iflt[n] = (int )*p_iflt;
      p_iflt++;
      printf("%d,",iflt[n]);
      }
     printf("\n");
     printf("typvars=%s\n",typvars_0);
     printf("nomvars=%s\n",nomvars_0); 
     printf("grtyp=%s,",f_grtyp);
     printf("nxgrid=%d,nygrid=%d\n",(int )*f_nxgrid,(int )*f_nygrid);
     printf("gxtyp=%s\n",f_gxtyp);
     printf("ig1=%d,ig2=%d,ig3=%d,ig4=%d\n",
           (int) *f_ig1,(int) *f_ig2,(int) *f_ig3,(int) *f_ig4);
     printf("ixyg1=%d,ixyg2=%d,ixyg3=%d,ixyg4=%d\n",
           (int) *f_ixyg1, (int) *f_ixyg2, (int) *f_ixyg3, (int) *f_ixyg4);
     printf("\n");
 debug only*/ 

 slab_descrt.slb1 = 'SLB1';
 slab_descrt.nBytes = 4 * (15 + nrows * (8 + nxtra) + (int) *f_nxgrid + (int) *f_nygrid);
 slab_descrt.slab_id = (int ) *f_snum;
  
 l3 = (l3 < MAX_LEN) ? l3 : MAX_LEN-1;
 strncpy(grtyp,f_grtyp,l3);
 grtyp[l3] ='\0';
 
/* check gridtype requested */
 if (l3 == 0 || *grtyp==' '){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: invalid GRTYP='%s'\n",file_table[ix].file_name,snum,grtyp);
      return(slab_exit(-2));
      }
 n=0;
 while (*grtyp != GRTYPES[n] && n < MAX_GRTYP) n++;
 if (*grtyp != GRTYPES[n]){
    fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: invalid GRTYP='%s'\n",
            file_table[ix].file_name,snum,grtyp);
      return(slab_exit(-2));
      }
 slab_descrt.Igrtyp = *grtyp; /* endian proof method */
 slab_descrt.Igrtyp <<= 24;
 
 slab_descrt.ig1    = (int ) *f_ig1;
 slab_descrt.ig2    = (int ) *f_ig2;
 slab_descrt.ig3    = (int ) *f_ig3;
 slab_descrt.ig4    = (int ) *f_ig4;
 slab_descrt.Nrows  = (int ) *f_mtout;
 slab_descrt.Niout  = (int ) *f_nio;
 slab_descrt.Njout  = (int ) *f_njo;
 slab_descrt.nxgrid = (int ) *f_nxgrid;
 slab_descrt.nygrid = (int ) *f_nygrid;
 slab_descrt.Nextra = (int ) *f_np;
 
 slab_descrt.Igrtyp_ = ' ' << 24;  /* endian proof method */
 slab_descrt.ig1_   = -2;
 slab_descrt.ig2_   = -2;
 slab_descrt.ig3_   = -2;
 slab_descrt.ig4_   = -2;

 if(*f_grtyp == 'Z' || *f_grtyp == '#') /* check descriptor dimensions vs grid dimensions for Z and # grids */
    {
      if(njo != *f_nygrid) {
	  fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: nygrid should be equal to njo for Z grid\n",
                        file_table[ix].file_name,snum);
          fprintf(stderr," nygrid = %d njo = %d\n",(int )*f_nygrid, njo);
	  return(slab_exit(-2));
      }
      if( *f_nxgrid != nio && *f_nxgrid != (nio+1) ) {
          fprintf(stderr,
                 "\n***ERROR in SLABDSC(%s)slabid %d: nxgrid should be equal to nio or (nio+1) for Z grid\n",
                 file_table[ix].file_name,snum);
          fprintf(stderr," nxgrid = %d nio = %d\n",(int )*f_nxgrid, nio);
          return(slab_exit(-2));
      }
    }

 if(*f_grtyp == 'Y') {
    if( (nio*njo) != *f_nxgrid) {
	  fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: nxgrid should be equal to nio*njo for Y grid\n",
                         file_table[ix].file_name,snum);
          fprintf(stderr," nxgrid = %d nio = %d njo = %d\n",(int )*f_nxgrid, nio,njo);
	  return(slab_exit(-2));
    }
    if( (nio*njo) !=  *f_nygrid) {
	  fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: nygrid should be equal to nio*njo for Y grid\n",
                         file_table[ix].file_name,snum);
          fprintf(stderr," nygrid = %d nio = %d njo = %d\n",(int )*f_nygrid, nio,njo);
	  return(slab_exit(-2));
    }
 }
 if(*f_grtyp == 'Z' || *f_grtyp == 'Y' || *f_grtyp == '#'){  /* collect >> and ^^ projection descriptors */
      l1 = (l1 < MAX_LEN) ? l1 : MAX_LEN-1;
      strncpy(grtyp_,f_gxtyp,l1);
      grtyp_[l1] = '\0';
      slab_descrt.Igrtyp_ = *grtyp_; /* endian proof method */
      slab_descrt.Igrtyp_ <<= 24;

      if (l1 == 0 || *grtyp_==' '){
      fprintf(stderr,"\n***ERROR in SLABDSC(%s)slabid %d: invalid GXTYP='%s'\n",file_table[ix].file_name,snum,grtyp_);
      return(slab_exit(-2));
      }
      
      slab_descrt.ig1_ = (int ) *f_ixyg1;
      slab_descrt.ig2_ = (int ) *f_ixyg2;
      slab_descrt.ig3_ = (int ) *f_ixyg3;
      slab_descrt.ig4_ = (int ) *f_ixyg4;
 }

 
 taille = (sizeof(slab_descrt) / sizeof(int));
 iVal = (int *) &slab_descrt;
 put_in_buffer(fd,iVal,intBuffer,pos,taille); /* SLB1 block header */
 pVal = f_ip1;
 put_in_buffer(fd,pVal,intBuffer,pos,nrows);
 pVal = f_ip2;
 put_in_buffer(fd,pVal,intBuffer,pos,nrows);
 pVal = f_ip3;
 put_in_buffer(fd,pVal,intBuffer,pos,nrows);
 pVal = f_nbits;
 put_in_buffer(fd,pVal,intBuffer,pos,nrows);
 pVal = f_datyp;
 put_in_buffer(fd,pVal,intBuffer,pos,nrows);
 iVal = nomvars_0;
 put_in_buffer(fd,iVal,intBuffer,pos,nrows);
 iVal = typvars_0;
 put_in_buffer(fd,iVal,intBuffer,pos,nrows);

 fBuffer = (float *) intBuffer;  

 if( (*f_grtyp == 'Z')  || ( *f_grtyp == 'Y') || ( *f_grtyp == '#') )  /* write optional part for >> and ^^ descriptors */
 {
#ifdef DEBUG
  printf("f_xgrid= ");
  {int npts=*f_nxgrid ; float *ffVal = (float *)f_xgrid; while(npts--){printf("  %f ",ffVal[npts]);}}
  printf("\n");
  printf("f_xgrid= ");
  {int npts=*f_nygrid ; float *ffVal = (float *)f_ygrid; while(npts--){printf("  %f ",ffVal[npts]);}}
  printf("\n");
#endif
     fVal = f_xgrid;
     put_in_buffer(fd,fVal, fBuffer, pos, *f_nxgrid);
     fVal = f_ygrid;
     put_in_buffer(fd,fVal, fBuffer, pos, *f_nygrid);
 }

 pVal = f_iflt;
 put_in_buffer(fd,pVal,intBuffer,pos,nrows);  /* filtering information (deprecated) */
 
 if(nxtra != 0) /* write optional extras if supplied  */
 {
    taille = (int )(*f_mtout)  *  (int )(*f_np);
    fVal = f_xp;
    put_in_buffer(fd,fVal,fBuffer,pos,taille);
 }

 file_table[ix].pos = pos;
 free (nomvars_0);
 nomvars_0 = NULL;
 free (typvars_0);
 typvars_0 = NULL;
 return(0);
 }/* end slabdsc */

/****************************************************************************
 *                              S L A B X T R F                             *
 *                                                                          *
 *Object:                                                                   *
 *       Fast extraction of valid data field into the slab file             *
 *Arguments:                                                                *
 *     IN   f_hand   :    file handler of slab file                         *
 *     IN   f_snum   :    slab number                                       *
 *     IN   f_nx     :    f_xnio dimension                                  *
 *     IN   f_xnio   :    output indicator for grid position (dim : f_nx)   *
 *     IN   f_mt     :    mtas dimension                                    *
 *     IN   f_mtval  :    raw slab; value of elements (dim : f_nx, f_mt)    *
 *                                                                          *
 ****************************************************************************/
ftnword f77name(slabxtrf)(ftnword *f_hand, ftnword *f_snum, ftnword *f_nx,
		      ftnword *f_xnio,ftnword *f_mt, ftnfloat *f_mtval)
{                                               
  return f77name(slabxtr)(f_hand, f_snum, f_nx, f_xnio,f_mt,NULL,NULL,NULL, f_mtval);
}
/****************************************************************************
 *                              S L A B X T R                               *
 *                                                                          *
 *Object:                                                                   *
 *       Extraction of valid data field into the slab file                  *
 *Arguments:                                                                *
 *     IN   f_hand   :    file handler of slab file                         *
 *     IN   f_snum   :    slab number                                       *
 *     IN   f_nx     :    f_xnio dimension                                  *
 *     IN   f_xnio   :    output indicator for grid position (dim : f_nx)   *
 *     IN   f_mt     :    mtas dimension                                    *
 *     IN   f_mtas   :    output indicator for row in slab  (dim: f_mt)     *
 *     IN   f_mtadd  :    values to add for each row in slab                *
 *     IN   f_mtmult :    values to multiply for each row in slab           *
 *     IN   f_mtval  :    raw slab; value of elements (dim : f_nx, f_mt)    *
 *                                                                          *
 ****************************************************************************/

ftnword f77name (slabxtr)(ftnword *f_hand, ftnword *f_snum, ftnword *f_nx,
		      ftnword *f_xnio,ftnword *f_mt,ftnword *f_mtas,
		      ftnfloat *f_mtadd, ftnfloat *f_mtmult, ftnfloat *f_mtval)
{                                               
 int i, j, ix, ij = 0, 
     taille,
     nrows = 0,
     nX = 0,
     nBytes, Nrows, Nx, fd, snum;
 int n;

 Nrows = (int ) *f_mt;
 Nx    = (int ) *f_nx; 
 snum  = (int ) *f_snum;

 if ( (ix = get_file_index( (int ) *f_hand)) < 0 ) return(slab_exit(-3));
 fd = (int) *f_hand;

 intBuffer = (int *)file_table[ix].buffer; 
 pos = file_table[ix].pos;
 fBuffer = (float *) intBuffer;


 if(f_mtas) {  /* take all rows if no indicator (null pointer) */
   for(j = 0; j < Nrows; j++)
     if(f_mtas[j] != 0)
       nrows++;
 }else{
   nrows=Nrows;
 }
 
 if(nrows != file_table[ix].nrows[snum]) {
    fprintf(stderr,"***ERROR in SLABXTR(%s)slabid %d:\n",
             file_table[ix].file_name,snum);
    fprintf(stderr,"  nrows in mtas(=%d) must be equal to SLABDSC mtout(=%d)\n",
             nrows,file_table[ix].nrows[snum]);
    return(slab_exit(-2));
 }

 for(i = 0; i < Nx; i++)
    if(f_xnio[i] != 0)
       nX++;

 if((f_mtas==NULL) && (nX != Nx)) {
    fprintf(stderr,"***ERROR in SLABXTRF(%s)slabid %d:\n",
             file_table[ix].file_name,snum);
    fprintf(stderr," HOLES NOT ALLOWED along X\n");
    return(slab_exit(-2));
 }

 file_table[ix].count[snum] -= nX;
 
 data_block.slb2 = 'SLB2';
 data_block.nBytes = 4 * (3 + nX * (1 + nrows));
 data_block.slab_id = (int ) *f_snum;
 data_block.nX = nX; 
 data_block.Nrows = nrows;
 
 iVal = (int *) &data_block;
 taille = (sizeof(data_block) / sizeof(int));

#ifdef DEBUG
 fprintf(stderr,"SLB2 header, pos avant=%d\n",pos);
#endif
 put_in_buffer(fd,iVal,intBuffer,pos,taille);
#ifdef DEBUG
 fprintf(stderr,"SLB2 header, pos apres=%d\n",pos);
#endif

 put_in_buffer(fd,f_xnio,intBuffer,pos,Nx);  /* insert column position markers into buffer */

 if(f_mtas==NULL && f_mtmult==NULL && f_mtadd==NULL) {  /* bounced call from slabxtrf */
   put_in_buffer(fd,f_mtval,intBuffer,pos,Nx*Nrows);
   file_table[ix].pos = pos;
   return(0);
 }

 for(j = 0; j < Nrows; j++){ /* insert data into buffer after scaling */
     if(f_mtas[j] != 0) {
        if((pos + nX) <= BUFSIZE) {
            if(nX == Nx ) {             /* Si le nombre de colonnes a extraire est egal au nombre   */
               for(i = 0; i < Nx; i++){ /* de colonnes de chaque slab     */
                   fBuffer[pos] = (float)(f_mtval[ij]*f_mtmult[j] + f_mtadd[j]);
                   pos++;
                   ij++;
               }/* end for */
            }else{  /* Si le nombre de colonnes a extraire est different de celui de la slab */
                 for(i = 0; i < Nx; i++){ 
                     if(f_xnio[i] != 0) { /* je cherche les valeurs valides une par une  */
                        fBuffer[pos] = (float)(f_mtval[ij]*f_mtmult[j] + f_mtadd[j]);
                        pos++;
                     }
                     ij++;
                 }/* end for */

            }
        }else{   /* Si la taille qui reste dans le buffer est insuffisante pour le reste des donnees */
             for(i = 0; i < Nx; i++) {
                 if(pos >= BUFSIZE) {     
                    taille = sizeof(float) * pos;  /* alors, on ecrit le buffer dans */
                    nBytes =  write_buf(fd, (unsigned INT_32 *)fBuffer, pos);

		    if(nBytes != (sizeof(float) *  pos)) {
		        fprintf(stderr,"\n***ERROR in SLABXTR(%s)slabid %d: WRITE ERROR in slab file\n",
                                       file_table[ix].file_name,snum);
                        fprintf(stderr,"tried to write %d bytes, wrote %d, pos*sizeof(float)=%d\n",
                                       taille,nBytes,BUFSIZE*sizeof(float));
                        fprintf(stderr,"Trying to write slab data\n");
			return(slab_exit(-2));
		    }
	            pos = 0;             /* le fichier slab et on met le   */
                                         /* buffer a zero                  */
		 }

                 if(f_xnio[i] != 0) {    /* a extraire,                    */
                    fBuffer[pos] =(float)(f_mtval[ij]*f_mtmult[j] + f_mtadd[j]);
                    pos++;
                 }
                 
                 ij++;

             }  /* end for */
             
        }
     
     }else{  /* f_mtas[j] == 0 */
         ij += Nx;
     }

 }/* end for */

 file_table[ix].pos = pos;
 return(0);
 }/* end slabxtr */

/*****************************************************************************
 *                              S L A B E N D                                *
 * Object :                                                                  *
 *         Put an end slab indicator into a file                             *
 * Arguments :                                                               *
 *            IN   f_hand : file handler of slab file                        *
 *            IN   sf_end : section or file end indicator                    *
 *                                                                           *
 *****************************************************************************/

ftnword f77name(slabend)(ftnword *f_hand, char *f_sf_hand, F2Cl l1)
 {
  int end, taille, i, ix, fd; 
  int nBytes;
  int n;

  if ((ix = get_file_index( (int ) *f_hand)) < 0 ) return(slab_exit(-3));
  fd = (int) *f_hand;

  intBuffer = (int *)file_table[ix].buffer;
  pos = file_table[ix].pos;
  fBuffer = (float *) intBuffer;

  /* end marker (SLB9) in endian proof way */
  end = (f_sf_hand[0] << 24) | (f_sf_hand[1]  << 16) | 
        (f_sf_hand[2] << 8) | (f_sf_hand[3]);

  /* check to see if number of values written is equal to what was requested */
  /* the count field originally contains the number of columns to write and  */
  /* is decremented by slabxtr for each column written. must be zero now     */
  for (i=0; i < MAX_SLAB_TYPES; i++)
    if(0 != file_table[ix].count[i] )
      {
      fprintf(stderr,"\n***ERROR in SLABEND(%s)slabid %d\n", file_table[ix].file_name,i);
      fprintf(stderr,"   Value of nio*njo must be equal to number of valid values in xnio\n");
      fprintf(stderr,"   No. of selected elements in xnio = %d, nio*njo=%d\n",
                      file_table[ix].count[i],file_table[ix].count[i]);
      return(slab_exit(-2));   
      }

  slab_end.id_end = end;
  slab_end.nBytes = 0;

  if(slab_end.id_end != 'SLB9')
     {
      fprintf(stderr,"\n***ERROR in SLABEND(%s):end indicator of slabend must be -> SLB9\n\n",file_table[ix].file_name);
      return(slab_exit(-2));
      }

  iVal = (int *) &slab_end;

  taille = sizeof(slab_end) / sizeof(INT_32);
  /* add SLB9 block to buffer */
  put_in_buffer(fd,iVal,intBuffer,pos,taille);

  if(pos == 0)
     {
      taille=sizeof(slab_end);
      nBytes = write_buf(fd,(unsigned INT_32 *)intBuffer,taille/sizeof(INT_32));
      if(nBytes != taille)
	 {
          fprintf(stderr,"\n***ERROR in SLABEND(%s): WRITE ERROR in slab file, errno=%d\n",file_table[ix].file_name,errno);
	  return(slab_exit(-2));
	  }
      }
  
  else{
       taille = sizeof(int) * pos;
       nBytes = write_buf(fd,(unsigned INT_32 *)intBuffer,pos);
       if(nBytes != taille)
	  {
	  fprintf(stderr,"\n***ERROR in SLABEND(%s): WRITE ERROR in slab file, errno=%d\n",file_table[ix].file_name,errno);
	  return(slab_exit(-2));
	  }
       
       }

  free(intBuffer);
  file_table[ix].buffer = NULL;
  file_table[ix].pos = 0;
   for (i=0;i < MAX_SLAB_TYPES; i++){
       file_table[ix].count[i] = 0;
       file_table[ix].nrows[i] = 0;
       file_table[ix].nio[i] = 0;
       file_table[ix].i1[i] = 0;
       file_table[ix].ni[i] = 0;
       file_table[ix].njo[i] = 0;
       file_table[ix].j1[i] = 0;
       file_table[ix].nj[i] = 0;
       }
  close (fd) ;
  f_index[ix] = -1; /* reset to -1 */
  return(0); /* indicates that SLB9 was added */
  }/* end slabend */
