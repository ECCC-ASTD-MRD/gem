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

#include<rpnmacros.h>
#include<fnom.h>
#include<stdio.h>
#include<stdlib.h>
#include<fcntl.h>

#ifdef WIN32	/*CHC/NRC*/
#include <io.h>
#endif

#define Max_File 10


ftnword f77name(open_db_file) (ftnword * );
ftnword f77name(close_db_file) (ftnword * );
ftnword f77name(read_db_file) (ftnword *, ftnword *, ftnword *);
ftnword f77name(write_db_file) (ftnword *, ftnword *, ftnword  *);
ftnword f77name(rewind_db_file) (ftnword *);
static int file_index(int );
static void init_TB_package(int);

static int BUF_SIZE=2*1024*1024 ;
static int BUF_SIZE_W ;
static int done=0;

typedef struct {
  int iun, fd;
  int pos_in, pos_out;
  int nb_elmt_lu, nb_elmt_ecrit, nb_elmt_Buf_IN;
  ftnword *pBuf_IN, *pBuf_OUT;
} Twin_Buffer;

static Twin_Buffer TB[Max_File];


/***************************************************************************** 
 *                            F N O M _ I N D E X                            *
 *                                                                           * 
 *Object:                                                                    *
 *       Find file position in file table.                                   *
 *                                                                           * 
 *Argument:                                                                  * 
 *          IN /  IUN : fortran unit                                         *
 *          OUT/  I   : position of the file                                 *
 *                                                                           *
 *****************************************************************************/

static int fnom_index(int iun)
{
   int i;
   
   for (i=0; i < MAXFILES; i++)           /* find file index in file table */
      if (FGFDT[i].iun == iun) return(i);
   return(-1);
   
   }

/****************************************************************************
 *                     I N I T _ T B _ P A C K A G E                        *
 *                                                                          *
 * Object :                                                                 *
 *    Initialize Twin_Buffer structure table.                               *
 *                                                                          *
 *                                                                          *
 ****************************************************************************/
static void init_TB_package(int inutile)
{
  int i;
  char *envar;

  envar=getenv("DB_FILE_CONFIG");
  if ( envar != NULL ) BUF_SIZE = atoi (envar)*sizeof(ftnword);


  BUF_SIZE_W = (BUF_SIZE / sizeof(ftnword));
  for (i=0; i < Max_File; i++) {
    TB[i].fd = -1;
    TB[i].iun = -1;
    TB[i].pos_in = 0;
    TB[i].pos_out = 0;
    TB[i].nb_elmt_lu = 0;
    TB[i].nb_elmt_ecrit = 0;
    TB[i].nb_elmt_Buf_IN = 0;
    TB[i].pBuf_IN = NULL;
    TB[i].pBuf_OUT = NULL;
  }
  done = 1;
}


/****************************************************************************
 *                          O P E N _ D B _ F I L E                         *
 *                                                                          *
 * Object :                                                                 *
 *    Return file descriptor. Allocate space for buffers.                   *
 *                                                                          *
 * Arguments :                                                              *
 *             IN/  IUN : fortran unit number                               *
 *                                                                          *
 ****************************************************************************/

ftnword f77name(open_db_file)(ftnword *iun)                                    
{
  int indx, f, i;

  if (! done) {
    fprintf(stderr,"Initializing  Dual Buffer Package\n");
    init_TB_package(123);
    fprintf(stderr,"Dual Buffer Package Initialized\n");
    fprintf(stderr,"Buffer size =%u Words\n",BUF_SIZE_W);
  }

  indx = fnom_index(*iun);
  if (indx < 0) {
    fprintf(stderr,
	    "open_db_file error: file (unit=%d) not connected with fnom\n",
	    *iun);
    exit(4);
  }

  f = -1;
  for (i=0; i < Max_File; i++)
    if (TB[i].iun == -1) {
      f = i;
      break;
    }
  if (f == -1) {
    fprintf(stderr,"open_db_file error: too may files\n");
    exit(5);
  }

  TB[f].fd = FGFDT[indx].fd;
  TB[f].iun = *iun;
  
  if((TB[f].pBuf_IN = malloc (BUF_SIZE)) == NULL)
    {
      fprintf(stderr,"CAN'T ALLOCATE MEMORY FOR BUFFER IN\n");
      exit(1);
    }
  
  if((TB[f].pBuf_OUT = malloc (BUF_SIZE)) == NULL)
    {
      fprintf(stderr,"CAN'T ALLOCATE MEMORY FOR BUFFER OUT\n");
      exit(1);
    }
   
  return(TB[f].fd);
  
} /* end open_db_file */


/****************************************************************************
 *                       R E A D _ D B _ F I L E                            *
 *                                                                          *
 * OBJECT :                                                                 *
 *         read data in the bucket.                                         *
 *                                                                          *
 * ARGUMENTS:                                                               *
 *            IN/ NB     : number of elements to read                       *
 *            IN/ fd     : file descripter                                  *
 *        OUT-IN/ bucket : in which data are returned                       *
 *                                                                          *
 ****************************************************************************/

ftnword f77name(read_db_file)( ftnword *iun, ftnword *bucket, ftnword *NB) 
{
   int i, nbytes,  data_a_lire;
   int fd, f;
   ftnword *p_bucket, *p_bufin;
 
   f = -1;
   for (i=0; i < Max_File; i++)
     if (TB[i].iun == *iun) {
       f = i;
       break;
     }
   if (f == -1) {
     fprintf(stderr,"read_db_file error: file (unit=%d) not open\n",*iun);
     exit(7);
   }

   fd = TB[f].fd;
   p_bucket = &(bucket[0]);
   data_a_lire = *NB;

   while (data_a_lire) {

     p_bufin = &(TB[f].pBuf_IN[TB[f].pos_in]);

     if (data_a_lire <= TB[f].nb_elmt_Buf_IN) {
       for (i=0; i < data_a_lire; i++)
	 p_bucket[i] = p_bufin[i];
       TB[f].pos_in += data_a_lire;
       TB[f].nb_elmt_Buf_IN -= data_a_lire;
       data_a_lire = 0;
     }
     
     else {
       for (i=0; i < TB[f].nb_elmt_Buf_IN; i++)
	 p_bucket[i] = p_bufin[i];
       p_bucket += TB[f].nb_elmt_Buf_IN;
       p_bufin += TB[f].nb_elmt_Buf_IN;
       data_a_lire -= TB[f].nb_elmt_Buf_IN;
       lseek(fd,TB[f].nb_elmt_lu * sizeof(ftnword),SEEK_SET);
       nbytes = read(fd,TB[f].pBuf_IN,BUF_SIZE);
       if (nbytes <= 0) {
	 fprintf(stderr,"read_db_file error: try to read past end of file\n");
	 return(-1);
       }
       TB[f].nb_elmt_lu += nbytes/sizeof(ftnword);
       TB[f].nb_elmt_Buf_IN = nbytes/sizeof(ftnword);
       TB[f].pos_in = 0;
     }

   }
   
   return (TB[f].nb_elmt_lu) ;

}/* end read_db_file */



/****************************************************************************
 *                        W R I T E _ D B _ F I L E                         *
 *                                                                          *
 * OBJECT:                                                                  *
 *        Write bucket in the file                                          *
 *                                                                          *
 * ARGUMENTS:                                                               *
 *           IN/ NB       : number of elements to write from bucket.        *
 *           IN/ bucket   : data to write.                                  *
 *       IN-OUT/ fd       : file descriptor.                                *
 *                                                                          *
 ****************************************************************************/
 
ftnword f77name(write_db_file)(ftnword *iun,ftnword *bucket,ftnword *NB)
{
   int data_a_ecrire = *NB;
   int i, nbytes, room_left, f, fd;
   ftnword *p_bucket, *p_bufout;

   f = -1;
   for (i=0; i < Max_File; i++)
     if (TB[i].iun == *iun) {
       f = i;
       break;
     }
   if (f == -1) {
     fprintf(stderr,"write_db_file error: file (unit=%d) not open\n",*iun);
     exit(7);
   }

   fd = TB[f].fd;

   if(((TB[f].nb_elmt_ecrit + data_a_ecrire + TB[f].pos_out) >=
       (TB[f].nb_elmt_lu + TB[f].pos_in)) && (TB[f].nb_elmt_lu != 0)) {
     fprintf(stderr,"write_db_file error: rewriting over unread data\n");
     exit(1);
   }

   room_left = BUF_SIZE_W - TB[f].pos_out;   /* room left in ouput buffer */
   p_bucket = &(bucket[0]);            /* position into bucket */

   while (data_a_ecrire) {
     
     p_bufout = &(TB[f].pBuf_OUT[TB[f].pos_out]);
     
     if (data_a_ecrire <= room_left) {
       for (i=0; i < data_a_ecrire; i++)
	 p_bufout[i] = p_bucket[i];
       p_bufout += data_a_ecrire;
       TB[f].pos_out += data_a_ecrire;
       room_left -= data_a_ecrire;
       data_a_ecrire = 0;
     }

     else {
       for (i=0; i < room_left; i++) 
	 p_bufout[i] = p_bucket[i];
       p_bufout += room_left;
       p_bucket += room_left;
       lseek(fd,TB[f].nb_elmt_ecrit * sizeof(ftnword),SEEK_SET);
       nbytes = write(fd,TB[f].pBuf_OUT,BUF_SIZE);
       if (nbytes != BUF_SIZE) {
	 fprintf(stderr,"write_db_file error: can't write of file\n");
	 exit(2);
       }
       TB[f].nb_elmt_ecrit += BUF_SIZE_W;
       data_a_ecrire -= room_left;
       TB[f].pos_out = 0;
       room_left = BUF_SIZE_W;
     }

   }

   return (TB[f].nb_elmt_ecrit);
}


/****************************************************************************
 *                     R E W I N D _ D B _ F I L E                          *
 *                                                                          *
 * OBJECT:                                                                  *
 *        Clear buffers and return in the initial position of the file and  *
 *        the two buffer.                                                   *
 * ARGUMENTS:                                                               *
 *           IN/ fd               : file descripter.                        *
 *                                                                          *
 ****************************************************************************/

ftnword f77name(rewind_db_file)(ftnword *iun)
{
   int i, fd, nbytes, f;

   f = -1;
   for (i=0; i < Max_File; i++)
     if (TB[i].iun == *iun) {
       f = i;
       break;
     }
   if (f == -1) {
     fprintf(stderr,"rewind_db_file error: file (unit=%d) not open\n",*iun);
     exit(7);
   }

   fd = TB[f].fd;

   TB[f].nb_elmt_lu = 0;

   if (TB[f].pos_out > 0) {
     lseek(fd,TB[f].nb_elmt_ecrit * sizeof(ftnword),SEEK_SET);
     nbytes = write(fd,TB[f].pBuf_OUT,TB[f].pos_out*sizeof(ftnword));
     if (nbytes != TB[f].pos_out*sizeof(ftnword)) {
       fprintf(stderr,"rewind_db_file error: can't write of file\n");
       exit(2);
     }
   }
   lseek(fd,0L,SEEK_SET);                   /* rewind the file.              */
   
   for(i=0; i < BUF_SIZE_W; i++)/* Initialiser les buffer a zero.*/
     {
       TB[f].pBuf_IN[i] = 0;
       TB[f].pBuf_OUT[i] = 0;
     }
   TB[f].nb_elmt_ecrit = 0;
   TB[f].nb_elmt_Buf_IN = 0;
   TB[f].pos_out = 0;
   TB[f].pos_in = 0;
   return(0);
}

/****************************************************************************
 *                        C L O S E _ D B _ F I L E                         *
 *                                                                          *
 * Object :                                                                 *
 *    Close file and flush output buffer.                                   *
 *                                                                          *
 * Arguments :                                                              *
 *             IN/  IUN : fortran unit number                               *
 *                                                                          *
 ****************************************************************************/

ftnword f77name(close_db_file)(ftnword *iun)                                    
{
  int ier, f, i;
 
  f = -1;
  for (i=0; i < Max_File; i++)
    if (TB[i].iun == *iun) {
      f = i;
      break;
    }
  if (f == -1) {
    fprintf(stderr,"close_db_file error: file (unit=%d) not open\n",*iun);
    exit(7);
  }

  ier = f77name(rewind_db_file)(iun);

  TB[f].fd = -1;
  TB[f].iun = -1;
  TB[f].pos_in = 0;
  TB[f].pos_out = 0;
  TB[f].nb_elmt_lu = 0;
  TB[f].nb_elmt_ecrit = 0;
  free(TB[f].pBuf_IN);
  free(TB[f].pBuf_OUT);
  TB[f].pBuf_IN = NULL;
  TB[f].pBuf_OUT = NULL;
  return(0);
}
