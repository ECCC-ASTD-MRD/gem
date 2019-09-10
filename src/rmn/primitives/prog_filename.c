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
 /* prog_filename.c
   Function for FORTRAN to create a filename given some parameters
   to give the form of "ppYYYYMMDDhh[mnss-XX-YY]_ddd[U]"

   original author: Vivian Lee
   revision 001 : M.Valin 2001/01/31 minor bugfix

   How to use:

      integer prog_filename
      external prog_filename
      character*50 filename
      character*2 prefix,unit
      integer ierr,date,hour,min,sec,npex,npey,num,numlen

   ierr=prog_filename(filename,prefix,date,hour,min,sec,npex,npey,num,numlen,unit)

   Arguments:

   ierr     - O - return code of function; 0=no error, -1=error

   prefix   - I - 2 alphabetic characters (A-Z,a-z)
   date     - I - 8 digit number from newdate, ie:20010120
   hour     - I - number between 00 - 24
   min      - I - minutes; 00 - 59, if -1, omitted
   sec      - I - seconds; 00 - 59, if -1, omitted
   npex     - I - PE row;  00 - 99, if -1, omitted
   npey     - I - PE col;  00 - 99, if -1, omitted
   num      - I - multiple of unit;  >= 0
   numlen   - I - number of digits for num, optional
   unit     - I - 1 character to represent type of units, optional

   filename - O - string of characters in this format:
                  ppYYYYMMDDhh[mnss-XX-YY]_ddd[U]
                  Anything surrounded by square brackets are optional
                  Some examples:
                       md20010104123059-05-01_00000000
                       pr2001010412_00000000
                       pr2000010500_024H

                  pp       - 2 lowercase alphabetic characters for prefix
                  YYYYMMDD - year+month+day (from newdate)
                  hh       - hour
                  mn       - minute
                  ss       - second
                  XX       - PE on row XX
                  YY       - PE on column YY
                  ddd      - is a count/multiple of units in the integration
                  U        - type of unit (1 character)
                             suggested representations for unit:
                             P - timesteps
                             D - days
                             H - hours
                             M - months
                             S - seconds
*/

#include <rpnmacros.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef WIN32	/*CHC/NRC*/
#else
#include <unistd.h>
#endif

#include <ctype.h>

#define FNAME_LEN 50

ftnword f77name(prog_filename)(char *f_nom, unsigned char *f_prefix, 
       ftnword *f_date, ftnword *f_hour, ftnword *f_min, ftnword *f_sec, 
       ftnword *f_npex, ftnword *f_npey, ftnword *f_num, ftnword *f_numlen, 
       char *f_unit, F2Cl l1, F2Cl l2, F2Cl l3)

{                                               

 int i,len,err,numlen,c1,c2;
 char format[20], units[2], fname[FNAME_LEN];
 unsigned char prefix[3];
 char *s;

 for (i=0;i<l1;i++) f_nom[i]=' ';
 s="Bad_Filename";
 /* return a bad file name in case of error */
 for (i=0;i<l1 && *s;i++) f_nom[i]=*s++;

 if (l1 < 16) { /* length of f_nom is less than 16 */
     fprintf(stderr,"prog_filename: length of output filename is less than 16\n");
     return(-1);
     }

 if (l2 < 2) {  /* length of prefix is less than 2 */
     fprintf(stderr,"prog_filename: length of prefix is less than 2\n");
     return(-1);
     }
 c1=f_prefix[0];
 c2=f_prefix[1];
 if (isalpha(c1) && isalpha(c2)) {
     prefix[0] = tolower(f_prefix[0]);
     prefix[1] = tolower(f_prefix[1]);
     prefix[2] = '\0';
     }
 else {
     fprintf(stderr,"prog_filename: prefix contains improper characters\n");
     return(-1);
     }
 if (*f_date<0 || *f_date>99999999){
     fprintf(stderr,"prog_filename: date<0 or date>99999999\n");
     return(-1);
     }
 if (*f_hour<0 || *f_hour>23) {
     fprintf(stderr,"prog_filename: hour<0 or hour>23\n");
     return(-1);
     }

 numlen=*f_numlen;
 if (*f_numlen<3 || *f_numlen>9) numlen=3;

 if (*f_num<0) {
     fprintf(stderr,"prog_filename: num<0\n");
     return(-1);
     }

 units[0] = (l3>=1 && *f_unit != ' ') ? tolower(*f_unit) : '\0' ;
 units[1] = '\0';


 snprintf(fname,FNAME_LEN,"%s%08d%02d",prefix,*f_date,*f_hour);
 for (s=fname ; *s ; s++);
 len=s-fname;

 if(*f_min!=-1 && *f_sec!=-1) {
    if (*f_min <0 || *f_min >59) {
         fprintf(stderr,"prog_filename: minutes<0 or minutes>59\n");
         return(-1);
         }
    if (*f_sec <0 || *f_sec >59) {
         fprintf(stderr,"prog_filename: seconds<0 or seconds>59\n");
         return(-1);
         }
    snprintf(s,FNAME_LEN-len,"%02d%02d",*f_min,*f_sec);
    while(*s) s++;
    len=s-fname;
    }

 if(*f_npex!=-1 && *f_npey!=-1) {
    if (*f_npex <0 || *f_npex >99) {
         fprintf(stderr,"prog_filename: npex<0 or npex>99\n");
         return(-1);
         }
    if (*f_npey <0 || *f_npey >99) {
         fprintf(stderr,"prog_filename: npey<0 or npey>99\n");
         return(-1);
         }
    snprintf(s,FNAME_LEN-len,"-%02d-%02d",*f_npex,*f_npey);
    while(*s) s++;
    len=s-fname;
    }

 strncpy(format,"_%03d%s",19);
 format[3] = '0' + numlen; /* replace the 3 with the proper length */

 snprintf(s,FNAME_LEN-len,format,*f_num,units);
 while(*s) s++;
 len=s-fname;

 for(i=0; i<len && i<l1 ; i++) f_nom[i]=fname[i];
 while (i<l1) f_nom[i++] = ' ';
    /*printf("f_nom=%s length=%d\n",fname,strlen(f_nom));*/
 return(0);
}
