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
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>



ftnword f77name(bmf_char2i)(char *f_char, ftnword *f_length,
			    ftnword *f_outint, ftnword *f_outlen, F2Cl l1)

{
  int length = *f_length, outlen = *f_outlen ;
  if (length == 0) {
    fprintf(stderr,"bmf_char2i WARNING: char length = 0, exiting \n");
    return(0);
  }
  /*  if (length > l1) {
    fprintf(stderr,"bmf_char2i: length too big for char variable \n");
    return(-1);
    }*/

   if (outlen*sizeof(outlen) < length*sizeof(f_char[0])+sizeof(outlen)) {
    fprintf(stderr,"bmf_char2i: integer array size too small \n");
    return(-1);
  }
   f_outint[0]= length;
   f_outint++;
   strncpy((char *)f_outint,f_char,length);
   f_outint--;
   return((length-1)/sizeof(length)+2);
}
ftnword f77name(bmf_i2char)(char *f_char, ftnword *f_length,
			    ftnword *f_outint, ftnword *f_outlen, F2Cl l1)

{
  int length= *f_length , outlen=*f_outlen;
 
  /* if (length > l1) {
    fprintf(stderr,"bmf_i2char: length too big for char variable \n");
    return(-1);
    }*/
 
   if (outlen*sizeof(outlen) < length*sizeof(f_char[0])+sizeof(outlen)) {
    fprintf(stderr,"bmf_i2char: integer array size too small \n");
    return(-1);
  }
   f_length = (long *)f_outint[0];
   f_outint++;
   strncpy(f_char,(char *)f_outint,length);
   f_outint--;
   return(length);
}
