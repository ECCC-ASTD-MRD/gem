/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2007  Environnement Canada
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
#include <stdlib.h>

//! FORTRAN style trimmed source to arbitrarily padded destination copy
//! @param src Source string
//! @param dest Destination string
//! @param lsrc Length of source string
//! @param ldest Length of destination string
//! @param pad Character used for padding
//!
//! @return Trimmed source length
//!
//! If destination pointer is NULL, no copy takes place.
int ftn2c_string_copy(unsigned char* src, unsigned char* dest, int lsrc, int ldest, unsigned char pad)
{
   int i;
   /* if there is a null before lsrc characters in the source string, act as if it was terminated at null */
   for (i = 0 ; src[i] != 0 && i < lsrc; i++);
   lsrc = i;
   while (src[lsrc - 1] == ' ') lsrc--;  /* ignore trailing blanks in source string */
   if (dest == NULL) return(lsrc);     /* no destination, just return trimmed source length */
   if(lsrc > ldest) return(-1);  /* OOPS, trimmed source longer than destination */
   if (pad == 0 && lsrc == ldest) return(-1);  /* OOPS, not enough space for padding */
   for (i = 0; i < lsrc; i++) dest[i] = src[i] ;  /* copy src to dest */
   if (pad)
      while(i < ldest) dest[i++] = pad;          /* pad destination */
   else
      dest[i] = pad;
   return(lsrc) ; /* return number of significant characters copied */
}

//! C String array to Fortran String Array
int ftn2c_cstra_fstra(unsigned char** src, unsigned char* dest, int lsrc, int ldest, int nitems, unsigned char pad) {
   int ii;
   if (nitems <= 0) return(-1);
   for (ii = 0; ii < nitems; ii++) {
      if (ftn2c_string_copy(src[ii], dest, lsrc, ldest, pad) < 0) return(-1);
      dest += ldest;
   }
   return(0);
}

//! Fortran String array to C String Array
int ftn2c_fstra_cstra(unsigned char *src, unsigned char **dest, int lsrc, int ldest, int nitems, unsigned char pad) {
   int ii;
   if (nitems <= 0) return(-1);
   for (ii = 0; ii < nitems; ii++) {
      if (ftn2c_string_copy(src, dest[ii], lsrc, ldest, pad) < 0) return(-1);
      src += lsrc;
   }
   return(0);
}
