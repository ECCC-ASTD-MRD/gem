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

#include <stdio.h>
#include <rpnmacros.h>


/*****************************************************************************
 *                                                                           * 
 *               F R E A D / F W T I T E  1 6 - 3 2 - 6 4                    * 
 *                                                                           *
 *             FREAD and FWRITE interface with bytes swap                    *
 *                                                                           *
 *Author                                                                     *
 *  M. Lepine - April 1999                                                   *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/



static int endian_int=1;
static char *little_endian=(char *)&endian_int;


/***************************************************************************** 
 *                            F R E A D 1 6                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads nitems elements of data, each size bytes long                     *
 *   and swap 8 bits by 8 bits of each 2 bytes elements                      *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  ptr     pointer to array to receive  data                            * 
 *  IN  size    size in bytes of elements of data                            * 
 *  IN  nitems  number of items to read                                      * 
 *  IN  stream  file pointer                                                 * 
 *                                                                           * 
 *****************************************************************************/
size_t fread16(void *ptr, size_t size, size_t nitems, FILE *stream)
{
  size_t nr;
  int i, n2=(size*nitems)/2;    /* number of 2 bytes */
  unsigned short *pt2 = (unsigned short *) ptr;

  if (*little_endian) {
    if ((size & 1) != 0) {
      fprintf(stderr,"fread16 error: size=%d must be a multiple of 2\n",size);
      return(-1);
    }
    
    nr = fread(ptr,size,nitems,stream);
    
    for (i=0; i < n2; i++) {
      *pt2 = (*pt2 >> 8) | (*pt2 << 8);
      pt2++;
    }
  }
  else
    nr = fread(ptr,size,nitems,stream);
  return((size_t) nr);
}

/***************************************************************************** 
 *                            F R E A D 3 2                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads nitems elements of data, each size bytes long                     *
 *   and swap each bytes for each 4 bytes elements                           *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  ptr     pointer to array to receive  data                            * 
 *  IN  size    size in bytes of elements of data                            * 
 *  IN  nitems  number of items to read                                      * 
 *  IN  stream  file pointer                                                 * 
 *                                                                           * 
 *****************************************************************************/
size_t fread32(void *ptr, size_t size, size_t nitems, FILE *stream)
{
  size_t nr;
  int i, n4=(size*nitems)/4;    /* number of 4 bytes */
  unsigned INT_32 *pt4 = (unsigned INT_32 *) ptr;

  if (*little_endian) {
    if ((size & 3) != 0) {
      fprintf(stderr,"fread64 error: size=%d must be a multiple of 4\n",size);
      return(-1);
    }
    
    nr = fread(ptr,size,nitems,stream);
    
    for (i=0; i < n4; i++) {
      *pt4 = (*pt4>>24) | (*pt4<<24) | ((*pt4>>8)&0xFF00) | ((*pt4&0xFF00)<<8);
      pt4++;
    }
  }
  else
    nr = fread(ptr,size,nitems,stream);

  return((size_t) nr);
}

/***************************************************************************** 
 *                            F R E A D 6 4                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Reads nitems elements of data, each size bytes long                     *
 *   and swap each bytes of each 32 bits elements and then swap every        *
 *   2 32 bits elements                                                      * 
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  ptr     pointer to array to receive  data                            * 
 *  IN  size    size in bytes of elements of data                            * 
 *  IN  nitems  number of items to read                                      * 
 *  IN  stream  file pointer                                                 * 
 *                                                                           * 
 *****************************************************************************/
size_t fread64(void *ptr, size_t size, size_t nitems, FILE *stream)
{
  size_t nr;
  int i, n4=(size*nitems)/4;    /* number of 4 bytes */
  unsigned INT_32 *pt4 = (unsigned INT_32 *) ptr;
  INT_32 temp;

  if (*little_endian) {
    if ((size & 3) != 0) {
      fprintf(stderr,"fread64 error: size=%d must be a multiple of 4\n",size);
      return(-1);
    }
    
    nr = fread(ptr,size,nitems,stream);
    
    for (i=0; i < n4; i++) {
      *pt4 = (*pt4>>24) | (*pt4<<24) | ((*pt4>>8)&0xFF00) | ((*pt4&0xFF00)<<8);
      pt4++;
    }
    
    pt4 = (unsigned INT_32 *) ptr;
    for (i=0; i < n4/2; i++) {
      temp = *pt4;
      *pt4 = *(pt4+1);
      pt4++;
      *pt4 = temp;
      pt4++;
    }
  }
  else
    nr = fread(ptr,size,nitems,stream);

  return((size_t) nr);
}

/***************************************************************************** 
 *                          F W R I T E 1 6                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Writes nitems elements of data, each size bytes long                    *
 *   and swap 8 bits by 8 bits of each 2 bytes elements                      *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  ptr     pointer to array to receive  data                            * 
 *  IN  size    size in bytes of elements of data                            * 
 *  IN  nitems  number of items to write                                     * 
 *  IN  stream  file pointer                                                 * 
 *                                                                           * 
 *****************************************************************************/
size_t fwrite16(void *ptr, size_t size, size_t nitems, FILE *stream)
{
  size_t nr;
  int i, n2=(size*nitems)/2;    /* number of 2 bytes */
  unsigned short *pt2 = (unsigned short *) ptr;

  if (*little_endian) {
    if ((size & 1) != 0) {
      fprintf(stderr,"fwrite16 error: size=%d must be a multiple of 2\n",size);
      return(-1);
    }
    
    for (i=0; i < n2; i++) {
      *pt2 = (*pt2 >> 8) | (*pt2 << 8);
      pt2++;
    }
    
    nr = fwrite(ptr,size,nitems,stream);
    
    pt2 = (unsigned short *) ptr;  
    for (i=0; i < n2; i++) {
      *pt2 = (*pt2 >> 8) | (*pt2 << 8);
      pt2++;
    }
  }
  else
    nr = fwrite(ptr,size,nitems,stream);

  return((size_t) nr);
}

/***************************************************************************** 
 *                          F W R I T E 3 2                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Writes nitems elements of data, each size bytes long                    *
 *   and swap each bytes for each 4 bytes elements                           *
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  ptr     pointer to array to receive  data                            * 
 *  IN  size    size in bytes of elements of data                            * 
 *  IN  nitems  number of items to write                                     * 
 *  IN  stream  file pointer                                                 * 
 *                                                                           * 
 *****************************************************************************/
size_t fwrite32(void *ptr, size_t size, size_t nitems, FILE *stream)
{
  size_t nr;
  int i, n4=(size*nitems)/4;    /* number of 4 bytes */
  unsigned INT_32 *pt4 = (unsigned INT_32 *) ptr;

  if (*little_endian) {
    if ((size & 3) != 0) {
      fprintf(stderr,"fwrite64 error: size=%d must be a multiple of 4\n",size);
      return(-1);
    }
    
    for (i=0; i < n4; i++) {
      *pt4 = (*pt4>>24) | (*pt4<<24) | ((*pt4>>8)&0xFF00) | ((*pt4&0xFF00)<<8);
      pt4++;
    }
    
    nr = fwrite(ptr,size,nitems,stream);
    
    pt4 = (unsigned INT_32 *) ptr;    
    for (i=0; i < n4; i++) {
      *pt4 = (*pt4>>24) | (*pt4<<24) | ((*pt4>>8)&0xFF00) | ((*pt4&0xFF00)<<8);
      pt4++;
    }
  }
  else
    nr = fwrite(ptr,size,nitems,stream);

  return((size_t) nr);
}

/***************************************************************************** 
 *                          F W R I T E 6 4                                  *
 *                                                                           * 
 *Object                                                                     * 
 *   Writes nitems elements of data, each size bytes long                    *
 *   and swap each bytes of each 32 bits elements and then swap every        *
 *   2 32 bits elements                                                      * 
 *                                                                           * 
 *Arguments                                                                  * 
 *                                                                           * 
 *  IN  ptr     pointer to array to receive  data                            * 
 *  IN  size    size in bytes of elements of data                            * 
 *  IN  nitems  number of items to write                                     * 
 *  IN  stream  file pointer                                                 * 
 *                                                                           * 
 *****************************************************************************/
size_t fwrite64(void *ptr, size_t size, size_t nitems, FILE *stream)
{
  size_t nr;
  int i, n4=(size*nitems)/4;    /* number of 4 bytes */
  unsigned INT_32 *pt4 = (unsigned INT_32 *) ptr;
  INT_32 temp;

  if (*little_endian) {
    if ((size & 3) != 0) {
      fprintf(stderr,"fwrite64 error: size=%d must be a multiple of 4\n",size);
      return(-1);
    }
    
    for (i=0; i < n4; i++) {
      *pt4 = (*pt4>>24) | (*pt4<<24) | ((*pt4>>8)&0xFF00) | ((*pt4&0xFF00)<<8);
      pt4++;
    }
    
    pt4 = (unsigned INT_32 *) ptr;
    for (i=0; i < n4/2; i++) {
      temp = *pt4;
      *pt4 = *(pt4+1);
      pt4++;
      *pt4 = temp;
      pt4++;
    }
    
    nr = fwrite(ptr,size,nitems,stream);
    
    pt4 = (unsigned INT_32 *) ptr;  
    for (i=0; i < n4; i++) {
      *pt4 = (*pt4>>24) | (*pt4<<24) | ((*pt4>>8)&0xFF00) | ((*pt4&0xFF00)<<8);
      pt4++;
    }
    
    pt4 = (unsigned INT_32 *) ptr;
    for (i=0; i < n4/2; i++) {
      temp = *pt4;
      *pt4 = *(pt4+1);
      pt4++;
      *pt4 = temp;
      pt4++;
    }
  }
  else
    nr = fwrite(ptr,size,nitems,stream);

  return((size_t) nr);
}
