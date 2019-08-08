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
static int endian_int=1;
static char *little_endian=(char *)&endian_int;

void f77name(movr8r4)(double *r8, float *r4, ftnword *f_lng)
{
  int lng = *f_lng;

  while (lng--)
    *r4++ = *r8 ++;

}
void f77name(movr4r8)(float *r4, double *r8, ftnword *f_lng)
{
  int lng = *f_lng;

  while (lng--)
    *r8++ = *r4 ++;

}

void f77name(move6432)(unsigned INT_32 *src, unsigned INT_32 *dst, ftnword *f_lng)
{
  int lng = *f_lng;
  register unsigned INT_32 t1,t2;

  if (*little_endian) {
    while (lng--) {
      t1 = *src++;
      t2 = *src++;
      *dst++ = t2;
      *dst++ = t1;
    }
  }
  else {
    while (lng--) {
      t1 = *src++;
      t2 = *src++;
      *dst++ = t1;
      *dst++ = t2;
    }
  }
}
