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
#include <string.h>
#include <rpnmacros.h>
#define string_copy(dest,src,l) while(--l >= 0) dest[l]=src[l]

ftnword f77name(mkdir_c)(char *filename, F2Cl lng1)
{
  int rcode;
  char fname[4097];

  if (lng1 > 4096) {
    printf("mkdir_c error: file name > 4096 char\n");
    return((ftnword) -1);
  }
  while (filename[lng1-1] == ' ' && lng1 > 0) lng1--;
  strncpy(fname,filename,lng1);
  fname[lng1] = '\0';
  rcode = mkdir(fname,0777);
  if (rcode == -1) perror("mkdir_c error");
  return((ftnword)rcode);
}
