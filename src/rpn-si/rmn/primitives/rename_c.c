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

ftnword f77name(rename_c)(char *oldname, char *newname, F2Cl lng1, F2Cl lng2)
{
  int rcode;
  char old[256], new[256];

  if (lng1 > 256 || lng2 > 256) {
    printf("rename_c error: oldname or newname > 256 char\n");
    return((ftnword) -1);
  }
  while (oldname[lng1-1] == ' ' && lng1 > 0) lng1--;
  while (newname[lng2-1] == ' ' && lng2 > 0) lng2--;
  strncpy(old,oldname,lng1);
  old[lng1] = '\0';
  strncpy(new,newname,lng2);
  new[lng2] = '\0';
  rcode = rename(&old[0],&new[0]);
  if (rcode == -1) perror("rename_c error");
  return((ftnword)rcode);
}
