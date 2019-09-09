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
main ()
{
   int handle=445466625;
   int index, addr, cluster, sign;
   
   index = handle & 0x7f;
   addr = (handle >> 7) & 0x3fffff;
   cluster = (handle >> 29) & 0x3;
   sign = (handle >> 31);
   fprintf(stdout,"[SEQ]\t sign=%d \t cluster=%d \t addr=%d \t index=%d\n",
	   sign,cluster,addr,index);
   addr = (handle >> 7) & 0x3fffff;
   cluster = (handle >> 29) & 0x3;
   fprintf(stdout,"[RND]\t sign=%d \t cluster=%d \t addr=%d \t index=%d\n",
	   sign,cluster,addr,index);
}
