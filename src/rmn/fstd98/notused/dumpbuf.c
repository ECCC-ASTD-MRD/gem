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

#include "qstdir.h"
#include <stdio.h>
void c_dumpbuf(buffer_interface_ptr buf, int n)
{
int i;      

      fprintf(stdout,"\n");
      fprintf(stdout," <<<  D U M P B U F  >>>\n");
      fprintf(stdout,"  buf->nwords =      %d\n",buf->nwords);
      fprintf(stdout,"  buf->nbits  =      %d\n",buf->nbits);
      fprintf(stdout,"  buf->data_index =  %d\n",buf->data_index);
      fprintf(stdout,"  buf->record_index =%d\n",buf->record_index);
      fprintf(stdout,"  buf->iun    =      %d\n",buf->iun);
      fprintf(stdout,"  buf->aux_index =   %d\n",buf->aux_index);
      fprintf(stdout,"  buf->buf9   =      %d\n",buf->buf9);
      fprintf(stdout,"\n");

      for (i=0; i < n; i+=5)
        fprintf(stdout," %08X  %08X  %08X  %08X  %08X\n",buf->data[i],
        buf->data[i+1],buf->data[i+2],buf->data[i+3],buf->data[i+4]);
      return;
}

void f77name(dumpbuf)(word *buf, ftnword *n)
{
int nn = *n;

#if defined(NEC64)
   BUF_C;
     c_dumpbuf(buf+1,nn);
   BUF_F;
#else
   c_dumpbuf(buf,nn);
#endif
return;
}
