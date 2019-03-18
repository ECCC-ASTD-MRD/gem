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

#include "ezscint.h"
#include "ez_funcdef.h"


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezget_subgridids)(wordint *gdid, wordint *subgrid)
{
   return c_ezget_subgridids(*gdid, subgrid);
}

wordint c_ezget_subgridids(wordint gdid, wordint *subgrid)
{
  wordint i,gdrow_id, gdcol_id;

  c_gdkey2rowcol(gdid, &gdrow_id, &gdcol_id);
  if (Grille[gdrow_id][gdcol_id].nsubgrids == 0) 
    {
    *subgrid=gdid;
    return 1;  
    }
  for (i=0; i<Grille[gdrow_id][gdcol_id].nsubgrids; i++)
   {
   subgrid[i]=Grille[gdrow_id][gdcol_id].subgrid[i];
   }
  return Grille[gdrow_id][gdcol_id].nsubgrids;
}
