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

void EliminerGrille(wordint gdid)
{
  wordint i, index;
  wordint gdrow_id, gdcol_id;
    
  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
  
   if (Grille[gdrow_id][gdcol_id].access_count > 0)
    {
    Grille[gdrow_id][gdcol_id].access_count--;
    }
   
   if (Grille[gdrow_id][gdcol_id].access_count == 0)
    {
    if (Grille[gdrow_id][gdcol_id].flags & LAT)
        {
        free(Grille[gdrow_id][gdcol_id].lat);
        free(Grille[gdrow_id][gdcol_id].lon);
        Grille[gdrow_id][gdcol_id].lat = NULL;
        Grille[gdrow_id][gdcol_id].lon = NULL;
        }

    if (Grille[gdrow_id][gdcol_id].flags & AX)
        {
        free(Grille[gdrow_id][gdcol_id].ax);
        free(Grille[gdrow_id][gdcol_id].ay);
        Grille[gdrow_id][gdcol_id].ax = NULL;
        Grille[gdrow_id][gdcol_id].ay = NULL;
        }

    if (Grille[gdrow_id][gdcol_id].ncx != NULL)
        {
        free(Grille[gdrow_id][gdcol_id].ncx);
        free(Grille[gdrow_id][gdcol_id].ncy);
        Grille[gdrow_id][gdcol_id].ncx = NULL;
        Grille[gdrow_id][gdcol_id].ncy = NULL;
        }
    Grille[gdrow_id][gdcol_id].flags = (int)0;
    }
   

   for (i=0; i < Grille[gdrow_id][gdcol_id].n_gdin_for; i++)
      {
      index = ez_find_gdin_in_gset(gdid, Grille[gdrow_id][gdcol_id].gdin_for[i]);
      c_ezfreegridset(Grille[gdrow_id][gdcol_id].gdin_for[i], index);
      }
   }
