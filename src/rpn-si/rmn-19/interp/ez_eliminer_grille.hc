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

#ifdef MUTEX
//JP
extern pthread_mutex_t EZ_MTX;
#endif

void EliminerGrille(wordint gdid)
{
  wordint i, index;
  wordint gdrow_id, gdcol_id;
    
#ifdef MUTEX
// JP
   pthread_mutex_lock(&EZ_MTX);
#endif

  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
  
   if (Grille[gdrow_id][gdcol_id].access_count > 0)
    {
    Grille[gdrow_id][gdcol_id].access_count--;
    }
   
   if (Grille[gdrow_id][gdcol_id].access_count == 0)
    {
    if (Grille[gdrow_id][gdcol_id].flags & LAT)
        {
        if (Grille[gdrow_id][gdcol_id].lat != NULL)
           {
           free(Grille[gdrow_id][gdcol_id].lat);
           free(Grille[gdrow_id][gdcol_id].lon);
           }
        }

    if (Grille[gdrow_id][gdcol_id].flags & AX)
        {
        if (Grille[gdrow_id][gdcol_id].ax != NULL)
           {
           free(Grille[gdrow_id][gdcol_id].ax);
           free(Grille[gdrow_id][gdcol_id].ay);
           }
        }        

    if (Grille[gdrow_id][gdcol_id].ncx != NULL)
        {
        if (Grille[gdrow_id][gdcol_id].ncx != NULL)
           {
           free(Grille[gdrow_id][gdcol_id].ncx);
           free(Grille[gdrow_id][gdcol_id].ncy);
           }
        }
 
    for (i=0; i < Grille[gdrow_id][gdcol_id].n_gdin_for; i++)
        {
        index = ez_find_gdin_in_gset(gdid, Grille[gdrow_id][gdcol_id].gdin_for[i]);
        c_ezfreegridset(Grille[gdrow_id][gdcol_id].gdin_for[i], index);
        }
        
      gr_list[Grille[gdrow_id][gdcol_id].grid_index] = (_Grille *) NULL;
      memset(&(Grille[gdrow_id][gdcol_id]),(int)NULL, sizeof(_Grille));      
      nGrilles--;
    }
   

#ifdef MUTEX
// JP
   pthread_mutex_unlock(&EZ_MTX);
#endif

   }
