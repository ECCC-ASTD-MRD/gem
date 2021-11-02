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
wordint ez_defzones(wordint gdin, wordint gdout)
{
wordint i;
wordint extrap;
int lcl_ngdin;

wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out, npts, idx_gdin;
   
c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
idx_gdin = c_find_gdin(gdin, gdout);

if (Grille[gdrow_out][gdcol_out].gset[idx_gdin].flags & ZONES)
      {
      return 0;
      }
   

npts = Grille[gdrow_out][gdcol_out].ni * Grille[gdrow_out][gdcol_out].nj;
extrap = EZ_NO_EXTRAP;
switch (Grille[gdrow_in][gdcol_in].grtyp[0])
   {
   case 'N':
   case 'S':
   case '!':
      extrap = EZ_EXTRAP;
   break;
      
   case 'L':
      if (Grille[gdrow_out][gdcol_out].extension == 0)
         {
         extrap = EZ_EXTRAP;   
         }
      else
         {
         extrap = EZ_NO_EXTRAP;
         }

   case '#':
   case 'Z':
   case 'Y':
      switch(Grille[gdrow_in][gdcol_in].grref[0])
         {
         case 'N':
         case 'S':
         extrap = EZ_EXTRAP;
         break;
         
         case 'E':
         case 'L':
         if (358.0 > (Grille[gdrow_in][gdcol_in].ax[Grille[gdrow_in][gdcol_in].ni-1] - Grille[gdrow_in][gdcol_in].ax[0]))
            {
            extrap = EZ_EXTRAP;
            }
         break;
         }
      break;
   }

   for (i=0; i < NZONES; i++)
      {
      Grille[gdrow_out][gdcol_out].gset[idx_gdin].zones[i].npts = 0;
      }

   switch (extrap)
      {
      case EZ_EXTRAP:
         ez_defzone_dehors(gdin, Grille[gdrow_out][gdcol_out].gset[idx_gdin].x, 
               Grille[gdrow_out][gdcol_out].gset[idx_gdin].y, npts, 
               &(Grille[gdrow_out][gdcol_out].gset[idx_gdin].zones[DEHORS]));
         break;
         
      case EZ_NO_EXTRAP:
         ez_defzone_polenord(gdin, Grille[gdrow_out][gdcol_out].gset[idx_gdin].x, 
               Grille[gdrow_out][gdcol_out].gset[idx_gdin].y, npts, 
               &(Grille[gdrow_out][gdcol_out].gset[idx_gdin].zones[POLE_NORD]));
         ez_defzone_polesud(gdin, Grille[gdrow_out][gdcol_out].gset[idx_gdin].x, 
               Grille[gdrow_out][gdcol_out].gset[idx_gdin].y, npts, 
               &(Grille[gdrow_out][gdcol_out].gset[idx_gdin].zones[POLE_SUD]));
         ez_defzone_sud(gdin, Grille[gdrow_out][gdcol_out].gset[idx_gdin].x, 
               Grille[gdrow_out][gdcol_out].gset[idx_gdin].y, npts, 
               &(Grille[gdrow_out][gdcol_out].gset[idx_gdin].zones[AU_SUD]));
         ez_defzone_nord(gdin, Grille[gdrow_out][gdcol_out].gset[idx_gdin].x, 
               Grille[gdrow_out][gdcol_out].gset[idx_gdin].y, npts, 
               &(Grille[gdrow_out][gdcol_out].gset[idx_gdin].zones[AU_NORD]));
      }

   Grille[gdrow_out][gdcol_out].gset[idx_gdin].flags |= ZONES;
   return 0;
   }
