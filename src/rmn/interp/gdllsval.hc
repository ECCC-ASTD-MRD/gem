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
wordint f77name(gdllsval)(wordint *gdid, ftnfloat *zout, ftnfloat *zin, ftnfloat *lat, ftnfloat *lon, wordint *n)
{
   wordint icode;
   
   icode = c_gdllsval(*gdid, zout, zin, lat, lon, *n);
   return icode;
}

wordint c_gdllsval(wordint gdid, ftnfloat *zout, ftnfloat *zin, ftnfloat *lat, ftnfloat *lon, wordint n)
{
   ftnfloat *x, *y;
   wordint ier,gdrow_id,gdcol_id;
   
  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);

   x = (ftnfloat *)malloc(n * sizeof(float));
   y = (ftnfloat *)malloc(n * sizeof(float));
   
   if (Grille[gdrow_id][gdcol_id].nsubgrids > 0 )
      {
         ier = c_gdxyfll(gdid, x, y, lat, lon, n);
      }
   else
      {
         ier = c_gdxyfll_orig(gdid, x, y, lat, lon, n);
      }
   ier = c_gdxysval(gdid, zout, zin, x, y, n);
   
   free(x);
   free(y);
   return 0;
}

