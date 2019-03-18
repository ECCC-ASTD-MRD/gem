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
wordint f77name(ezgdef_fll)(wordint *ni, wordint *nj, ftnfloat *lat, ftnfloat *lon)
{
  wordint gdid;

  gdid = c_ezgdef_fll(*ni, *nj, lat, lon);
  return gdid;
}

wordint c_ezgdef_fll(wordint ni, wordint nj,ftnfloat *lat, ftnfloat *lon)
   {
   wordint gdid;
   wordint gdrow_id, gdcol_id;

   ftnfloat swlat, swlon, dlat, dlon;
   wordint ig1, ig2, ig3, ig4;
   char grtyp[2], grref[2];

   swlat = 0.0;
   swlon = 0.0;
   dlat  = 1.0;
   dlon  = 1.0;
   strcpy(grtyp, "Y");
   strcpy(grref, "L");

   f77name(cxgaig)(grref, &ig1, &ig2, &ig3, &ig4, &swlat, &swlon, &dlat, &dlon);
   gdid = c_ezgdef_fmem(ni, nj, grtyp, grref, ig1, ig2, ig3, ig4, lon, lat);

   return gdid;
}

