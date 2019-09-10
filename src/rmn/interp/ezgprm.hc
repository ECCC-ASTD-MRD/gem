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
wordint f77name(ezgprm)(wordint *gdid, char *grtyp, wordint *ni, wordint *nj,
             wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, F2Cl lengrtyp)
{
  int i;
  wordint icode;
  char lgrtyp[2];
  lgrtyp[0]= ' ';
  lgrtyp[1]= '\0';
  icode = c_ezgprm(*gdid, lgrtyp, ni, nj, ig1, ig2, ig3, ig4);
  grtyp[0]=lgrtyp[0];
  return icode;
}

wordint   c_ezgprm(wordint gdid, char *grtyp, wordint *ni, wordint *nj, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4)
{
  wordint gdrow_id, gdcol_id;
    
  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
   
   *grtyp = Grille[gdrow_id][gdcol_id].grtyp[0];
   *ni    = Grille[gdrow_id][gdcol_id].ni;
   *nj    = Grille[gdrow_id][gdcol_id].nj;
   *ig1   = Grille[gdrow_id][gdcol_id].fst.ig[IG1];
   *ig2   = Grille[gdrow_id][gdcol_id].fst.ig[IG2];
   *ig3   = Grille[gdrow_id][gdcol_id].fst.ig[IG3];
   *ig4   = Grille[gdrow_id][gdcol_id].fst.ig[IG4];

   return 0;
}
