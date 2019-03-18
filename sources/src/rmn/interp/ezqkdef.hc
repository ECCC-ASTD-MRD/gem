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
wordint f77name(ezqkdef)(wordint *ni, wordint *nj, char *grtyp, 
                    wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, wordint *iunit, F2Cl lengrtyp)
{
   wordint icode;
   char cgrtyp[2];

   cgrtyp[0] = grtyp[0];
   cgrtyp[1] = '\0';
   icode = c_ezqkdef(*ni, *nj, cgrtyp, *ig1, *ig2, *ig3, *ig4, *iunit);
   return icode;
}

wordint c_ezqkdef(wordint ni, wordint nj, char *grtyp,
             wordint ig1, wordint ig2, wordint ig3, wordint ig4, wordint iunit)
{
  wordint icode;

  icode = c_ezgdef_ffile(ni, nj, grtyp, ig1, ig2, ig3, ig4, iunit);
  
  return icode;
}
