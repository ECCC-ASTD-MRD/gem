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
wordint f77name(ezgdef)(wordint *ni, wordint *nj, char *grtyp, char *grref, 
			wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, 
			ftnfloat *ax, ftnfloat *ay, F2Cl lengrtyp, F2Cl lengrref)
{
  wordint icode;
  char lgrtyp[2],lgrref[2];
  
  lgrtyp[0] = grtyp[0];
  lgrtyp[1] = '\0';
  
  lgrref[0] = grref[0];
  lgrref[1] = '\0';
  
  icode = c_ezgdef(*ni, *nj, lgrtyp, lgrref, *ig1, *ig2, *ig3, *ig4, ax, ay);
  return icode;
}

wordint c_ezgdef(wordint ni, wordint nj, char *grtyp, char *grref,
		 wordint ig1, wordint ig2, wordint ig3, wordint ig4, ftnfloat *ax, ftnfloat *ay)
{
  wordint found,source;
  char typeGrille;
  

  found = -1;
  typeGrille = grtyp[0];
  
  if (grtyp[0] == '#') 
    {
    fprintf(stderr, "The '#' grid type is not supported with ezgdef.\nPlease use ezgdef_ffile or ezgdef_fmem\n");
    return -1;
    }

  switch(typeGrille)
    {
    case 'Y':
    case 'Z':
    case '#':
      if ((0 == strcmp(grref, "FILE")) || (0 == strcmp(grref, "file")))
	{
	source = FICHIER;
	}
      else
	{
	source = MEMOIRE;
	}       
      break;
      
    default:
      source = MEMOIRE;
      break;
    }
  
  switch (source)
    {
    case MEMOIRE:
      found = c_ezgdef_fmem(ni, nj, grtyp, grref, ig1, ig2, ig3, ig4, ax, ay);
      break;
      
    case FICHIER:
      found = c_ezgdef_ffile(ni, nj, grtyp, ig1, ig2, ig3, ig4, ig4);
      break;
    }
  
  return found;
}

