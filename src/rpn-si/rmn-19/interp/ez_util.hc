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


wordint f77name(ezsetgdout)(int *gdout)
  {
  return c_ezsetgdout(*gdout);
  }

wordint f77name(ezgetgdout)()
  {
  return c_ezgetgdout();
  }

wordint c_ezgetgdin()
  {
  return iset_gdin; 
  }

wordint c_ezgetgdout()
  {
  if (iset_gdout == UNDEFINED)
     {
     if (iset_gdin != UNDEFINED)
        {
	iset_gdout = iset_gdin;
	}
     else
       {
       iset_gdout = UNDEFINED;
       }
     }
  return iset_gdout;
  }

wordint c_ezsetgdout(gdout)
wordint gdout;
  {
  iset_gdout = gdout;
  return 0;
  }
