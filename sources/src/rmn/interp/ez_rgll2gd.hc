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

#include <stdio.h>
#include <rpnmacros.h>

f77name(ez_rgll2gd)(ftnfloat *z1, ftnfloat *z2, ftnfloat *xlon, wordint *ni, wordint *nj, 
                    char *grtyp, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4,F2Cl lengrtyp)
{
  ftnfloat *tmplat;
  wordint gdid, npts,n;
  char lgrtyp[2];
  
  npts = *ni * *nj;
  tmplat = (float *)malloc(npts*sizeof(ftnfloat));
  for (n=0; n < npts; n++)
    {
    tmplat[n] = 0.0;
    }

  ftnstrclean(grtyp,lengrtyp);
  strcpy(lgrtyp, grtyp);

  gdid  = c_ezqkdef(*ni,*nj, lgrtyp, *ig1,  *ig2,  *ig3,  *ig4, 0);

  c_gduvfwd(gdid, z1, z2, z1, z2, tmplat, xlon, npts);
  
  
  free(tmplat);
  return 0;
}
