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


wordint f77name(ezgenpole)(ftnfloat *vpolnor, ftnfloat *vpolsud, ftnfloat *fld,
                           wordint *ni, wordint *nj, wordint *vecteur, 
                           char *grtyp, wordint *hem, F2Cl lengrtyp)
{
   return c_ezgenpole(vpolnor, vpolsud, fld, *ni, *nj, *vecteur, grtyp, *hem);

}

wordint c_ezgenpole(ftnfloat *vpolnor, ftnfloat *vpolsud, ftnfloat *fld,
                           wordint ni, wordint nj, wordint vecteur, 
                           char *grtyp, wordint hem)
{
   wordint lni, lnj, lvecteur, lhem;
   ftnfloat *x, *y, *z, *lat, *lon, *gausslat;

   lni = ni; 
   lnj = nj; 
   lvecteur = vecteur; 
   lhem = hem;

   x =    (ftnfloat *) malloc(2*lni*sizeof(ftnfloat));
   y =    (ftnfloat *) malloc(2*lni*sizeof(ftnfloat));
   z =    (ftnfloat *) malloc(2*lni*sizeof(ftnfloat));
   lat =  (ftnfloat *) malloc(2*lni*sizeof(ftnfloat));
   lon =  (ftnfloat *) malloc(2*lni*sizeof(ftnfloat));
   gausslat = (ftnfloat *) malloc(2*lnj*sizeof(ftnfloat));
   
   f77name(ez_genpole)(vpolnor, vpolsud, fld, &lni, &lnj, &lvecteur, grtyp, &lhem,
                       x,y,z,lat,lon,gausslat, &groptions.degre_interp,1);

   free(x);
   free(y);
   free(z);
   free(lat);
   free(lon);
   free(gausslat);

   return 0;

}
