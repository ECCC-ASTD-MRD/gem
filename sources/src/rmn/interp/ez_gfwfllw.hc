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

void f77name(ez_gfwfllw)(ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *latin, ftnfloat *lonin,
			 ftnfloat *xlatingf, ftnfloat *xloningf, 
			 wordint *ni, wordint *nj,
			 char *grtyp, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4)
{
  c_ezgfwfllw(uullout, vvllout, latin, lonin, xlatingf, xloningf, 
	      ni, nj, grtyp, ig1, ig2, ig3, ig4);
    }

void c_ezgfwfllw(ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *latin, ftnfloat *lonin,
                  ftnfloat *xlatingf, ftnfloat *xloningf, 
                  wordint *ni, wordint *nj,
                  char *grtyp, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4)
{

  /*
    gfwfllw -> GeF Winds From LatLon Winds
    Ce sous-programme effectue la rotation des vents d'un systeme de coordonne
    non tourne a un systeme de coordonnee tourne.
    latin, lonin sont les latlons vraies
    xlatingf, xloningf sont les latlons sur la grille tournee
  */

  wordint zero = 0;
  wordint npts = *ni * *nj;
  wordint trois = 3;
  ftnfloat r[9], ri[9], xlon1, xlat1, xlon2, xlat2;
  ftnfloat *uvcart, *xyz;
  char grtypl[2];

  uvcart = (ftnfloat *) malloc(3*npts*sizeof(ftnfloat));
  xyz    = (ftnfloat *) malloc(3*npts*sizeof(ftnfloat));
  
  
  f77name(cigaxg)(grtyp, &xlat1, &xlon1, &xlat2, &xlon2, ig1, ig2, ig3, ig4, 1);
  f77name(ez_crot)(r, ri, &xlon1, &xlat1, &xlon2, &xlat2);
  
  grtypl[0] = 'L';
  f77name(ez_gdwfllw)(uullout,vvllout,lonin,ni,nj,grtypl, &zero, &zero, &zero, &zero, 1);
  
  f77name(ez_uvacart)(xyz, uullout, vvllout, lonin, latin, ni, nj);
  f77name(mxm)(r, &trois, xyz, &trois, uvcart, &npts);
  f77name(ez_cartauv)(uullout, vvllout, uvcart, xloningf, xlatingf, ni, nj);

  free(uvcart);
  free(xyz);
}


