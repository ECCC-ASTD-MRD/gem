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
wordint ez_corrval_aunord(ftnfloat *zout, ftnfloat *zin, wordint gdin, wordint gdout)
{
  wordint i;
  wordint npts;
  ftnfloat *temp, *vals, *temp_y;
  ftnfloat ay[4], poleval;
  wordint ni, nj, i1, i2, j1, j2;
  wordint un = 1;
  wordint quatre = 4;
  wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out, idx_gdin;
  _gridset *gset;
  _Grille *lgd;

  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
  idx_gdin = c_find_gdin(gdin, gdout);

  gset = &(Grille[gdrow_out][gdcol_out].gset[idx_gdin]);
  npts = gset->zones[AU_NORD].npts;
  lgd = (_Grille *) &(Grille[gdrow_in][gdcol_in]);
  if (npts > 0)
    {
    ni = lgd->ni;
    nj = lgd->j2 - lgd->j1 + 1;

    i1 = 1;
    i2 = ni;

    j1 = lgd->j2-2;
    j2 = j1 + 3;

    temp = (ftnfloat *) malloc(4 * ni * sizeof(ftnfloat));
    vals = (ftnfloat *) malloc(npts * sizeof(ftnfloat));
    f77name(ez_calcpoleval)(&poleval, &zin[(nj-1)*ni], &ni, lgd->ax, &(lgd->grtyp),
&(lgd->grref),1,1);
    f77name(ez_fillnpole)(temp, zin, &ni, &(lgd->j1), &(lgd->j2), &poleval);

    switch (groptions.degre_interp)
      {
      case CUBIQUE:
	   switch (lgd->grtyp[0])
	     {
	     case 'Z':
	     case 'E':
	     case 'G':
    if  (lgd->ay[lgd->j2-1] == 90.0)
       {
          ay[0] = lgd->ay[lgd->j2-4];
          ay[1] = lgd->ay[lgd->j2-3];
          ay[2] = lgd->ay[lgd->j2-2];
          ay[3] = lgd->ay[lgd->j2-1];
          f77name(ez_irgdint_3_wnnc)(vals,gset->zones[AU_NORD].x, gset->zones[AU_NORD].y,&npts,
                      lgd->ax, ay, temp,&ni, &j1, &j2, &(lgd->extension));
       }
    else
       {
          ay[0] = lgd->ay[lgd->j2-3];
          ay[1] = lgd->ay[lgd->j2-2];
          ay[2] = lgd->ay[lgd->j2-1];
          ay[3] = 90.0;
          f77name(ez_irgdint_3_wnnc)(vals,gset->zones[AU_NORD].x,
                      gset->zones[AU_NORD].y,&npts,
                      lgd->ax, ay, temp,
                      &ni, &j1, &j2, &(lgd->extension));
       }
	       break;

	     default:
	       f77name(ez_rgdint_3_wnnc)(vals,gset->zones[AU_NORD].x,
				         gset->zones[AU_NORD].y,&npts,
				         temp,&ni, &j1, &j2, &(lgd->extension));
	       break;
	     }
   	break;

      case LINEAIRE:
	   temp_y = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
	   for (i=0; i < npts; i++)
	     {
        temp_y[i] = gset->zones[AU_NORD].y[i] - (1.0 * (lgd->j2-3));
	     }
	   f77name(ez_rgdint_1_w)(vals,gset->zones[AU_NORD].x,temp_y,&npts,temp,&ni, &un, &quatre, &(lgd->extension));
	   free(temp_y);
	   break;

      case VOISIN:
	   temp_y = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
	   for (i=0; i < npts; i++)
	     {
	     temp_y[i] = gset->zones[AU_NORD].y[i] - (1.0 * (lgd->j2-3));
	     }
	   f77name(ez_rgdint_0)(vals,gset->zones[AU_NORD].x,temp_y,&npts,temp,&ni, &un, &quatre);
	   free(temp_y);
	   break;
      }


    for (i=0; i < gset->zones[AU_NORD].npts; i++)
      {
      zout[gset->zones[AU_NORD].idx[i]] = vals[i];
      }

    free(vals);
    free(temp);
    }
  return 0;
}

