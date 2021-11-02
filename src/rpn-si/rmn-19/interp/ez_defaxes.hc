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

void c_ezdefaxes(wordint gdid, ftnfloat *ax, ftnfloat *ay)
{
  wordint i,j;
  ftnfloat *temp, dlon;
  wordint zero, deuxnj;

  _Grille *gr;

  wordint gdrow_id, gdcol_id;

  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);

  gr = &Grille[gdrow_id][gdcol_id];
  switch (gr->grtyp[0])
    {
    case '#':
    case 'Z':
      f77name(cigaxg)(&gr->grref,&gr->fst.xgref[XLAT1], &gr->fst.xgref[XLON1], &gr->fst.xgref[XLAT2], &gr->fst.xgref[XLON2],
		      &gr->fst.igref[IG1], &gr->fst.igref[IG2], &gr->fst.igref[IG3], &gr->fst.igref[IG4],1);

      Grille[gdrow_id][gdcol_id].ax = (ftnfloat *) malloc(gr->ni*sizeof(ftnfloat));
      Grille[gdrow_id][gdcol_id].ay = (ftnfloat *) malloc(gr->nj*sizeof(ftnfloat));

      memcpy(Grille[gdrow_id][gdcol_id].ax,ax,gr->ni*sizeof(ftnfloat));
      memcpy(Grille[gdrow_id][gdcol_id].ay,ay,gr->nj*sizeof(ftnfloat));
      ez_calcxpncof(gdid);
      ez_calcntncof(gdid);
      break;

    case 'Y':
      Grille[gdrow_id][gdcol_id].ax = (ftnfloat *) malloc(gr->ni*gr->nj*sizeof(ftnfloat));
      Grille[gdrow_id][gdcol_id].ay = (ftnfloat *) malloc(gr->ni*gr->nj*sizeof(ftnfloat));
      memcpy(Grille[gdrow_id][gdcol_id].ax,ax,gr->ni*gr->nj*sizeof(ftnfloat));
      memcpy(Grille[gdrow_id][gdcol_id].ay,ay,gr->ni*gr->nj*sizeof(ftnfloat));

      ez_calcxpncof(gdid);
      break;

    case 'G':
      gr->grref[0] = 'L';
      gr->fst.xgref[SWLAT] = 0.0;
      gr->fst.xgref[SWLON] = 0.0;
      gr->fst.xgref[DLAT] = 1.0;
      gr->fst.xgref[DLON] = 1.0;
      f77name(cxgaig)(&gr->grref,&gr->fst.igref[IG1], &gr->fst.igref[IG2], &gr->fst.igref[IG3], &gr->fst.igref[IG4],
		      &gr->fst.xgref[SWLAT], &gr->fst.xgref[SWLON], &gr->fst.xgref[DLAT], &gr->fst.xgref[DLON],1);

      Grille[gdrow_id][gdcol_id].ax = (ftnfloat *) malloc(gr->ni*sizeof(ftnfloat));
      dlon = 360. / (ftnfloat) gr->ni;
      for (i=0; i < gr->ni; i++)
	      {
	      Grille[gdrow_id][gdcol_id].ax[i] = (ftnfloat)i * dlon;
	      }

      zero = 0;
      ez_calcxpncof(gdid);

      switch (Grille[gdrow_id][gdcol_id].fst.ig[IG1])
	      {
	      case GLOBAL:
	        Grille[gdrow_id][gdcol_id].ay = (ftnfloat *) malloc(gr->nj*sizeof(ftnfloat));
	        temp    = (ftnfloat *) malloc(gr->nj*sizeof(ftnfloat));
	        f77name(ez_glat)(Grille[gdrow_id][gdcol_id].ay,temp,&gr->nj,&zero);
	        free(temp);
	        break;

	      case NORD:
	      case SUD:
	        deuxnj = 2 * gr->nj;
	        Grille[gdrow_id][gdcol_id].ay = (ftnfloat *) malloc(deuxnj*sizeof(ftnfloat));
	        temp    = (ftnfloat *) malloc(deuxnj*sizeof(ftnfloat));
	        f77name(ez_glat)(Grille[gdrow_id][gdcol_id].ay,temp,&deuxnj,&zero);
	        free(temp);
	        break;
	      }


      ez_calcntncof(gdid);
      Grille[gdrow_id][gdcol_id].flags |= AX;
      break;

    default:
      ez_calcxpncof(gdid);
      break;
    }


  Grille[gdrow_id][gdcol_id].flags |= AX;

}
