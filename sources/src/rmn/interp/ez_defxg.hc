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

void c_ezdefxg(wordint gdid)
{

  _Grille *gr;

  wordint gdrow_id, gdcol_id;

  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);

  gr = &Grille[gdrow_id][gdcol_id];

  switch (gr->grtyp[0])
    {
    case 'A':
    case 'G':
      gr->fst.xg[DLON]  = 360. /gr->ni;
      gr->fst.xg[SWLON] = 0.0;
      switch (gr->fst.ig[IG1])
	{
	case 0:
	  gr->fst.xg[DLAT] = 180./gr->nj;
	  gr->fst.xg[SWLAT] = -90. + 0.5*gr->fst.xg[DLAT];
	  break;

	case 1:
	  gr->fst.xg[DLAT] = 90./gr->nj;
	  gr->fst.xg[SWLAT] = 0.5*gr->fst.xg[DLAT];
	  gr->needs_expansion = OUI;
	  break;

	case 2:
	  gr->fst.xg[DLAT] = 90./gr->nj;
	  gr->fst.xg[SWLAT] = -90. + 0.5*gr->fst.xg[DLAT];
	  gr->needs_expansion = OUI;
	  break;

	default:
	  fprintf(stderr, "<ez_gdef_fmem> 'A' grid has to be Global/North/South\n");
	  break;
	}

   switch(gr->fst.ig[IG2])
	   {
	   case 1:
	     gr->fst.axe_y_inverse = OUI;
	     break;

	   default:
	     break;
	   }

      break;

    case 'B':
      gr->fst.xg[DLON] = 360. /(gr->ni-1);
      gr->fst.xg[SWLON] = 0.0;
      switch (gr->fst.ig[IG1])
	      {
	      case 0:
	        gr->fst.xg[DLAT] = 180./(gr->nj-1);
	        gr->fst.xg[SWLAT] = -90.;
	        break;

	      case 1:
	        gr->fst.xg[DLAT] = 90./(gr->nj-1);
	        gr->fst.xg[SWLAT] = 0.;
	        gr->needs_expansion = OUI;
	        break;

	      case 2:
	        gr->fst.xg[DLAT] = 90./(gr->nj-1);
	        gr->fst.xg[SWLAT] = -90.;
	        gr->needs_expansion = OUI;
	        break;

	      default:
	        fprintf(stderr, "<ezgdef_fmem> 'B' grid has to be Global/North/South\n");
	        break;
	      }

      switch(gr->fst.ig[IG2])
	      {
	      case 1:
	        gr->fst.axe_y_inverse = OUI;
	        break;

	      default:
	        break;
	      }
      break;

    case 'E':
      f77name(cigaxg)(&gr->grtyp,&gr->fst.xg[XLAT1],&gr->fst.xg[XLON1],&gr->fst.xg[XLAT2],&gr->fst.xg[XLON2],
		      &gr->fst.ig[IG1],&gr->fst.ig[IG2],&gr->fst.ig[IG3],&gr->fst.ig[IG4],1);
      /*      gr->fst.xg[DLAT] = 180./gr->nj;
	      gr->fst.xg[DLON] = 360./(gr->ni-1);
	      gr->fst.xg[SWLON] = 0.0;
	      gr->fst.xg[SWLAT] = -90. + 0.5*gr->fst.xg[DLAT];
      */
      break;

    case 'H':
    case 'Y':
    case '!':
      break;

    case '#':
    case 'Z':
      if (gr->grref[0] == 'N') gr->fst.hemisphere = 1;
      if (gr->grref[0] == 'S') gr->fst.hemisphere = 2;
      if (gr->grref[0] == 'E')
         {
         f77name(cigaxg)(&gr->grref,&gr->fst.xgref[XLAT1], &gr->fst.xgref[XLON1], &gr->fst.xgref[XLAT2], &gr->fst.xgref[XLON2],
            &gr->fst.igref[IG1], &gr->fst.igref[IG2], &gr->fst.igref[IG3], &gr->fst.igref[IG4],1);
         }

    break;

    case 'L':
      f77name(cigaxg)(&gr->grtyp,&gr->fst.xg[SWLAT], &gr->fst.xg[SWLON], &gr->fst.xg[DLAT], &gr->fst.xg[DLON],
		      &gr->fst.ig[IG1], &gr->fst.ig[IG2], &gr->fst.ig[IG3], &gr->fst.ig[IG4],1);
      break;

    case 'N':
      f77name(cigaxg)(&gr->grtyp,&gr->fst.xg[PI], &gr->fst.xg[PJ], &gr->fst.xg[D60], &gr->fst.xg[DGRW],
		      &gr->fst.ig[IG1], &gr->fst.ig[IG2], &gr->fst.ig[IG3], &gr->fst.ig[IG4],1);
      gr->fst.hemisphere = 1;
      break;

    case 'S':
      f77name(cigaxg)(&gr->grtyp,&gr->fst.xg[PI], &gr->fst.xg[PJ], &gr->fst.xg[D60], &gr->fst.xg[DGRW],
		      &gr->fst.ig[IG1], &gr->fst.ig[IG2], &gr->fst.ig[IG3], &gr->fst.ig[IG4],1);
      gr->fst.hemisphere = 2;
      break;

    case 'T':
      f77name(cigaxg)(&gr->grtyp,&gr->fst.xg[TD60], &gr->fst.xg[TDGRW], &gr->fst.xg[CLAT], &gr->fst.xg[CLON],
		      &gr->fst.ig[IG1], &gr->fst.ig[IG2], &gr->fst.ig[IG3], &gr->fst.ig[IG4],1);
      break;

    case 'X':
      fprintf(stderr,"<c_ezgdef> There is no support for grid type 'X'\n");
      return;

    default:
      fprintf(stderr,"<c_ezgdef> Grid type not supported\n");
      return;
    }

  return;

}

