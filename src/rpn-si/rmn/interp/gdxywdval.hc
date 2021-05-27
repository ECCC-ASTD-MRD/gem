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
wordint f77name(gdxywdval)(wordint *gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint *n)
{
   wordint icode;
   
   icode = c_gdxywdval(*gdin, uuout, vvout, uuin, vvin, x, y, *n);
   return icode;

}
wordint c_gdxywdval(wordint gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint n)
{
  wordint ier,j, icode, yin_gdid, yan_gdid,lni,lnj;

  wordint gdrow_id, gdcol_id,yin_gdrow_id,yin_gdcol_id;
  ftnfloat *tmplat, *tmplon, *tmpy;
  ftnfloat *uuyin, *vvyin, *uuyan, *vvyan;
  ftnfloat *tmpuu, *tmpvv;
  
  tmplat = (ftnfloat *) malloc(n * sizeof(ftnfloat));
  tmplon = (ftnfloat *) malloc(n * sizeof(ftnfloat));
  tmpuu = (ftnfloat *) malloc(n * sizeof(ftnfloat));
  tmpvv = (ftnfloat *) malloc(n * sizeof(ftnfloat));
  
  c_gdkey2rowcol(gdin,  &gdrow_id,  &gdcol_id);
  if (Grille[gdrow_id][gdcol_id].nsubgrids > 0)
      {
      yin_gdid=Grille[gdrow_id][gdcol_id].subgrid[0];
      yan_gdid=Grille[gdrow_id][gdcol_id].subgrid[1];
      c_gdkey2rowcol(yin_gdid,  &yin_gdrow_id,  &yin_gdcol_id);
      lni = Grille[yin_gdrow_id][yin_gdcol_id].ni;
      lnj = Grille[yin_gdrow_id][yin_gdcol_id].nj;
      tmpy = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      uuyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      vvyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      uuyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      vvyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      for (j=0; j< n; j++)
        {
          if (y[j] > Grille[yin_gdrow_id][yin_gdcol_id].nj)
             {
             tmpy[j]=y[j]-Grille[yin_gdrow_id][yin_gdcol_id].nj;
             }
          else
             {
             tmpy[j]=y[j];
             }
        }
      icode = c_gdxyvval_orig(yin_gdid, tmpuu, tmpvv, uuin, vvin, x, tmpy, n);
      icode = c_gdllfxy_orig (yin_gdid, tmplat, tmplon, x, tmpy, n);
      icode = c_gdwdfuv_orig (yin_gdid, uuyin,vvyin,tmpuu,tmpvv,tmplat,tmplon,n);

      icode = c_gdxyvval_orig(yan_gdid, tmpuu, tmpvv, &uuin[(lni*lnj)], &vvin[(lni*lnj)], x, tmpy, n);
      icode = c_gdllfxy_orig (yan_gdid, tmplat, tmplon, x, tmpy, n);
      icode = c_gdwdfuv_orig (yan_gdid, uuyan,vvyan,tmpuu,tmpvv,tmplat,tmplon,n);
      for (j=0; j< n; j++)
        {
          if (y[j] > Grille[yin_gdrow_id][yin_gdcol_id].nj)
             {
             uuout[j]=uuyan[j];
             vvout[j]=vvyan[j];
             }
          else
             {
             uuout[j]=uuyin[j];
             vvout[j]=vvyin[j];
             }
        }
      free(uuyin); free(vvyin);
      free(uuyan); free(vvyan);
      free(tmpy);

      }
  else
      {
      ier = c_gdxyvval(gdin, tmpuu, tmpvv, uuin, vvin, x, y, n);
      ier = c_gdllfxy_orig(gdin, tmplat, tmplon, x, y, n);
      ier = c_gdwdfuv(gdin, uuout, vvout, tmpuu, tmpvv, tmplat, tmplon, n);
      }

  free(tmplat);
  free(tmplon);
  free(tmpuu);
  free(tmpvv);

  return 0;
}

