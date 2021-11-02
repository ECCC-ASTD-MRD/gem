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
wordint f77name(gdxyvval)(wordint *gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint *n)
{
   wordint icode;

   icode = c_gdxyvval(*gdin, uuout, vvout, uuin, vvin, x, y, *n);
   return icode;
}

wordint c_gdxyvval(wordint gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint n)
{
  wordint j, icode, yin_gdid, yan_gdid, ni, nj;
  ftnfloat *uuyin, *vvyin, *uuyan, *vvyan;
  ftnfloat *tmpy;

  wordint gdrow_id, gdcol_id,yin_gdrow_id,yin_gdcol_id;

  c_gdkey2rowcol(gdin,  &gdrow_id,  &gdcol_id);
  if (Grille[gdrow_id][gdcol_id].nsubgrids > 0)
      {
      yin_gdid=Grille[gdrow_id][gdcol_id].subgrid[0];
      yan_gdid=Grille[gdrow_id][gdcol_id].subgrid[1];
      c_gdkey2rowcol(yin_gdid,  &yin_gdrow_id,  &yin_gdcol_id);
      ni = Grille[yin_gdrow_id][yin_gdcol_id].ni;
      nj = Grille[yin_gdrow_id][yin_gdcol_id].nj;
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
      icode = c_gdxyvval_orig(yin_gdid,uuyin,vvyin,uuin,vvin,x,tmpy,n);
      icode = c_gdxyvval_orig(yan_gdid,uuyan,vvyan,&uuin[ni*nj],&vvin[ni*nj],x,tmpy,n);
 
      for (j=0; j < n; j++)
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
      free(tmpy);
      free(uuyan); free(vvyan); 
      free(uuyin); free(vvyin);
      return icode;
    }
  else
   {
   icode = c_gdxyvval_orig(gdin,uuout,vvout,uuin,vvin,x,y,n);
/*  printf("gdxyvval reg val %f,%f uuout=%f vvout=%f\n",x[j],y[j],uuout[j],vvout[j]); */
   return icode;
   }
}

wordint c_gdxyvval_orig(wordint gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint n)
{
  groptions.vecteur = VECTEUR;


  groptions.symmetrie = SYM;
  c_gdxysint(uuout,uuin, gdin, x, y, n);
  groptions.symmetrie = ANTISYM;
  c_gdxysint(vvout,vvin, gdin, x, y, n);
  groptions.symmetrie = SYM;

 groptions.vecteur = SCALAIRE;

  return 0;

}
