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

wordint c_ezgdef_yymask(_Grille *subgd)
{
  wordint ni,nj,yni,ynj,i,j,k,ii;
  wordint i0,i1,j0,j1;
  wordint ig1ref,ig2ref,ig3ref,ig4ref;
  wordint mask_gdrow_id, mask_gdcol_id, mask_gdid;
  
  ftnfloat *ax,*ay;

  k=0;
  for (i=0; i < subgd->ni; i++)
       {
       if (subgd->ax[i] >= 45.0 && subgd->ax[i] <= 315.0)
           {
            k++;
            if (k == 1) i0=i;
            i1=i;
           }
       }
     yni=k;
     k=0;
     for (i=0; i < subgd->nj; i++)
       {
       if (subgd->ay[i] >= -45.0 && subgd->ay[i] <= 45.0)
           {
            k++;
            if (k == 1) j0=i;
            j1=i;
           }
       }
     ynj=k;
     mask_gdid=c_ezgdef_fmem(yni,ynj,subgd->grtyp,subgd->grref,subgd->fst.igref[IG1],subgd->fst.igref[IG2],subgd->fst.igref[IG3],subgd->fst.igref[IG4],&subgd->ax[i0],&subgd->ay[j0]);
     subgd->mymaskgrid = mask_gdid;
     subgd->mymaskgridi0=i0;
     subgd->mymaskgridi1=i1;
     subgd->mymaskgridj0=j0;
     subgd->mymaskgridj1=j1;

    if (groptions.verbose > 0)
      {
       c_gdkey2rowcol(mask_gdid,  &mask_gdrow_id,  &mask_gdcol_id);
       printf("Subgd.mymaskgrid   = %d\n", subgd->mymaskgrid);
       printf("Subgd.mymaskgridi0 = %d pt=%f\n", subgd->mymaskgridi0, subgd->ax[i0]);
       printf("Subgd.mymaskgridi1 = %d pt=%f\n", subgd->mymaskgridi1, subgd->ax[i1]);
       printf("Subgd.mymaskgridj0 = %d pt=%f\n", subgd->mymaskgridj0, subgd->ay[j0]);
       printf("Subgd.mymaskgridj1 = %d pt=%f\n", subgd->mymaskgridj1, subgd->ay[j1]);
      }
    return 0;

}
 
