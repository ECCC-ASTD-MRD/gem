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
wordint ez_corrvec_ausud(ftnfloat *uuout, ftnfloat *vvout,
         ftnfloat *uuin, ftnfloat *vvin,
         wordint gdin, wordint gdout)
{
ftnfloat *polar_uu_in, *polar_vv_in, *corr_uus, *corr_vvs, *temp_y, ay[4];
wordint ni, nj, i1, i2, j1, j2, degree,npts,i;
wordint quatre = 4;
wordint un = 1;

wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out, idx_gdin;
_gridset *gset;
_Grille laGrille;


c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
idx_gdin = c_find_gdin(gdin, gdout);

gset = &(Grille[gdrow_out][gdcol_out].gset[idx_gdin]);
laGrille = Grille[gdrow_in][gdcol_in];

npts = gset->zones[AU_SUD].npts;
ni = laGrille.ni;
nj = laGrille.j2 - laGrille.j1 + 1;

i1 = 1;
i2 = ni;
j1 = laGrille.j1 - 1;
j2 = j1 + 3;
degree = 3;

polar_uu_in = (ftnfloat *) malloc(4 * ni * sizeof(ftnfloat));
polar_vv_in = (ftnfloat *) malloc(4 * ni * sizeof(ftnfloat));
corr_uus = (ftnfloat *) malloc(npts * sizeof(ftnfloat));
corr_vvs = (ftnfloat *) malloc(npts * sizeof(ftnfloat));

ez_calcspolarwind(polar_uu_in, polar_vv_in, uuin, vvin, ni, nj, gdin);

switch (groptions.degre_interp)
   {
   case CUBIQUE:
      switch (laGrille.grtyp[0])
         {
         case 'Z':
         case 'E':
         case 'G':
         if (laGrille.ay[0] == -90.0)
            {
            ay[0] = laGrille.ay[0];
            ay[1] = laGrille.ay[1];
            ay[2] = laGrille.ay[2];
            ay[3] = laGrille.ay[3];
            }
         else
            {
            ay[0] = -90.;
            ay[1] = laGrille.ay[0];
            ay[2] = laGrille.ay[1];
            ay[3] = laGrille.ay[2];
               
            }
         f77name(ez_irgdint_3_wnnc)(corr_uus,gset->zones[AU_SUD].x,
                     gset->zones[AU_SUD].y,&npts,
                     laGrille.ax, ay, polar_uu_in,
                     &ni, &j1, &j2, &laGrille.extension);
         f77name(ez_irgdint_3_wnnc)(corr_vvs,gset->zones[AU_SUD].x,
                     gset->zones[AU_SUD].y,&npts,
                     laGrille.ax, ay, polar_vv_in,
                     &ni, &j1, &j2, &laGrille.extension);
         break;

         default:
         f77name(ez_rgdint_3_wnnc)(corr_uus,gset->zones[AU_SUD].x,
                     gset->zones[AU_SUD].y,&npts,
                     polar_uu_in,&ni, &j1, &j2, &laGrille.extension);
         f77name(ez_rgdint_3_wnnc)(corr_vvs,gset->zones[AU_SUD].x,
                     gset->zones[AU_SUD].y,&npts,
                     polar_vv_in,&ni, &j1, &j2, &laGrille.extension);
         break;
         }
      break;

   case LINEAIRE:
      f77name(ez_rgdint_1_w)(corr_uus,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,polar_uu_in,&ni, &j1, &j2, &laGrille.extension);
      f77name(ez_rgdint_1_w)(corr_vvs,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,polar_vv_in,&ni, &j1, &j2, &laGrille.extension);
      break;

   case VOISIN:
      f77name(ez_rgdint_0)(corr_uus,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,polar_uu_in,&ni, &j1, &j2);
      f77name(ez_rgdint_0)(corr_vvs,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,polar_vv_in,&ni, &j1, &j2);
      break;
   }


for (i=0; i < gset->zones[AU_SUD].npts; i++)
   {
   uuout[gset->zones[AU_SUD].idx[i]] = corr_uus[i];
   vvout[gset->zones[AU_SUD].idx[i]] = corr_vvs[i];
   }

free(polar_uu_in);
free(polar_vv_in);
free(corr_uus);
free(corr_vvs);

return 0;
}
