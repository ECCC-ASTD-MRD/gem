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
wordint ez_corrval(ftnfloat *zout, ftnfloat *zin, wordint gdin, wordint gdout)
{
  wordint i,ierc;
  ftnfloat valmax, valmin,fudgeval;
  wordint fudgeval_set;
  wordint degIntCourant;
  wordint npts,nj;
  ftnfloat vpolnor, vpolsud;
  ftnfloat *temp;
  wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out, idx_gdin;
  _Grille *lgdin, *lgdout;

  extern ftnfloat f77name(amax)();
  extern ftnfloat f77name(amin)();

  gdin = iset_gdin;
  gdout= iset_gdout;
  fudgeval_set = 0;
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);

  lgdin = &(Grille[gdrow_in][gdcol_in]);
  lgdout = &(Grille[gdrow_out][gdcol_out]);
  
  idx_gdin = c_find_gdin(gdin, gdout);

  nj = lgdin->j2 - lgdin->j1 +1;
  ierc = 0; /* no extrapolation */
  
  if (lgdout->gset[idx_gdin].zones[DEHORS].npts > 0)
    {
    ierc = 2; /* code to indicate extrapolation */
    if (groptions.degre_extrap == ABORT)
      {
      fprintf(stderr, "<ez_corrval> There are points on the destination grid that lie outside the source grid\n");
      fprintf(stderr, "<ez_corrval> aborting at your request!\n\n\n");
      return -1;
      }

    f77name(ez_aminmax)(&valmin,&valmax,zin,&(lgdin->ni), &nj);
    if (groptions.degre_extrap >= MAXIMUM)
      {
      if (groptions.vecteur == VECTEUR)
	      {
	      fudgeval = 0.0;
         fudgeval_set = 1;
	      }
      else
         {
         switch (groptions.degre_extrap)
           {
           case MAXIMUM:
             fudgeval = valmax + 0.05 * (valmax - valmin);
                  fudgeval_set = 1;
             if (groptions.verbose > 0)
               {
               fprintf(stderr, "<ez_corrval>: maximum: %f \n", fudgeval);
               }
             break;

           case MINIMUM:
             fudgeval = valmin - 0.05 * (valmax - valmin);
                  fudgeval_set = 1;
             if (groptions.verbose > 0)
               {
               fprintf(stderr, "<ez_corrval>: minimum: %f \n", fudgeval);
               }
             break;

           case VALEUR:
             fudgeval = groptions.valeur_extrap;
                  fudgeval_set = 1;
             if (groptions.verbose > 0)
               {
               fprintf(stderr, "<ez_corrval>: valeur: %f \n", fudgeval);
               }
             break;
           }
         }

      if (fudgeval_set == 0) 
         {
        fprintf(stderr, "Error : ezcorrval : fudgeval not set \n");
         }
      for (i=0; i < lgdout->gset[idx_gdin].zones[DEHORS].npts; i++)
	      {
	      zout[lgdout->gset[idx_gdin].zones[DEHORS].idx[i]] = fudgeval;
	      }
      }
    else
      {
      degIntCourant = groptions.degre_interp;
      groptions.degre_interp = groptions.degre_extrap;
      temp = (ftnfloat *) malloc(lgdout->gset[idx_gdin].zones[DEHORS].npts*sizeof(ftnfloat));

      c_gdinterp(temp, zin, gdin, lgdout->gset[idx_gdin].zones[DEHORS].x,
		 lgdout->gset[idx_gdin].zones[DEHORS].y, lgdout->gset[idx_gdin].zones[DEHORS].npts);

     for (i=0; i < lgdout->gset[idx_gdin].zones[DEHORS].npts; i++)
         {
         zout[lgdout->gset[idx_gdin].zones[DEHORS].idx[i]] = temp[i];
         }
      free(temp);
      groptions.degre_interp = degIntCourant;
      }
   }

  if (groptions.vecteur == VECTEUR)
    {
    return ierc;
    }


  if (lgdout->gset[idx_gdin].zones[AU_NORD].npts > 0)
    {
    ez_corrval_aunord(zout, zin, gdin, gdout);
    }

  if (lgdout->gset[idx_gdin].zones[AU_SUD].npts > 0)
    {
    ez_corrval_ausud(zout, zin, gdin, gdout);
    }


  if (lgdout->gset[idx_gdin].zones[POLE_NORD].npts > 0 || lgdout->gset[idx_gdin].zones[POLE_SUD].npts > 0)
    {
    if (lgdin->grtyp[0] != 'w')
      {
      npts = lgdin->ni * lgdin->j2;
      f77name(ez_calcpoleval)(&vpolnor, &(zin[(nj-1)*lgdin->ni]), &(lgdin->ni),
lgdin->ax, lgdin->grtyp, lgdin->grref,1,1);
      for (i=0; i < lgdout->gset[idx_gdin].zones[POLE_NORD].npts; i++)
	      {
	      zout[lgdout->gset[idx_gdin].zones[POLE_NORD].idx[i]] = vpolnor;
	      }

      f77name(ez_calcpoleval)(&vpolsud, zin, &(lgdin->ni), lgdin->ax,
lgdin->grtyp, lgdin->grref,1,1);
      for (i=0; i < lgdout->gset[idx_gdin].zones[POLE_SUD].npts; i++)
	      {
	      zout[lgdout->gset[idx_gdin].zones[POLE_SUD].idx[i]] = vpolsud;
	      }
      }
    }
  if ((lgdin->grtyp[0] == 'Z' || lgdin->grtyp[0] == '#') && lgdin->grref[0] == 'E' && lgdout->grtyp[0] == 'B')
    {
    f77name(ez_corrbgd)(zout, &(lgdout->ni), &(lgdout->nj), &(lgdout->fst.ig[IG1]));
    }

  return ierc;
}
