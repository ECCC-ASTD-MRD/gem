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
wordint ez_calcxy(wordint gdin, wordint gdout)
{
   wordint coordonnee, ni_in, nj_in, ni_out, nj_out, ninj_in, ninj_out;
   wordint i,j,ier;
   wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out, npts, cur_gdin, previous_val_polar_correction;
   int lcl_ngdin, idx_gdin;
   static wordint found = -1;
   static wordint ncalls = 0;

   _ygrid *ygrid;
   float *gdout_lat, *gdout_lon;

  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
  idx_gdin = c_find_gdin(gdin, gdout);
  /* Mettre du code au cas ou gdx_gdin == -1 */

  if (Grille[gdrow_out][gdcol_out].gset[idx_gdin].flags & XXX)
      {
      return 0;
      }

   /* Dans un premier temps on calcule la position x-y de tous les points sur la grille */

   ni_in =  Grille[gdrow_in][gdcol_in].ni;
   nj_in =  Grille[gdrow_in][gdcol_in].nj;
   ninj_in = ni_in * nj_in;

   ni_out = Grille[gdrow_out][gdcol_out].ni;
   nj_out = Grille[gdrow_out][gdcol_out].nj;
   ninj_out = ni_out * nj_out;

   Grille[gdrow_out][gdcol_out].gset[idx_gdin].x = (ftnfloat *) malloc(ninj_out*sizeof(ftnfloat));
   Grille[gdrow_out][gdcol_out].gset[idx_gdin].y = (ftnfloat *) malloc(ninj_out*sizeof(ftnfloat));

   switch(Grille[gdrow_in][gdcol_in].grtyp[0])
      {
      case 'A':
      case 'B':
      case 'E':
      case 'L':
      case 'N':
      case 'S':
      case 'T':
      case '!':
        f77name(ez_ll2rgd)(Grille[gdrow_out][gdcol_out].gset[idx_gdin].x,
                           Grille[gdrow_out][gdcol_out].gset[idx_gdin].y,
                           Grille[gdrow_out][gdcol_out].lat, Grille[gdrow_out][gdcol_out].lon, &ninj_out,
                           &ni_in, &nj_in, &Grille[gdrow_in][gdcol_in].grtyp,
                           &Grille[gdrow_in][gdcol_in].fst.ig[IG1], &Grille[gdrow_in][gdcol_in].fst.ig[IG2],
                           &Grille[gdrow_in][gdcol_in].fst.ig[IG3], &Grille[gdrow_in][gdcol_in].fst.ig[IG4],
                           &groptions.symmetrie, Grille[gdrow_in][gdcol_in].ay);
        break;


      case '#':
      case 'Z':
      case 'G':
         coordonnee = RELATIF;
         f77name(ez_ll2igd)(Grille[gdrow_out][gdcol_out].gset[idx_gdin].x,
                            Grille[gdrow_out][gdcol_out].gset[idx_gdin].y,
                            Grille[gdrow_out][gdcol_out].lat, Grille[gdrow_out][gdcol_out].lon, &ninj_out,
                            &ni_in,&nj_in,&Grille[gdrow_in][gdcol_in].grtyp, &Grille[gdrow_in][gdcol_in].grref,
                            &Grille[gdrow_in][gdcol_in].fst.igref[IG1], &Grille[gdrow_in][gdcol_in].fst.igref[IG2],
                            &Grille[gdrow_in][gdcol_in].fst.igref[IG3], &Grille[gdrow_in][gdcol_in].fst.igref[IG4],
                            Grille[gdrow_in][gdcol_in].ax, Grille[gdrow_in][gdcol_in].ay,
                            &coordonnee);
         if (Grille[gdrow_in][gdcol_in].grtyp[0] == 'G')
            {
            if (Grille[gdrow_in][gdcol_in].fst.ig[IG1] == NORD)
              {
              for (j=0; j < ni_out*nj_out; j++)
                {
                Grille[gdrow_out][gdcol_out].gset[idx_gdin].y[j] -= nj_in;
                }
               }
            }
         break;

      case 'Y':
         previous_val_polar_correction = groptions.polar_correction;
         groptions.polar_correction = NON;
         Grille[gdrow_out][gdcol_out].gset[idx_gdin].ygrid.n_wts = groptions.wgt_num;
/*         fprintf(stderr, "(ez_calcxy) %d\n", groptions.wgt_num);*/
         ygrid = &(Grille[gdrow_out][gdcol_out].gset[idx_gdin].ygrid);
         ygrid->lat =  (float *) malloc(ninj_in*sizeof(ftnfloat));
         ygrid->lon =  (float *) malloc(ninj_in*sizeof(ftnfloat));
         gdout_lat =  (float *) malloc(ninj_out*sizeof(ftnfloat));
         gdout_lon =  (float *) malloc(ninj_out*sizeof(ftnfloat));
         ygrid->wts =  (float *) malloc(ninj_out * groptions.wgt_num*sizeof(ftnfloat));
         ygrid->idx =  (int *) malloc(ninj_out * groptions.wgt_num*sizeof(wordint));
         ygrid->mask = (int *) malloc(ninj_out*sizeof(wordint));
         ier = c_gdll(gdin, ygrid->lat, ygrid->lon);
         ier = c_gdll(gdout, gdout_lat, gdout_lon);

         if (Grille[gdrow_in][gdcol_in].mask == NULL)
            {
            f77name(ez_calcxy_y)(ygrid->wts, ygrid->idx,
               Grille[gdrow_out][gdcol_out].gset[idx_gdin].x, Grille[gdrow_out][gdcol_out].gset[idx_gdin].y, gdout_lat, gdout_lon, ygrid->lat, ygrid->lon, ygrid->mask,
               &ni_in, &nj_in, &ni_out, &nj_out, &(groptions.wgt_num));
            }
         else
            {
            f77name(ez_calcxy_y_m)(ygrid->wts, ygrid->idx,
               Grille[gdrow_out][gdcol_out].gset[idx_gdin].x, Grille[gdrow_out][gdcol_out].gset[idx_gdin].y, gdout_lat, gdout_lon, ygrid->mask, ygrid->lat, ygrid->lon, Grille[gdrow_in][gdcol_in].mask,
               &ni_in, &nj_in, &ni_out, &nj_out, &(groptions.wgt_num));
            }

         groptions.polar_correction = previous_val_polar_correction;
         free(gdout_lat);
         free(gdout_lon);
      break;

      default:
        break;
      }

   Grille[gdrow_out][gdcol_out].gset[idx_gdin].flags |= XXX;

   return 0;
}

