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
wordint ez_defzone_dehors(wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zone)
{
  ftnfloat *tmpx, *tmpy;
  wordint *tmpidx;
  wordint nhits;
  wordint i;
  
  wordint offsetleft, offsetright, ix, iy;
  
   wordint gdrow_in, gdcol_in;
   int lcl_ngdin;
   
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  
  tmpx =   (ftnfloat *) malloc(npts*sizeof(ftnfloat));
  tmpy =   (ftnfloat *) malloc(npts*sizeof(ftnfloat));
  tmpidx = (wordint   *) malloc(npts*sizeof(wordint));
  
/*
  if (groptions.degre_interp == CUBIQUE)
    {
    offsetright = 2;
    offsetleft = 1;
    }
  else
    {
    offsetright = 0;
    offsetleft = 0;
    }
*/
  
  offsetright = 0;
  offsetleft = 0;
  if (groptions.verbose > 0)
    {
    fprintf(stderr, "degre_extrap: %d offset left: %d offset right: %d\n", groptions.degre_extrap, offsetleft, offsetright);
    }
  nhits = 0;
  for (i=0; i < npts; i++)
    {
    ix = (wordint)(x[i]+0.5);
    iy = (wordint)(y[i]+0.5);
    if (ix < (1+offsetleft) || iy < (1+offsetleft) || ix > (Grille[gdrow_in][gdcol_in].ni-offsetright) || iy > (Grille[gdrow_in][gdcol_in].nj-offsetright))
      {
      tmpx[nhits]  = x[i];
      tmpy[nhits]  = y[i];
      tmpidx[nhits]=i;
      nhits++;
      }
    }

  if (nhits > 0)
    {
    zone->npts = nhits;
    zone->x =   (ftnfloat *) malloc(zone->npts*sizeof(ftnfloat));
    zone->y =   (ftnfloat *) malloc(zone->npts*sizeof(ftnfloat));
    zone->idx = (wordint *) malloc(zone->npts*sizeof(wordint));
    if (groptions.verbose > 0)
      {
      fprintf(stderr, "Nombre de points dehors: %d\n", zone->npts); 
      }
    for (i=0; i < zone->npts; i++)
      {
      zone->x[i] = tmpx[i];      
      zone->y[i] = tmpy[i];     
      zone->idx[i] = tmpidx[i];
      }
    }
  
  free(tmpx);
  free(tmpy);
  free(tmpidx);
  
  return 0;
}
