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

wordint ez_defzone_sud(wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zone)
{
  ftnfloat *tmpx, *tmpy;
  wordint *tmpidx;
  wordint nhits, i;
  wordint jmin;


  wordint gdrow_in, gdcol_in;
    
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  
  tmpx =   (ftnfloat *) malloc(npts*sizeof(ftnfloat));
  tmpy =   (ftnfloat *) malloc(npts*sizeof(ftnfloat));
  tmpidx = (wordint  *) malloc(npts*sizeof(wordint));
  
  nhits = 0;
  jmin = Grille[gdrow_in][gdcol_in].j1+1;
  for (i=0; i < npts; i++)
    {
    if ((int)y[i] < jmin)
      {
      tmpx[nhits] = x[i];
      tmpy[nhits] = y[i];
      tmpidx[nhits]=i;
      nhits++;
      }
    }
  
  zone->npts = nhits;
  if (nhits > 0)
    {
    zone->x = (ftnfloat *) malloc(nhits*sizeof(ftnfloat));
    zone->y = (ftnfloat *) malloc(nhits*sizeof(ftnfloat));
    zone->idx = (wordint *) malloc(nhits*sizeof(wordint));
    if (groptions.verbose > 0)
      {
      fprintf(stderr, "Nombre de points entre le pole et nj=2 : %d\n", nhits); 
      }
    
    for (i=0; i < zone->npts; i++)
      {
      zone->x[i]   = tmpx[i];      
      zone->y[i]   = tmpy[i];     
      zone->idx[i] = tmpidx[i];
      }
    }
  
  free(tmpx);
  free(tmpy);
  free(tmpidx);
  return 0;
}
