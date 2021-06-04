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

wordint ez_defzone_polesud(wordint gdin, ftnfloat *x, ftnfloat *y, wordint npts, _zone *zone)
{

  ftnfloat *tmpx, *tmpy;
  ftnfloat latpolesud, lonpolesud, xpolesud, ypolesud;
  wordint nhits, i;
  wordint *tmpidx;


  wordint gdrow_in, gdcol_in;
    
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  
  tmpx =   (ftnfloat *) malloc(npts*sizeof(ftnfloat));
  tmpy =   (ftnfloat *) malloc(npts*sizeof(ftnfloat));
  tmpidx = (wordint   *) malloc(npts*sizeof(wordint));
  
  nhits = 0;
  
  if (Grille[gdrow_in][gdcol_in].grtyp[0] == 'Z' && Grille[gdrow_in][gdcol_in].grref[0] == 'E')
    {
    xpolesud = 0.5 * Grille[gdrow_in][gdcol_in].ni;
    ypolesud = 0.5;
    }
  else
    {
    latpolesud = -90.0;
    lonpolesud = 0.0;
    c_gdxyfll_orig(gdin, &xpolesud, &ypolesud,  &latpolesud, &lonpolesud, 1);
    }


  for (i=0; i < npts; i++)
    {
    if ((fabs(y[i]-ypolesud) < 1.0e-3))
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
      fprintf(stderr, "Nombre de points au pole sud: %d\n", nhits); 
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
