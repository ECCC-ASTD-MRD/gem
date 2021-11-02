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

#include <stdio.h>
#include <rpnmacros.h>

static ftnfloat *eziglat = NULL;
static ftnfloat *eziglon = NULL;

f77name(ez_igscint)(ftnfloat *zo, wordint *li, wordint *lj, ftnfloat *xlat, ftnfloat *xlon, 
                    ftnfloat *zi, wordint *ni, wordint *nj, 
                    char* grtyp, char *grref, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, 
                    ftnfloat *sym, ftnfloat *ax, ftnfloat *ay, F2Cl lengrtyp, F2Cl lengrref)
{
  wordint gdin, gdout, npts, i, ier;
  char ogrtyp[2], ogrref[2];
  wordint oig1, oig2, oig3, oig4;
  ftnfloat xg1, xg2, xg3, xg4;
  ftnfloat *tmplon;
  
  ftnstrclean(grtyp,lengrtyp);
  ftnstrclean(grref,lengrref);
  
  npts = *li * *lj;
  tmplon = (ftnfloat *) malloc(npts * sizeof(ftnfloat));
  for (i=0; i < npts; i++)
    {
    tmplon[i] = xlon[i] < 0.0 ? xlon[i] + 360.0 : xlon[i];
    }
  
  gdin  = c_ezgdef_fmem(*ni, *nj, grtyp, grref,  *ig1,  *ig2,  *ig3,  *ig4, ax,   ay);
  ier = c_gdllsval(gdin, zo, zi, xlat, tmplon, npts);
  free(tmplon);
  return 0;
}

f77name(ez_rgscint)(ftnfloat *zo, wordint *li, wordint *lj, ftnfloat *xlat, ftnfloat *xlon, 
                    ftnfloat *zi, wordint *ni, wordint *nj, 
                    char* grtyp, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, ftnfloat *sym,
                    F2Cl lengrtyp)
{
  wordint gdin, gdout;
  char ogrtyp[2], ogrref[2], igrtyp[2],igrref[2];
  wordint oig1, oig2, oig3, oig4, npts, i, ier;
  ftnfloat xg1, xg2, xg3, xg4;
  ftnfloat *tmplon;
  
  npts = *li * *lj;
  tmplon = (ftnfloat *) malloc(npts * sizeof(ftnfloat));
  for (i=0; i < npts; i++)
    {
    tmplon[i] = xlon[i] < 0.0 ? xlon[i] + 360.0 : xlon[i];
    }
  
  gdin  = c_ezgdef_fmem(*ni, *nj, grtyp, NULL,  *ig1,  *ig2,  *ig3,  *ig4, NULL, NULL);
  ier = c_gdllsval(gdin, zo, zi, xlat, tmplon, npts);
  free(tmplon);
  return 0;
}

f77name(ez_iguvint)(ftnfloat *spdo, ftnfloat *psio, wordint *li, wordint *lj, ftnfloat *xlat, ftnfloat *xlon, 
                    ftnfloat *ui, ftnfloat *vi, wordint *ni, wordint *nj, 
                    char* grtyp, char *grref, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, 
                    ftnfloat *sws, ftnfloat *ax, ftnfloat *ay, F2Cl lengrtyp, F2Cl lengrref)
{
  ftnfloat *tmplon;
  wordint gdin, gdout, i, ier;
  char ogrtyp[2], ogrref[2];
  ftnfloat xg1, xg2, xg3, xg4;
  wordint npts, oig1, oig2, oig3, oig4;
  
  npts = *li * *lj;
  tmplon = (ftnfloat *) malloc(npts * sizeof(ftnfloat));
  for (i=0; i < npts; i++)
    {
    tmplon[i] = xlon[i] < 0.0 ? xlon[i] + 360.0 : xlon[i];
    }
  
  ftnstrclean(grtyp,lengrtyp);
  ftnstrclean(grref,lengrref);
  gdin  = c_ezgdef_fmem(*ni,*nj, grtyp, grref,  *ig1,  *ig2,  *ig3,  *ig4, ax, ay);
  ier = c_gdllwdval(gdin, spdo, psio, ui, vi, xlat, tmplon, npts);
  free(tmplon);
}

f77name(ez_rguvint)(ftnfloat *spdo, ftnfloat *psio, wordint *li, wordint *lj, ftnfloat *xlat, ftnfloat *xlon, 
                    ftnfloat *ui, ftnfloat *vi, wordint *ni, wordint *nj, 
                    char* grtyp, wordint *ig1, wordint *ig2, wordint *ig3, wordint *ig4, ftnfloat *sws, F2Cl lengrtyp)
{
  ftnfloat *tmplon;
  wordint gdin, gdout, i, ier;
  char ogrtyp[2], ogrref[2];
  ftnfloat xg1, xg2, xg3, xg4;
  wordint npts, oig1, oig2, oig3, oig4;
  
  npts = *li * *lj;
  tmplon = (ftnfloat *) malloc(npts * sizeof(ftnfloat));
  for (i=0; i < npts; i++)
    {
    tmplon[i] = xlon[i] < 0.0 ? xlon[i] + 360.0 : xlon[i];
    }
  
  ftnstrclean(grtyp,lengrtyp);
  gdin  = c_ezgdef(*ni,*nj, grtyp, NULL,  *ig1,  *ig2,  *ig3,  *ig4, NULL, NULL);
  ier = c_gdllwdval(gdin, spdo, psio, ui, vi, xlat, tmplon, npts);
  free(tmplon);
}

