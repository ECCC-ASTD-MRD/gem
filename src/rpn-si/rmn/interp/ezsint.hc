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
wordint f77name(ezsint)(ftnfloat *zout, ftnfloat *zin)
{
   wordint icode;
   
   icode = c_ezsint(zout, zin);
   return icode;
}
wordint c_ezsint(ftnfloat *zout, ftnfloat *zin)
{
  wordint icode,gdin,gdout;
  wordint gdrow_in,gdcol_in, gdrow_out,gdcol_out;
   
  if (iset_gdin == UNDEFINED || iset_gdout == UNDEFINED)
    {
    fprintf(stderr,"<c_ezsint> Source or target grid undefined! Aborting...\n");
    return -1;
    }
  
  
  gdin = iset_gdin;
  gdout= iset_gdout;
  
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   
  if (iset_gdin == iset_gdout)
    {
    memcpy(zout, zin, Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    return 1;
    }


  if (Grille[gdrow_in][gdcol_in].nsubgrids > 0 || Grille[gdrow_out][gdcol_out].nsubgrids > 0)
      {
/* get the subgrids and interpolate accordingly */
      icode = c_ezyysint(zout,zin,gdout,gdin);
      iset_gdin=gdin;
      iset_gdout=gdout;
      return icode;
      }
  icode = c_ezsint_orig(zout, zin);
  return icode;
}

wordint c_ezsint_orig(ftnfloat *zout, ftnfloat *zin)
{
  wordint gdin, gdout;
  wordint ier,ierc;
  wordint npts;
  wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out;
  ftnfloat *lzin, *lxzin;
  
  lzin  = NULL;
  lxzin = NULL;
  ierc  = 0;
  
  if (iset_gdin == UNDEFINED || iset_gdout == UNDEFINED)
    {
    fprintf(stderr,"<c_ezsint_orig> Source or target grid undefined! Aborting...\n");
    return -1;
    }
  
  
  gdin = iset_gdin;
  gdout= iset_gdout;
  
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   
  if (iset_gdin == iset_gdout)
    {
    memcpy(zout, zin, Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    return 1;
    }
  
  if (Grille[gdrow_in][gdcol_in].fst.axe_y_inverse == 1)
    {
    lzin = (ftnfloat *) malloc(Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    memcpy(lzin, zin, Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    f77name(permut)(lzin, &Grille[gdrow_in][gdcol_in].ni, &Grille[gdrow_in][gdcol_in].nj);
    }
  else
    {
    lzin = zin;
    }
  
  if (Grille[gdrow_in][gdcol_in].needs_expansion == OUI)
    {
    lxzin = (ftnfloat *) malloc(2*Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    ez_xpnsrcgd(gdin, lxzin, lzin);
    }
  else
    {
    lxzin = lzin;
    }
  
  ier = ez_calclatlon(gdout);
  ier = ez_calcxy(gdin, gdout);
  npts = Grille[gdrow_out][gdcol_out].ni*Grille[gdrow_out][gdcol_out].nj;
  
  ier = ez_interp(zout, lxzin, gdin, gdout);
  
  if (groptions.polar_correction == OUI)
    {
    ier = ez_defzones(gdin, gdout);
    ierc= ez_corrval(zout, lxzin, gdin, gdout);
    }
  
  if (lzin != zin && lzin != NULL)
    {
    free(lzin);
    }
  
  if (lxzin != lzin && lxzin != zin && lxzin != NULL)
    {
    free(lxzin);
    }
  
  return ierc;
}
