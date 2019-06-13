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
wordint f77name(ezuvint)(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin)
{
   wordint icode;
   
   icode = c_ezuvint(uuout, vvout, uuin, vvin);
   return icode;
}

wordint c_ezuvint(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin)
{
  wordint icode,gdin,gdout;
   
  wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out;
   
  gdin = iset_gdin;
  gdout= iset_gdout;
  
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);

  if (Grille[gdrow_in][gdcol_in].nsubgrids > 0 || Grille[gdrow_out][gdcol_out].nsubgrids > 0)
      {
      icode = c_ezyyuvint(uuout,vvout,uuin,vvin,gdout,gdin);
      iset_gdin=gdin;
      iset_gdout=gdout;
      return icode;
      }
  icode = c_ezuvint_orig(uuout, vvout, uuin,vvin);
  return icode;

}

wordint c_ezuvint_orig(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin)
{
  wordint gdin,gdout;
  wordint npts, ier, ierc,ierc1,ierc2;
  ftnfloat *uullout = NULL;
  ftnfloat *vvllout = NULL;
   
  wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out;
   
  gdin = iset_gdin;
  gdout= iset_gdout;
  ierc = 0;
  ierc1 = 0;
  ierc2 = 0;
  
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
  
  npts = Grille[gdrow_out][gdcol_out].ni*Grille[gdrow_out][gdcol_out].nj;
  ier = ez_calclatlon(gdout);
  
  groptions.vecteur = VECTEUR;
  
  groptions.symmetrie = SYM;
  ierc1=c_ezsint(uuout,uuin);
  groptions.symmetrie = ANTISYM;
  ierc2=c_ezsint(vvout,vvin);
  groptions.symmetrie = SYM;
  if (ierc1 == 2 || ierc2 == 2) 
      {
      ierc = 2;
      }
  
  if (groptions.polar_correction == OUI)
    {
    ier=ez_corrvec(uuout, vvout, uuin, vvin, gdin, gdout);
    }
  
  uullout = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
  vvllout = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
  
  c_gdwdfuv(gdin, uullout, vvllout, uuout, vvout,
            Grille[gdrow_out][gdcol_out].lat, Grille[gdrow_out][gdcol_out].lon, npts);
  c_gduvfwd(gdout, uuout, vvout, uullout, vvllout,
            Grille[gdrow_out][gdcol_out].lat, Grille[gdrow_out][gdcol_out].lon, npts);
  
  groptions.vecteur = SCALAIRE;
  free(uullout);
  free(vvllout);
  
  return ierc;
}

