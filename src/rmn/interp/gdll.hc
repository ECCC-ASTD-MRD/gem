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
wordint f77name(ezll)(wordint *gdid, ftnfloat *lat, ftnfloat *lon)
{
   wordint icode;
   
   icode = c_gdll(*gdid, lat, lon);
   return icode;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdll)(wordint *gdid, ftnfloat *lat, ftnfloat *lon)
{
   wordint icode;
   
   icode = c_gdll(*gdid, lat, lon);
   return icode;
}

wordint c_gdll(wordint gdid, ftnfloat *lat, ftnfloat *lon)
{
   wordint icode;
  wordint gdrow_id, gdcol_id, ni, nj;
  wordint yin_gdid,yan_gdid;
  wordint yin_gdrow_id, yin_gdcol_id;
  wordint yan_gdrow_id, yan_gdcol_id;
    
  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
  if (Grille[gdrow_id][gdcol_id].nsubgrids > 0 )
    {
    yin_gdid = Grille[gdrow_id][gdcol_id].subgrid[0];
    yan_gdid = Grille[gdrow_id][gdcol_id].subgrid[1];
/*    printf("gdll: gdid for yin=%d,gdid for yan=%d\n",yin_gdid,yan_gdid); */
    c_gdkey2rowcol(yin_gdid, &yin_gdrow_id, &yin_gdcol_id);
    c_gdkey2rowcol(yan_gdid, &yan_gdrow_id, &yan_gdcol_id);
    ni = Grille[yin_gdrow_id][yin_gdcol_id].ni;
    nj = Grille[yin_gdrow_id][yin_gdcol_id].nj;
/*  printf("gdll: ni=%d, nj=%d\n",ni,nj); */
    icode=c_gdll_orig(yin_gdid,lat,lon);
    icode=c_gdll_orig(yan_gdid,&lat[ni*nj],&lon[ni*nj]);
    }
  else
    {
    icode=c_gdll_orig(gdid,lat,lon);
    }
  return icode;
}

wordint c_gdll_orig(wordint gdid, ftnfloat *lat, ftnfloat *lon)
{
  wordint gdrow_id, gdcol_id;
    
  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
   
   ez_calclatlon(gdid);
   if (Grille[gdrow_id][gdcol_id].flags & LAT)
      {
      memcpy(lon, Grille[gdrow_id][gdcol_id].lon, Grille[gdrow_id][gdcol_id].ni*Grille[gdrow_id][gdcol_id].nj*sizeof(ftnfloat));
      if (Grille[gdrow_id][gdcol_id].fst.axe_y_inverse == 0)
         {
         memcpy(lat, Grille[gdrow_id][gdcol_id].lat, Grille[gdrow_id][gdcol_id].ni*Grille[gdrow_id][gdcol_id].nj*sizeof(ftnfloat));
         }
      else
         {
         memcpy(lat, Grille[gdrow_id][gdcol_id].lat, Grille[gdrow_id][gdcol_id].ni*Grille[gdrow_id][gdcol_id].nj*sizeof(ftnfloat));
         }
      }
   else
      {
      fprintf(stderr, "Erreur! A l'aide! Descripteurs manquants!\n");
      return -1;
      }
   return 0;
}
