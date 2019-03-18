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

void ez_xpnsrcgd(wordint gdid, ftnfloat *zout, ftnfloat *zin)
{
   _Grille gr;
  wordint gdrow_id, gdcol_id;
    
  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
   gr = Grille[gdrow_id][gdcol_id];
   
   switch (gr.grtyp[0])
     {
     case 'A':
     case 'G':
       f77name(ez_xpngdag2)(zout,zin,&gr.ni,&gr.nj,&gr.j1,&gr.j2,&gr.fst.ig[IG1],&groptions.symmetrie);
       break;

     case 'B':
       f77name(ez_xpngdb2)(zout,zin,&gr.ni,&gr.nj,&gr.j1,&gr.j2,&gr.fst.ig[IG1],&groptions.symmetrie);
       break;

     default:
       break;
     }
}
