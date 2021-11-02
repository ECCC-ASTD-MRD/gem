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
wordint f77name(ezgfstp)(wordint *gdid, 
                         char *nomvarx, char *typvarx, char *etiketx,
                         char *nomvary, char *typvary, char *etikety,
                         wordint *ip1, wordint *ip2, wordint *ip3, wordint *dateo, 
                         wordint *deet, wordint *npas, wordint *nbits,
                         F2Cl lennomvarx, F2Cl lentypvarx, F2Cl lenetiketx,
                         F2Cl lennomvary, F2Cl lentypvary, F2Cl lenetikety)
{
   wordint i;
   wordint icode;
   
   char lnomvarx[5],lnomvary[5],letiketx[13],letikety[13],ltypvarx[2],ltypvary[2];

   c_ezgfstp(*gdid, lnomvarx, ltypvarx, letiketx, lnomvary, ltypvary,letikety,
             ip1, ip2, ip3, dateo, deet, npas, nbits);

   for (i=0; i < lennomvarx; i++)
      {
      nomvarx[i] = ' ';
      nomvary[i] = ' ';
      }
   
   for (i=0; i < lentypvarx; i++)
      {
      typvarx[i] = ' ';
      typvary[i] = ' ';
      }
   
   for (i=0; i < lenetiketx; i++)
      {
      etiketx[i] = ' ';
      etikety[i] = ' ';
      }
   
   strncpy(nomvarx, lnomvarx, 4);
   strncpy(nomvary, lnomvary, 4);
   strncpy(typvarx, ltypvarx, 2);
   strncpy(typvary, ltypvary, 2);
   strncpy(etiketx, letiketx, 12);
   strncpy(etikety, letikety, 12);

   icode = 0;
   return icode;
}

wordint c_ezgfstp(wordint gdid, char *nomvarx, char *typvarx, char *etiketx, 
              char *nomvary, char *typvary, char *etikety,
              wordint *ip1, wordint *ip2, wordint *ip3, wordint *dateo, wordint *deet, wordint *npas, wordint *nbits)
{
  wordint gdrow_id, gdcol_id;
    
  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
   
   strcpy(nomvarx, Grille[gdrow_id][gdcol_id].fst.nomvarx);
   strcpy(typvarx, Grille[gdrow_id][gdcol_id].fst.typvarx);
   strcpy(etiketx, Grille[gdrow_id][gdcol_id].fst.etiketx);

   strcpy(nomvary, Grille[gdrow_id][gdcol_id].fst.nomvary);
   strcpy(typvary, Grille[gdrow_id][gdcol_id].fst.typvary);
   strcpy(etikety, Grille[gdrow_id][gdcol_id].fst.etikety);

   *ip1 = Grille[gdrow_id][gdcol_id].fst.ip1;
   *ip2 = Grille[gdrow_id][gdcol_id].fst.ip2;
   *ip3 = Grille[gdrow_id][gdcol_id].fst.ip3;

   *dateo = Grille[gdrow_id][gdcol_id].fst.date;
   *deet = Grille[gdrow_id][gdcol_id].fst.deet;
   *npas = Grille[gdrow_id][gdcol_id].fst.npas;
   *nbits = Grille[gdrow_id][gdcol_id].fst.nbits;

   return 0;
}

