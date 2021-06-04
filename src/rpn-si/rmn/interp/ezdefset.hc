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

void reallocate_gridset_table(int gdid);
void   allocate_gridset_table(int gdid);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezdefset)(wordint *gdout, wordint *gdin)
{
   wordint icode;

   icode = c_ezdefset(*gdout, *gdin);
   return icode;
}

wordint c_ezdefset(wordint gdout, wordint gdin)
{
  /* d'abord trouver si l'ensemble est deja defini */

   wordint i;
   wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out, npts, cur_gdin;
   int lcl_ngdin, idx_gdin;
   static wordint found = -1;
   static wordint ncalls = 0;
   wordint nsets = 0;
   _Grille *gr;

   if (gdout == UNDEFINED)
   {
      if (gdin != UNDEFINED)
      {
         gdout = gdin;
      }
   }

   if (gdout == gdin)
   {
      iset_gdin = gdin;
      iset_gdout = gdout;
   }


   c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
   c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   gr = &(Grille[gdrow_out][gdcol_out]);

   nsets = Grille[gdrow_out][gdcol_out].n_gdin;

   if (nsets == 0)
   {
      allocate_gridset_table(gdout);
      gr->log_chunk_gdin = cur_log_chunk;
   }

   if (nsets >= primes[MAX_LOG_CHUNK])
   {
      for (i=0; i < nsets; i++)
      {
         if (gr->gset[i].gdin != -1)
         {
            c_ezfreegridset(gdin, i);
         }
      }
      nsets = 0;
      allocate_gridset_table(gdout);
      gr->log_chunk_gdin = cur_log_chunk;
   }


   found = -1;

   idx_gdin = gdin % primes[gr->log_chunk_gdin];
   if (gr->gset[idx_gdin].gdin == gdin)
   {
      found = 1;
      iset_gdin = gdin;
      iset_gdout = gdout;
      return 1;
   }

   i = idx_gdin;
   while ((found == -1) && (i != idx_gdin-1) && (gr->gset[i].gdin != -1))
   {
      if (gdin == gr->gset[i].gdin)
      {
         found = i;
         iset_gdin = gdin;
         iset_gdout = gdout;
         Grille[gdrow_out][gdcol_out].idx_last_gdin = Grille[gdrow_out][gdcol_out].gset[i].gdin;
         return 1;
      }
      else
      {
         i++;
         if (0 == (i % (primes[gr->log_chunk_gdin])))
         {
            i = 0;
         }
      }
   }

   /* si on se rend jusqu'ici alors c'est que le set n'a pas ete trouve */

   /* On initialise le vecteur de gdout pour lequel gdin est utilise en entree
      Ceci sera utile si le vecteur de grilles deborde */

   gr->gset[i].gdin = gdin;
   cur_gdin = gdin;
   gr->n_gdin++;

   npts = gr->ni * gr->nj;
   gr->gset[i].x = malloc (sizeof(ftnfloat)*npts);
   gr->gset[i].y = malloc (sizeof(ftnfloat)*npts);
   gr->gset[i].use_sincos_cache = NON;

   if (gr->n_gdin >= (primes[gr->log_chunk_gdin]/2))
   {
      reallocate_gridset_table(gdout);
   }

   if (Grille[gdrow_in][gdcol_in].n_gdin_for == 0)
   {
      Grille[gdrow_in][gdcol_in].gdin_for = malloc(CHUNK *sizeof(int));
      for (i=0; i < CHUNK; i++)
         {
            Grille[gdrow_in][gdcol_in].gdin_for[i] = -1;
         }
      Grille[gdrow_in][gdcol_in].gdin_for[0] = gdout;
      Grille[gdrow_in][gdcol_in].n_gdin_for++;
   }
   else
   {
      if (0 == (Grille[gdrow_in][gdcol_in].n_gdin_for % CHUNK))
      {
         Grille[gdrow_in][gdcol_in].gdin_for = (int *) realloc(Grille[gdrow_in][gdcol_in].gdin_for, (Grille[gdrow_in][gdcol_in].n_gdin_for+CHUNK)*sizeof(int));
      }
      Grille[gdrow_in][gdcol_in].gdin_for[Grille[gdrow_in][gdcol_in].n_gdin_for] = gdout;
      Grille[gdrow_in][gdcol_in].n_gdin_for++;
   }

   if (groptions.verbose > 0)
   {
      printf("gdin : %d gdout: %d\n", gdin, gdout);
      printf("cur_gdin                           = %03d\n", cur_gdin);
      printf("n_gdin                             = %03d\n", Grille[gdrow_out][gdcol_out].n_gdin);
      printf("Grille[%03d][%03d].gset[%03d].gdin = %d\n", gdrow_out, gdcol_out, cur_gdin, gdin);
   }
   iset_gdin = gdin;
   iset_gdout = gdout;
   return 1;
}

void reallocate_gridset_table(int gdid)
{
   wordint i;
   wordint gdrow_id, gdrow_out, gdcol_id, gdcol_out, npts, cur_gdin, log_chunk_gdin;
   wordint newIndex, inserted, curIndex;
   int lcl_ngdin, oldChunkSize, newChunkSize;
   int cur_chunk;
   static wordint found = -1;
   static wordint ncalls = 0;
   _Grille *gr;
   _gridset *gset, *newTable;

   c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
   gr = &(Grille[gdrow_id][gdcol_id]);
   cur_chunk = gr->log_chunk_gdin;
   oldChunkSize = primes[cur_chunk];
   newChunkSize = primes[cur_chunk + 1];
   newTable = (_gridset *) calloc(sizeof(_gridset), newChunkSize);

   for (i=0 ; i < newChunkSize; i++)
   {
      newTable[i].gdin = -1;
   }

   for (i=0 ; i < oldChunkSize; i++)
   {
      if (gr->gset[i].gdin != -1)
      {
         newIndex = gr->gset[i].gdin % newChunkSize;
         if (newTable[newIndex].gdin == -1)
         {
            memcpy(&(newTable[newIndex]), &(gr->gset[i]), sizeof(_gridset));
         }
         else
         {
            newIndex++;
            inserted = -1;
            while (inserted == -1) /** && curIndex != (newIndex-1)) (a reverifier-- Yves-20120228)**/
            {
               fprintf(stderr, "reallocate_gridset_table -- should not be here\n ");
               if (newTable[newIndex].gdin == -1)
               {
                  memcpy(&(newTable[newIndex]), &(gr->gset[i]), sizeof(_gridset));
                  inserted = 1;
               }
               else
               {
                  newIndex++;
                  if (0 == (newIndex % newChunkSize))
                  {
                     newIndex = 0;
                  }
               }
            }
         }
      }
   }

   free(gr->gset);
   gr->gset = newTable;
   gr->log_chunk_gdin++;
}

void allocate_gridset_table(int gdid)
{
   wordint i;
   wordint gdrow_id, gdrow_out, gdcol_id, gdcol_out, npts, cur_gdin, log_chunk_gdin;
   wordint newIndex, inserted, curIndex;
   int lcl_ngdin, chunkSize;
   int cur_chunk;
   _Grille *gr;
   _gridset *gset, *newTable;

   c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
   gr = &(Grille[gdrow_id][gdcol_id]);
   chunkSize = primes[cur_log_chunk];
   gr->gset = (_gridset *) calloc(sizeof(_gridset), chunkSize);

   for (i=0 ; i < chunkSize; i++)
   {
      gr->gset[i].gdin = -1;
   }
}
