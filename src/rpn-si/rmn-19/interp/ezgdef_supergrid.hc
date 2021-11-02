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
wordint f77name(ezgdef_supergrid)(wordint *ni, wordint *nj, char *grtyp, char *grref, wordint *vercode, wordint *nsubgrids, wordint *subgrid, F2Cl lengrtyp, F2Cl lengrref)
{
  wordint gdid,i;
  char lgrtyp[2],lgrref[2];

  lgrtyp[0] = grtyp[0];
  lgrtyp[1] = '\0';
   
  lgrref[0] = grref[0];
  lgrref[1] = '\0';
  gdid = c_ezgdef_supergrid(*ni, *nj, lgrtyp, lgrref, *vercode, *nsubgrids,subgrid);
  return gdid;
}

wordint c_ezgdef_supergrid(wordint ni, wordint nj, char *grtyp, char *grref, wordint vercode,wordint nsubgrids, wordint *subgrid)
{
  wordint sub_gdrow_id, sub_gdcol_id, sub_gdid;
  wordint mask_gdrow_id, mask_gdcol_id, mask_gdid;
  wordint  ii,i,j,k,i0,j0,i1,j1,yni,ynj,gdid;
  wordint gdrow_in, gdcol_in, grid_index,newgrsize;
  ftnfloat *ax,*ay;
  unsigned int grid_crc;
   _Grille newgr,*maskgd;
    
  if (nsubgrids <= 1)
    {
    fprintf(stderr,"<c_ezgdef_supergrid> nsubgrids given is less than 2! Aborting...\n");
    return -1;
    }
  if (vercode != 1)
    {
    fprintf(stderr,"<c_ezgdef_supergrid> invalid vercode! Aborting...\n");
    return -1;
    }
  memset(&newgr, (int)0, sizeof(_Grille));
  strcpy(newgr.fst.etiketx, "            ");
  strcpy(newgr.fst.etikety, "            ");
  strcpy(newgr.fst.typvarx, "  ");
  strcpy(newgr.fst.typvary, "  ");
  strcpy(newgr.fst.nomvarx, "^>  ");
  strcpy(newgr.fst.nomvary, "^>  ");
  if (vercode == 1)
    {
    sub_gdid=subgrid[0];
    c_gdkey2rowcol(sub_gdid,  &sub_gdrow_id,  &sub_gdcol_id);

/*    strcpy(newgr.grtyp, grtyp);
    strcpy(newgr.grref, grref);
*/
    newgr.grtyp[0]=grtyp[0];
    newgr.grref[0]=grref[0];
    RemplirDeBlancs(newgr.fst.nomvarx, 5);
    RemplirDeBlancs(newgr.fst.typvarx, 3);
    RemplirDeBlancs(newgr.fst.etiketx, 13);
    RemplirDeBlancs(newgr.fst.nomvary, 5);
    RemplirDeBlancs(newgr.fst.typvary, 3);
    RemplirDeBlancs(newgr.fst.etikety, 13);
    newgr.ni = ni;
    newgr.nj = nj;
    newgr.idx_last_gdin = -1;
    /* create tictac arrays to add uniqueness in supergrid*/
    ax = (ftnfloat *) malloc(newgr.ni*sizeof(ftnfloat));
    ay = (ftnfloat *) malloc(newgr.nj*sizeof(ftnfloat));
    memcpy(ax,Grille[sub_gdrow_id][sub_gdcol_id].ax,newgr.ni*sizeof(ftnfloat));
    memcpy(ay,Grille[sub_gdrow_id][sub_gdcol_id].ay,newgr.nj*sizeof(ftnfloat));
    newgr.fst.ip1      = Grille[sub_gdrow_id][sub_gdcol_id].fst.ip1;
    newgr.fst.ip2      = Grille[sub_gdrow_id][sub_gdcol_id].fst.ip2;
    newgr.fst.ip3      = Grille[sub_gdrow_id][sub_gdcol_id].fst.ip3;
/*   to add more uniqueness to the super-grid index Yin-Yang grid, we also
add   the rotation of YIN */
    newgr.fst.ig[IG1]  = Grille[sub_gdrow_id][sub_gdcol_id].fst.ig[IG1];
    newgr.fst.ig[IG2]  = Grille[sub_gdrow_id][sub_gdcol_id].fst.ig[IG2];
    newgr.fst.ig[IG3]  = Grille[sub_gdrow_id][sub_gdcol_id].fst.ig[IG3];
    newgr.fst.ig[IG4]  = Grille[sub_gdrow_id][sub_gdcol_id].fst.ig[IG4];
    newgr.fst.xg[IG1]  = 0.0;
    newgr.fst.xg[IG2]  = 0.0;
    newgr.fst.xg[IG3]  = 0.0;
    newgr.fst.xg[IG4]  = 0.0;
    newgr.fst.igref[IG1]=vercode;
    newgr.fst.igref[IG2]=0;
    newgr.fst.igref[IG3]=0;
    newgr.fst.igref[IG4]=0;
    newgr.nsubgrids= nsubgrids;
     }
    newgrsize = sizeof(_Grille);
    grid_crc = ez_calc_crc((int *)&newgr,&newgrsize,ax,ay,newgr.ni,newgr.nj);

    free(ax); free(ay);
    grid_index = grid_crc % primes_sq[cur_log_chunk];
    if (gr_list[grid_index] == NULL)
      {
      gdid = c_ez_addgrid(grid_index, &newgr);
      }
    else
      {
      gdid = c_ez_findgrid(grid_index, &newgr);
      if (gdid == -1)
        {
        gdid = c_ez_addgrid(grid_index, &newgr);
        }
      else
        {
        return gdid;
        }
      }
    c_gdkey2rowcol(gdid, &gdrow_in, &gdcol_in);
    strcpy(Grille[gdrow_in][gdcol_in].fst.nomvarx,newgr.fst.nomvarx);
    strcpy(Grille[gdrow_in][gdcol_in].grtyp,newgr.grtyp);
    strcpy(Grille[gdrow_in][gdcol_in].grref,newgr.grref);
    Grille[gdrow_in][gdcol_in].ni = newgr.ni;
    Grille[gdrow_in][gdcol_in].nj = newgr.nj;
    Grille[gdrow_in][gdcol_in].idx_last_gdin=-1;
    Grille[gdrow_in][gdcol_in].fst.ip1 = newgr.fst.ip1;
    Grille[gdrow_in][gdcol_in].fst.ip2 = newgr.fst.ip2;
    Grille[gdrow_in][gdcol_in].fst.ip3 = newgr.fst.ip3;
    Grille[gdrow_in][gdcol_in].fst.ig[IG1] = newgr.fst.ig[IG1];
    Grille[gdrow_in][gdcol_in].fst.ig[IG2] = newgr.fst.ig[IG2];
    Grille[gdrow_in][gdcol_in].fst.ig[IG3] = newgr.fst.ig[IG3];
    Grille[gdrow_in][gdcol_in].fst.ig[IG4] = newgr.fst.ig[IG4];
    Grille[gdrow_in][gdcol_in].fst.xg[IG1] = newgr.fst.xg[IG1];
    Grille[gdrow_in][gdcol_in].fst.xg[IG2] = newgr.fst.xg[IG2];
    Grille[gdrow_in][gdcol_in].fst.xg[IG3] = newgr.fst.xg[IG3];
    Grille[gdrow_in][gdcol_in].fst.xg[IG4] = newgr.fst.xg[IG4];
    Grille[gdrow_in][gdcol_in].fst.igref[IG1] = newgr.fst.igref[IG1];
    Grille[gdrow_in][gdcol_in].fst.igref[IG2] = newgr.fst.igref[IG2];
    Grille[gdrow_in][gdcol_in].fst.igref[IG3] = newgr.fst.igref[IG3];
    Grille[gdrow_in][gdcol_in].fst.igref[IG4] = newgr.fst.igref[IG4];
    Grille[gdrow_in][gdcol_in].nsubgrids = nsubgrids;
    Grille[gdrow_in][gdcol_in].subgrid = (wordint *) malloc(nsubgrids*sizeof(wordint));

    for (ii=0; ii < nsubgrids; ii++)
       {
       Grille[gdrow_in][gdcol_in].subgrid[ii] = subgrid[ii];
       c_gdkey2rowcol(subgrid[ii],  &sub_gdrow_id,  &sub_gdcol_id);
       c_ezgdef_yymask(&(Grille[sub_gdrow_id][sub_gdcol_id]));
       if (groptions.verbose > 0)
          {
          printf("Grille[%02d].subgrid[%d] has maskgrid=%d\n",gdid,subgrid[ii],Grille[sub_gdrow_id][sub_gdcol_id].mymaskgrid);
          }
       }

    if (groptions.verbose > 0)
    {
    printf("Grille[%02d].nomvarx=%s\n",gdid,Grille[gdrow_in][gdcol_in].fst.nomvarx);
    printf("Grille[%02d].nomvary=%s\n",gdid,Grille[gdrow_in][gdcol_in].fst.nomvary);
    printf("Grille[%02d].etikx=%s\n",gdid,Grille[gdrow_in][gdcol_in].fst.etiketx);
    printf("Grille[%02d].etiky=%s\n",gdid,Grille[gdrow_in][gdcol_in].fst.etikety);
    printf("Grille[%02d].grtyp = '%c'\n", gdid, Grille[gdrow_in][gdcol_in].grtyp[0]);
    printf("Grille[%02d].grref = '%c'\n", gdid, Grille[gdrow_in][gdcol_in].grref[0]);
    printf("Grille[%02d].ni    = %d\n",   gdid, Grille[gdrow_in][gdcol_in].ni);
    printf("Grille[%02d].nj    = %d\n",   gdid, Grille[gdrow_in][gdcol_in].nj);
    printf("Grille[%02d].ip1   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.ip1);
    printf("Grille[%02d].ip2   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.ip2);
    printf("Grille[%02d].ip3   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.ip3);
    printf("Grille[%02d].ig1   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.ig[IG1]);
    printf("Grille[%02d].ig2   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.ig[IG2]);
    printf("Grille[%02d].ig3   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.ig[IG3]);
    printf("Grille[%02d].ig4   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.ig[IG4]);
    printf("Grille[%02d].ig1ref = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.igref[IG1]);
    printf("Grille[%02d].ig2ref = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.igref[IG2]);
    printf("Grille[%02d].ig3ref = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.igref[IG3]);
    printf("Grille[%02d].ig4ref = %d\n",   gdid, Grille[gdrow_in][gdcol_in].fst.igref[IG4]);
    printf("Grille[%02d].nsubgrids = %d\n", gdid, Grille[gdrow_in][gdcol_in].nsubgrids);
    printf("Grille[%02d].subgrid[0]   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].subgrid[0]);
    printf("Grille[%02d].subgrid[1]   = %d\n",   gdid, Grille[gdrow_in][gdcol_in].subgrid[1]);

    printf("Grille[%02d].fst.xg[1]   = %f\n",   gdid, Grille[gdrow_in][gdcol_in].fst.xg[1]);
    printf("Grille[%02d].fst.xg[2]   = %f\n",   gdid, Grille[gdrow_in][gdcol_in].fst.xg[2]);
    printf("Grille[%02d].fst.xg[3]   = %f\n",   gdid, Grille[gdrow_in][gdcol_in].fst.xg[3]);
    printf("Grille[%02d].fst.xg[4]   = %f\n",   gdid, Grille[gdrow_in][gdcol_in].fst.xg[4]);
    }
    strcpy(Grille[gdrow_in][gdcol_in].fst.nomvarx, newgr.fst.nomvarx);
    strcpy(Grille[gdrow_in][gdcol_in].fst.typvarx, newgr.fst.typvarx);
    strcpy(Grille[gdrow_in][gdcol_in].fst.etiketx, newgr.fst.etiketx);
    strcpy(Grille[gdrow_in][gdcol_in].fst.nomvary, newgr.fst.nomvary);
    strcpy(Grille[gdrow_in][gdcol_in].fst.typvary, newgr.fst.typvary);
    strcpy(Grille[gdrow_in][gdcol_in].fst.etikety, newgr.fst.etikety);
    return gdid;
}

