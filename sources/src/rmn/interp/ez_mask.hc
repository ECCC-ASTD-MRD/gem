/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "ezscint.h"
#include "ez_funcdef.h"

#define BITPOS(i) (i - ((i >> 5) << 5))
#define GETMSK(fld,i) ((fld[i >> 5]  & (1 << BITPOS(i))) >> BITPOS(i))
#define SETMSK(fld,i) (fld[i >> 5] | (fld[i] << BITPOS(i)))


int f77name(gdsetmask)(int *gdid, int *mask)
   {
   return c_gdsetmask(*gdid, mask);
   }

int f77name(gdgetmask)(int *gdid, int *mask)
   {
   return c_gdgetmask(*gdid, mask);
   }

int f77name(ezsint_m)(float *zout, float *zin)
   {
   return c_ezsint_m(zout, zin);
   }

int f77name(ezuvint_m)(float *uuout, float *vvout, float *uuin, float *vvin)
   {
   return c_ezuvint_m(uuout, vvout, uuin, vvin);
   }

int f77name(ezsint_mdm)(float *zout, int *mask_out, float *zin, int *mask_in)
   {
   return c_ezsint_mdm(zout, mask_out, zin, mask_in);
   }

int f77name(ezuvint_mdm)(float *uuout, float *vvout, int *mask_out, float *uuin, float *vvin, int *mask_in)
   {
   return c_ezuvint_mdm(uuout, vvout, mask_out, uuin, vvin, mask_in);
   }

int f77name(ezsint_mask)(int *mask_out, int *mask_in)
   {
   return c_ezsint_mask(mask_out, mask_in);
   }

/* ------------------------------------------------------------------------- */

int c_gdsetmask(int gdid, int *mask)
   {
   int gdrow, gdcol;
   int ni, nj;
   c_gdkey2rowcol(gdid, &gdrow, &gdcol);
   if (Grille[gdrow][gdcol].nsubgrids > 0)
      {
       fprintf(stderr, "<gdsetmask> This operation is not supported for 'U' grids.\n");
       return -1;
      }
   ni = Grille[gdrow][gdcol].ni;
   nj = Grille[gdrow][gdcol].nj;

   if (Grille[gdrow][gdcol].mask != NULL)
      {
      free(Grille[gdrow][gdcol].mask);
      }
   Grille[gdrow][gdcol].mask = (int *)malloc(ni*nj*sizeof(int));
   memcpy(Grille[gdrow][gdcol].mask, mask, ni*nj*sizeof(int));
   return 0;
   }


/* ------------------------------------------------------------------------- */

int c_gdgetmask(int gdid, int *mask)
   {
   int gdrow, gdcol;
   int ni, nj;

   c_gdkey2rowcol(gdid, &gdrow, &gdcol);
   if (Grille[gdrow][gdcol].nsubgrids > 0)
      {
       fprintf(stderr, "<gdgetmask> This operation is not supported for 'U' grids.\n");
       return -1;
      }
   ni = Grille[gdrow][gdcol].ni;
   nj = Grille[gdrow][gdcol].nj;

   if (Grille[gdrow][gdcol].mask != NULL)
      {
      memcpy(mask, Grille[gdrow][gdcol].mask, ni*nj*sizeof(int));
      return 0;
      }
   else
      {
      mask = NULL;
      return -1;
      }
   }


/* ------------------------------------------------------------------------- */

int c_ezsint_m(float *zout, float *zin)
   {
       fprintf(stderr, "<ezsint_m> This operation is currently not implemented.\n");
    return 0;
   } 


/* ------------------------------------------------------------------------- */

int c_ezuvint_m(float *uuout, float *vvout, float *uuin, float *vvin)
   {
       fprintf(stderr, "<ezuvint_m> This operation is currently not implemented.\n");
    return 0;
   }


/* ------------------------------------------------------------------------- */

int c_ezsint_mdm(float *zout, int *mask_out, float *zin, int *mask_in)
   {
   wordint gdin, gdout, gdrow_out, gdcol_out;
   wordint              gdrow_in,  gdcol_in;
   wordint methode = 2;
   wordint ni_out, nj_out;

   gdin = c_ezgetgdin();
   gdout = c_ezgetgdout();

   c_ezdefset(gdout, gdin);

   c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   c_gdkey2rowcol(gdin, &gdrow_in, &gdcol_in);
   if (Grille[gdrow_out][gdcol_out].nsubgrids > 0 || 
       Grille[gdrow_in][gdcol_in].nsubgrids > 0)
      {
       fprintf(stderr, "<ezsint_mdm> This operation is not supported for 'U' grids.\n");
       return -1;
      }
   ni_out = Grille[gdrow_out][gdcol_out].ni;
   nj_out = Grille[gdrow_out][gdcol_out].nj;
   c_ezsint(zout, zin);
   c_ezsint_mask(mask_out, mask_in);
   f77name(lorenzo_mask_fill)(zout, mask_out, &ni_out, &nj_out, &methode);
   return 0;

   }


/* ------------------------------------------------------------------------- */

int c_ezuvint_mdm(float *uuout, float *vvout, int *mask_out, float *uuin, float *vvin, int *mask_in)
   {
   wordint gdin, gdout, gdrow_out, gdcol_out;
   wordint              gdrow_in,  gdcol_in;
   wordint methode = 2;
   wordint ni_out, nj_out;

   gdin = c_ezgetgdin();
   gdout = c_ezgetgdout();

   c_ezdefset(gdout, gdin);

   c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   c_gdkey2rowcol(gdin, &gdrow_in, &gdcol_in);
   if (Grille[gdrow_out][gdcol_out].nsubgrids > 0 || 
       Grille[gdrow_in][gdcol_in].nsubgrids > 0)
      {
       fprintf(stderr, "<ezuvint_mdm> This operation is not supported for 'U' grids.\n");
       return -1;
      }
   ni_out = Grille[gdrow_out][gdcol_out].ni;
   nj_out = Grille[gdrow_out][gdcol_out].nj;
   c_ezsint_mask(mask_out, mask_in);
   c_ezuvint(uuout, vvout, uuin, vvin);
   f77name(lorenzo_mask_fill)(uuout, mask_out, &ni_out, &nj_out, &methode);
   f77name(lorenzo_mask_fill)(vvout, mask_out, &ni_out, &nj_out, &methode);
   return 0;
   }


/* ------------------------------------------------------------------------- */

int c_ezsint_mask(int *mask_out, int *mask_in)
   {
   char grtyp_in[2], grtyp_out[2];
   int ni_gdin, ni_gdout, nj_gdin, nj_gdout;
   int ig1_gdin, ig2_gdin, ig3_gdin, ig4_gdin, ig1_gdout, ig2_gdout, ig3_gdout, ig4_gdout;
   int i, j, k, ier,npts_in, npts_out, idx_gdin, gdrow_out, gdcol_out;
   wordint              gdrow_in,  gdcol_in;
   unsigned int bitpos;
   float *fmask_in, *fmask_out, *x, *y;
   char current_option[32], interp_degree[32];
   _ygrid *ygrid;

   wordint gdin, gdout;

   gdin = c_ezgetgdin();
   gdout = c_ezgetgdout();
   c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   c_gdkey2rowcol(gdin, &gdrow_in, &gdcol_in);
   if (Grille[gdrow_out][gdcol_out].nsubgrids > 0 || 
       Grille[gdrow_in][gdcol_in].nsubgrids > 0)
      {
       fprintf(stderr, "<ezsint_mask> This operation is not supported for 'U' grids.\n");
       return -1;
      }

   c_ezdefset(gdout, gdin);
   idx_gdin = c_find_gdin(gdin, gdout);
   ier = c_ezgprm(gdin, grtyp_in, &ni_gdin, &nj_gdin, &ig1_gdin, &ig2_gdin, &ig3_gdin, &ig4_gdin);
   ier = c_ezgprm(gdout, grtyp_out, &ni_gdout, &nj_gdout, &ig1_gdout, &ig2_gdout, &ig3_gdout, &ig4_gdout);


   npts_in  = ni_gdin*nj_gdin;
   npts_out = ni_gdout*nj_gdout;

   if (grtyp_in[0] == 'Y')
      {
      ygrid = &(Grille[gdrow_out][gdcol_out].gset[idx_gdin].ygrid);
      memcpy(mask_out, ygrid->mask, ni_gdout*nj_gdout*sizeof(int));
      }
   else
      {
      x = (float *) Grille[gdrow_out][gdcol_out].gset[idx_gdin].x;
      y = (float *) Grille[gdrow_out][gdcol_out].gset[idx_gdin].y;
      f77name(qqq_ezsint_mask)(mask_out, x, y, &ni_gdout, &nj_gdout, mask_in, &ni_gdin, &nj_gdin);
      }
   return 0;
   }


/* ------------------------------------------------------------------------- */

int f77name(ezget_mask_zones)(int *mask_out, int *mask_in)
   {
    return c_ezget_mask_zones(mask_out, mask_in);
   }


/* ------------------------------------------------------------------------- */

int c_ezget_mask_zones(int *mask_out, int *mask_in)
   {
   char grtyp_in[2], grtyp_out[2];
   int ni_gdin, ni_gdout, nj_gdin, nj_gdout;
   int ig1_gdin, ig2_gdin, ig3_gdin, ig4_gdin, ig1_gdout, ig2_gdout, ig3_gdout, ig4_gdout;
   int i, j, k, ier,npts_in, npts_out, idx_gdin, gdrow_out, gdcol_out;
   wordint              gdrow_in,  gdcol_in;
   unsigned int bitpos;
   float *x, *y;
   char current_option[32], interp_degree[32];

   wordint gdin, gdout;

   strcpy(interp_degree,"interp_degree");
   gdin = c_ezgetgdin();
   gdout = c_ezgetgdout();
   c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   c_gdkey2rowcol(gdin, &gdrow_in, &gdcol_in);
   if (Grille[gdrow_out][gdcol_out].nsubgrids > 0 || 
       Grille[gdrow_in][gdcol_in].nsubgrids > 0)
      {
       fprintf(stderr, "<ezget_mask_zones> This operation is not supported for 'U' grids.\n");
       return -1;
      }

   c_ezdefset(gdout, gdin);
   idx_gdin = c_find_gdin(gdin, gdout);
   ier = c_ezgprm(gdin, grtyp_in, &ni_gdin, &nj_gdin, &ig1_gdin, &ig2_gdin, &ig3_gdin, &ig4_gdin);
   ier = c_ezgprm(gdout, grtyp_out, &ni_gdout, &nj_gdout, &ig1_gdout, &ig2_gdout, &ig3_gdout, &ig4_gdout);

   npts_in  = ni_gdin*nj_gdin;
   npts_out = ni_gdout*nj_gdout;

    x = (float *) Grille[gdrow_out][gdcol_out].gset[idx_gdin].x;
    y = (float *) Grille[gdrow_out][gdcol_out].gset[idx_gdin].y;

   f77name(qqq_ezget_mask_zones)(mask_out, x, y, &ni_gdout, &nj_gdout, mask_in, &ni_gdin, &nj_gdin);
    return 0;
   }


/* ------------------------------------------------------------------------- */


