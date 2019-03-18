/*
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
wordint c_ezyy_use1subgrid(wordint gdout,wordint gdin,wordint *yin2yin,wordint *yan2yin,wordint *yin2yan,wordint *yan2yan)
{
  wordint icode,i,j,k,ni,nj,ninj,inyin,inyan;
  wordint yin_gdin,yan_gdin,yin_gdout,yan_gdout,yyin,yyout;
  wordint yin_gdrow_in, yin_gdcol_in, yin_gdrow_out, yin_gdcol_out;
  wordint yan_gdrow_in, yan_gdcol_in, yan_gdrow_out, yan_gdcol_out;
  wordint     gdrow_in,     gdcol_in,     gdrow_out,     gdcol_out;
  ftnfloat *dlat,*dlon,*xpx,*ypx;
  
  _Grille *lgdin, *lgdin_yan, *lgdout;
 /*  need only access to either Yin or Yang info for the lat and lon val */
   
  yyin=0; yyout=0;
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
  yyin=1;
  yin_gdin = Grille[gdrow_in][gdcol_in].subgrid[0];
  yan_gdin = Grille[gdrow_in][gdcol_in].subgrid[1];
  c_gdkey2rowcol(yin_gdin,  &yin_gdrow_in,  &yin_gdcol_in);
  c_gdkey2rowcol(yan_gdin,  &yan_gdrow_in,  &yan_gdcol_in);
  if (Grille[gdrow_out][gdcol_out].nsubgrids > 0)
     {/* will define yin2yin,yan2yin,yin2yan,yan2yan flags */
     yyout=1;
     yin_gdout = Grille[gdrow_out][gdcol_out].subgrid[0];
     yan_gdout = Grille[gdrow_out][gdcol_out].subgrid[1];
     c_gdkey2rowcol(yin_gdout,  &yin_gdrow_out,  &yin_gdcol_out);
     c_gdkey2rowcol(yan_gdout,  &yan_gdrow_out,  &yan_gdcol_out);
     }
  else
     { /* will define yin2yin,yan2yin */
     yin_gdout = gdout;
     yin_gdrow_out = gdrow_out;
     yin_gdcol_out = gdcol_out;
     }

  lgdin = &(Grille[yin_gdrow_in ][yin_gdcol_in ]);
  lgdin_yan = &(Grille[yan_gdrow_in ][yan_gdcol_in ]);
  lgdout= &(Grille[gdrow_out][gdcol_out]);
  if (yyout == 1)
     {
     lgdout = &(Grille[yin_gdrow_out ][yin_gdcol_out ]);
     }
  k=0;
  ni= lgdout->ni;
  nj= lgdout->nj;   
  ninj= ni*nj;
  dlat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  dlon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  xpx  = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  ypx  = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));

  if (yyout == 0)
     {
     inyin= 1;
     inyan= 0;
     icode = c_gdll(yin_gdout,dlat,dlon);
     icode = c_gdxyfll_orig(yin_gdin,xpx,ypx,dlat,dlon,ninj);
     for (j=0; j< nj; j++)
         {
         for (i=0; i<ni; i++)
             {
             k=(j*ni)+i;
             if (xpx[k] < 1.0 || xpx[k] > ni || ypx[k] < 1.0 || ypx[k] > nj)
                {
                inyin = 0;
                inyan = 1;
                break;
                }
             }
         }
     if (inyin == 0) /* check with YANG grid */
         {
         icode = c_gdxyfll_orig(yan_gdin,xpx,ypx,dlat,dlon,ninj);
         for (j=0; j< nj; j++)
             {
             for (i=0; i<ni; i++)
                 {
                 k=(j*ni)+i;
                 if (xpx[k] < 1.0 || xpx[k] > ni || ypx[k] < 1.0 || ypx[k] > nj)
                    {
                    inyan = 0;
                    break;
                    }
                 }
             }
         }
     free(dlat);
     free(dlon);
     free(xpx);
     free(ypx);
     if (inyin == 1) *yin2yin=1;
     if (inyan == 1) *yan2yin=1;
     if (inyin == 0 && inyan == 0)
        {
        fprintf(stderr,"<ezyy1subgrid>: Option Use1Subgrid selected. Cannot find coverage within one subgrid,Exiting...\n\n");
        return -1;
        }
     }
  if (yyout == 1) /* output is a YIN-YANG grid */
     {
     dlat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
     dlon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
     xpx  = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
     ypx  = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
     inyin= 1;
     inyan= 0;
     icode = c_gdll(yan_gdout,dlat,dlon);
     icode = c_gdxyfll_orig(yin_gdin,xpx,ypx,dlat,dlon,ninj);
     for (j=0; j< nj; j++)
         {
         for (i=0; i<ni; i++)
             {
             k=(j*ni)+i;
             if (xpx[k] < 1.0 || xpx[k] > ni || ypx[k] < 1.0 || ypx[k] > nj)
                {
                inyin = 0;
                inyan = 1;
                break;
                }
             }
         }
     if (inyin == 0) /* check other grid */
         {
         icode = c_gdxyfll_orig(yan_gdin,xpx,ypx,dlat,dlon,ninj);
         for (j=0; j< nj; j++)
             {
             for (i=0; i<ni; i++)
                 {
                 k=(j*ni)+i;
                 if (xpx[k] < 1.0 || xpx[k] > ni || ypx[k] < 1.0 || ypx[k] > nj)
                    {
                    inyan = 0;
                    break;
                    }
                 }
             }
         }
     free(dlat);
     free(dlon);
     free(xpx);
     free(ypx);
     if (inyin == 1) *yin2yan=1;
     if (inyan == 1) *yan2yan=1;
     if (inyin == 0 && inyan == 0)
        {
        fprintf(stderr,"<ezyy1subgrid>: Option Use1Subgrid selected. Cannot find coverage within one subgrid,Exiting...\n\n");
        return -1;
        }
     }
}

