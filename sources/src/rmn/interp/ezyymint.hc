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

wordint c_ezyymint(wordint gdout, wordint gdin, wordint ni, wordint nj, ftnfloat *maskout, ftnfloat *dlat, ftnfloat *dlon, ftnfloat *yinlat, ftnfloat *yinlon, wordint *yyincount, ftnfloat *yanlat, ftnfloat *yanlon, wordint *yyancount)
{
  wordint ivalue,icode,i,j,k,yin_mgid;
  wordint gdrow_in, gdcol_in, mask_gdrow, mask_gdcol;
  wordint yincount,yancount,yni,ynj;
  ftnfloat *yin_fld, global_extrap_value, local_extrap_value;
  char interp_degree[32],extrap_degree[32],extrap_value[32],local_val[32];
  char global_interp_degree[32],global_extrap_degree[32];
  
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  yin_mgid=Grille[gdrow_in][gdcol_in].mymaskgrid;
  c_gdkey2rowcol(yin_mgid,  &mask_gdrow,  &mask_gdcol);
  yni=Grille[mask_gdrow][mask_gdcol].ni;
  ynj=Grille[mask_gdrow][mask_gdcol].nj;

  yin_fld = (ftnfloat *) malloc(yni*ynj*sizeof(ftnfloat));
  memset(yin_fld,0.0,yni*ynj*sizeof(ftnfloat));
  /*get original options*/
  strcpy(interp_degree,"interp_degree");
  icode = c_ezgetopt(interp_degree,global_interp_degree);
  strcpy(extrap_degree,"extrap_degree");
  strcpy(extrap_value,"extrap_value");
  icode = c_ezgetopt(extrap_degree,global_extrap_degree);
  ivalue = 0;
  if (0 == strcmp(global_extrap_degree,"value"))
    {
    icode = c_ezgetval(extrap_value,&global_extrap_value);
    ivalue = 1;
    }
  strcpy(local_val,"nearest");
  icode = c_ezsetopt(interp_degree, local_val);

  local_extrap_value = 1.0;
  icode = c_ezsetval(extrap_value,local_extrap_value);
  strcpy(local_val,"value");
  icode = c_ezsetopt(extrap_degree, local_val);
  icode = c_ezdefset(gdout,yin_mgid);
  icode = c_ezsint_orig(maskout,yin_fld);
  /*masking is done,reset original interp options*/
  icode = c_ezsetopt(interp_degree, global_interp_degree);
    if (ivalue == 1)
    {
    icode = c_ezsetval(extrap_value, global_extrap_value);
    }
  icode = c_ezsetopt(extrap_degree, global_extrap_degree);
  free(yin_fld);

/* now create the destination grids */
  yancount=0;
  yincount=0;
  for (j=0; j<nj; j++)
      {
      for (i=0;i<ni; i++)
         {
         k=(j*ni)+i; 
         if (maskout[k] == 1.0)
            {
            yanlat[yancount]=dlat[k];
            yanlon[yancount]=dlon[k];
            yancount++;
            }
         else
            {
            yinlat[yincount]=dlat[k];
            yinlon[yincount]=dlon[k];
            yincount++;
            }
         }
      }
  *yyincount = yincount;
  *yyancount = yancount;
  return icode;
}

