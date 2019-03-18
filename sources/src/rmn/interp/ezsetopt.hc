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
#include <ctype.h>

extern int f77name(longueur)(char *string, int stringlength);


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezsetopt)(char *option, char *value, F2Cl llenoption, F2Cl llenvalue)
{
   wordint i, icode;
   wordint longueur_option, longueur_value;
   F2Cl lenoption=llenoption, lenvalue=llenvalue;

   char local_opt[32], local_val[32];
   longueur_option = f77name(longueur)(option, lenoption);
   longueur_value = f77name(longueur)(value, lenvalue);

   longueur_option = longueur_option < 32 ? longueur_option : 31;
   longueur_value  = longueur_value  < 32 ? longueur_value  : 31;

   for (i=0; i < longueur_option; i++)
     {
     local_opt[i] = option[i];
     }
   for (i=0; i < longueur_value; i++)
     {
     local_val[i] = value[i];
     }

   local_opt[longueur_option] = '\0';
   local_val[longueur_value] = '\0';
   icode = c_ezsetopt(local_opt, local_val);

   return icode;
}

wordint c_ezsetopt(char *option, char *value)
   {
   char local_opt[32], local_val[32];
   wordint i, option_ok, value_ok;


   option_ok = 0;
   value_ok = 0;

   memset(local_opt, (int) '\0', 32);
   memset(local_val, (int) '\0', 32);

   strcpy(local_opt, option);
   strcpy(local_val, value);

   for (i=0; i < strlen(local_opt); i++)
     local_opt[i] = (char) tolower((int)local_opt[i]);

   for (i=0; i < strlen(local_val); i++)
     local_val[i] = (char) tolower((int)local_val[i]);

   if (0 == strcmp(local_opt, "correction_polaire")) strcpy(local_opt, "polar_correction");
   if (0 == strcmp(local_opt, "degre_interp")) strcpy(local_opt, "interp_degree");
   if (0 == strcmp(local_opt, "degre_extrap")) strcpy(local_opt, "extrap_degree");
   if (0 == strcmp(local_opt, "use_1sousgrille")) strcpy(local_opt, "use_1subgrid");

   if (0 == strcmp(local_val, "oui")) strcpy(local_val, "yes");
   if (0 == strcmp(local_val, "ouiouioui")) strcpy(local_val, "yesyesyes");
   if (0 == strcmp(local_val, "non")) strcpy(local_val, "no");
   if (0 == strcmp(local_val, "voisin")) strcpy(local_val, "nearest");
   if (0 == strcmp(local_val, "lineair")) strcpy(local_val, "linear");
   if (0 == strcmp(local_val, "lineaire")) strcpy(local_val, "linear");
   if (0 == strcmp(local_val, "cubique")) strcpy(local_val, "cubic");
   if (0 == strcmp(local_val, "neutre")) strcpy(local_val, "neutral");
   if (0 == strcmp(local_val, "valeur")) strcpy(local_val, "value");

   if (0 == strcmp(local_opt, "use_1subgrid"))
      {
      option_ok = 1;
      value_ok = 1;
      if (0 == strcmp(local_val, "yes"))
         {
         groptions.use_1subgrid= 1;
         }
      else if (0 == strcmp(local_val, "no"))
         {
         groptions.use_1subgrid= 0;
         }
      else
	      {
	      value_ok = 0;
	      }
      }

 
   if (0 == strcmp(local_opt, "verbose"))
      {
      option_ok = 1;
      value_ok = 1;
      if (0 == strcmp(local_val, "yes"))
         {
         groptions.verbose = 1;
         }
      else if (0 == strcmp(local_val, "yesyesyes"))
         {
         groptions.verbose = 2;
         }
      else if (0 == strcmp(local_val, "no"))
	      {
	      groptions.verbose = 0;
	      }
            else
	      {
	      value_ok = 0;
	      }
      }

   if (0 == strcmp(local_opt, "polar_correction"))
      {
      option_ok = 1;
      value_ok = 1;
      if (0 == strcmp(local_val, "yes"))
         {
         groptions.degre_interp = 1;
         }
      else if (0 == strcmp(local_val, "no"))
         {
         groptions.degre_interp = 0;
         }
      else
	      {
	      value_ok = 0;
	      }
      }

   if (0 == strcmp(local_opt, "interp_degree"))
      {
      option_ok = 1;
      value_ok = 1;
      if (0 == strcmp(local_val, "nearest"))
         {
         groptions.degre_interp = 0;
         }
      else if (0 == strcmp(local_val, "linear"))
         {
         groptions.degre_interp = 1;
         }
      else if (0 == strcmp(local_val, "cubic"))
         {
         groptions.degre_interp = 3;
         }
      else if (0 == strcmp(local_val, "average"))
         {
         groptions.degre_interp = 4;
         }
      else if (0 == strcmp(local_val, "sph_average"))
         {
         groptions.degre_interp = 5;
         }
      else
         {
         value_ok = 0;
         }
      }

   if (0 == strcmp(local_opt, "extrap_degree"))
      {
      option_ok = 1;
      value_ok = 1;
      if (0 == strcmp(local_val, "neutral"))
         {
         groptions.degre_extrap = groptions.degre_interp;
         }
      else if (0 == strcmp(local_val, "nearest"))
         {
         groptions.degre_extrap = 0;
         }
      else if (0 == strcmp(local_val, "linear"))
         {
         groptions.degre_extrap = 1;
         }
      else if (0 == strcmp(local_val, "cubic"))
         {
         groptions.degre_extrap = 3;
         }
      else if (0 == strcmp(local_val, "maximum"))
         {
         groptions.degre_extrap = MAXIMUM;
         if (groptions.verbose == 1)
            {
            fprintf(stderr, "Extrapolation set to maximum value\n");
            }
         }
      else if (0 == strcmp(local_val, "minimum"))
         {
         groptions.degre_extrap = MINIMUM;
         if (groptions.verbose == 1)
            {
            fprintf(stderr, "Extrapolation set to minimum value\n");
            }
         }
      else if (0 == strcmp(local_val, "value"))
         {
         groptions.degre_extrap = VALEUR;
         if (groptions.verbose == 1)
            {
            fprintf(stderr, "Extrapolation set to value: %f\n", groptions.valeur_extrap);
            }
         }
      else if (0 == strcmp(local_val, "abort"))
         {
         groptions.degre_extrap = ABORT;
         if (groptions.verbose == 1)
            {
            fprintf(stderr, "Extrapolation set to ABORT\n");
            }
         }
      else
	      {
	      value_ok = 0;
	      }
      }

   if (0 == strcmp(local_opt, "missing_interp_alg"))
      {
      }
   else
      {
      }

   if (0 == strcmp(local_opt, "cloud_interp_alg"))
      {
      value_ok = 0;
      option_ok = 0;
      if (0 == strcmp(local_val, "linear"))
         {
         value_ok = 1;
         option_ok = 1;
         groptions.cld_interp_alg = LINEAIRE;
         }
      if (0 == strcmp(local_val, "distance"))
         {
         value_ok = 1;
         option_ok = 1;
         groptions.cld_interp_alg = DISTANCE;
         }
      }
   else
      {
      }

   if ((value_ok + option_ok) != 2)
      {
      if (option_ok == 0)
         {
         fprintf(stderr, "ezsetopt : option not recognized : %s\n", option);
         }
      if (value_ok == 0)
         {
         fprintf(stderr, "ezsetopt : value not recognized : %s\n", value);
         }
      return -1;
      }
   else
      {
     return 0;
      }

   }

