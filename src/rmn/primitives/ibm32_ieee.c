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

#include <stdio.h>
#include <rpnmacros.h>
void c_ibm32_ieee (unsigned long *tab_data_IBM, int nb_data)
{
   int signbit;
   int expos;
   int long Mantis;
   int i;
   
   for(i=0; i<nb_data; i++)
      {
      signbit = tab_data_IBM[i] >> 31; 
      
      expos = (tab_data_IBM[i] >> 24) & 0x7f;
      
      Mantis = (tab_data_IBM[i] & 0xffffff);
      
      expos = ((expos - 64) << 2) + 128 - 2;
      
      if(Mantis != 0)
	 {
	 while((Mantis & 0x800000) == 0)
	    { 
	    Mantis <<= 1;
	    expos--;
	    }
	 Mantis &= 0x7fffff;
	 
	 tab_data_IBM[i] = Mantis | (expos << 23) | (signbit << 31);
	 }
      
      if(expos <= 0)                               /* || (Mantis == 0))*/
	 {
	 tab_data_IBM[i] = 0.0;
         expos = 0;
         }
      else{
          if(expos >= 255)
	     {
	     fprintf(stderr,"c_ibm32_ieee ERROR: Overflow in data field\n");
	     exit(1);
	     }
	  }

      }/* end for */
   
}/* end transfert_IBM_IEEE */

void f77name(ibm32_ieee)(unsigned long *tab_data_IBM, int *f_nb_data)
{
  int nb_data = *f_nb_data;
  c_ibm32_ieee(tab_data_IBM,nb_data);
}
