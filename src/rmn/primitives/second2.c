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

#if defined (HP) || defined(Alpha)
#include <rpnmacros.h>


#include <sys/times.h>
#include <unistd.h>

/*

   Fonction retournant la somme des temps usager et systeme d'un
   travail sur une machine HP.   ( BD, RPN - 08 octobre 1993 )

*/

wordfloat
f77name(second) ( )

{

   struct tms buffer;
   clock_t elapsed;
   wordint ticks;
   wordfloat hold;

   ticks = sysconf(_SC_CLK_TCK) ;
   elapsed = times(&buffer) ;
   hold = buffer.tms_utime + buffer.tms_stime ;

   return hold / ticks ;
   
}

#endif
