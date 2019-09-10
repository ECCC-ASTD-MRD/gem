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

#if defined(NEC)
#endif
#include <stdlib.h>
#include <rpnmacros.h>
/*
   La routine getenvc fouille l'environnement pour une chaine
   de type "name=xxx", retournant "xxx" ou une chaine vide 
   dans value. ( BD, RPN - 04 mars 1992 )

*/

void
f77name(getenvc) ( name, value, len1, len2 )
F2Cl  len1, len2;
char name[1], value[1];

{

   wordint i;
   unsigned size;
   char *temp, *hold;

/* Transfert de name dans une chaine C */

   size = len1+len2+1 ;
   temp = (char *) malloc( size ) ;

   for ( i=0 ; 
         i < len1 && name[i] != ' ' ; 
         i++ ) 
         *(temp+i) = name[i] ;

   *(temp+i) = '\0' ;

/* Appel a la fonction C getenv */

   hold = (char *) getenv( temp ) ;

/*
   Si la chaine n'a pas ete trouvee, getenv retourne un
   pointeur null. value sera alors entierement vide.
   Sinon, on copie le contenu de hold dans value.
 
*/

   for ( i=0 ; i < len2 ; i++ ) value[i] = ' ' ;

   if ( hold != 0 && size != 1 ) {
        size = strlen( hold ) ;
        for ( i=0 ; i < size ; i++ ) value[i] = *(hold+i) ;
        }

   free( temp );

}
