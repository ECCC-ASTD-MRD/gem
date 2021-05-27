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

#ifdef WIN32	/*CHC/NRC*/
#include <string.h> 
#endif

#include <rpnmacros.h>
#define SHRINKFACTOR 1.3
/* 

   Fonction Combsort11 extraite de BYTE, volume 16 numero 4
   (avril 1991), pp. 315-320,  par Richard Box et Stephen Lacey. 
   Modifiee par B. Dugas, RPN - 24 avril 1991.

      CALL IPSORTC( index,list,size )  ...  version caracteres
      CALL  CSORTC( list,index,size )

      CALL IPSORT8( index,list,size )  ...  version nombres reels "double"
      CALL  CSORTD( list,index,size )

      CALL IPSORTI( index,list,size )  ...  version nombres entiers
      CALL  CSORTE( list,index,size )

      CALL IPSORT ( index,list,size )  ...  version nombres reels "single"
      CALL  CSORTR( list,index,size )

   Parametres ...

      list  -   Valeurs a trier (respectivement, de type caractere, 
                entier ou reel (simple ou double), selon que l'on fasse
                appel a csortc, csorte, csortr ou csortd).
      index - (*Dans les routines de type csort*) Valeurs entieres de
                1 a "size". On ne deplace pas les elements de "list", 
                mais plutot leurs indices. Si index[0] est nul ou negatif,
                tous les elements de ce champs sont initialises (juste-
                ment, de 1 a "size" ...).
              (*Dans les routines de type ipsort*) Valeurs de retour de
                1 a "size" denotant l'ordre croissant des valeurs des
                elements de "list".
      size  -   Nombre entier d'elements dans "list" et "index".


   Notes ...

      Les routines CSORTE, CSORTD et CSORTR sont a toutes fins identiques.
      La seule difference reside dans la declaration des types contenus
      dans list. Les routines IPSORTx appellent ces dernieres.

*/
void 
f77name(csortr) ( list, index, size )
wordfloat *list ;
wordint *index, *size ;

{

   int switches = 0 ;
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faut-il initialiser les index's */

   if ( *index <= 0 ) 
      for ( i = 0 ; i < *size ; *( index + i ) = i+1 , ++i );
   
/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( *(list + *(index + i) - 1) > *(list + *(index + j) - 1) )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );
/* 
   Like the bubble sort and shell, we check for
   a clean pass, i.e. switches = 0 and gap = 1.
*/

}

void 
f77name(ipsort) ( index, list, size )
wordint *index, *size ;
wordfloat *list ;

{
   *index = -1 ;
   f77name(csortr) ( list, index, size ) ;
}


void 
f77name(csortd) ( list, index, size )
double *list ;
wordint *index, *size ;

{

   int switches = 0 ; 
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faut-il initialiser les index's */

   if ( *index <= 0 ) 
      for ( i = 0 ; i < *size ; *( index + i ) = i+1 , ++i );
   
/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( *(list + *(index + i) - 1) > *(list + *(index + j) - 1) )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );
/* 
   Like the bubble sort and shell, we check for
   a clean pass, i.e. switches = 0 and gap = 1.
*/

}

void 
f77name(ipsort8) ( index, list, size )
wordint *index, *size ;
double  *list ;

{
   *index = -1 ;
   f77name(csortd) ( list, index, size ) ;
}

void 
f77name(csorte) ( list, index, size )
wordint *list, *index, *size ;

{

   int switches = 0 ; 
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faut-il initialiser les index's */

   if ( *index <= 0 ) 
      for ( i = 0 ; i < *size ; *( index + i ) = i+1 , ++i );
   
/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( *(list + *(index + i) - 1) > *(list + *(index + j) - 1) )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );
/* 
   Like the bubble sort and shell, we check for
   a clean pass, i.e. switches = 0 and gap = 1.
*/

}

void 
f77name(ipsorti) ( index, list, size )
wordint *index, *list, *size ;

{
   *index = -1 ;
   f77name(csorte) ( list, index, size ) ;
}

void 
f77name(csortc) ( list, index, size, len )
char *list ;
wordint *index, *size, len ;

{

   int switches = 0 ;
   wordint gap = *size, i = 0, j = 0, top = 0 ;

/* Faut-il initialiser les index's */

   if ( *index <= 0 ) 
      for ( i = 0 ; i < *size ; *( index + i ) = i+1 , ++i );
   
/* Faire le tri en comparant les elements de "list" */

   do
      {

      gap = (wordint)((wordfloat)gap/SHRINKFACTOR);
      switch (gap)
	 {
	 case 0:  /* The smallest gap is 1 - bubble sort */
	    gap = 1;
	    break;
	 case 9:  /* This is what makes this Combsort11 */
	 case 10:
 	    gap = 11;
	    break;
	 default:
	    break;
	 }

      switches = 0; /* Dirty Pass flag */
      top = *size-gap;
	    
      for ( i=0 ; i < top ; ++i )
	 {
	 j = i+gap;
	 if ( strncmp ( list + (*(index + i) - 1)*len,
                        list + (*(index + j) - 1)*len, len ) > 0 )
	    { /* Swap indices */
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    *(index+j) = *(index+i) ^ *(index+j) ;
	    *(index+i) = *(index+i) ^ *(index+j) ;
	    ++switches;
	    } /* End of Swap */
	 }    /* End of Pass */

      } while ( switches || (gap > 1) );

/* Like the bubble sort and shell, we check for a clean pass */

} 

void 
f77name(ipsortc) ( index, list, size, len )
wordint *index, *size, len ;
char *list ;

{
   *index = -1 ;
   f77name(csortc) ( list, index, size, len ) ;
}
