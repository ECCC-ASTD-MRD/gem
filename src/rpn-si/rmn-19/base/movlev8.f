*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
**s/p movlev8 - deplacer un bloc memoire par tranches de 8 octets
*
           subroutine movlev8(src,dest,n)
           integer n
           real*8 src(*),dest(*)

*auteur:  J. Caveen - septembre 1993
*
*objet(movlev8)  - sous programme servant a deplacer le contenu d'un bloc
*                  memoire par tranches de 8 octets a la fois.
*                  movlev8 verifie avant de faire le transfer si on 
*                  est sur une machie a 64 bits ou a 32 bits.
*
*arguments
*                src   -  tableau contenant les valeurs a transferer
*                dest  -  tableau destination
*                n     -  longueur des tableaux en mots machine
*
***

            real z1(2)
            real*8 z2(2)

            integer difz1,difz2, inn, i

*
*           verifier si on est sur une machine a 64 ou 32 bits
*
            difz1 = loc(z1(2)) - loc(z1(1))
            difz2 = loc(z2(2)) - loc(z2(1))
*
*           si difz2 = difz1, machine a 64 bits
*
            if(difz2 .gt. difz1) then
                   inn = n /2
            else
                   inn = n
            endif

            do 10 i = 1, inn
                  dest(i) = src(i)
 10         continue

            return
            end

