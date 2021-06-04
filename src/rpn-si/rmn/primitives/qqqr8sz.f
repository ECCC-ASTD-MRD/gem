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
      INTEGER FUNCTION qqqr8sz( )

      IMPLICIT   none

***    AUTEUR : B.Dugas - 14 juillet 1993.

***    DESCRIPTION : Retourne la taille d'un mot REAL*8 en
***                  unites de mots REALs.

      REAL*8     ZD(2)
      REAL       ZS(2)

      qqqr8sz = ( LOC( ZD(2) ) - LOC( ZD(1) ) )
     +        / ( LOC( ZS(2) ) - LOC( ZS(1) ) )

      RETURN
      END
