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

***S/P MOVLEV EQUIVALENT CRAY DE MOVLEV CDC
*
      SUBROUTINE MOVLEV(A,B,N)
      INTEGER A(N),B(N)
*
*AUTEUR   M.VALIN DRPN 1983
*
*LANGAGE fortran
*
*OBJECT(MOVLEV)
*      COPIER A DANS B
*
*ARGUMENTS
*        A     TABLEAU DE DEPART
*        B     TABLEAU DE DESTINATION
*        N     NOMBRE DE MOTS A COPIER DE A A B
**
      INTEGER LITEND
      DATA LITEND /1/
      INTEGER *2 LITTLE(2)
      EQUIVALENCE (LITTLE(1),LITEND)
      

 01   CONTINUE
      DO 10 I=1,N
10    B(I) = A(I)
      RETURN
*
*** ENTRY MOVE3216
*
      ENTRY MOVE3216(A,B,N)
*
*AUTEUR   M.Lepine DRPN 1999
*
*LANGAGE fortran
*
*OBJECT(MOVLEV)
*      COPIER A DANS B
*
*ARGUMENTS
*        A     TABLEAU DE DEPART
*        B     TABLEAU DE DESTINATION
*        N     NOMBRE DE MOTS A COPIER DE A A B
**
      IF (LITTLE(1) .EQ. 1) THEN
*         PRINT *,'DEBUG THIS COMPUTER IS LITTLE ENDIAN'
      ELSE
*         PRINT *,'DEBUG THIS COMPUTER IS BIG ENDIAN'
         GOTO 01
      ENDIF
      DO 20 I=1,N
         B(I) = IOR(ISHFT(A(I),16), IAND(ISHFT(A(I),-16),65535))
 20   CONTINUE
      RETURN

*** ENTRY MOVE832
*
      ENTRY MOVE832(A,B,N)
*
*AUTEUR   M.Lepine DRPN 1999
*
*LANGAGE fortran
*
*OBJECT(MOVLEV)
*      COPIER A DANS B
*
*ARGUMENTS
*        A     TABLEAU DE DEPART
*        B     TABLEAU DE DESTINATION
*        N     NOMBRE DE MOTS A COPIER DE A A B
**
      IF (LITTLE(1) .EQ. 1) THEN
*         PRINT *,'DEBUG THIS COMPUTER IS LITTLE ENDIAN'
      ELSE
*         PRINT *,'DEBUG THIS COMPUTER IS BIG ENDIAN'
         GOTO 01
      ENDIF
      DO 30 I=1,N
         B(I) = IOR(IOR(IOR(ISHFT(IAND(A(I),255),24),
     %                      IAND(ISHFT(A(I),8),16711680)),
     %                      IAND(ISHFT(A(I),-8),65280)),
     %                      IAND(ISHFT(A(I),-24),255))
 30   CONTINUE
      RETURN
      END

