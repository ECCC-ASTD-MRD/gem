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
        SUBROUTINE UP2LOW(S1,S2)
        CHARACTER * (*) S1,S2
        INTEGER I,L,L2
        L = MIN(LEN(S1),LEN(S2))
        L2 = LEN(S2)
        DO 10 I = 1,L
          IF ((ICHAR(S1(I:I)).GE.65) .AND. (ICHAR(S1(I:I)).LE.90)) THEN
            S2(I:I) = CHAR(ICHAR(S1(I:I)) + 32)
          ELSE
            S2(I:I) = S1(I:I)
          ENDIF
 10     CONTINUE
	DO 15 I = L+1,L2
	  S2(I:I) = ' '
 15     CONTINUE
        RETURN
*
        ENTRY LOW2UP(S1,S2)
        L = MIN(LEN(S1),LEN(S2))
        L2 = LEN(S2)
        DO 20 I = 1,L
          IF ((ICHAR(S1(I:I)).GE.97) .AND. (ICHAR(S1(I:I)).LE.122)) THEN
            S2(I:I) = CHAR(ICHAR(S1(I:I)) - 32)
          ELSE
            S2(I:I) = S1(I:I)
          ENDIF
 20     CONTINUE
        DO 30 I = L+1,L2
          S2(I:I) = ' '
 30     CONTINUE
        RETURN
        END

