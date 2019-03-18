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
      SUBROUTINE SCOF2(SC,LR,LM,KIND) 
C 
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****
C     * ZERO OR DOUBLE THE COMPLEX COEFFICIENTS IN SC.
C     * DEPENDING ON THE VALUE OF KIND. 
C 
      COMPLEX SC(1) 
      DIMENSION LR(1) 
C-----------------------------------------------------------------------
C 
C     * IF KIND=0 SET ALL OF SC TO (0.,0.). 
C 
      IF(KIND.NE.0) GO TO 30
      MR = 0
      LMA = LM
      IF(LR(LM).EQ.0) LMA = LM-1
      DO 20 M=1,LMA 
      JR = LR(M)
      DO 25 N=1,JR
      MN = MR + N 
   25 SC(MN) = (0.,0.)
   20 MR = MR + JR
C 
C     * IF KIND=2 DOUBLE ALL OF SC. 
C 
   30 IF(KIND.NE.2) GO TO 99
      LMA = LM
      IF(LR(LM).EQ.0) LMA = LM-1
      MR = 0
      DO 40 M=1,LMA 
      JR = LR(M)
      DO 45 N=1,JR
      MN = MR + N 
   45 SC(MN) = SC(MN) + SC(MN)
   40 MR = MR + JR
C 
   99 RETURN
      END 
