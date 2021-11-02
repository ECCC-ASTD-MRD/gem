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
      SUBROUTINE ALPDL8 (DELALP,ALP,LALP,LM,KHEM)
C 
C     ***** DEC 1976 - CLAUDE GIRARD **** 
*     * DEC 23/93 - D.Dugas, RPN: Version REAL/REAL*8.
*
C     * CALCUL DE L,OPERATEUR  DELALP = -N*(N+1)*ALP
C     * 
C     * N.B. - VOIR CALCUL DE ALP DANS ALPNM2 
C 
      REAL*8    DELALP(1),ALP(1)
      INTEGER   LM,LALP(LM),KHEM
*
      INTEGER   MULT,MR,M,N,JR,MN,NS
      REAL*8    FNS1
C     ------------------------------------------------------------------
                      MULT = 1
      IF (KHEM.NE.0)  MULT = 2 
C
      MR = 0
C 
      DO 50 M=1,LM
         JR = LALP(M)*MULT 
         DO 40 N=1,JR
            MN   = MR+N 
C 
            NS   = M+N-2
            FNS1 = NS*(NS+1)
C 
            DELALP(MN) = -FNS1*ALP(MN)
C 
   40    CONTINUE
         MR = MR+JR
   50 CONTINUE
C 
      RETURN
      END 
