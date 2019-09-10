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
       SUBROUTINE EPSIL2(EPSI,LALP,LM,KHEM) 
C 
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****
C     * CALCULATES EPSILON(N,M) = SQRT((N**2-M**2)/(4*N**2 - 1))
C     * FOR N FROM 0 TO LALP(M)-1, AND M FROM 0 TO LM-1.
C     * EPSILON IS A FIELD OF CONSTANTS USED IN THE SPECTRAL MODELS.
C 
      DIMENSION EPSI(1) 
      DIMENSION LALP(1) 
C-----------------------------------------------------------------------
C 
      MULT = 1
      IF(KHEM.NE.0) MULT=2
      MR = 0
      DO 20 M=1,LM
      MS=M-1
      N1=1
      IF(MS.EQ.0) N1=2
      LR = LALP(M)*MULT 
      DO 30 N=N1,LR 
      MN = MR+N 
      NS=MS+N-1 
      FNUM=FLOAT(NS**2 - MS**2) 
      FDEN=FLOAT(4*NS**2 - 1) 
   30 EPSI(MN) = SQRT(FNUM/FDEN)
      MR = MR+LR
   20 CONTINUE
      EPSI(1) = 0.
C 
      RETURN
      END 
