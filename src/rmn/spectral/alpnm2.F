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
      SUBROUTINE ALPNM2(ALP,LALP,LM,SINLAT,EPSI,KHEM) 
C 
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****
C     * PUTS LEGENDRE POLYNOMIALS IN ALP FOR ONE LATITUDE.
C     * SINLAT IS THE SINE OF THE REQUIRED LATITUDE.
C     * EPSI IS A FIELD OF CONSTANTS THE SAME SIZE AS ALP.
C 
C     * THE SYMMETRIC AND ANTISYMMETRIC VALUES IN ALP ARE INTERLEAVED 
C     *  IN EACH ROW.   E.G. ROW 1 = ( 0,1,2,3....LALP ). 
C 
      DIMENSION ALP(1),EPSI(1)
      DIMENSION LALP(1) 
C-----------------------------------------------------------------------
      MULT = 1
      IF(KHEM.NE.0) MULT = 2
      COS2=1.-SINLAT**2 
      PROD=1. 
      A=1.
      B=0.
C 
C     * LOOP 30 COVERS LONGITUDINAL WAVE NUMBERS 0 TO LM-1. 
C 
      MR = 0
      DO 30 M=1,LM
C 
C     * COMPUTE THE FIRST TWO ELEMEMTS OF THE ROW.
C 
      MN=MR+1 
      ALP(MN) = SQRT(.5*PROD) 
C 
C     * NOW COMPUTE ELEMENTS 2 TO LR
C 
      LR = LALP(M)*MULT 
      IF(LR.EQ.1)  GO TO 30 
      ALPILM=0. 
      DO 20 N=2,LR
      MN = MR+N 
      ALP(MN)=(SINLAT*ALP(MN-1)-EPSI(MN-1)*ALPILM)/EPSI(MN) 
      ALPILM=ALP(MN-1)
   20 CONTINUE
C 
      A=A+2.
      B=B+2.
      PROD=PROD*COS2*A/B
C 
      MR = MR + LR
   30 CONTINUE
C 
      RETURN
      END 
