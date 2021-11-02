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
      SUBROUTINE ALPDR2(DALP,ALP,LALP,LM,SINLAT,EPSI,KHEM)
C 
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****
C    * CALCULATES N-S DERIVATIVES OF EACH ASSOCIATED LEGENDRE POLYNOMIAL
C     * DALP MUST CONTAIN N-S DERIVATIVE OF ALP.
C     * ALP WILL CONTAIN LEGENDRE POLYNOMIALS FOR ONE LATITUDE. 
C     * EPSI CONTAINS PREVIOUSLY CALCULATED CONSTANTS.
C 
      DIMENSION DALP(1),ALP(1),EPSI(1)
      DIMENSION LALP(1) 
C-----------------------------------------------------------------------
      MULT=1
      IF(KHEM.NE.0) MULT=2
      MR = 0
C 
      DO 30 M=1,LM
C 
      LR = LALP(M)*MULT 
      ALPILM=0. 
      DO 20 N=1,LR
      MN = MR+N 
      FNS=FLOAT(M+N-2)
      DALP(MN)=(FNS+FNS+1.)*EPSI(MN)*ALPILM-FNS*SINLAT*ALP(MN)
      ALPILM=ALP(MN)
   20 CONTINUE
C 
      MR = MR+LR
   30 CONTINUE
C 
      RETURN
      END 
