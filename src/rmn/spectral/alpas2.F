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
      SUBROUTINE ALPAS2( ALP, LA, LRS, LRA, LM, WRKL )
C 
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****
C     *  ALP CAN CONTAIN THE LEGENDRE POLYNOMIALS CALCULATED
C     *  BY SUBROUTINE ALPMN2, OR THEIR N-S DERIVATIVES CALCULATED
C     * BY SUBROUTINE ALPDR2,OR -N*(N+1)*ALP CALCULATED 
C     * BY SUBROUTINE ALPDL2. 
C     * THE SYMMETRIC AND ANTISYMMETRIC VALUES IN EACH ROW OF ALP 
C     * ARE ORIGINALLY INTERLEAVED. E.G. ROW 1 = (0,1,2,3,...LALP). 
C     * THIS SUBROUTINE SEPERATES ALP INTO SYMMETRIC AND ANTI-
C     * SYMMETRIC PARTS.  WRKS IS AN SCM WORK FIELD WITH SAME LENGTH
C     *  AS ALP.
C 
      DIMENSION ALP(1), WRKL(1) 
      DIMENSION LRS(1),LRA(1) 
C-----------------------------------------------------------------------
C 
C     * TRANSFERS ALP TO WRKS 
C 
      LA2=2*LA
      DO 60 MN=1,LA2
   60 WRKL(MN) = ALP(MN)
C 
C     * PUT SYMMETRIC VALUES FROM ALP(1) TO ALP(LAX)
C     *  PUT ANTI-SYMMETRIC VALUES FROM ALP(LAX+1) TO ALP(2*LAX)
C 
      MHS=0 
      MHA=LA
      MF = 0
      DO 70 M=1,LM
C 
      LHS=LRS(M)
      LHA=LRA(M)
      DO 80 N=1,LHS 
      MNS=MHS+N 
      MNA=MHA+N 
      NNS = MF+2*N-1
      NNA = NNS+1 
      ALP(MNS) = WRKL(NNS)
      ALP(MNA) = WRKL(NNA)
   80 CONTINUE
C 
      MF=MF+2*LHS 
      MHS=MHS+LHS 
      MHA=MHA+LHA 
   70 CONTINUE
C 
      RETURN
      END 
