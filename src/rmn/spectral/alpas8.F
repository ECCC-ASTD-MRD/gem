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
      SUBROUTINE ALPAS8( ALP, LA, LRS, LRA, LM, OLDW )

***   ****   OCT 1975  -  JOHN D. HENDERSON  ****
***    DEC 23/93 - D.Dugas, RPN: Version REAL/REAL*8.

***     ALP CAN CONTAIN THE LEGENDRE POLYNOMIALS (ALPNM82)
***                  OR THEIR N-S DERIVATIVES    (ALPDR8)
***                  OR -N*(N+1)*ALP             (ALPDL8)

***    THE SYMMETRIC AND ANTISYMMETRIC VALUES IN EACH ROW OF ALP 
***    ARE ORIGINALLY INTERLEAVED. E.G. ROW 1 = (0,1,2,3,...LALP). 
***    THIS SUBROUTINE SEPERATES ALP INTO SYMMETRIC AND ANTI-SYMMETRIC
***    PARTS.  WRKL IS NOW ALLOCATED ON DEMAND. THE FIELD OLDW
***    IS ONLY MEANT TO PRESERVE THE CALL SEQUENCE.

      REAL      OLDW
      REAL*8    ALP(*)
      INTEGER   LM,LRS(LM),LRA(LM),LA

      REAL*8    WRKL
      POINTER ( IW,  WRKL(1) )

      INTEGER   LA2,MN,MHS,MHA,MF,M,N,IER
      INTEGER   LHS,LHA,MNS,MNA,NNS,NNA

      EXTERNAL  HPALLOC,HPDEALLC
*-----------------------------------------------------------------------

      CALL HPALLOC( IW,LA*2,IER,8 )

***    TRANSFERS ALP TO WRKL

      LA2=2*LA
      DO 60 MN=1,LA2
         WRKL(MN) = ALP(MN)
   60 CONTINUE

***    PUT SYMMETRIC VALUES FROM ALP(1) TO ALP(LAX)
***    PUT ANTI-SYMMETRIC VALUES FROM ALP(LAX+1) TO ALP(2*LAX)

      MHS = 0 
      MHA = LA
      MF  = 0
      DO 70 M=1,LM

         LHS = LRS(M)
         LHA = LRA(M)
         DO 80 N=1,LHS 

            MNS = MHS+N 
            MNA = MHA+N 
            NNS = MF+2*N-1
            NNA = NNS+1 

            ALP(MNS) = WRKL(NNS)
            ALP(MNA) = WRKL(NNA)

   80    CONTINUE

         MF  = MF+2*LHS 
         MHS = MHS+LHS 
         MHA = MHA+LHA 
   70 CONTINUE

      CALL HPDEALLC( IW,IER,1 )

      RETURN
      END 
