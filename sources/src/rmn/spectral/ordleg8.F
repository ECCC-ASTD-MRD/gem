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

***S/P ORDLEG8 - THIS ROUTINE IS A SUBSET OF BELOUSOVS ALGORITHM
*                USED TO CALCULATE ORDINARY LEGENDRE POLYNOMIALS.


      SUBROUTINE ORDLEG8(SX,COA,IR)

*AUTEUR   - D. ROBERTSON

*LANGAGE  - fortran

*OBJET(ORDLEG)
*         - THIS ROUTINE IS A SUBSET OF BELOUSOVS ALGORITHM
*           USED TO CALCULATE ORDINARY LEGENDRE POLYNOMIALS.

*LIBRAIRIES
*         - SOURCE  RMBSOURCELIB,ID=RMNP     DECK=ORDLEG
*         - OBJET   RMNLIB,ID=RMNP

*APPEL    - CALL ORDLEG8(SX,COA,IR)

*ARGUMENTS
*   OUT   - SX  - LEGENDRE POLYNOMIAL EVALUATED AT COA
*   IN    - COA - COSINE OF COLATITUDE
*   IN    - IR  - WAVE NUMBER

*     Version codee en REAL*8 par Bernard Dugas, 7 janvier 1994.

*     -----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
*     ----------------------------------------------------------------- 

*   RDTODG = 180/PIE ; DGTORD = PIE/180
*     DATA PIE    /3.1415926535898/
*     DATA RDTODG /57.295779513082/
*     DATA DGTORD /1.7453292519943E-2/

*-------------------------------------------------------------------------------
      PI     = 4.0*ATAN(1.D0)
      RDTODG = 180./PI
      DGTORD = PI/180.

      SQR2   = SQRT(2.D0)
      IRPP   = IR + 1
      IRPPM  = IRPP - 1
      DELTA  = ACOS(COA)
      SIA    = SIN(DELTA)

      THETA  = DELTA
      C1     = SQR2

      DO 20 N=1,IRPPM
         FN    = N
         FN2   = 2.*FN
         FN2SQ = FN2*FN2
         C1    = C1* SQRT(1.0-1.0/FN2SQ)
   20 CONTINUE

      N   = IRPPM
      ANG = FN*THETA
      S1  = 0.0
      C4  = 1.0
      A   =-1.0
      B   = 0.0
      N1  = N+1

      DO 27 KK=1,N1,2
         K   = KK-1
         IF (K.EQ.N) C4 = 0.5*C4
         S1  = S1+C4* COS(ANG)
         A   = A+2.0
         B   = B+1.0
         FK  = K
         ANG = THETA*(FN-FK-2.0)
         C4  = (A*(FN-B+1.0)/(B*(FN2-A)))*C4
   27 CONTINUE

      SX = S1*C1

      RETURN
      END
