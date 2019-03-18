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
      SUBROUTINE ALPNM8(ALP,LALP,LM,SINLAT,OLDE,KHEM) 

***    AUTHOR: B.DUGAS, SEPT 1993
***    (THIS IS A RECODED VERSION BASED ON 1) A FIRST VERSION 
***     BY J.HENDERSON AND ON 2) COMMENTS REGARDING BELOUSOV'S
***     FORMULA AS FOUND IN MWR[118], PP. 2248-2251)

***    PUTS LEGENDRE POLYNOMIALS IN ALP FOR ONE LATITUDE.
***    SINLAT IS THE SINE OF THE REQUIRED LATITUDE.
***    EPSI IS NOW CALCULATED ON DEMAND BY ROUTINE EPSIL8. 
***    OLDE IS USED TO PRESERVE THE OLD CALLING SEQUENCE.

***    THE SYMMETRIC AND ANTISYMMETRIC VALUES IN ALP ARE INTERLEAVED 
***     IN EACH ROW.   E.G. ROW 1 = ( 0,1,2,3....LALP ). 

      IMPLICIT  none

      REAL      OLDE
      REAL*8    ALP(*),SINLAT
      INTEGER   LM,LALP(LM),KHEM

      LOGICAL   RHOMBOID
      INTEGER   NM0,NM,NM1,NM2,LRLM
      INTEGER   MULT,M,N,IR,LR,MR,LASTN
      REAL*8    ALPNM,ALPNMM1,ALPNMM2,AC,AD,AE
      REAL*8    COS2, PROD,A,B, EPSI

      POINTER ( IC,AC(0:LM-1,0:1) )
      POINTER ( ID,AD(0:LM-1,0:1) )
      POINTER ( IE,AE(0:LM-1,0:1) )

      POINTER ( IP,EPSI(1) )

      COMMON   /COMLRLM/ LRLM
      COMMON   /ALPNRT/ IC,ID,IE,IP

      EXTERNAL  EPSIL8

*-----------------------------------------------------------------------
      IF ((LALP(1)*1000+LM) .NE. LRLM)
     +    CALL EPSIL8( LALP,LM,KHEM )

                    MULT = 1
      IF(KHEM.NE.0) MULT = 2

      COS2   = 1.-SINLAT**2 

      PROD = 1. 
      A    = 1.
      B    = 0.

      IR    = LALP(1)*MULT-1

***    COMPUTE THE FIRST TWO COLUMNS.

      MR = 0 
      DO M=1,MIN(2,LM)

      NM      = MR+1 

      ALPNMM2 = 0.0
      ALPNMM1 = SQRT(.5*PROD) 

      ALP(NM) = ALPNMM1

***    NOW COMPUTE ELEMENTS 2 TO LR

      LR = LALP(M)*MULT 

      IF (LR.GT.1)                                             THEN

          DO  100 NM=MR+2,MR+LR

              ALPNM   = ( SINLAT*ALPNMM1 - EPSI(NM-1)*ALPNMM2 )
     +                /              EPSI(NM) 

              ALP(NM) = ALPNM

              ALPNMM2 = ALPNMM1
              ALPNMM1 = ALPNM

  100     CONTINUE

          A    = A+2.
          B    = B+2.
          PROD = PROD*COS2*A/B

          MR = MR + LR

      END IF
c     WRITE(6,'(I5,2E25.15)') NM-1,ALP(NM-1),ALP(NM-2)
      END DO

      IF (LM.LT.3) RETURN

***    DETERMINE THE TYPE OT TRUNCATION.

      IF (LALP(1).EQ.LALP(3))                                  THEN
          RHOMBOID = .TRUE.
      ELSE
          RHOMBOID = .FALSE.
      END IF

***    CALCULATE THE REST USING BELOUSOV'S FORMULA.

      NM0 = 1

      DO  M=2,LM-1

          NM2   = NM0
          NM1   = NM0+1
          NM    = NM0+(LALP(M-1)+LALP(M))*MULT
          N     = 0

          ALP(NM) = AC(M,N)*ALP(NM2)-AD(M,N)*ALP(NM1)*SINLAT

          IF (RHOMBOID)                                        THEN
              LASTN = IR-1
          ELSE
              LASTN = LALP(M+1)*MULT-1
          END IF

          DO  200 N=1,LASTN

              NM2     = NM2+1
              NM1     = NM1+1
              NM      = NM +1

              ALP(NM) = AC(M,N)*ALP(NM2 )
     +                - AD(M,N)*ALP(NM1 )*SINLAT
     +                + AE(M,N)*ALP(NM-1)*SINLAT

  200     CONTINUE

          IF (RHOMBOID)                                        THEN

***        CALCULATE THE LAST ROW OF A RHOMBOIDAL TRUNCATION.

              NM      =   NM+1
              ALP(NM) = ( SINLAT*ALP(NM-1) - EPSI(NM-1)*ALP(NM-2) )
     +                /                 EPSI(NM)

          END IF

          A    = A+2.
          B    = B+2.
          PROD = PROD*COS2*A/B

          NM0 = NM0+LALP(M-1)*MULT

      END DO

      RETURN
*-----------------------------------------------------------------------

      END 
