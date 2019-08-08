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
      SUBROUTINE SPAGG8(GG,ILG,ILATH,COA,KHEM,P,LR,LM,KSM,
     +                  ALP,EPSI,LALP,LA,WRKS,WRKL) 

***   ****   OCT 1975  -  JOHN D. HENDERSON  ****
***    PRODUCES GAUSSIAN GRID FROM SPECTRAL COEFF IN P(LR,LM). 
***      THERE IS LM EAST-WEST WAVES FROM M=0 TO LM-1
***      FOR EACH M THERE IS LR(M) COMPLEX COEFFICIENTS. 
***    IF GLOBAL GG=(ILG,ILATH*2), OTHERWISE GG=(ILG,ILATH). 
***    LONGITUDE IS ZERO AT LEFT OF GG AND INCREASES EASTWARD
***    LATITUDE IS SOUTHERN MOST AT BOTTOM OF GG AND INCREASES NORTHWAR
***    LEFT COLUMN IS COPIED INTO THE RIGHT COLUMN.

***    KHEM  0 = GLOBAL,  1 = N HEM ONLY,  2 = S HEM ONLY. 
***    KSM   +1=SYMMETRIC,  0=GLOBAL,  -1=ANTISYMMETRIC
***    COA(ILATH) CONTAINS THE COSINES OF THE COLATITUDE (N TO S). 
***    ALP IS A WORK FIELD FOR LEGENDRE POLYNOMIALS. 
***    EPSI IS A FIELD OF CONSTANTS THE SAME SIZE AS ALP.
***    WRKS DOIT AVOIR UNE DIMENSION D'AU MOINS (ILG+1)
***    WRKL IS A WORK ARRAY OF LENGTH  LA  COMPLEX WORDS 
***    ON DOIT S'ASSURER QUE WRKL A UNE DIMENSION D'AU MOINS (ILG+2) 
***    FAST FOURIER TRANSFORM REQUIRES THAT ILG BE FACTORIZABLE
***                                      P    Q     R
***    IN THE FORM :          ILG   =  2  * 3   * 5
***    ILG MUST BE > 8 

***    REVISION 1 : C.BEAUDOIN NOV. 85 
***               : UTILISATION DE FFT91A AU LIEU DE FFGFW2
***    Version codee en Real/REAL*8 par Bernard Dugas, 4 janvier 1994.

      INTEGER   ILG,ILATH,KHEM
      INTEGER   LM,KSM,LA,LALP(LM),LR(LM) 
      REAL*8    ALP(LA),EPSI(LA),COA(ILATH)
      REAL      GG(ILG,1),WRKL(*),WRKS(*) 
      COMPLEX   P(LA)
*
      INTEGER   L,ILGH,MAXF,ILP,NH,IH,IHGG
      REAL*8    SINLAT
*
CCCC  EXTERNAL  FFGFW2
      EXTERNAL  ALPNM8,ALPAS8,MOVLEV,FFT91A,SPAF8
*-----------------------------------------------------------------------
***    CALCULATE CONSTANTS.

      ILGH = ILG/2
      MAXF = LM-1 
      ILP  = 1 
      IF (KSM.LT.0) 
     +ILP  = LA+1 

***    THE NORTHERN HEMISPHERE IS DONE FIRST, FOLLOWED BY THE
***    SOUTHERN HEMISPHERE. EITHER CAN BE OMITTED IF NOT NEEDED. 
***    NH=1 FOR NORTHERN HEMISPHERE, NH=2 FOR SOUTHERN HEMISPHERE. 

      DO 100 NH=1,2 

         IF(NH.EQ.1.AND.KHEM.EQ.2)  GOTO 100
         IF(NH.EQ.2.AND.KHEM.EQ.1)  GOTO 100

***       LATITUDE LOOP IN ONE HEMISPHERE.
***       IH  COUNTS FROM THE POLE TO THE EQUATOR. 

         DO 50 IH=1,ILATH

***          IHGG  COUNTS FROM BOTTOM TO TOP OF GAUSSIAN GRID.

                                        IHGG=IH
            IF (NH.EQ.1.AND.KHEM.EQ.1)  IHGG=ILATH+1-IH
            IF (NH.EQ.1.AND.KHEM.EQ.0)  IHGG=ILATH+1-IH+ILATH

            IF (NH.EQ.1)  SINLAT= COA(IH)
            IF (NH.EQ.2)  SINLAT=-COA(IH)

            CALL ALPNM8(ALP,LALP,LM,SINLAT,EPSI,KSM)
            IF (KHEM.NE.0) CALL ALPAS8( ALP, LA, LALP, LR, LM, WRKL )

***          ON REMPLIT WRKL AVEC DES ZEROS EN PREVISION DU "CALL FFT91A"

            DO 10 L=2*LM+1,ILG+2
               WRKL(L) = 0.0 
   10       CONTINUE

            CALL SPAF8( WRKL,P,LR,LM,ALP(ILP),LA )
CCCC        CALL FFGFW2( GG(1,IHGG),ILG,WRKL,ILGH,MAXF,ILG,WRKS,1 )
            CALL FFT91A( WRKL,WRKS,1,ILG,+1 )
            CALL MOVLEV( WRKL,GG(1,IHGG),ILG ) 

   50    CONTINUE

  100 CONTINUE

      RETURN
      END 
