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
      SUBROUTINE GWAQD8(Q,D,LRQ,LRD,LM,UGG,VGG,ILG,ILATH,COA,SIA, 
     +                  WOCS,KHEM,ALP,DALP,EPSI,LALP,LA,WRKS,WRKL,A)

***   ****   OCT 1975  -  JOHN D. HENDERSON  ****
***    CONVERTS GAUSSIAN GRIDS OF U,V IN UGG,VGG(ILG,NLAT).
***    TO SPECTRAL COEFFICIENTS OF VORTICITY Q(LRQ,LM) 
***     AND DIVERGENCE D(LRD,LM).
***    IF GLOBAL, NLAT=ILAT.  IF HEMISPHERIC, NLAT=ILATH.

***    U ET V A L ENTREE SONT EN NOEUDS. 

***    KHEM  0 = GLOBAL,  1 = N HEM ONLY,  2 = S HEM ONLY. 
***    SIA(ILATH) CONTAINS THE SINES OF THE COLATITUDE (N TO S). 
***    COA(ILATH) CONTAINS THE COSINES OF THE COLATITUDE (N TO S). 
***    WOCS(ILATH) CONTAINS (GAUSSIAN WEIGHTS)/(SIN(COLAT)**2) (N TO S).
***    ALP IS A WORK FIELD FOR LEGENDRE POLYNOMIALS. 
***    DALP IS FOR THE N-S DERIVATIVES OF ALP. 
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

***    DECLARATION DES PARAMETRES

      INTEGER  LM,LA,ILG,ILATH,KHEM
      INTEGER  LRQ(LM),LRD(LM),LALP(LM)
      COMPLEX  Q(LA),D(LA),WRKL(*) 
      REAL     UGG(ILG,1),VGG(ILG,1),WRKS(1),A
      REAL*8   ALP(LA),DALP(LA),EPSI(LA)
      REAL*8   COA(ILATH),WOCS(ILATH),SIA(ILATH)

***    DECLARATION DES VARIABLES LOCALES

      REAL*8   WEIGHT,SINLAT
      INTEGER  ILGH,IHGG,ILP,NH,IH,I
      REAL     CON,CON1

***    DECLARATION DES REFERENCES EXTERNES

CCCC  EXTERNAL FFWFG2
      EXTERNAL SCOF2,ALPNM8,ALPDR8,ALPAS8
      EXTERNAL MOVLEV,FFT91A,FASP8,EWFDR2
*-----------------------------------------------------------------------
***    INITIAL CONSTANTS. SET Q,D TO ZERO. 

      CON1 = 0.514791/A
      ILGH = ILG/2
      ILP  = 1 
      IF (KHEM.NE.0) 
     +ILP  = LA+1

      CALL SCOF2( Q,LRQ,LM,0 )
      CALL SCOF2( D,LRD,LM,0 )

***    THE NORTHERN HEMISPHERE IS DONE FIRST, FOLLOWED BY THE
***    SOUTHERN HEMISPHERE. EITHER CAN BE OMITTED IF NOT NEEDED. 
***    NH=1 FOR NORTHERN HEMISPHERE, NH=2 FOR SOUTHERN HEMISPHERE. 

      DO 100 NH=1,2 

         IF (NH.EQ.1.AND.KHEM.EQ.2)  GOTO 100
         IF (NH.EQ.2.AND.KHEM.EQ.1)  GOTO 100

***       LATITUDE LOOP IN ONE HEMISPHERE.
***       IH COUNTS FROM THE POLE TO THE EQUATOR. 

         DO 50 IH=1,ILATH

***          IHGG COUNTS FROM BOTTOM TO TOP OF GAUSSIAN GRID.

                                        IHGG = IH
            IF (NH.EQ.1.AND.KHEM.EQ.1)  IHGG = ILATH+1-IH
            IF (NH.EQ.1.AND.KHEM.EQ.0)  IHGG = ILATH+1-IH+ILATH

            IF (NH.EQ.1)  SINLAT = COA(IH)
            IF (NH.EQ.2)  SINLAT =-COA(IH)

            WEIGHT = WOCS(IH) 
            CON    = SIA(IH)*CON1

***          COMPUTE ALP,DALP AND REORDER ROWS IF GRIDS ARE HEMISPHERIC. 

            CALL ALPNM8( ALP,LALP,LM,SINLAT,EPSI,KHEM ) 
            CALL ALPDR8( DALP,ALP,LALP,LM,SINLAT,EPSI,KHEM )
            IF (KHEM.NE.0)  CALL ALPAS8(  ALP,LA,LRD,LRQ,LM,WRKL ) 
            IF (KHEM.NE.0)  CALL ALPAS8( DALP,LA,LRD,LRQ,LM,WRKL ) 

***          TERMS DEPENDING ON U. 

            DO 10 I=1,ILG 
               UGG(I,IHGG) = UGG(I,IHGG)*CON 
   10       CONTINUE

CCCC        CALL FFWFG2( WRKL,ILGH,UGG(1,IHGG),ILG,LM-1,ILG,WRKS,1 )
            CALL MOVLEV( UGG(1,IHGG),WRKL,ILG )
            CALL FFT91A( WRKL,WRKS,1,ILG,-1 )
            CALL FASP8( Q,LRQ,LM,WRKL,DALP(ILP),LALP, WEIGHT )
            CALL EWFDR2( WRKL,WRKL,WRKL,ILGH,1,LM,1., 0 )
            CALL FASP8( D,LRD,LM,WRKL,ALP,LALP, WEIGHT )

***          TERMS DEPENDING ON V. 

            DO 20 I=1,ILG 
               VGG(I,IHGG) = VGG(I,IHGG)*CON 
   20       CONTINUE

CCCC        CALL FFWFG2( WRKL,ILGH,VGG(1,IHGG),ILG,LM-1,ILG,WRKS,1 )
            CALL MOVLEV( VGG(1,IHGG),WRKL,ILG )
            CALL FFT91A( WRKL,WRKS,1,ILG,-1 )
            CALL FASP8( D,LRD,LM,WRKL,DALP,LALP,-WEIGHT )
            CALL EWFDR2( WRKL,WRKL     ,WRKL,ILGH,1,LM,1., 0 )
            CALL FASP8( Q,LRQ,LM,WRKL,ALP(ILP),LALP, WEIGHT ) 

   50    CONTINUE

  100 CONTINUE

***    IF GRIDS ARE HEMISPHERIC, DOUBLE ALL COEFF IN Q AND D.

      IF (KHEM.NE.0) CALL SCOF2( Q,LRQ,LM,2 )
      IF (KHEM.NE.0) CALL SCOF2( D,LRD,LM,2 )

      RETURN
      END 
