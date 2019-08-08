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
      SUBROUTINE SPAGG0(GG,ILG,ILATH,COA,KHEM,P,LR,LM,KSM,
     1                  ALP,EPSI,LALP,LA,WRKS,WRKL) 
C 
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****
C     * PRODUCES GAUSSIAN GRID FROM SPECTRAL COEFF IN P(LR,LM). 
C     *   THERE IS LM EAST-WEST WAVES FROM M=0 TO LM-1
C     *   FOR EACH M THERE IS LR(M) COMPLEX COEFFICIENTS. 
C     * IF GLOBAL GG=(ILG,ILATH*2), OTHERWISE GG=(ILG,ILATH). 
C     * LONGITUDE IS ZERO AT LEFT OF GG AND INCREASES EASTWARD
C     * LATITUDE IS SOUTHERN MOST AT BOTTOM OF GG AND INCREASES NORTHWAR
C     * LEFT COLUMN IS COPIED INTO THE RIGHT COLUMN.
C 
C     * KHEM  0 = GLOBAL,  1 = N HEM ONLY,  2 = S HEM ONLY. 
C     * KSM   +1=SYMMETRIC,  0=GLOBAL,  -1=ANTISYMMETRIC
C     * COA(ILATH) CONTAINS THE COSINES OF THE COLATITUDE (N TO S). 
C     * ALP IS A WORK FIELD FOR LEGENDRE POLYNOMIALS. 
C     * EPSI IS A FIELD OF CONSTANTS THE SAME SIZE AS ALP.
C     * WRKS DOIT AVOIR UNE DIMENSION D'AU MOINS (ILG+1)
C     * WRKL IS A WORK ARRAY OF LENGTH  LA  COMPLEX WORDS 
C     * ON DOIT S'ASSURER QUE WRKL A UNE DIMENSION D'AU MOINS (ILG+2) 
C     * FAST FOURIER TRANSFORM REQUIRES THAT ILG BE FACTORIZABLE
C     *                                   P    Q     R
C     * IN THE FORM :          ILG   =  2  * 3   * 5
C     * ILG MUST BE > 8 
C 
C     * REVISION 1 : C.BEAUDOIN NOV. 85 
C     *            : UTILISATION DE FFT91A AU LIEU DE FFGFW2
C 
      COMPLEX P(1)
      DIMENSION GG(ILG,1),COA(1)
      DIMENSION ALP(1),EPSI(1),WRKS(1),WRKL(1)
      DIMENSION LR(1),LALP(1) 
C-----------------------------------------------------------------------
C     * CALCULATE CONSTANTS.
C 
      ILGH=ILG/2
      MAXF=LM-1 
      ILP=1 
      IF(KSM.LT.0) ILP=LA+1 
C 
C     * THE NORTHERN HEMISPHERE IS DONE FIRST, FOLLOWED BY THE
C     * SOUTHERN HEMISPHERE. EITHER CAN BE OMITTED IF NOT NEEDED. 
C     * NH=1 FOR NORTHERN HEMISPHERE, NH=2 FOR SOUTHERN HEMISPHERE. 
C 
      DO 100 NH=1,2 
      IF(NH.EQ.1.AND.KHEM.EQ.2)  GO TO 100
      IF(NH.EQ.2.AND.KHEM.EQ.1)  GO TO 100
C 
C     * LATITUDE LOOP IN ONE HEMISPHERE.
C     *  IH  COUNTS FROM THE POLE TO THE EQUATOR. 
C 
      DO 50 IH=1,ILATH
C 
C     *  IHGG  COUNTS FROM BOTTOM TO TOP OF GAUSSIAN GRID.
C 
                                 IHGG=IH
      IF(NH.EQ.1.AND.KHEM.EQ.1)  IHGG=ILATH+1-IH
      IF(NH.EQ.1.AND.KHEM.EQ.0)  IHGG=ILATH+1-IH+ILATH
C 
      IF(NH.EQ.1)  SINLAT= COA(IH)
      IF(NH.EQ.2)  SINLAT=-COA(IH)
C 
      CALL ALPNM2(ALP,LALP,LM,SINLAT,EPSI,KSM)
      IF(KHEM.NE.0) CALL ALPAS2( ALP, LA, LALP, LR, LM, WRKL )
C 
C     ON REMPLIT WRKL AVEC DES ZEROS EN PREVISION DU "CALL FFT91A"
C 
      DO 10 L=2*LM+1,ILG+2
      WRKL(L) = 0.0 
   10 CONTINUE
C 
      CALL SPAF2(WRKL,P,LR,LM,ALP(ILP),LALP)
C     CALL FFGFW2(GG(1,IHGG),ILG,WRKL,ILGH,MAXF,ILG,WRKS,1) 
      CALL FFT91A (WRKL,WRKS,1,ILG,+1)
      CALL MOVLEV (WRKL,GG(1,IHGG),ILG) 
C 
   50 CONTINUE
C 
  100 CONTINUE
C 
      RETURN
      END 
