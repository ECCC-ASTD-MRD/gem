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
***S/P CXGAIG - PASSE DES PARAMETRES (REELS) DESCRIPTEURS DE GRILLE
*              AUX PARAMETRES ENTIERS.
      
      SUBROUTINE CXGAIG(CGTYP,IG1,IG2,IG3,IG4,XG1,XG2,XG3,XG4)
      CHARACTER * 1 CGTYP
*     
*AUTEUR- M. VALIN  -  FEV 82
*     
*REVISION001  C. THIBEAULT  -  MARS 83  CONVERSION AU CODE CRAY
*     002  M. Lepine     -  fev  94  bug fix nint pour grille L
*     003  M. Lepine     -  fev  94  introduction du type de grille E
*     004  M. Lepine     -  nov  94  traduction de ratfor a fortran
*     005  M. Valin      -  fev 2013 permutation de 2 enonces (ig2<0)
*     006  M. Valin      -  mars 2018 grille +
*     
*LANGAGE- RATFOR
*     
*OBJET(XGAIG)
*     - PASSE DES PARAMETRES (REELS) DESCRIPTEURS DE GRILLE
*     AUX PARAMETRES ENTIERS.
*     
*LIBRAIRIES
*     - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=XGAIG
*     - OBJET   RMNLIB,ID=RMNP
*     
*APPEL- CALL XGAIG(CGTYP,IG1,IG2,IG3,IG4,XG1,XG2,XG3,XG4)
*     
*ARGUMENTS
*     IN    - CGTYP - TYPE DE GRILLE (VOIR OUVRIR)
*     OUT   - IG1   - DESCRIPTEUR DE GRILLE (ENTIER) VOIR OUVRIR
*     OUT   - IG2   - DESCRIPTEUR DE GRILLE (ENTIER) VOIR OUVRIR
*     OUT   - IG3   - DESCRIPTEUR DE GRILLE (ENTIER) VOIR OUVRIR
*     OUT   - IG4   - DESCRIPTEUR DE GRILLE (ENTIER) VOIR OUVRIR
*     IN    - XG1   - ** DESCRIPTEUR DE GRILLE (REEL),
*     IN    - XG2   -    CGTYP = 'N', PI, PJ, D60, DGRW
*     IN    - XG3   -    CGTYP = 'L', LAT0, LON0, DLAT, DLON,
*     IN    - XG4   -    CGTYP = 'A', 'B', 'G', XG1 = 0, GLOBAL
*     = 1, NORD
*     = 2, SUD **
*     CGTYP = 'E', LAT1, LON1, LAT2, LON2
*     CGTYP = '+', LAT, LON, dummy, dummy
*     
*MESSAGES- "ERREUR DANS LA DESCRIPTION DE LA GRILLE (IG1) (XGAIG)"
*     "ERREUR, MAUVAISE SPECIFICATION (LAT0) (XGAIG)"
*     "ERREUR, GRILLE INCONNUE (TYPE) (XGAIG)"
*     
*------------------------------------------------------------------
*     
*     
      REAL XXG2,XXG4,XLON
      REAL*8 :: XLON8, XLAT8
      INTEGER I2B
      LOGICAL, EXTERNAL :: VALIDE
      LOGICAL STATUS

      IF (CGTYP .EQ. 'N' .OR. CGTYP.EQ.'S') THEN
         IG1 = NINT(XG2 * 10.)
         IG2 = NINT(XG1 * 10.)
         IG3 = NINT(XG4 * 100.)
         IG4 = NINT(XG3 * 0.01)
 100     CONTINUE
         IF (IG3 .LT. 0) THEN
            IG3 = IG3 + 36000
            GOTO 100
         ENDIF
         IF(IG1.LT.0.OR.IG2.LT.0.OR.IG1.GT.2047.OR.
     %        IG2.GT.2047.OR.IG4.GT.32000) THEN
            IG1 = 0
            IG2 = 0
            IG3 = 0
            IG4 = 32768
            IF(XG3 .GT. 204700) THEN ! CA NE PASSE PAS EN HECTOMETRES
               IG3 = 32768      ! ON LE MET EN KILOMETRES
               IG1 = NINT(XG3*.001)
            ELSE
               IG3 = 0
               IG1 = NINT(XG3*.01)
            ENDIF
            IG2 = NINT(XG4*10)  ! DGRW EN DECIDEGRES
            IF(IG2.LT.0) THEN
               IG2 = ABS(IG2)
               IG4 = IG4 + 16384
            ENDIF
            IF(IG2.GT.1800) THEN
               IG2 = ABS(IG2 - 3600)
               IG4 = IG4 + 16384
            ENDIF
            IHEM = 1
            IF('S'.EQ.CGTYP) IHEM = 2
            CALL LLFXY(DLAT,DLON,1.-XG1,1.-XG2,XG3,XG4,IHEM)
            DLAT = 90. - DLAT   ! ANGLE PAR RAPPORT AU POLE NORD
            IF(DLON.LT.0) DLON = DLON + 360.
            IG3 = IG3 + NINT(DLON*32767./360.)
            IG4 = IG4 + NINT(DLAT*16383./180.)
         ENDIF
         
      ELSE IF (CGTYP .EQ. 'A' .OR. CGTYP .EQ. 'B' .OR.
     %        CGTYP .EQ. 'G')  THEN
         IG1 = XG1
         IG2 = XG2
         IG3 = 0
         IG4 = 0
         STATUS = VALIDE("IG1",IG1,0,2) ! VERIFIER SI IG1=0,1,OU 2
         STATUS = VALIDE("IG2",IG2,0,1) ! VERIFIER SI IG2=0 OU 1
         
      ELSE IF(CGTYP .EQ. 'C') THEN   ! C TYPE LAT LON GRID
        IG1 = NINT(180. / XG3)
        IG2 = NINT(360. / XG4)
        IG3 = NINT((90. + XG1) * 100.)
        IG4 = NINT(XG2 * 100.)
 200    CONTINUE
        IF (IG4 .LT. 0) THEN
           IG4 = IG4 + 36000
           GOTO 200
        ENDIF
        IF (IG3 .LT. 0) WRITE(6,601)

      ELSE IF (CGTYP .EQ. 'H')  THEN   ! LAMBERT CONFORME CENTREE
        IG1 = NINT(5.*XG4)             ! LAMBDA0 EN 5EME DE DEGRE
 300    CONTINUE
        IF (IG1 .LT. 0) THEN
           IG1 = IG1 + 1800
           GOTO 300
        ENDIF
        IG2 = NINT(.002*XG3)           ! DELTA S EN DEMI KILOMETRES
        IG3 = NINT(XG1)                ! 200*PHI1+PHI2 (DEMI DEGRES)
        IG4 = NINT(100.*(90.+XG2))     ! PHI0+90  EN CENTIDEGRES

      ELSE IF (CGTYP .EQ. 'L') THEN
        IG1 = NINT(XG3 * 100.)
        IG2 = NINT(XG4 * 100.)
        IG3 = NINT((90. + XG1) * 100.)
        IG4 = NINT(XG2 * 100.)
 400    CONTINUE
        IF (IG4 .LT. 0) THEN
           IG4 = IG4 + 36000
           GOTO 400
        ENDIF
        IF (IG3 .LT. 0) WRITE(6,601)

      ELSE IF (CGTYP .EQ. 'E')  THEN          !  GRILLE LAT,LON (GEF)
        STATUS = VALIDE("XG1",NINT(XG1),-90,90)
        STATUS = VALIDE("XG3",NINT(XG3),-90,90)
        XXG2 = XG2
        XXG4 = XG4
 500    CONTINUE
        IF (XXG2 .LT. 0) THEN
           XXG2 = XXG2 + 360.
           GOTO 500
        ENDIF
 600    CONTINUE
        IF (XXG4 .LT. 0) THEN
           XXG4 = XXG4 + 360.
           GOTO 600
        ENDIF
        IG1 = NINT((XG1+90.) * 40.)
        IG2 = NINT(XG3 * 40.)
        IG3 = NINT((XXG2+90.) * 40.)
C
C       bug de code, le +90 est de trop, ce qui peut causer un debordement
C       pour ig3
C
        if(ig3 .ge. 16384) ig3=ig3-16384
        IG4 = NINT(XXG4 * 40.)
        I2B = IAND(IG1,3)
        IG1 = ISHFT(IG1,-2)
        IG3 = IOR(ISHFT(IG3,2),I2B)
        if (ig2.lt.0) ig2 = ig2 + 7201
        I2B = IAND(IG2,3)
        IG2 = ISHFT(IG2,-2)
        IG4 = IOR(ISHFT(IG4,2),I2B)

      ELSE IF (CGTYP .EQ. '+')  THEN            !  point LAT,LON
        XLAT8 = XG1
        STATUS = VALIDE("XG1",NINT(XLAT8),-90,90)   ! -90, +90
        XLON8 = XG2
        if(XLON8 < 0) XLON8 = XLON8 + 360.0     ! -180, +180 -> 0, 360
        STATUS = VALIDE("XG2",NINT(XLON8),0,360)    ! 0, 360
        IG3  = nint( (XLAT8+100.)*100. )        ! compatibilite arriere, centidegres (10 -> 19000)
        IG4  = nint( XLON8*100. )               ! compatibilite arriere, centidegres (0 -> 36000)
        IG1  = nint( (XLAT8+100.)*100000. ) - IG3*1000  ! en 1/100000 de degre
        IG1  = IG1 + 1000                       ! correction, IG1 pourrait etre < 0  (500 -> 1500)
        IG2  = nint( XLON8*100000. ) - IG4*1000         ! en 1/100000 de degre
        IG2  = IG2 + 1000                       ! correction, IG2 pourrait etre < 0  (500 -> 1500)

      ELSE
        WRITE(6,602)
      ENDIF


  601 FORMAT(1H0,' ERREUR, MAUVAISE SPECIFICATION (LAT0) (XGAIG)')
  602 FORMAT(1H0,' ERREUR, GRILLE INCONNUE (TYPE) (XGAIG)')
      RETURN
      END

