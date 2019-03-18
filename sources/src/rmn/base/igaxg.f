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
***S/P CIGAXG - PASSE DES PARAMETRES (ENTIERS) DESCRIPTEURS DE GRILLE
*              AUX PARAMETRES REELS.


      SUBROUTINE CIGAXG(CGTYP,XG1,XG2,XG3,XG4,IG1,IG2,IG3,IG4)
      CHARACTER * 1 CGTYP
*
*AUTEUR   - M. VALIN  -  FEV 82
*
*REVISION 001  C. THIBEAULT  -  MARS 83  CONVERSION AU CODE CRAY
*         002  M. LEPINE  -  AVRIL 90 GRTYP DE TYPE CARACTERE POUR
*                            COMPATIBILITE AVEC LES FICHIERS STD 89
*         003  M. LEPINE  -  FEVRIER 94 AJOUT DU GRTYP E
*         004  M. LEPINE  -  NOVEMBRE 94 TRADUCTION DE RATFOR A FORTRAN
*         005  M. Valin   -  Fev 2013 Bug fix, enlever le then de trop
*
*LANGAGE  - FORTRAN
*
*OBJET(CIGAXG)
*         - PASSE DES PARAMETRES (ENTIERS) DESCRIPTEURS DE GRILLE
*           AUX PARAMETRES REELS.
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=IGAXG
*         - OBJET   RMNLIB,ID=RMNP
*
*APPEL    - CALL IGAXG(CGTYP,XG1,XG2,XG3,XG4,IG1,IG2,IG3,IG4)
*
*ARGUMENTS
*   IN    - CGTYP - TYPE DE GRILLE (VOIR OUVRIR)
*   OUT   - XG1   - ** DESCRIPTEUR DE GRILLE (REEL),
*   OUT   - XG2   -    IGTYP = 'N', PI, PJ, D60, DGRW
*   OUT   - XG3   -    IGTYP = 'L', LAT0, LON0, DLAT, DLON,
*   OUT   - XG4   -    IGTYP = 'A', 'B', 'G', XG1 = 0. GLOBAL,
*                                                 = 1. NORD
*                                                 = 2. SUD **
*   IN    - IG1   - DESCRIPTEUR DE GRILLE (ENTIER) VOIR OUVRIR
*   IN    - IG2   - DESCRIPTEUR DE GRILLE (ENTIER) VOIR OUVRIR
*   IN    - IG3   - DESCRIPTEUR DE GRILLE (ENTIER) VOIR OUVRIR
*   IN    - IG4   - DESCRIPTEUR DE GRILLE (ENTIER) VOIR OUVRIR
*
*MESSAGES - "ERREUR, MAUVAISE SPECIFICATION DE GRILLE, (TYPE) (IGAXG)"
*
*-------------------------------------------------------------------
*

      IF ((CGTYP .EQ. 'N') .OR. (CGTYP .EQ.'S')) THEN
        IF(IG4 .LT. 32768) THEN    ! ANCIEN STYLE DE CODAGE
          XG1 = IG2 * 0.1
          XG2 = IG1 * 0.1       ! PJ
          XG3 = IG4 * 100.      ! D60
          XG4 = IG3 * 0.01      ! DGRW
        ELSE                    ! NOUVEAU STYLE PLUS GENERAL
          JG3 = IG3
          JG4 = IG4
          JG4 = JG4 - 32768
          XG3 = IG1 * 100.      ! D60 NORMALEMENT EN HECTOMETRES
          IF(IG3 .GT. 32767) THEN   ! C'EST EN KILOMETRES
            XG3 = XG3 * 10.
            JG3 = JG3 - 32768
          ENDIF
          XG4 = IG2 * .1        ! ABS(DGRW)
          IF(JG4 .GT. 16383) THEN  ! DGRW NEGATIF
            XG4 = 360. - XG4
            JG4 = JG4 - 16384
          ENDIF
          DLAT = 90. -(JG4*180./16383.)
          DLON = (JG3*360./32767.)
          IHEM = 1
          IF('S'.EQ.CGTYP) IHEM = 2
          CALL XYFLL(XG1,XG2,DLAT,DLON,XG3,XG4,IHEM)
          XG1 = 1.0 - XG1
          XG2 = 1.0 - XG2
        ENDIF

      ELSE IF(CGTYP .EQ. 'C') THEN   ! C TYPE LAT LON GRID
        XG1 = IG3 * 0.01 - 90.
        XG2 = IG4 * 0.01
        XG3 = 180. / IG1
        XG4 = 360. / IG2

      ELSE IF ((CGTYP .EQ. 'A') .OR. (CGTYP .EQ. 'B') .OR.
     %        (CGTYP .EQ. 'G')) THEN
        XG1 = IG1
        XG2 = IG2
        XG3 = 0.
        XG4 = 0.

      ELSE IF(CGTYP .EQ. 'L') THEN
        XG1 = IG3 * 0.01 - 90.
        XG2 = IG4 * 0.01
        XG3 = IG1 * 0.01
        XG4 = IG2 * 0.01

      ELSE IF(CGTYP .EQ. 'H') THEN
        XG1 = IG3                 !  PHI1 / PHI2
        XG2 = .01*IG4 - 90.       !  PHI0
        XG3 = 500*IG2             !  DELTA S  (DEMI KILOMETRES)
        XG4 = IG1*.2              !  LAMBDA0
*                                 !  GRILLE LAMBERT CONFORME CENTREE NORD

      ELSE IF(CGTYP .EQ. 'E') THEN
        I2B = IAND(IG3,3)
        LG3 = ISHFT(IG3,-2)
        LG1 = IOR(ISHFT(IG1,2),I2B)
        I2B = IAND(IG4,3)
        LG4 = ISHFT(IG4,-2)
        LG2 = IOR(ISHFT(IG2,2),I2B)
        if (lg2.gt.3600)  lg2=lg2-7201
        XG1 = (LG1 -3600.0D0) / 40.0
C
C       bug de code, le +90 est de trop, ce qui peut causer un debordement
C       pour ig3
C
        if(lg3 .lt. 3559) lg3=lg3+16384
        XG2 = (LG3 -3600.0D0) / 40.0
        XG3 = LG2 / 40.0D0
        XG4 = LG4 / 40.0D0
*                                 !  GRILLE LAT,LON (GEF)
      ELSE
        WRITE(6,600)
      ENDIF

  600 FORMAT(1H0,' ERREUR, MAUVAISE SPECIFICATION DE GRILLE, (TYPE)',
     % '(IGAXG)')
      RETURN
      END
