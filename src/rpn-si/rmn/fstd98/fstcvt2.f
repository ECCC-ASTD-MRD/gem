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
***FONCTION Fstcvt2 HOLLERITH A CARACTERE OU L'INVERSE
      INTEGER FUNCTION Fstcvt2( NOM, TYP, ETIK, GRTP, CNOM, CTYP,CETIK,
     % CGRTP, HOLACAR)
      IMPLICIT NONE
      INTEGER Fstcvt
      INTEGER NOM, TYP, ETIK(2), GRTP
      CHARACTER *(*) CNOM
      CHARACTER *(*) CTYP, CGRTP
      CHARACTER *(*) CETIK
      LOGICAL HOLACAR

*AUTEUR  P. SARRAZIN   NOV  1989
*
*REVISION   000 NOUVELLE VERSION NOVEMBRE 1989
*         1.0.1 -  AJOUTER 2 TABLES POUR FICHIER RANDOM/SEQ STANDARD
*                   MODIFICATIONS 22 JANV 89
*           002 - M. Lepine - Fevrier 1998, eliminer le commno
*           003 - M. Lepine - Mars 1998, extensions pour le nombre de caract.
*
*LANGAGE ratfor
*
*OBJET(Fstcvt)
*         VARIABLE UTILISEE COMME HOLLERITH SERA TRANSFORME EN CARACTERE
*         OU L'INVERSE POUR NOMVAR,TYPVAR,GRID TYPE, ET ETIKET LE MOT
*         ETIKET AURA 2 LOCATIONS SUR LE SUN
*
*ARGUMENTS
*  IN OUT    NOM       HOLLERITH *4                       [NOMVAR]
*  IN OUT    TYP       HOLLERITH *2                       [TYPVAR]
*  IN OUT    ETIK      HOLLERITH *12   2 MOTS A4          [ETIKET]
*  IN OUT    GRTP      HOLLERITH *1                       [GRTYP]
*  OUT IN    CNOM      CARACTERE *4
*  OUT IN    CTYP      CARACTERE *2
*  OUT IN    CETIK     CARACTERE *12
*  OUT IN    CGRTP     CARACTRE *1
*  IN        HOLACAR   LOGICAL .TRUE.  HOLLERITH A CARATERE
*                      LOGICAL .FALSE. CARACTERE A HOLLERITH

      INTEGER I

      Fstcvt2 = 0
      IF((HOLACAR))THEN

*       TRANSFER STRING D'UNE LOCATION HOLLERITH EN CARACTERE

         IF((NOM.EQ. -1))THEN
            CNOM = ' '
         ELSE
            IF (LEN(CNOM) .GT. 2) THEN
              WRITE(CNOM,'(A4)') NOM
            ELSE
              WRITE(CNOM,'(A2)') NOM
            ENDIF
         ENDIF
         IF((TYP.EQ. -1))THEN
            CTYP = ' '
         ELSE
            IF (LEN(CTYP) .GT. 1) THEN
               WRITE(CTYP,'(A2)') TYP
            ELSE
               WRITE(CTYP,'(A1)') TYP
            ENDIF
         ENDIF
         IF((GRTP.EQ. -1))THEN
            CGRTP = ' '
         ELSE
            WRITE(CGRTP,'(A1)') GRTP
         ENDIF
         CETIK = ' '
         IF((ETIK(1).EQ. -1))THEN
            CETIK = ' '
         ELSE
            IF (LEN(CETIK) .GT. 8) THEN
               WRITE(CETIK,600)(ETIK(I),I=1,3)
            ELSE
               WRITE(CETIK,601)(ETIK(I),I=1,2)
            ENDIF
         ENDIF
      ELSE

*       TRANSFER STRING D'UNE LOCATION CARACTERE EN HOLLERITH*
         READ(CNOM,'(A4)') NOM
         IF((CNOM.EQ. ' '))THEN
            NOM = -1

         ENDIF
         READ(CTYP,'(A2)') TYP
         IF((CTYP.EQ. ' '))THEN
            TYP = -1

         ENDIF
         READ(CGRTP,'(A1)') GRTP
         IF((CGRTP.EQ. ' '))THEN
            GRTP = -1

         ENDIF
         IF (LEN(CETIK) .GT. 8) THEN
            READ(CETIK,600) (ETIK(I),I=1,3)
         ELSE
            READ(CETIK,601) (ETIK(I),I=1,2)
         ENDIF
         IF((CETIK.EQ. ' '))THEN
            ETIK(1) = -1
         ENDIF

      ENDIF
600   FORMAT(3A4)
601   FORMAT(2A4)
666   FORMAT(' HOLLERITH NOM=',A4,' TYP=',A2,' GRTP=',A1,' ETIK= ',
     %A12)
668   FORMAT(' HOLLERITH NOM=',A4,' TYP=',A2,' GRTP=',A1,' ETIK= ',
     %A4,A4,a4)
667   FORMAT(' CARACTERE NOM=',A4,' TYP=',A2,' GRTP=',A10,'ETIK= ',
     %A12)
      RETURN
      END
