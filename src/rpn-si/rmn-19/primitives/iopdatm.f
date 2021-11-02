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
***FONCTION IOPDATM  ALLER CHERCHER LE DATE TIME STAMP OPERATIONNEL
*
      INTEGER FUNCTION IOPDATM(FLAG)
      Implicit NONE
      CHARACTER*(*)FLAG
*
*AUTEUR  M.VALIN RPN MARS 1983
*
*LANGAGE RATFOR
*        CE MODULE CONTIENT DES DEPENDANCES CDC
*        ON Y PROCEDE A UN ATTACH (OU GET) INTERNE
*
*OBJET(IOPDATM)
*        ALLER CHERCHER A L'ENDROIT APPROPRIE LE DATE TIME STAMP OPERATIONNEL.
*        CE MODULE MARCHE SOUS SCOPE 2 ET SOUS NOS
*
*CALL STDLIB
*
*ARGUMENTS
*  E     FLAG      CHAINE DE CARACTERES (7), FLAG PEUT ETRE
*                  - 'NON' AUQUEL CAS IOPDATM RENVOIE 2010101011.
*                  - UNE CLE DE 7 CARACTERES OU MOINS CORRESPONDANT
*                    A D'AUTRES "DATE TIME STAMP".
*                  - UN NOMBRE DE 7 CHIFFRES YYJJJZZ OU YY EST L'ANNEE,
*                    JJJ LE JOUR DANS L'ANNEE (1 A 366) ET ZZ L'HEURE
*                    (00 - 24).
*                  - IOPDATM RENVOIE ALORS LE DATE TIME STAMP QUI
*                    CORRESPOND A CETTE DATE ET CETTE HEURE.
*                  - UN NOM DE FICHIER. IOPDATM LIRA LE PREMIER MOT DE CE
*                    FICHIER EN LE CONSIDERANT COMME UN DATE TIME STAMP.
*                    SI IOPDATM RENCONTRE UNE ERREUR, LE "STAMP" RENVOYE
*                    SERA 1010101011.
*
*MODULES
*        PF        POUR ATTACHER LE BON FICHIER CONTENANT LA DATE
*
*NOTES
*        IOPDATM SE SERT DU FICHIER FORTRAN 88 ET LE RETOURNE (RETURN) A
*        APRES USAGE.
*
**
      INTEGER ISTAMP,IST,JD,IYR,IMON,IWK,IDAY,ISTMP,L,LONGUEUR,IUN
      INTEGER ier,fnom,fclos
      CHARACTER* 10 IFFLG,IQUOI,IDNT
      CHARACTER * 128 DATAREP
      EXTERNAL system_time,fnom,fclos,jdatec,datec,newdate,longueur
      INTEGER i1,i2
      character *26 upper, lower

      upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      lower='abcdefghijklmnopqrstuvwxyz'
      IUN = 88
      IFFLG=FLAG
      do i1=1,len(IFFLG)
        i2=index(lower,IFFLG(i1:i1))
        if(i2.gt.0) IFFLG(i1:i1)=upper(i2:i2)
      enddo
      READ(IFFLG,'(A7)')IQUOI
      IF(IFFLG.EQ.'OUI') IQUOI='OPRUN'
      ISTAMP = 010101011
      IF(IFFLG.EQ.'NON') THEN
        ISTAMP = 010101011
      ELSE IF(IFFLG.EQ.'NOW') THEN
        call system_time(i1,i2)
        call newdate(ISTAMP,i1,i2,3)
      ELSE
         READ(IFFLG,'(I10)',ERR=100)ISTAMP
         if(ISTAMP.le.9936624) then  ! it is in YYJJJHH format
           IYR=1900+ISTAMP/100000
           if(iyr.lt.1950) iyr=iyr+100
           CALL JDATEC(JD,iyr,1,1)
           CALL DATEC(JD+MOD(ISTAMP/100,1000)-1,IYR,IMON,IDAY)
           call newdate(ISTAMP,IYR*10000+IMON*100+IDAY,
     %               MOD(ISTAMP,100)*1000000,3)
C          ISTAMP = IMON*10000000 + IDAY*100000 + (IYR-1900)*1000 +
C    %            MOD(ISTAMP,100)*10 + 1
         endif
         GOTO 101
 100     CONTINUE
         ISTAMP = 010101011
         CALL GETENV('CMC_OCMPATH',DATAREP)
         L = LONGUEUR(DATAREP)
         IF (L .GT. 0) THEN 
           IUN = 0
           IER = FNOM(IUN,DATAREP(1:L)//'/datafiles/data/uspmadt',
     %                'SEQ+FTN+FMT',0)                     
         ELSE
           CALL GETENV('AFSISIO',DATAREP)
           L = LONGUEUR(DATAREP)
           IUN = 0
           IER = FNOM(IUN,DATAREP(1:L)//'/datafiles/data/uspmadt',
     %               'SEQ+FTN+FMT',0)
         ENDIF
 200     READ(IUN,'(A7,1X,I9)',END=300)IDNT,ISTMP
         IF(IDNT .EQ. IQUOI) THEN
            ISTAMP = ISTMP
         ELSE
            GOTO 200
         ENDIF
 300     CLOSE(IUN)
 101     CONTINUE
         IF(ISTAMP.EQ. 010101011) THEN
           OPEN(IUN,FILE=FLAG,FORM='FORMATTED')
           READ(IUN,'(I9)',END=400)ISTAMP
 400       CLOSE(IUN)
           IER = FCLOS(IUN)
         ENDIF
      ENDIF

      IOPDATM = ISTAMP

      RETURN
      END
