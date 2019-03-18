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
***S/R CONSTNT  -  LIT LE FICHIER DE CONSTANTES COMMUNES CMC-RPN
*
      SUBROUTINE CONSTNT(VALEUR,FLAG,NOM,MODE0)
      IMPLICIT NONE
      REAL VALEUR
      INTEGER FLAG,MODE0
      CHARACTER *(*) NOM
      INTEGER IERROR
      external dummy_do_nothing_routine
      CALL CONSTNT_X(VALEUR,FLAG,NOM,MODE0,
     %               dummy_do_nothing_routine,-1,-1,-1,IERROR)
      RETURN
      END
      subroutine dummy_do_nothing_routine()  ! to support a mode 4 call
      return
      end subroutine dummy_do_nothing_routine

      SUBROUTINE CONSTNT_X(VALEUR,FLAG,NOM,MODE0,
     %                     BCAST_MPI,DATATYPE,ROOT,COMM,IERROR)
*
      IMPLICIT NONE
      REAL VALEUR
      INTEGER FLAG,MODE0,DATATYPE,ROOT,COMM,IERROR
      CHARACTER *(*) NOM
      EXTERNAL BCAST_MPI
*
*AUTEURS  N. BRUNET ET M. VALIN - FEV 91
*
*REVISION 001 G. PELLERIN - OCTOBRE 1991 - ADAPTATION A unix
*         002 M.VALIN aller chercher le fichier de constantes
*             operationnel par defaut suf si fichier "constantes"
*             present localement  AVRIL 93
*         003 M.VALIN janvier 2000 ajout d'un mode MPI avec
*             broadcast + quelques corrections
*         004 M. DESGAGNE - Decembre 2006 - Introduce BLOCK DATA
*             DATA_CONSTNT_X and return a valid flag for section
*             (NAME(1).EQ.' ')
*         005 M. Valin - Novembre 2014 corrrection de l'appel "broadcast"
*             l'argument DATATYPE ne sert plus
*
*OBJET(CONSTNT)
*     LIT LE FICHIER DE CONSTANTES COMMUNES CMC-RPN.
*     LA VALEUR DE L'ARGUMENT "MODE" INDIQUE LE MODE D'UTILISATION
*     DE LA ROUTINE. L'UTILISATEUR PEUT AVOIR UNE LISTE IMPRIMEE
*     COMPLETE DU FICHIER, OU PEUT OBTENIR LA VALEUR D'UNE
*     CONSTANTE OU PEUT AJOUTER OU MODIFIER ,LOCALEMENT AU TRAVAIL,
*     UNE CONSTANTE. Ce module peut etre utilise en mode MPI
*     et propager a ses partenaires les tables.
*
*PARAMETRES
*     MAXCNS - NOMBRE MAXIMUM DE CONSTANTES DANS LE FICHIER
      INTEGER MAXCNS
      PARAMETER(MAXCNS=200)
*
*FICHIERS
*     ACCES AU FICHIER DE CONSTANTES SE FAIT A L'INTERIEUR MEME
*     DE LA SOUS-ROUTINE.
*
*ARGUMENTS
*     VALEUR (E/S) - VALEUR DE LA CONSTANTE DEMANDEE, AJOUTEE
*                    OU MODIFIEE.
*     FLAG (S) - SI FLAG=0, LA DEMANDE FUT UN INSUCCES, I.E.
*                LA CONSTANTE N'EXISTE PAS POUR LE CAS OU ON VEUT
*                OBTENIR OU MODIFIER UNE VALEUR, OU LA CONSTANTE
*                QU'ON DESIRE AJOUTER EXISTE DEJA.
*              - SI FLAG=1, LA DEMANDE FUT UN SUCCES.
*     NOM (E)  - NOM DE LA CONSTANTE DONT ON DESIRE OBTENIR,
*                AJOUTER OU MODIFIER LA VALEUR.
*                MAXIMUM DE 8 CARACTERES.  EX. 'OMEGA'
*     MODE (E) - SI MODE=0, ON VEUT OBTENIR LA VALEUR D'UNE CONSTANTE
*              - SI MODE=1, ON VEUT UN IMPRIME COMPLET DU FICHIER
*              - SI MODE=2, ON VEUT AJOUTER (LOCALEMENT AU TRAVAIL)
*                           UNE CONSTANTE
*              - SI MODE=3, ON VEUT MODIFIER (LOCALEMENT AU TRAVAIL)
*                           UNE CONSTANTE.
*              - SI MODE=4, ON VEUT PROPAGER LE CONTENU DES TABLES
*                           VIA L'ARGUMENT BCAST_MPI. 
*              par MODE, on entend mod(MODE0,100)
*              MODE0 >  100 : VALEUR est un REAL*8
*              MODE0 <= 100 : VALEUR est un REAL*4
*     BCAST_MPI (E)   L'ARGUMENT BCAST_MPI DOIT ETRE RPN_COMM_BCST_X
*     DATATYPE,ROOT,COMM (E) :  ne servent plus 
*     IERROR (S)  code d'erreur renvoye par RPN_COMM_BCST_X
*

      INTEGER, EXTERNAL :: FNOM,FCLOS
      INTEGER MODE
*
*NOTE
*     SI ON AJOUTE UNE CONSTANTE (MODE=2), NE PAS OUBLIER DE
*     MODIFIER EN CONSEQUENCE LE COMDECK CONSDYN OU CONSPHY
*      QUI RESIDE DANS LA LIBRAIRIE 'PHYSIQUE' ET/OU 'DYNAMIQUE',
*     AINSI QUE LE PARAMETER (NBRE=..) DANS LA ROUTINE INCTDYN
*     ET DANS LA ROUTINE INCTPHY (ROUTINE QUI INITIALISE
*     LE COMMON 'CTESDYN' OU 'CTESPHY'.
*     avant de proceder a un appel en MODE 4, il faut que le "root pe" ait fait 
*     un appel en mode 0/1/2/3 pour faire initialiser les tables
*
*     ces tours de passe-passe sont necessaires pour eviter une reference 
*     externe a la librairie MPI, ce qui rendrait sa presence obligatoire
*     pour tout appel a constnt/constnt_x et constituerait un ennui serieux
*     pour tout programme non MPI qui se servirait de ces routines
*
**
*-------------------------------------------------------------------
      real *8 valeur8
      pointer (pv8,valeur8)
      INTEGER I, NCNS, IER, IUNREAD, istrt(2), iend(2),COUNT, dsize
      LOGICAL FEXIST
*
      REAL*8 VAL(MAXCNS)
*
      CHARACTER *8 NAME(MAXCNS), TNAME
      CHARACTER *42 TEXT
*
      COMMON/CONSTNT_PDATA/ istrt, NAME, VAL, NCNS, iend
      DATA NAME /MAXCNS * ' '/
*
*-----------------------------------------------------------------
      MODE = mod(mode0,100)
      pv8  = loc(valeur)
      IF (MODE .EQ. 4) THEN
        dsize = ( loc(istrt(2)) - loc(istrt(1)) )        ! taille d'un "integer"
        COUNT = ( loc(iend(1)) - loc(istrt(1)) ) / dsize ! nombre d' "integers" dans /CONSTNT_PDATA/
        CALL BCAST_MPI(istrt,COUNT,dsize,ROOT,COMM,IERROR)
!       ca donne  (root et comm seront ignores)
!       call RPN_COMM_bcst_x(istrt,COUNT,dsize,ROOT,COMM,IERROR)
!       et par ricochet
!       call RPN_COMM_bcast (istrt, count*dsize/4, "MPI_INTEGER",0, "grid", err)
        RETURN
      ENDIF
      IF (NAME(1).EQ.' ') THEN   ! lire les tables en provenance du fichier approprie

         FLAG=0
         IUNREAD=0
         INQUIRE(FILE='./constantes',EXIST=FEXIST)
         IF (FEXIST) THEN
           IER=FNOM(IUNREAD,'constantes','FTN+SEQ+FMT',0)
         ELSE
           IER=FNOM(IUNREAD,'@thermoconsts','FTN+SEQ+FMT',0)
         ENDIF
         if (IER.ne.0) return  ! pas trouve le fichier de constantes
*
         IF (MODE.EQ.1) THEN
            WRITE(6,600)
            WRITE(6,602)
         END IF

         DO I=1,MAXCNS
            READ(IUNREAD,'(2X,A8,2X,E20.13,2X,A42)',END=2,err=3)
     $          NAME(I),VAL(I),TEXT
            IF (MODE.EQ.1) WRITE(6,605)NAME(I),VAL(I),TEXT
            NCNS = I
         ENDDO
 2       FLAG = 1

 3       IER = FCLOS(IUNREAD)
*
         IF (MODE.EQ.1) RETURN
      END IF
*
*     ----  FIN DE L'INITIALISATION DE NAME ET VAL ----
*
      TNAME = NOM
*
*     OPTION  MODE=0, ALLER CHERCHER LA VALEUR D'UNE CONSTANTE
*
      IF (MODE.EQ.0) THEN
         FLAG = 0
         if (mode0.lt.100) then
           VALEUR = 0
         else
           valeur8=0
         endif
         DO I=1,NCNS
            IF (TNAME.EQ.NAME(I)) THEN
               FLAG = 1
               if (mode0.lt.100) then
                VALEUR = VAL(I)
               else
                valeur8=val(i)
               endif
               RETURN
            END IF
         enddo
*
*     OPTION  MODE=2, AJOUTER UNE CONSTANTE
*
      ELSE IF(MODE.EQ.2)THEN
         FLAG = 0
*        VERIFIER SI CONSTANTE EXISTE DEJA
         DO I=1,NCNS
            IF (TNAME.EQ.NAME(I)) THEN
               WRITE(6,620)TNAME
620            FORMAT(/,5X,'**** LA CTE',1X,A8,1X,'EXISTE DEJA',/)
               STOP   ! return avec flag = 0 ?
            END IF
         enddo
         IF (NCNS.LT.MAXCNS) THEN
            FLAG = 1
            NCNS = NCNS + 1
            if (mode0 .lt.100) then
              VAL(NCNS) = VALEUR
            else
              val(ncns)=valeur8
            endif
            NAME(NCNS) = TNAME
         ELSE
            WRITE(6,625)
625         FORMAT(/,5X,'**** NB DE CONSTANTES DEPASSE LA LIMITE',/)
            STOP   ! return avec flag = 0 ?
         ENDIF
*
*     OPTION  MODE=3,  MODIFIER LA VALEUR D'UNE CONSTANTE
*
      ELSE IF (MODE.EQ.3) THEN
         FLAG = 0
*        VERIFIE SI CONSTANTE EXISTE
         DO I=1,NCNS
            IF (TNAME.EQ.NAME(I)) THEN
               FLAG = 1
               if (mode0 .lt.100) then
                 VAL(I) = VALEUR
               else
                 val(i)=valeur8
               endif
               exit
            ENDIF
         enddo
         IF (FLAG.EQ.0) THEN
            WRITE(6,630)TNAME
            STOP
         ENDIF
*
      ENDIF
*
600   FORMAT(1H1,10X,'LISTE DES CONSTANTES COMMUNES CMC-RPN',///)
602   FORMAT(2X,'NOM',17X,'VALEUR',8X,'DESCRIPTION',10X,'UNITE',//)
605   FORMAT(2X,A8,2X,E20.13,2X,A42)
630   FORMAT(/,5X,'**** LA CTE A MODIFIER',1X,A8,1X,"N'EXISTE PAS",/)
      RETURN
      END
