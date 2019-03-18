!/* EDITFST - Collection of useful routines in C and FORTRAN
! * Copyright (C) 1975-2014  Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!**PROGRAMME EDITFST - COPIE UNE PARTIE D'UN FICHIER STANDARD DANS UN 
!                      AUTRE FICHIER DTANDARD.
      subroutine EDITFST
!****************************************
!********    VERSION  'UNIX'    *********
!****************************************
      use configuration
      IMPLICIT NONE 
!
!         AUTEURS                                         DATE   VERSION
!         VERSION ORIGINALE (COPYSTD) C. THIBEAULT  -     FEV. 83
!         Y.BOURASSA/ M.VALIN (EDITSTD)                   JUL. 86 1.3 
!         BOURASSA         CORRECTION DANS SEQCOPY        FEV. 87 1.4 
!         VALIN            AJOUT DE FASTIO                FEV. 87 1.5 
!         BOURASSA         DIRECTIVE "CRITSUP"            OCT. 89 1.6 
!            "             DIRECTIVE "EXCLURE"            OCT. 89 1.6 
!            "             ACCEPTE "OPDATE"               JAN. 90 1.7 
!            "             NEWEDIT VERSION UNIX/CRAY      JUL. 90 1.8 
!            "             VERSION EDITFST (CRAY/UNIX)    NOV. 90 2.0 
!            "             QLXINS, FICHIERS FST..ETC. 
!            "             CORRECTION DANS COPISEQ        DEC. 90 2.1 
!            "             OPTION IMAGE                   JAN. 91 2.2 
!            "             ARGDOPE DES ETIKETS            MAR. 91 2.3 
!            "             LEGERE MODIF DANS PRINT        AVR. 91 2.31
!            "                "     "   SNOM-DNOM         MAI. 91 2.4  
!            "             DIRECTIVE ZAP                  JUL. 91 2.5
!            "             LINQUAGE DES FICHIERS SOURCES  OCT. 91 2.6
!            "             ACCELERATION DU TRIAGE         NOV. 91 2.7
!            "             QLXINX POUR NEC, OTE FRETOUR,
!                          CHANHE APPEL A CCARD & HOLACAR JAN. 92 2.8
!            "             FIX BUG 1e  ALLOCATION MEMOIRE JAN. 92 2.9
!            "             AJOUTE +OLD AU FICHIER SOURCE  MAR. 92 3.0
!                          ACCES A $IN VERSION CRAY       
!                          APPEL A LOW2UP APRES CCARD
!                          OUVRE SOURCES AVEC FSTNOL 
!            "             PERMET FICHIERS FTN SUR UNIX
!            "             NEUTRALISE LA CLE NPD          MAR. 92 3.1
!            "             FIX BUG DANS EXDES             AVR. 92 3.2
!            "             (SAUVDEZ=0) ZAP(-1) APRES COPIEMAI. 93 3.3
!                          QLXINX POUR TYPE=2    
!                          MODIF DANS EXDES (RANGE INVERSE) 
!                          FIX BUG DANS STDCIPI/WEOFILE 
!                          OUVRE LES 1 @ 15 INPUTS QUI EXTSTENT  
!                          ENDIQUE L'ETAT DU PGM. A LA SORTIE
!                          EVITE LES ABORTS EN MODE INTERACTIF
!            "             CORRECTION FASTIO              JAN. 94 3.4
!                          PASSAGE DE 15 FICHIERS D'ENTREE A 35    
!            "             BUG FIX DATA NUMERO D'UNITES   JUIL 94 3.5
!            "             BUG FIX FICHIER SEQ EXISTANT   JUIL 94 3.6
!            "             BUG FIX LIMITE FICHIER SEQ     AVR. 95 3.7
!         Y. Bourassa      ENLEVE REFERENCES A 1950       MAI  95 3.8
!         M. Lepine        Utilisation des waio "cmcarc"  mai  95 3.9
!         M. Lepine        Reload avec waio modifie tres legerement sept 95 4.0
!         M. Lepine        "work around" pour bug -i 0 nov. 95 4.1
!         M. Lepine        Reload avec librmnx32stack.a avril 97 v4.2
!         M. Lepine        Reload avec c_baseio.o (2 write)  avril 97 v4.3
!         M. Lepine        modif message dans c_baseio (2 write)  avril 97 v4.4
!         M. Lepine        Reload avec fstd89.o (bugfix fstecr) mai 97 v4.5
!         M. Lepine        Extensions pour fstd98 - mars 98 - v5.0
!                          ouverture des fichiers sources en R/O.
!         M. Lepine        Bug fix date stamp superieur exdes - sept 98 - v5.1
!         M. Lepine        Rechargement avec dernier release - v5.2
!         M. Lepine        Ajout de la clef nrecmin  - v5.3
!         M. Lepine        Reload, bug fix c_fnom  - mars 1999 - v5.4
!         M. Lepine        Menage des if defined - mars 2000 - v5.5
!         M. Lepine        Reload pour datev (fstd89) - oct 2000 - v5.6
!         M. Lepine        Possibilite d'appel a convip pour les ip1 
!                          + fichier source read only - juil 2001 - v5.7
!         M. Lepine        Julhr en julsec - oct 2001 - v5.8
!         M. Lepine        fichier de sortie R/W obligatoire - Dec 2001 - v5.81
!         M. Lepine        verification de la limite de 10 elements - Fev 2002 - v5.82
!         M. Lepine        reload, bug fix fstecr mode reecriture - Mars 2002 - v5.83
!         M. Lepine        Codification des ip1 en integer*8 (kind,level) - Mai 2002 - v5.84
!         M. Lepine        Bug fix buftemp dans copystx - Oct 2002 - v5.85
!         M. Lepine        Bug fix arguments fstcvt dans critsup  - Nov 2002 - v5.86
!         M. Lepine        Tolerance d'erreur pour les IP1 (ip1equiv)  - Mars 2003 - v5.87
!         M. Lepine        Correction pour code 5 (hybrid)  - Avril 2003 - v5.88
!         M. Lepine        Correction ip1equiv cas ipa ipb = 0  - Juin 2003 - v5.89
!         M. Lepine        Reload a partir de rossby - Avril 2004- v5.90
!         M. Lepine        Eviter de traiter des fichiers endommages - Mai 2004- v5.91
!         M. Lepine        Version avec compresseurs et bugfix cmcarc Linux - Jan 2005- v5.92
!         M. Lepine        Reload avec librmn_x, bug fix requete ip1=0.0 - Fev 2005 - v5.93
!         M. Lepine        Bug fix reconnaissance des ip1 reel - Mars 2005 - v5.94
!         M. Lepine        Utilisation optionnelle des fichiers remote - Mars 2005 - v5.95
!         M. Lepine        Reload pour bugfix datev de excdes, librmnbeta - Mars 2005 v5.96 
!         M. Lepine        Reload pour bugfix longueur fstecr mode image, librmnbeta - Avril 2005 v5.97
!         M. Lepine        De 35 fichiers d'entree a 120, Juillet 2005 - v5.98 
!         M. Lepine        Remplacement de tous les fstabt par qqexit - v5.99
!         M. Lepine        Reload avec librmn_rc008.a - v6.00
!         M. Lepine        Fichier source read only dans stdcopi - v6.01 - Fev 2006 
!         M. Lepine        Bug fix longeur fstluk,fstecr datyp=6 mode image - v6.02 - Avril 2006
!         M. Lepine        Reload avec librmn_008 - v6.03 - Aout 2006
!         M. Lepine        Reload avec librmnbeta,nouveau compresseurs IEEE - v6.04 - Oct 2006
!         M. Lepine        Reload avec librmnbeta,bug fix valid record length = 0 - v6.05 - Fev 2007
!         M. Lepine        Taille des noms de fichiers augemntee a 1024 car. - v6.06 - Mars 2007
!         M. Lepine        Reload avec librmn_009 - v6.07 - Juin 2007
!         M. Lepine        Reload avec librmnbeta, correction fichiers > 2G - v6.08 - Sept 2007
!         M. Lepine        Reload avec librmnbeta, correction fichiers > 2G - v6.09 - Mars 2008
!         M. Lepine        Reload avec librmnbeta, correction lng mode image pour datyp=6 - v6.10 - Juillet 2008
!         M. Lepine        Reload avec librm_rc010, correction pour fichier cmcarc remote - v6.11 - Sept 2008
!         M. Lepine        Reload avec librm_rc010, correction pour alias dans gossip_sock - v6.12 - Oct 2008
!         M. Lepine        Correction pour empecher le remapping du datatype  - v6.13 - Fev 2009
!         M. Lepine        Reload avec librmn_012 et codebeta moduledate  - v6.14 - Dec 2010
!         M. Lepine        Reload avec librmn_012 et codebeta c_baseio_714 (remote)  - v6.15 - Oct 2011
!         M. Lepine        La correction pour empecher le remapping est faite dans fstecr.
!                          Reload avec librmn_013_rc2 - v6.16 - Mars 2012
!         M. Lepine        Reload avec librmn_013 - v6.17 - Sept 2012
!         B. Dugas         Correction dates etendues dans zap et copystx - v6.18 - Sept 2012
!         M. Valin         Ajustements pour traduction IP1/ip2/IP3, librmn_014, dryrun      - v6.19 - Fev 2014
!         M. Valin         refactoring, fold comdecks into a single module, rename to .F90  - v6.20 - Fev/Mar 2014
!                          replacement of record selection logic. version 015 or newer is NEEDED for librmn
!                          subroutine exdes no longer used
!         M. Valin         new code using extended selection features from standard file package
!                          rmnlib alpha 16 minimum - v7.0a - april 2015
!         M. Lepine        Remplacement ou elimination des variables a 128 car. pour les noms de fichiers
!                          v7.1a - oct 2015
!         M. Valin         Correction d'un bug de logique dans sauvdez - v7.2a - aout 2016
!         M. Lepine        Elimination des espaces blancs a l'impression - v7.3b - sept 2106
!         M. Lepine        Remettre l'initialisation du package convip en mode newstyle - v7.4b - sept 2016
!         M. Valin         Correction du traitement des desire/exclure dans excdes_new - v7.5b - sept 2016
!         M. Valin         Correction du traitement des selections avec delta dans excdes_new - v7.6b - oct 2016
!         M. Lepine        Correction dans sauvdez, remettre le compteur NREQ a zero - v7.7 - nov 2016
!         M. Valin         Bug fix pour le cas desire avec tous les arguments a -1 - v7.8 - sept 2017
!
!LANGAGE  - FTN77
!
!OBJET(EDITFST)
!         - COPIE UN FICHIER STANDARD (RANDOM OU SEQUENTIEL) DANS
!           UN AUTRE FICHIER STANDARD (RANDOM OU SEQUENTIEL).
!           UTILISE CCARD POUR RAMASSER LES PARAMETRES SUR L'ENONCE
!           D'EXECUTION DU PROGRAMME (DOIT AVOIR LA FORME SUIVANTE)
!           EDITFST -s(noms des fichiers sources)
!                   -d (nom  du  fichier  destination)
!                   -i (nom  du  fichier  stdinp)
!                   -l (nom  du  fichier  stdout)
!                   -ss               ( s=sequentiel sqi)         
!                   -sf               ( s=sequentiel fortran)
!                   -ds               ( d=sequentiel sqi)
!                   -df               ( d=sequentiel fortran)
!                   -n                ( pas de boite)
!                   -e                ( reecrire un enregistrement dans d)
!                   -f                ( fast IO)
!                   -eof 0<entier<15  ( marquer d si sequentiel)
!                   -m inform         ( diagnostiques)
!                   -m debug          ( diagnostiques et debug)
!                   -dryrun           ( dryrun, messages seulement, pas d'ecriture dans fichier de sortie)
!                   -t                ( changer le niveau de tolerence)
!                   -v                ( un voir de d)
!                   -nrecmin          ( nombre minimal d'enregistrement)
!           FICHIERS (25@39 = SOURCE) (3 = DESTINATION)
!#include "maxprms.cdk"
!               NMR = MAXIMUM DE REGIONS
!               NMS =    "     " SCORES 
!               NME =    "     " ETAPES 
!               NMN =    "     " NIVEAUX
!               NMM =    "     " MODELES
!               NMD =    "     " DESIRES/EXCLURES 
!#include "tapes.cdk"
!     - MEOF       - LEVEL D'EOF LOGIQUE A NE PAS PASSER
!     - COPIES     - NOMBRE D'AJOUTS AU FICHIER DESTINATION EN USAGE
!     - NDS        - NOMBRE D'ENREGISTREMENT FICHIER SOURCE EN USAGE
!     - NDD        -    "           "           "    DESTYINATION EN USAGE
!     - EOF        - TERMINE COPIE SEQUENTIELLE PAR UNE MARQUE=EOF
!     - CEOF       - CONTROLE LES MARQUES DE FIN DE FICHIER (DESTINATION)
!                    =-1 COPIE TOUS   "    "  "   "    "
!                    = I ECRIT DES    "    "  "   "    "    DE NIVEAU I
!                    = 0 AUCUN TRANSFER DE MARQUE DE FIN DE FICHIER
!     - LEOF       - VALEUR DE LA DERNIERE MARQUE DE FIN DE FICHIER LOGIQUE
!                    RENCONTREE DANS SOURCE
!     - LIMITE     - NOMBRE MAXIMUM DE COPIES PERMIS PAR L'USAGER
!     - NFS        -    "   DE FICHIERS SOURCES SOURCES
!     - NFSO       -    "    "    "        "       "    OUVERTS
!     - SOURCES    - TABLEAU DES NUMERO DE FICHIER SOURCE
!#include "fiches.cdk"
!     - IST        - INDICATEUR DU PREMIER POINT DE LA MEMOIRE FLOTTANTE
!     - NP         - NOMBRE DE PARAMETRES DANS SUB. APPELEE PAR DIRECTIVE
!     - FIXD       - INDIQUE QUE LES ENREGISTREMENTS DU FICHIER SOURCE
!                    SONT TOUS VALIDE EN MEME TEMPS.
!     - ECR        - =.TRUE. POUR RECRIRE.
!     - SSEQ       -    "    FICHIER SOURCE EST SEQUENTIEL (.FALSE.RANDOM)
!     - DSEQ       -    "       "     DESTIN "       "             "
!     - VS         -    "    DESIRE UN VOIR DU FICHIER SOURCE
!     - VD         -    "       "    "    "   "   "    "    DESTINATION
!     - OUVS       -    "    FICHIER SOURCE EST OUVERT
!     - OUVD       -    "       "    DESTIN  "     "
!#include "logiq.cdk"
!     - SCRI       - .TRUE. SI DES CRITERES SUPLEMENTAIRE EN FORCE
!     - XPRES      - .TRUE. POUR COPIE COMPLETE DE SOURCE A DESTINATION
!     - ESAIS      -    "   SI UNE TENTATIVE DE COPIE A ETE FAITE.
!     - DM1        -    "    " PAS DE DATES DANS LES DESIRES
!     - DEBUG      -    "    " EN MODE DEBUG
!     - DRYRUN     -    "    " EN MODE DRYRUN
!     - SELEC      -    "    " FICHIER DE DIRECTIVE PRESENT 
!     - BOX        -    "    " LA CLE NOBOX PAS DANS LA SEQUENSE D'APPEL
!     - DIAG       -    "    " DIAGNOSTIQUES SERONT A IMPRIMER
!     - INTERC     -    "    " EN MODE INTERACTIF
!     - ZA         -    "    " ON VEUT MODIFIER UN ETIQUETTE A LA SORTIE
!#include "desrs.cdk"
!     - JOURS      - PERIODE A UTILISER PAR DESIRE/EXCLURE
!     - NREQ       - NOMBRE DE DIRECTIVES DESIRE/EXCLURE RENCONTREES
!     - SAUV       -    "    "     "      A CONSERVER APRES COPIE
!     - NEXC       -    "    "     "      EXCLURE 
!     - DESEXC(I)  - =0 (POUR EXCLURE),    =-1 (POUR DESIRE)
!     - SATISF(I)  - 0 = DIRECTIVE INSATISFAITE
!     - REQ        - TABLEAU DES DESIRES/EXCLURES DE L'USAGER.
!     - SUP        - TABLEAU DES CLES SUPLEMENTAIRES.
!     - NIS        - CRITERE SUPLEMENTAIRE DE CELECTION # 1 
!     - NJS        -    "          "        "     "     # 2 
!     - NKS        -    "          "        "     "     # 3 
!     - IG1S       -    "          "        "     "     # 4 
!     - IG2S       -    "          "        "     "     # 5 
!     - IG3S       -    "          "        "     "     # 6 
!     - IG4S       -    "          "        "     "     # 7 
!     - REQN       - NOMBRE DE NOMVAR/REQUETE
!     - REQT       - NOMBRE DE TYPVAR/REQUETE
!     - REQE       - NOMBRE D'ETIKET/REQUETE
!     - Z1         - IP1  A OsCHANGER SI DIFFERENT DE -1
!     - Z2         - IP2  "    "     "     "      "  "
!     - Z3         - IP3  "    "     "     "      "  "
!     - ZD         - DATE "    "     "     "      "  "
!     - ZG1@ZG4    - IG1 A IG4 A ZAPPER
!#include "key.cdk"
!     - KLE        - CLES DE LA SEQUENCE D'APPEL
!     - DEF1       - VALEURE DES CLES PAR DEFAUT
!     - DEF2       -     "    "    "  SI PRESENTES
!     - PRINTR     - DN DU FICHIER DE SORTIE D'IMPRIMENTE
!#include "char.cdk"
!     - NS         - DN DU FICHIER SOURCE
!     - ND         -  "  "    "    DESTINATION
!     - SNOM       - TYPE  PASSE A FNOM FICHIER SOURCE
!     - DNOM       -   "     "   "   "     "    DESTINATION 
!     - ZE         - CARACTERES DE L'ETIQUETTE A ZAPPER SI DIFFERENT DE '????????'
!     - ZG         - TYPE DE GRILLE A ZAPPER
!     - ETI        - ETIQUETTE TEMPORAIRE
!     - ETIS       - ETIQUETTES DES DESIRES/EXCLURES
!     - ZT         - CARACTERE DU TYPEVAR A ZAPPER SI DIFFERENT DE '?'.
!     - TYP        - TYPE DE GRILLE
!     - TYPS       - TYPE DE GRILLEDES DES DESIRES/EXCLURES 
!     - GTY        - TYPE DE GRILLE (SOURCE)
!     - GTYS       - TYPE DE GRILLE DES CRITERES SUPLEMENTAIRES
!     - GTYPS      - TYPE DE GRILLE DES CRITERES SUPLEMENTAIRES
!     - ZN         - CARACTERES DU NOMVAR A ZAPPER SI DIFFERENT '??'
!     - NOM        - NOM DE VARIABLE
!     - NOMS       - NOMS DE VARIABLE DES DESIRES/EXCLURES
!     - ETAT       - INDIQUE L'ETAT DU PGM. DANS LA BOITE A LA FIN
      EXTERNAL    FERMED, SELECT, FNOM, CCARD, OUVRES, SAUVDEZ, FERMES
      EXTERNAL    FSTOPC, FSTOPL, EXDB, EXFIN, OUVRED, STDCOPI, MEMOIRH
      EXTERNAL    QQEXIT
      INTEGER     FSTOPC, FSTOPL, EXDB, EXFIN, OUVRED, FNOM, I
      integer junk
      LOGICAL     FASTIO
!      character(len=*), parameter :: current_version="v 1.19"
      character(len=16) :: RELEASE
      character(len=32) :: SUB_RELEASE
      character(len=4096) ,dimension(:), pointer, save:: def1,     def2
      character(len=4096), save :: PRINTR

      include 'version.inc'

      call config_init   ! initialize values in module "configuration"
      allocate(def1(NCCARDKEYS),def2(NCCARDKEYS))
      def1 = def1b
      def2 = def2b

      do i = 1,120
        sources(i) = 24+i
      enddo
      SAUV = 0
      CALL SAUVDEZ
  
!     EXTRACTION DES CLES DE LA SEQUENCE D'APPEL. 
      I    = -111
      CALL CCARD(KLE, DEF2, DEF1, NCCARDKEYS, I)
      def1b = def1
      READ(DEF1(3), '(I2)') EOF                                 ! -eof
      READ(DEF1(15),'(I8)') LIMITE                              ! -c
      READ(DEF1(25),'(I8)') NRECMIN                             ! -nrecmin
      PRINTR = DEF1(11)                                         ! -l  fichier pour listing
      ND     = DEF1(2)                                          ! -d  fichier destination
      VD     = (DEF1(6) .EQ.'OUI')  .OR.  (DEF1(18).EQ.'OUI')   ! -vd , -v  voir destination
      VS     = (DEF1(20).EQ.'OUI')                              ! -vs   voir source
      BOX    = (DEF1(7) .EQ.'NON')  .AND. (DEF1(19).EQ.'NON')   ! -nobox , -n
      ECR    = (DEF1(9) .EQ.'OUI')  .OR.  (DEF1(21).EQ.'OUI')   ! -ecr , -e
      SELEC  = (DEF1(10).NE.'NON')  .AND. (DEF1(10).NE.'NIL') .AND. (DEF1(10).NE.'0')  ! -i 
!
!     Contourner le bug du -i 0 en ouvrant l'unite 5 sur /dev/null
!
      IF (.NOT. SELEC) THEN  ! def1(10) = "$IN" si cle non specifiee, selec = .true. par defaut
        SELEC = .TRUE.
        DEF1(10) = '/dev/null'
      ENDIF

      DEBUG  = (DEF1(13)  .EQ.'DEBUG')       ! -m  (non par defaut)
      DRYRUN  = (DEF1(146)  .EQ.'DRYRUN')    ! -dryrun  (non par defaut)
      FASTIO = (DEF1(22).EQ.'OUI')           ! -f  (oui par defaut)
      DIAG   = (DEF1(8) .EQ.'OUI')  .OR.  (DEF1(13).EQ.'INFORM') 
      IF((DEF1(16).NE.'NON') .OR.    &  ! -ss , source is sequential sqi
         (DEF1(4).NE.'NON')  .OR.    &  ! -sseq , source is sequential
         (DEF1(23).NE.'NON') ) THEN    ! -ss -sseq -sf
         IF(DEF1(23) .NE.'NON') THEN   ! -sf , source is sequential fortran
            SNOM = 'STD+SEQ+R/O+OLD'
         ELSE
            SNOM = 'STD+SEQ+OLD+R/O'
         ENDIF
      ELSE                             ! standard random file
         SNOM = 'STD+RND+OLD+R/O'
      ENDIF

      IF((DEF1(17).NE.'NON') .OR.    &  ! -ds , destination is sequential sqi
         (DEF1(5).NE.'NON')  .OR.    &  ! -dseq , destination is sequential
         (DEF1(24).NE.'NON') ) THEN    ! -df , destination is sequential fortran
         IF(DEF1(24) .NE.'NON') THEN
            DNOM = 'STD+SEQ+FTN'
         ELSE
            DNOM = 'STD+SEQ'
         ENDIF
      ELSE
         DNOM = 'STD+RND'
      ENDIF
      I = FNOM(6, PRINTR, 'SEQ', 0)    ! open listing file

!     IMPRIME L'INDICATIF DE DEBUT DU PROGRAMME.
      IF( BOX ) THEN
         I = EXDB('EDITFST',trim(RELEASE),'NON')
      ELSE
         WRITE(6,*)'***   E D I T F S T   '//trim(RELEASE)//'   ***'
      ENDIF


      IF( DIAG ) THEN
         I = FSTOPC('MSGLVL', 'INFORM', .FALSE.)
      ELSE
         I = FSTOPC('MSGLVL', DEF1(13), .FALSE.)   ! -m
      ENDIF
      IF(DEF1(14) .NE. 'FATALE') THEN
         I = FSTOPC('TOLRNC', DEF1(14), .FALSE.)   ! -t
      ELSE
         I = FSTOPC('TOLRNC', DEF1(12), .FALSE.)   ! -k
      ENDIF
      I = FSTOPL('FASTIO', FASTIO, .FALSE.)
      I = FSTOPL('IMAGE',  .TRUE., .FALSE.)

!     ON EMPECHE LE REMAPPING DU DATATYPE EN MODE IMAGE
!      call c_no_datyp_remap()
!     Verification du remapping faite directement dans fstecr pour le mode image

!     COMPTER LES FICHIERS SOURCES
      DO I=1,120
         IF(DEF1(I+25) .NE. ' ') NFS = NFS+1  ! -s
      ENDDO
      IF(NFS .GT. 0) THEN
         CALL OUVRES( DEF1(26) )    ! ouvrir le premier fichier
         IF(NFSO .EQ. 0) THEN
            PRINT*,' *********************************'
            PRINT*,'***         PAS DE COPIE        ***'
            PRINT*,'***    FICHIER SOURCE INCONNU   ***'
            PRINT*,' *********************************'
            ETAT = 'ABORT'
            GO TO 30
         ENDIF
      ENDIF
   
!     OUVRE LE FICHIER DESTINATION
      IF(DEF1(2) .NE. ' ') I = OUVRED( DEF1(2) )   ! -d 

!     LIRE UN JEU DE DIRECTIVES

      IF( SELEC ) THEN
         I = FNOM(5, DEF1(10), 'SEQ', 0)
         IF(DEF1(10) .EQ. '$IN') THEN   !  (-i ) directives from stdin, prompt for directives
            INTERAC = .TRUE.
            PRINT*,'DIRECTIVES ?'
            PRINT*,'TAPER  END A LA FIN DES DIRECTIVES'
            PRINT*,'TYPE  END  AFTER LAST DIRECTIVE'
         ENDIF
         CALL SELECT                    ! process directives
      ENDIF

!     SI PAS DE DIRECTIVES STDCOPI OU SEQCOPI ALORS
      IF( .NOT.ESAIS ) THEN
         IF( OUVS .AND. OUVD ) THEN
            NP = 1
            CALL STDCOPI( -1,-1,-1,-1 )
         ELSE
            PRINT*,               ' *********************************'
            PRINT*,               '***         PAS DE COPIE        ***'
            IF(.NOT. OUVS) PRINT*,'***    FICHIER SOURCE INCONNU   ***'
            IF(.NOT. OUVD) PRINT*,'*** FICHIER DESTINATION INCONNU ***'
            PRINT*,               ' *********************************'
            ETAT = 'ABORT'
         ENDIF
      ENDIF

!     TOUT EST TERMINE , FERME LES FICHIERS
      CALL FERMES
      CALL FERMED
!     IMPRIMER L'INDICATIF DE FIN DU PGM.
   30 IF( BOX ) THEN
         I = EXFIN('EDITFST', ETAT, 'NON')
      ELSE
         IF(ETAT .EQ. 'ABORT') THEN
            WRITE(6,*)'***   E D I T F S T   A V O R T E   ***'
         ELSE
            WRITE(6,*)'***   E D I T F S T   T E R M I N E   ***'
         ENDIF
      ENDIF
      IF(ETAT .EQ. 'ABORT') CALL QQEXIT(50)  ! get error exit code back to shell
      STOP
      END 
      
      character *128 function product_id_tag()
      product_id_tag='$Id: editfst.F90 6.20 2015-10-07 18:53:41Z armnlib $'
      return
      end
