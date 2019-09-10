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
!** S/P COPYSTX COPIE UN FICHIER STANDARD EN TOUT OU EN PARTIE.
      SUBROUTINE COPYSTX
      use ISO_C_BINDING
!      use convert_ip123
!      use format_ip123_kind
       use configuration
      IMPLICIT NONE 
       include 'convert_ip123.inc'
       include 'excdes.inc'
  
!AUTEURS
!         - C. THIBEAULT  FEV 83
!         - Y. BOURASSA   FEV 86
!           "     "       NOV 89 INGORER LISTE OU GANGE DE RECORDS 
!           "     "              CORRIGE BUG QUAND EXPRESS = .TRUE.
!           "     "       OCT 90 ACCEPTE FICHIERS STD89. 
!           "     "       NOV 90 SI MIX DE DESIRE ET EXCLURE, EXCLURE
!                                LES ENREGISTREMENTS DESIRES SEULEMENT.
!           "     "       JUL 91 VERSION ZAPPER
!           "     "       NOV 91 TESTING DES PARAMETRES DE SELECTION DANS 
!                                UN ORDRE DIFFERENT
!           "     "       JAN 92 BUG PREMIERE ALLOCATION MEMOIRE
!Revision 009   M. Lepine - mars 98 - extensions pour fstd98
!Revision 010   M. Lepine - Oct  98 - Ajout du test de nrecmin
!Revision 010   M. Lepine - Dec  01 - Termine dans le cas d'une erreur d'ecriture
!Revision 011   M. Lepine - Mai  02 - code de traitement des IP1 en valeurs reelle
!Revision 012   M. Lepine - Oct  02 - save sur buftemp
!Revision 013   M. Lepine - Nov  05 - remplacement de fstabt par qqexit
!Revision 014   M. Valin  - Fev  14 - mode dryrun, nouveau traitement IP1/2/3
!                         - Mai  14 - remplacement de la logique de selection
!Revision 015   M. Lepine - Sep  16 - elimination des espaces blancs a l'impression
!Revision 016   M. Lepine - Sep  16 - dump_request_table en mode debug seulement
!
!LANGAGE  - FTN77
      integer, dimension ( : ), pointer, save :: buftemp => NULL() ! CHAMP pour lire ce qu'on doit copier, (LONGUEUR INITIALE = 0)

      EXTERNAL     FSTPRM, FSTSUI, FSTECR, FSTINF, FSTWEO
      EXTERNAL     CRITSUP, DESIRE, FSTLUK, FSTFRM, FSTEOF
      EXTERNAL     QQEXIT
      INTEGER      FSTINF, FSTEOF, IG1, VALID, NI, I, IREC, DATE, SWA
      INTEGER      FSTPRM, FSTFRM, IG2, XTRA2, NJ, J, DEET,       LNG
      INTEGER      FSTECR, FSTSUI, IG3, XTRA3, NK, K, DLFT, DTYP, UBC
      INTEGER      FSTLUK, FSTWEO, IG4, NBITS, NPAS
      integer      IP(4)
      integer      IP1,IP2,IP3 
      character(len=4) :: nomvar
      logical :: can_translate
      real :: p1, p2, p3
      integer :: kind1, kind2, kind3, matches
      character (len=2) :: strkind1, strkind2, strkind3
      integer :: nrecords
      integer :: date_1, date_2
  
      LOGICAL      FIRSTP, BONNE, OK
      interface ! can IP1/2/3 for variable 'name' be translated to value/kind ?
        function fstcantranslate(name) result (yesno) BIND(C,name='FstCanTranslateName')
        use ISO_C_BINDING
        integer(C_INT) :: yesno
        character(C_CHAR), dimension(*), intent(IN) :: name
        end function
      end interface
      integer, save :: NM=0
!      DATA         NM, IST / 0, 0/

      OK     = .NOT.FIXD .AND. .NOT.DM1 
!     OK     = .TRUE. UNE DATE DANS LES DESIRES ET ENREGISTREMENTS
!                     PAS NECESSAIREMENT VALIDES EN MEME TEMP
!     DM1    = .TRUE. PAS DE DATES DANS LES DESIRES
!     FIXD   = .TRUE. LES ENREGISTREMENTS DU FICHIER SOURCE SONT TOUS 
!                     VALIDES EN MEME TEMPS.
!     BONNE  = .TRUE. SI LA DATE DU PREMIER ENREGISTREMENT ACCEPTABLE 
!     DONC IF(FIXD .AND. .NOT.BONNE) INUTILE DE CHERCHER PLUS LOIN

      IF( DEBUG ) call Dump_Request_table()
      nrecords = 0

   10 BONNE  = .FALSE.
      FIRSTP = .TRUE.
      IF( DEBUG ) WRITE(6,*)' FIRSTP=',FIRSTP,'  OK=',OK,' BONNE=',BONNE
!     OBTENIR LA CLE DU PROCHAIN ENREGISTREMENT QUI NOUS INTERESSE
      IREC   = FSTINF(SOURCES(1), NI, NJ, NK, -1, ' ', -1, -1, -1,' ', ' ')
!
      do while(irec >= 0)
!
      I = FSTPRM(IREC, DATE, DEET, NPAS, NI, NJ, NK, NBITS, DTYP,    &
                 IP1, IP2, IP3, TYP, NOM, ETI, GTY, IG1, IG2,        &
                 IG3, IG4, SWA, LNG, DLFT, UBC, VALID, XTRA2, XTRA3)
!
      IF(NBITS.GT.48 .AND. DTYP.EQ.1) THEN
         WRITE(6,*)'IMPOSSIBLE DE COPIER ENREGISTREMENT NO.',IREC,' NBITS =',NBITS 
         GO TO 140
      ENDIF
      nrecords = nrecords + 1
      write(nomvar,'(A4)')NOM
      can_translate = (0 /= fstcantranslate(nomvar))   ! est-ce qu'on peut traduire ip1/2/3 pour cette variable ?

      IF(FIRSTP .OR. OK) THEN 
         IP(4) = VALID   ! date valid
         IF(DEBUG .AND. FIRSTP) WRITE(6,*)'ENRG. #1 DATE ORIG = ',DATE,' VALID =',IP(4)
         FIRSTP = .FALSE.
         IF( DEBUG ) WRITE(6,*)' FIRSTP=',FIRSTP,' OK=',OK,  ' BONNE=',BONNE
      ENDIF

!     SI ON DEMANDE TOUT LE FICHIER.
      IF( XPRES ) THEN
         BONNE  = .TRUE.
      ENDIF
!
!     CONTROLE DE LA MEMOIRE TAMPON AVANT LECTURE
!     si on se rend ici, c'est que l'enregistrement est a copier
!
      if(dryrun) then
        I = 0    ! fake status of succesful read if dry run
      else
        IF(LNG .GT. NM) THEN  ! make sure buffer is large enough to receive data
          IF (associated(buftemp)) deallocate(buftemp)
          NM = LNG
          allocate(BUFtemp(NM+10))
        ENDIF
        I = FSTLUK(BUFtemp, IREC, NI, NJ, NK)  ! on lit l'enregistrement
      endif
!
!     ==================   logique pour directive ZAP  ==================
!
      IF( ZA ) THEN    ! on "zappe ?"
         IF(Z1 .NE. -1)  IP1 = Z1  ! zap ip1
         IF(Z2 .NE. -1)  IP2 = Z2  ! zap ip2
         IF(Z3 .NE. -1)  IP3 = Z3  ! zap ip3
         IF(ZD .NE. -1)  DATE  = ZD  ! zap origin date
         IF(ZT .NE. '??') THEN   ! zap type
            IF(ZT(1:1) .NE. '?') TYP(1:1) = ZT(1:1)
            IF(ZT(2:2) .NE. '?') TYP(2:2) = ZT(2:2)
         ENDIF
         IF(ZN .NE. '??') THEN   ! zap name
            IF(ZN(1:1) .NE. '?') NOM(1:1) = ZN(1:1)
            IF(ZN(2:2) .NE. '?') NOM(2:2) = ZN(2:2)
            IF(ZN(3:3) .NE. '?') NOM(3:3) = ZN(3:3)
            IF(ZN(4:4) .NE. '?') NOM(4:4) = ZN(4:4)
         ENDIF
         IF(ZE .NE. '????????????') THEN   ! zap ETIKET
            DO 130 I=1,12
               IF(ZE(I:I) .NE. '?') ETI(I:I) = ZE(I:I)
  130          CONTINUE
         ENDIF
      ENDIF
!
      if(dryrun) then  ! dry run debug mode, tell what we would be copying
!
        i = decode_ip(p1,kind1,p2,kind2,p3,kind3,ip1,ip2,ip3)
        strkind1=kind_to_string(kind1)
        strkind2=kind_to_string(kind2)
        strkind3=kind_to_string(kind3)
        call newdate(date,date_1,date_2,-3)   ! translate date time stamp to printable format
        if(can_translate)then
          write(6,667)'DRYRUN-select: ',date_1,date_2,TYP,NOM,ETI,P1,strkind1,P2,strkind2,P3,strkind3,DATE,GTY,IG1,IG2,IG3,IG4
        else
          write(6,666)'DRYRUN-select: ',date_1,date_2,TYP,NOM,ETI,IP1,IP2,IP3,DATE,GTY,IG1,IG2,IG3,IG4
        endif
666     format(A,2(I8.8,1X),A3,A5,A13,3I14       ,I12,A2,4I10)
667     format(A,2(I8.8,1X),A3,A5,A13,3(G12.5,A2),I12,A2,4I10)
!
      else   !  real mode, write into the output file
!
        I = FSTECR(BUFtemp, BUFtemp, -NBITS, 3, DATE, DEET, NPAS, NI, NJ, NK,  &
                   IP1, IP2, IP3, TYP, NOM, ETI, GTY, IG1, IG2, IG3, IG4, DTYP, ECR)
        if (i .lt. 0) then
          write(6,*) 'ERROR: (copystx) write error, ABORTING'
          call qqexit(55)
        endif
      endif
!
      COPIES = COPIES + 1
      LIMITE = LIMITE - 1
      IF(LIMITE .EQ. 0) GO TO 180            ! nombre maximum d'enregistrements a copier atteint
140   IREC = FSTSUI(SOURCES(1), NI, NJ, NK)  ! prochain enregistrement
!
      END DO  ! (while(irec >= 0)
!
  160 IF(SSEQ) THEN                          ! le fichier source est sequentiel
         LEOF = FSTEOF(SOURCES(1))
         print*,'apres fsteof leof=',leof
         IF(DIAG .OR. DEBUG) WRITE(6,*)'RENCONTRE UN EOF',LEOF, ' DANS ', SOURCES(1),'...'
         IF(LEOF.GT.15 .OR. LEOF.LT.1) THEN
            WRITE(6,*) LEOF," N'EST PAS ACCEPTABLE COMME EOF LOGIQUE"
            call qqexit(30)
         ENDIF
         IF(DSEQ .AND. CEOF.NE.0) THEN
            K = CEOF
            IF(CEOF .LT. 0) K = LEOF
            IF(K .LT. 15) THEN
               I = FSTWEO(3, K)
               IF(I.EQ.0 .AND. (DIAG .OR. DEBUG)) THEN
                  WRITE(6,*)'EOF LOGIQUE ',K,' AJOUTEE AU FICHIER',TRIM(ND)
               ELSEIF(I .NE. 0) THEN
                  WRITE(6,*)'IMPOSSIBLE D''ECRIRE UNE MARQUE DE ', 'NIVEAU ',K,' DANS ', TRIM(ND)
                  call qqexit(31)
               ENDIF
            ENDIF
         ENDIF
!        DEVONS-NOUS CONTINUER PASSE LE EOF RENCONTRE DANS SOURCE?
         IF(LEOF .LT. MEOF) GO TO 10
      ENDIF
  
!     DOIT-ON ECRIRE UN EOF AVANT DE FERMER?
      IF(DSEQ .AND. EOF.GT.0) THEN           ! le fichier destination est sequentiel
         I = FSTWEO(3, EOF)
         IF(I.EQ.0 .AND. (DIAG .OR. DEBUG)) THEN
            WRITE(6,*)' MARQUE DE NIVEAU',K,' ECRITE DANS ', TRIM(ND)
         ELSEIF(I .NE. 0) THEN
            WRITE(6,*)' IMPOSSIBLE D''ECRIRE UNE MARQUE DE NIVEAU', K,' DANS ', TRIM(ND)
            call qqexit(32)
         ENDIF
      ENDIF
  
  180 WRITE(6,*) COPIES,' ENREGISTREMENT(S) COPIES DANS ', TRIM(ND)
!      WRITE(6,*) nrecords,' ENREGISTREMENT(S) LUS DANS ', NS

      IF (COPIES .LT. NRECMIN) THEN
         WRITE(6,*) ' NOMBRE MINIMAL D ENREGISTREMENT INSATISFAIT'
         WRITE(6,*) ' NRECMIN=',NRECMIN,' NOMBRE TROUVE = ',COPIES
         CALL QQEXIT(12)
      ENDIF

      RETURN
  
      END 
