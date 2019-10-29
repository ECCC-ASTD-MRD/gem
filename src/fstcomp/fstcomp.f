C*DECK FSTCOMP
***PROGRAM FSTCOMP - COMPARAISON DE 2 FICHIERS STANDARDS
*
      PROGRAM  FSTCOMP
      IMPLICIT NONE
*
*AUTEURS
*     V1.0 ORIGINALE (COMPSTD)   M. LEPINE   (MARS 87)
*REVISIONS
* 001 V1.? AJOUT DU MARQUAGE EN CAS D'ERREUR PLUS GRANDE QUE LA LIMITE
*          SPECIFIEE
* 002 V4.0 VERSION UNIX/FST). Y BOURASSA (JAN 91)
* 003 V4.1 POSSIBILITE D'IGNORER LES TYPVAR,NOMVAT,ETIKET,IP1,IP2 ET IP3
*          LORS DE LA RECHERCHE DANS B. Y BOURASSA (JAN 91)
* 004 V4.2 MODIFICATIONS AUX SORTIES. Y BOURASSA (FEV 91)
* 005 V4.3 ACCEPTE LES FICHIERS BINAIRES. Y BOURASSA (AVR 91)
* 006 V4.4 APPEL A CCARD (AUG 91)
* 007 V5.0 APPEL A LOW2UP AORES CCARD, BUG DANS 2IE PARAMS
* 008 V5.1 (JUIL 92)
* 009 V6.0 (M. Lepine, Jan. 93) tenir compte du format IEEE
* 010 V6.1 (V. Lee, Dec. 1997) infinite loop when NI/NJ/NK does not match.
* 011 V6.2 (M. Lepine, avr. 1998) reload avec librmnbeta
* 012 V6.3 (M. Lepine, avr. 1998) reload avec fstd98
* 013 V6.4 (M. Lepine, avr. 1998) reload avec fstd98 et librmn32stack_LR.a
* 014 V6.5 (M. Valin,  avr. 1999) certains calculs faits en real*8
* 015 V6.6 (M. Lepine, Mai  2000) appel a incdatr au lieu de incdat
* 016 V6.7 (M. Lepine, Oct  2000) reload pour datev (fstd89)
* 017 V6.8 (M. Lepine, Jan  2002) mods pour tenir compte des extensions fstd2000
* 018 V6.9 (M. Lepine, Fev  2003) appel a convip pour afficher les niveaux
* 019 V7.0 (M. Lepine, Jan  2004) option reduction32 pour type E64, correction de format
* 019 V7.1 (M. Lepine, Oct  2004) Ajout datatype compresse (>128) et datatype 6
* 020 V7.2 (M. Lepine, Fev  2005) Reload avec librmn_x
* 021 V7.3 (M. Lepine, Mars 2005) Ajout de la fonctionnalite des fichiers remote
* 022 V7.4 (M. Lepine, Fev  2006) Appel a ip1_all 
* 023 V7.5 (M. Lepine, Mai  2006) Reload, bug fix float_packer 
* 024 V7.6 (Y. Chartier, Oct. 2006) Reload pour compresseur point flottant
* 025 V7.7 (M. Lepine, Fev. 2007) Comparaison avec l'erreur du a l'algorithme de compaction 
* 026 V7.8 (M. Lepine, Juin 2007) Reload avec librmn_009
* 027 V7.9 (M. Lepine, Avr. 2010) Correction pour affichage de facteur de correlation negatif
* 028 V8.0 (M. Lepine, Juin 2011) Variables d'exception au convip du IP1
* 026 V8.1 (M. Lepine, Mai  2012) Reload avec librmn_013
* 027 V8.2 (M. Lepine, Nov  2012) En cas d'erreur, retourner un code d'erreur avec qqexit
* 028 V8.3 (M. Lepine, Mars 2014) Utilisation du fichier $ARMNLIB/data/exception_vars
* 029 V8.4 (M. Lepine, Juil 2014) Remettre ARMNLIB a la place de ARMNLIB_DATA
* 030 V8.5 (M. Lepine, Dec  2014) Reload avec librmn_015.1
* 031 V8.6 (M. Lepine, Fev  2015) Reload avec librmn_015.2
* 031 V8.7a (M. Lepine, Mars  2016) Reload avec librmn_Alpha_016
* 032 V8.7a (M. Lepine, Avril  2016) Ajout de l'option pour ignorer la verification de grille
* 033 V8.8 (M. Lepine, Dec  2016) Ajout de statistiques additionnelles
* 034 V8.9 (M. Lepine, Jan  2017) Eviter le traitement des enregistrements '!!' qui contiennent
*                                 un melange d'entiers, reels et caracteres
*
*OBJET(FSTCOMP)
*     ETABLIT DES STATISTIQUES DE COMPARAISON ENTRE DEUX FICHIERS
*     STANDARDS SEQUENTIEL OU ACCES DIRECT, PRODUIT UN RAPPORT.
*
*MODULES
      EXTERNAL FSTLUK, FSTOUV, FSTFRM, FSTVOI, CCARD, FNOM,
     X         FSTINF, INCDATR, RCMP1D, FSTSUI, EXFIN, FSTOPC,  EXDB,
     X         FSTPRM, FSTNBR, ICMP1D, FSTRWD, ABORT, LOW2UP, convip_plus,
     %         fstopl, ip1_all, qqexit, longueur
*
      CHARACTER*1  GRTYPA, GRTYPB
      CHARACTER*2  TYPVAR, TYPVAB
      CHARACTER*4  NOMVAR, NOMVAB
      CHARACTER*8  CLE(22)
      CHARACTER*12 ETIKET, ETIKB
      CHARACTER*12 NOMA, NOMB, NOMC
      CHARACTER*40  NA, NB
      CHARACTER*128 DEF1(22), DEF2(22), NOMD

      LOGICAL TD, TE, TT, AS, AF, BS, BF, VA, VB, DI, LN, ECRIT,
     X        P1, P2, P3, TN ,T, TG, EXCEPTION
      INTEGER KA, KB, N1, N2, X1, X2, LIMITE, I, J, K, L, N, MX,
     X        DATE, DEET, NPAS, NI, NJ, NK, IP1, IP2, IP3,
     X        IG1, IG2, IG3, IG4, IB1, IB2, IB3, IB4, SWA, LNG, DLTF,
     X        UBC, EX1, EX2, EX3, NBITS, DATYPA, IDATE, NBIT2,
     X        DATYPB, FSTRWD, IP1B, IP2B, IP3B,
     X        FSTLUK, FSTOUV, FSTFRM, FSTVOI, FNOM, FSTOPC, EXFIN,
     X        FSTINF, FSTSUI, FSTPRM, FSTNBR, EXDB, fstopl, longueur,
     X        TABLO(0:6,0:6)
      integer ier, kind, ip1_all, PACK_ERR, PACK_ERR2, ind
      integer lvar, iunexpv
      real Level
      character *30 string
      character *128 exception_vars
      character * 512 ARMNLIB_var
      character * 16 RELEASE
      character * 32 SUB_RELEASE
      REAL *8 NHOURS

      COMMON/BUFR/ BUF(1)
      INTEGER      BUF

      REAL, ALLOCATABLE, DIMENSION (:) :: XX1,XX2

*      DATA exception_vars /'^^  >>  !!  '/
      DATA exception_vars /'^^  >>  !!   ^>  '/

      DATA CLE  /'A:', 'B:', 'L',    'AS',  'BS ',  'AF', 'BF',  'LI',
     X           'ND',  'NE',  'D',      'N',   'VA',  'VB',  'NT',
     X           'N1',  'N2',  'N3',  'NN',  'X', 'PACKERR', 'NG' /

      DATA DEF1 /'A', 'B', '$OUT', 'NON', 'NON', 'NON', 'NON', '-7',
     X           'NON', 'NON', 'WARNIN', 'NON', 'NON', 'NON', 'NON',
     X           'NON', 'NON', 'NON', 'NON', 'X', '0', 'NON'/

      DATA DEF2 /'A', 'B', '$OUT', 'SQI', 'SQI', 'FTN', 'FTN', '-7',
     X           'OUI', 'OUI', 'INFORM', 'OUI', 'VA',  'VB',  'OUI',
     X           'OUI', 'OUI', 'OUI', 'OUI', 'X', '1', 'OUI'/

      DATA       N, NOMVAB, TYPVAB, ETIKB, IDATE, IP1B, IP2B, IP3B
     X         / 0, ' ',    ' ',    ' ',   4*-1/
      DATA        ECRIT/ .FALSE. /

*     VALIDE QUAND LA CLE 'X'.NE.'R'
      DATA TABLO/ 2, 1, 2, 3, 2, 1, 1,
     X            1, 1, 3, 3, 3, 1, 1,
     X            2, 3, 2, 3, 2, 3, 3,
     X            3, 3, 3, 3, 3, 3, 3,
     X            2, 3, 2, 3, 2, 3, 3,
     X            1, 1, 3, 3, 3, 1, 1,
     X            1, 1, 3, 3, 3, 1, 1 /
*     0=BINAIRE 1=REEL 2=ENTIER 3=CARACTERE 4=ENTIER SANS SIGNE 5=IEEE
*     NOTE:QUAND LA CLE('X'.EQ.'R') DATA TABLO(0,0)=1

C*IF DEF, NEC64
C      DEF1(8) = '-13'
C      DEF2(8) = '-13'
C*ENDIF
*     EXTRACTION DES CLES DE LA SEQUENCE D'APPEL

      include 'version.inc'

      I = -1
      CALL CCARD(CLE, DEF2, DEF1, 22, I)
      DO 1 I=4,11
         CALL LOW2UP(DEF1(I), DEF1(I))
    1    CONTINUE
      READ(DEF1(8), '(I8)') LIMITE
      READ(DEF1(21),'(I8)') PACK_ERR 
      
      VA = DEF1(13) .NE. 'NON'
      IF(VA .AND. (DEF1(13).NE.'VA')) DEF1(1) = DEF1(13)
      VB = DEF1(14) .NE. 'NON'
      IF(VB .AND. (DEF1(14).NE.'VB')) DEF1(2) = DEF1(14)

      AF = INDEX(DEF1(6),'FTN') .NE. 0
      IF( AF ) THEN
         NOMA = 'STD+SEQ+FTN'
         AS   = .FALSE.
      ELSE
         AS = INDEX(DEF1(4),'SQI') .NE. 0
         IF( AS ) THEN
            NOMA = 'STD+SEQ'
         ELSE
            NOMA = 'STD+RND'
         ENDIF
      ENDIF

      BF = INDEX(DEF1(7),'FTN') .NE. 0
      IF( BF ) THEN
         NOMB = 'STD+SEQ+FTN'
         BS   = .FALSE.
      ELSE
         BS = INDEX(DEF1(5),'SQI') .NE. 0
         IF( BS ) THEN
            NOMB = 'STD+SEQ'
         ELSE
            NOMB = 'STD+RND'
         ENDIF
      ENDIF
      I  = FNOM(6, DEF1(3), 'SEQ+FTN', 0)
      TD = DEF1( 9) .EQ. 'NON'
      TE = DEF1(10) .EQ. 'NON'
      TT = DEF1(15) .EQ. 'NON'
      P1 = DEF1(16) .EQ. 'NON'
      P2 = DEF1(17) .EQ. 'NON'
      P3 = DEF1(18) .EQ. 'NON'
      TN = DEF1(19) .EQ. 'NON'
      TG = DEF1(22) .EQ. 'NON'
      DI = DEF1(11) .EQ. 'INFORM'
      LN = DEF1(12) .EQ. 'OUI'
      IF(DEF1(20) .EQ. 'R') TABLO(0,0) = 1

      IF( LN ) THEN
         WRITE(6,*)'* * *  FSTCOMP '//trim(RELEASE)//'  * * *'
      ELSE
         L = EXDB('FSTCOMP', trim(RELEASE), 'NON')
      ENDIF
      L = FSTOPC('MSGLVL', DEF1(11), .FALSE.)
      ier = fstopl('REDUCTION32',.true.,.false.)

      CALL GETENV('ARMNLIB',ARMNLIB_var)
      lvar = LONGUEUR(ARMNLIB_var)
      IF (lvar .gt. 0) THEN
        iunexpv=0
	ier = fnom(iunexpv,
     %     ARMNLIB_var(1:lvar)//'/data/exception_vars_ok',
     %     'SEQ+FTN+FMT+OLD+R/O',0)
        IF (ier .lt. 0) THEN
          print *,'$ARMNLIB_DATA/exception_vars file not found;'//
     %            ' using internal exception list'
        ELSE
          READ(iunexpv,'(a)') exception_vars
        ENDIF
      ENDIF
*      print *,'Debug exception_vars=',exception_vars
      
*     SI A=RND & B=SEQ CHANGE [A POUR B] & [B POUR A]
      IF((BF.OR.BS) .AND. .NOT.(AF.OR.AS)) THEN
         NOMD    = DEF1(1)
         DEF1(1) = DEF1(2)
         DEF1(2) = NOMD
         NOMC    = NOMA
         NOMA    = NOMB
         NOMB    = NOMC
         T  = BS
         BS = AS
         AS = T
         T  = AF
         AF = BF
         BF = T
      ENDIF

      IF(DI .AND. (DEF1(1).EQ.DEF1(2)))
     X   WRITE(6,*)' ATTENTION ON COMPARE LE FICHIER "A" AVEC LUI-MEME'

*     OUVRE LE FICHIER 1
      NA = DEF1(1)
      L  = FNOM  (1, DEF1(1), NOMA//'+OLD+R/O+REMOTE', 0)
      IF (L .LT. 0) THEN
        CALL QQEXIT(1)
        STOP
      ENDIF
      L  = FSTOUV(1, NOMA)
      N1 = FSTNBR(1)
      IF(N1 .Le. 0) THEN
         IF( DI ) WRITE(6,*)' FICHIER SEQUENTIEL ', NA
         IF( .NOT.(AS.OR.AF) ) THEN
            WRITE(6,*)' FICHIER DECLARE RND ', NA
            GO TO 60
         ENDIF
      ELSE
         IF( DI ) WRITE(6,*)' ',N1,' ENREGISTREMENTS DANS ', NA
         IF(N1 .EQ. 0) GO TO 90
      ENDIF

*     OUVRE LE FICHIER 2
      NB = DEF1(2)
      L  = FNOM  (2, DEF1(2), NOMB//'+OLD+R/O+REMOTE', 0)
      IF (L .LT. 0) THEN
        CALL QQEXIT(2)
        STOP
      ENDIF
      L  = FSTOUV(2, NOMB)
      N2 = FSTNBR(2)
      IF(N2 .Le. 0) THEN
         IF( DI ) WRITE(6,*)' FICHIER SEQUENTIEL ', NB
         IF( .NOT.BS ) THEN
            WRITE(6,*)' FICHIER DECLARE RND ', NB
            GO TO 80
         ENDIF
      ELSE
         IF( DI ) WRITE(6,*)' ',N2,' ENREGISTREMENTS DANS ', NB
         IF(N2 .EQ. 0) GOTO 80
      ENDIF

*     ECRIRE L'ENTETE DE PAGE
      IF (PACK_ERR .eq. 0) THEN      
         WRITE(6,600)
      ELSE
         WRITE(6,700)
      ENDIF
*     RESERVE LA MEMOIRE
      IF(AS .OR. AF) THEN
         L  = FSTRWD(1)
      ENDIF
      KA = FSTINF(1, NI, NJ, NK, -1, ' ', -1, -1, -1,' ', ' ')
      N  = NI*NJ*NK
  10  MX = N
      ALLOCATE(XX1(MX))
      ALLOCATE(XX2(MX))
*      CALL MEMOIRH(BUF, X1, MX)
*      CALL MEMOIRH(BUF, X2, MX)

*     RAMASSE LES IDENTIFICATEURS DU RECORD KA
  20  L = FSTPRM(KA, DATE, DEET, NPAS, NI, NJ, NK, NBITS, DATYPA,
     X           IP1, IP2, IP3, TYPVAR, NOMVAR, ETIKET, GRTYPA,
     X           IG1, IG2, IG3, IG4, SWA, LNG, DLTF, UBC, EX1,
     X           EX2, EX3)
      IF (NOMVAR == '!!') GOTO 35
*     SI LA DATE EST A CONSIDERER
      IF( TD ) THEN
         IF(DATE .EQ. 0) THEN
            IDATE = 0
         ELSE
            IF ((DEET*NPAS) .EQ. 0) THEN
               IDATE = DATE
            ELSE
               nhours = (deet*npas)/3600.
               CALL INCDATR(IDATE, DATE, nhours)
            ENDIF
         ENDIF
      ENDIF
*     SI L'ETIKETTE EST A CONSIDERER
      IF( TE ) ETIKB = ETIKET

*     SI LE TYPEVAR EST A CONSIDERER
      IF( TT ) TYPVAB = TYPVAR

*     SI LE NOMVAR EST A CONSIDERER
      IF( TN ) NOMVAB = NOMVAR

*     SI LE IP1 EST A CONSIDERER
      IF( P1 ) then
        IF (INDEX(exception_vars,nomvar) .ne. 0) then
          IP1B = IP1
          EXCEPTION = .TRUE.
        ELSE
          call convip_plus(ip1,level,kind,-1,string,.false.)
          IP1B = IP1_ALL(level,kind)
          EXCEPTION = .FALSE.
        ENDIF
      ENDIF

*     SI LE IP2 EST A CONSIDERER
      IF( P2 ) IP2B = IP2

*     SI LE IP3 EST A CONSIDERER
      IF( P3 ) IP3B = IP3

      IF(BS .OR. BF) L = FSTRWD(2)
      KB = FSTINF(2, I, J, K, IDATE, ETIKB, IP1B, IP2B, IP3B, TYPVAB,
     X            NOMVAB)
      IF(KB .LT. 0)THEN
         WRITE(6,601) NOMVAB, TYPVAB, IP1B, IP2B, IP3B, IDATE, NB
         GOTO 60
      ENDIF

*     VERIFICATION DES DIMENSIONS DE LA GRILLE SI PRESENT
      IF((NI.NE.I) .OR. (NJ.NE.J) .OR. (NK.NE.K)) THEN
         WRITE(6,603)KA,KB,NOMVAB,I, J, K, NI, NJ, NK
         GOTO 60
      ENDIF

*     RAMASSE LES IDENTIFICATEURS DU RECORD KB
      L = FSTPRM(KB, DATE, DEET, NPAS, I, J, K, NBIT2, DATYPB,
     X           IP1, IP2, IP3, TYPVAR, NOMVAR, ETIKET, GRTYPB, IB1,
     X           IB2, IB3, IB4, SWA, LNG, DLTF, UBC, EX1, EX2, EX3)

*     VERIFICATION DES PARAMETRES DE LA GRILLE
      IF (TG) THEN
        IF((GRTYPA.NE. GRTYPB) .OR. (IG1.NE.IB1) .OR. (IG2.NE.IB2) .OR.
     X     (IG3.NE.IB3) .OR. (IG4.NE.IB4)) THEN
           WRITE(6,602) NA, GRTYPA, IG1, IG2, IG3, IG4,
     X                  NB, GRTYPB, IB1, IB2, IB3, IB4
           GOTO 60
        ENDIF
      ENDIF

      IF(NBIT2 .NE. NBITS) PRINT*,'NBITSA=',NBITS,' NBITSB=',NBIT2
      IF ( DI ) PRINT*,'COMPARE DATYPA=',DATYPA,'  @  DATYPB=',DATYPB

*     TOUT EST OK LIT ET COMPARE
      L = FSTLUK(XX1, KA, NI, NJ, NK)
*      L = FSTLUK(BUF(X1), KA, NI, NJ, NK)
*      L = FSTLUK(BUF(X2), KB, NI, NJ, NK)
      L = FSTLUK(XX2, KB, NI, NJ, NK)
      IF ((mod(DATYPA,128) .ne. 1) .and. (mod(DATYPA,128) .ne. 6)) THEN
        PACK_ERR2 = 0
      ELSE
        PACK_ERR2 = PACK_ERR
      ENDIF
      IF ((mod(DATYPA,128) .gt. 6).or.(mod(DATYPB,128) .gt. 6)) goto 30
      GO TO (40, 50, 30) TABLO(mod(DATYPA,128),mod(DATYPB,128))
  30  WRITE(6,*)' *  PAS DE COMPARAISON  *  DATYPA=',DATYPA,
     X                                    ' DATYPB=',DATYPB
      GO TO 60
  35  WRITE(6,*)' **   SKIPPING RECORD "!!", CAN''T COMPARE  **' 
      GO TO 60
  40  CALL RCMP1D(XX1, XX2, N, 6, KA, KB, NOMVAR, ETIKB,
     X            IP1, IP2, IP3, LIMITE, MIN(NBITS,NBIT2), PACK_ERR2,
     X            EXCEPTION, DATE, TYPVAR, NI, NJ, NK)
      GO TO 60
  50  CALL ICMP1D(XX1, XX2, N, 6, KA, KB, NOMVAR, ETIKB,
     X            IP1, IP2, IP3, EXCEPTION)
  60  KA = FSTSUI(1, NI, NJ, NK)
      IF(KA .GE. 0) THEN
         N = NI*NJ*NK
         IF(N .LE. MX) GO TO 20
         DEALLOCATE(XX1)
         DEALLOCATE(XX2)
*         CALL MEMOIRH(BUF, X1, 0)
*         CALL MEMOIRH(BUF, X2, 0)
         GOTO 10
      ENDIF

   70 IF(VA.AND..NOT.AF .OR. VB.AND..NOT.BF) THEN
         IF(.NOT.DI) L = FSTOPC('MSGLVL', 'INFORM', .FALSE.)
         IF( VA ) THEN
            IF( AS ) L = FSTRWD(1)
            WRITE(6,*)' DN=', NA
            L = FSTVOI(1, DEF1(1))
         ENDIF
         IF( VB ) THEN
            IF( BS ) L = FSTRWD(2)
            WRITE(6,*)' DN=', NB
            L = FSTVOI(2, DEF1(2))
         ENDIF
      ENDIF
   80 ier = FSTFRM(2)
   90 ier = FSTFRM(1)

      IF( LN ) THEN
         WRITE(6,*)'* * *  FSTCOMP TERMINE  * * *'
      ELSE
         L = EXFIN('FSTCOMP', 'NORMAL', 'NON')
      ENDIF

  600 FORMAT('  NOM    ETIKET           IP1',
     X       '            IP2       IP3  E-REL-MAX',
     X       '  E-REL-MOY    VAR-A      C-COR        MOY-A',
     X       '        BIAIS      E-MAX      E-MOY')


  601 FORMAT(' PAS TROUVE ',A4,' ',A2,' IP123=', 3I8, I10,' DANS ',A40)
  602 FORMAT(' ',A40,' GRTYP IG1@4=', A1,1X, 4I6,/
     X       ' ',A40,' GRTYP IG1@4=', A1,1X, 4I6)
  603 FORMAT(2I6,A4,' -LES DIMENSIONS TROUVEES SONT',3I5,
     %       ' CHERCHE',3I5)
     
  700 FORMAT('  NOM    ETIKET           IP1',
     X       '            IP2       IP3  E-REL-MAX',
     X       '  E-REL-MOY    VAR-A      C-COR        MOY-A',
     X       '        BIAIS      E-MAX      E-MOY     TOLERANCE')

     
      STOP
      END
C*DECK RCMP1D
***S/P RCMP1D  COMPARAISON DE DEUX CHAMPS REELS DE UNE DIMENSION
*
      SUBROUTINE RCMP1D(A, B, N, IUN, NUMA, NUMB, NOMVAR, ETIKET, IP1,
     X                  IP2, IP3, LIMITE, NBITS, PACK_ERR, EXCEPTION,
     x                  DATE, TYPVAR, NI, NJ, NK)

      IMPLICIT NONE
      INTEGER  N, IUN, LIMITE, NUMA, NUMB, IP1, IP2, IP3, NBITS
      INTEGER  PACK_ERR
      INTEGER  NI, NJ, NK, DATE
      LOGICAL EXCEPTION
      CHARACTER*12 ETIKET
      CHARACTER*4 NOMVAR
      CHARACTER*2 TYPVAR
      REAL     A(N), B(N), MAXABS, SUMABS, ERRABS
*
*AUTEURS  VERSION ORIGINALE (REALCMP)  M.VALIN DRPN 1987
*         VERSION (RCMP1D)  Y.BOURASSA DRPN JAN 1990
*         Ajout de l'argument exception - M. Lepine Mars 2014
*
*ARGUMENTS
* ENTRE  A,B     CHAMPS REELS A COMPARER
*   "    N       DIMENSION DE A ET B
*   "    IUN     UNITE SUR LAQUELLE ON ECRIT LES RESULTATS
*   "    NUMA    NUMERO DE L'ENREGISTREMENT SUR FICHIER A
*   "    NUMB    NUMERO DE L'ENREGISTREMENT SUR FICHIER B
*   "    NOMVAR  NOM DE VARIABLE
*   "    ETIKET  ETIKET DU CHAMP A COMPARER
*   "    IP1@3   PARAMETRES DE SELECTION
*   "    LIMITE  ERREUR MAXIMUM TOLOREE
*   "    NBITS   NOMBRE DE BITS UTILISE POUR LE PACKING
*   "    PACK_ERR NOMBRE D'UNITE D'ERREUR DU A L'ALGORITHME DE PACKING
*                 A UTILISER POUR DETERMINER SI "A" COMPARE A "B"
**
      INTEGER   I, kind, irange
      CHARACTER*15 Level
      REAL      rlevel
      REAL*8    SA, SB, SA2, SB2, ERR, DERR, ERRMAX, ABAR, BBAR,
     X          AA, BB, FN, ERRLIM, VARA, VARB, SAB
      REAL MIN_A, MAX_A, MIN_B, MAX_B, RANGE_A, RANGE_B, DEUX_EXP_NB
      REAL ratio_max, ratio
      REAL ERR_UNIT
      integer nbdiff
      EXTERNAL statfldx

      nbdiff = 0
      ERRLIM = 10.**LIMITE
      DEUX_EXP_NB = 2.0 ** NBITS
      SA     = 0.
      SB     = 0.
      SAB    = 0.
      SA2    = 0.
      SB2    = 0.
      ERRMAX = 0.
      ERR    = 0.
      SUMABS = 0.
      MAXABS = 0.
      ratio_max = 0.
      MIN_A = A(1)
      MAX_A = A(1)
      MIN_B = B(1)
      MAX_B = B(1)
      DO 10 I=1,N
         AA     = A(I)
         BB     = B(I)
         MIN_A = MIN(MIN_A,A(I))
         MAX_A = MAX(MAX_A,A(I))
         MIN_B = MIN(MIN_B,B(I))
         MAX_B = MAX(MAX_B,B(I))
         SA     = SA+AA
         SB     = SB+BB
         IF(AA .NE. BB) THEN
            if (aa .ne. 0.) ratio = (max(aa,bb) - min(aa,bb)) / aa * 100
            if (ratio > ratio_max) ratio_max = ratio
!            write(6,888) 'Debug difference A vs B au point I=',i,
!     %       'AA, AA-BB=',aa,aa-bb
            nbdiff = nbdiff +1
 888        format(a,i8,2x,a,e14.7,2x,e14.7)
            ERRABS = ABS(AA-BB)
            SUMABS = SUMABS+ERRABS
            MAXABS = MAX(ERRABS,MAXABS)
            derr=0.0
            IF(AA .NE. 0.) THEN
               DERR = ABS(1.-BB/AA)
            ELSEIF(BB .NE. 0.)THEN
               DERR = ABS(1.-AA/BB)
            ENDIF
            ERRMAX = MAX(ERRMAX,DERR)
            ERR    = ERR+DERR
         ENDIF
   10    CONTINUE
!      print *,'Debug nbdiff = ',nbdiff,' sur un total de ',n

      RANGE_A = MAX_A - MIN_A
      RANGE_B = MAX_B - MIN_B
      irange = TRANSFER(RANGE_A,1)
      irange = ISHFT(ISHFT(irange,-23) +1,23)
      RANGE_A = TRANSFER(irange,1.0)
      FN   = FLOAT(N)
      ERR  = ERR/FN
      ABAR = SA/FN
      BBAR = SB/FN
      DO 20 I=1,N
         AA  = A(I)-ABAR
         BB  = B(I)-BBAR
         SAB = SAB+AA*BB
         SA2 = SA2+AA*AA
   20    SB2 = SB2+BB*BB
      SUMABS = SUMABS/FN
      VARA   = SA2/FN
      VARB   = SB2/FN
!      print *,'Debug+ sab avant=',sab
      IF(SA2*SB2 .NE. 0.) THEN
         SAB    = SAB/SQRT(SA2*SB2)
!      print *,'Debug+ sa2,sb2,sqtr,sab=',sa2,sb2,SQRT(SA2*SB2),sab
      ELSEIF(SA2.EQ.0. .AND. SB2.EQ.0.) THEN
         SAB = 1.0
      ELSEIF(SA2 .EQ. 0.) THEN
         SAB = SQRT(VARB)
      ELSE
         SAB = SQRT(VARA)
      ENDIF
!      print *,'Debug+ sa2,sb2,vara,varb,sab=',sa2,sb2,vara,varb,sab
      IF (EXCEPTION) THEN
        WRITE(level,'(i5)') ip1
      ELSE
        CALL convip_plus(ip1,rlevel,kind,-1,level,.true.)
      ENDIF
      ERR_UNIT = RANGE_A / DEUX_EXP_NB
      
      IF ((ERRMAX .LE. ERRLIM) .and. (PACK_ERR .eq.0)) THEN
         WRITE(IUN,600) NOMVAR, ETIKET, level, IP2, IP3,
     X                  ERRMAX, ERR, VARA, SAB, ABAR, BBAR-ABAR,
     X                  MAXABS, SUMABS
      ELSE IF (PACK_ERR .gt. 0) THEN
         IF (MAXABS .le. (PACK_ERR*ERR_UNIT*1.001)) THEN
           WRITE(IUN,602) NOMVAR, ETIKET, level, IP2, IP3,
     X                  ERRMAX, ERR, VARA, SAB, ABAR, BBAR-ABAR,
     X                  MAXABS, SUMABS, PACK_ERR*ERR_UNIT
         ELSE       
            WRITE(IUN,603) NOMVAR, ETIKET, level, IP2, IP3,
     X                  ERRMAX, ERR, VARA, SAB, ABAR, BBAR-ABAR,
     X                  MAXABS, SUMABS, PACK_ERR*ERR_UNIT
         ENDIF
      ELSE
         WRITE(IUN,601) NOMVAR, ETIKET, level, IP2, IP3,
     X                  ERRMAX, ERR, VARA, SAB, ABAR, BBAR-ABAR,
     X                  MAXABS, SUMABS
      ENDIF
 
      if ((nbdiff .ne. 0) .and. (PACK_ERR .gt. 0)) then
        write(6,900) '  <Difference> Number of elements differing is',
     %    nbdiff, ' out of ',n,
     %    ' elements, Ratio = ',(float(nbdiff) / float(n)) * 100.0,'%'
 900    format(a,i7,a,i8,a,f10.4,a)
        call statfldx(nomvar,typvar,ip1,ip2,ip3,date,etiket,A,ni,nj,nk)
        call statfldx(nomvar,typvar,ip1,ip2,ip3,date,etiket,B,ni,nj,nk)
        write(6,*) ' '
!        write(6,901) 'Err_Max=',ratio_max
! 901    format(a,e8.2,'%')
!        write(6,889) 'Debug MIN_A, MAX_A, MIN_B - MIN_A = ',MIN_A,
!     %                    MAX_A, MIN_B - MIN_A, MAX_B - MAX_A
! 889    format(a,2x,e14.7,2x,e14.7,2x,e14.7,2x,e14.7)
      endif

*  600 FORMAT('  CLEA CLEB NOM  ETIKET    IP1 IP2 IP3  E-REL-MAX',
*     X       '  E-REL-MOY   VAR-A       C-COR         MOY-A',
*     X       '         BIAIS      E-MAX      E-MOY')
  600 FORMAT(' ', '  ', A4, '  ', A12, a15, 2I9, 4(1X,1PE10.4),
     X       2(1X,1PE12.4), 2(1X,1PE10.4) )
  601 FORMAT(' ', ' <', A4, '> ', A12, a15, 2I9, 3(1X,1PE10.4),
     X       3(1X,1PE12.4), 2(1X,1PE10.4) )
  602 FORMAT(' ', '  ', A4, '  ', A12, a15, 2I9, 4(1X,1PE10.4),
     X       2(1X,1PE12.4), 2(1X,1PE10.4), 2X, 1PE10.4 )
  603 FORMAT(/' ', ' <', A4, '> ', A12, a15, 2I9, 3(1X,1PE10.4),
     X       3(1X,1PE12.4), 2(1X,1PE10.4), 2X, 1PE10.4 )

      RETURN
      END
C*DECK ICMP1D
***S/P ICMP1D DIFFERENCES ENTRE DEUX CHAMPS ENTIERS DE UNE DIMENSION
*
      SUBROUTINE ICMP1D(A, B, N, IUN, NUMA, NUMB, NOMVAR, ETIKET, IP1,
     X                  IP2, IP3, EXCEPTION)

      IMPLICIT    NONE
      INTEGER     N, A(N), B(N), IUN, NUMA, NUMB, IP1, IP2, IP3
      CHARACTER*12 ETIKET
      CHARACTER*4 NOMVAR
*
*AUTEURS  VERSION ORIGINALE (INTEMP)  M.VALIN DRPN 1987
*         VERSION (ICMP1D)  Y.BOURASSA DRPN JAN 1990
*
*ARGUMENTS
* ENTRE  A,B     CHAMPS ENTIERS A COMPARER
*   "    N       DIMENSION DE A ET B
*   "    IUN     UNITE SUR LAQUELLE ON ECRIT LES RESULTATS
*   "    NUMA    NUMERO DE L'ENREGISTREMENT SUR FICHIER A
*   "    NUMB    NUMERO DE L'ENREGISTREMENT SUR FICHIER B
*   "    NOMVAR  NOM DE VARIABLE
*   "    ETIKET  ETIKET DE L'ENREGISTREMENT
*   "    IP1@3   PARAMETRES DE SELECTION
**
      INTEGER I, J, K, MD, NC, kind
      CHARACTER*15 Level
      REAL      rlevel
      LOGICAL EXCEPTION

      MD  = 0
      NC  = 0

      DO 10 I=1,N
         if(A(I) .NE. B(I)) THEN
            NC = NC+1
            K  = ABS(A(I)-B(I))
            IF(K .GT. MD) THEN
               J  = N
               MD = K
            ENDIF
         ENDIF
10    CONTINUE

      IF (EXCEPTION) THEN
        WRITE(level,'(i5)') ip1
      ELSE
        CALL convip_plus(ip1,rlevel,kind,-1,level,.true.)
      ENDIF
      IF (MD .EQ. 0) THEN
         WRITE(IUN,600) NOMVAR, ETIKET, level, IP2, IP3
      ELSE
         WRITE(IUN,601) NOMVAR, ETIKET, level, IP2, IP3,
     X                  NC, MD, J
      ENDIF

 600  FORMAT(' ',  '  ', A4, '  ', A12, a15, 2I5,' SONT EGAUX')
 601  FORMAT(' ',  ' <', A4, '> ', A12, a15, 2I5,' ONT',I6,' POINTS ',
     X       'INEGAUX, L''ERREUR MAX.=',I10,'  AU POINT',I6)

      RETURN
      END
      
      character *128 function product_id_tag()
      product_id_tag='$Id: fstcomp.f 178 2015-02-26 14:44:36Z armnlib $'
      return
      end
