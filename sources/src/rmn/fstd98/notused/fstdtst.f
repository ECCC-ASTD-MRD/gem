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




      PROGRAM TESTSUN



      IMPLICIT LOGICAL(A-Z)

*
*AUTEUR   P. CADIEUX  AVR 1983
*
*REVISION  002  P. SARRAZIN MAI 89 FICHIER STANDARD  89
*
*LANGAGE  FORTRAN
*
*OBJET(STDTEST)
*         CE PROGRAMME CONSTITUE UN TEST ASSEZ COMPLET DE L'ENSEMBLE
*         DU PROGICIEL DES FICHIERS STANDARDS.  IL A ETE DEVELOPPE
*         POUR LA VERIFICATION DE LA VERSION 89 (NOS + CRAY + MFR + SUN).
*
*MODULES
*
**
*
*REVISION
*         1.0.1 -  AJOUTER 2 TABLES POUR FICHIER RANDOM/SEQ STANDARD
*                   MODIFICATIONS 22 JANV 89
*

*
*  MACRO ICHAMP(MOT,FIN,LG)   EXTRAIT UN CHAMP D'UN MOT
*
*   MOT   MOT QUI CONTIENT LE CHAMP
*   FIN   NUMERO DU DERNIER BIT (A DROITE) DU CHAMP, EN NUMEROTATION
*         GAUCHE > DROITE (LE BIT 0 EST A DROITE DU MOT).
*   LG    LONGUEUR, EN BITS, DU CHAMP
*

*
*  MACRO IUNPAK(BASE, BITPOS, LG)  OBTENIR UN CHAMP D'UN TABLEAU
*
*  BASE    TABLEAU CONTENANT LE CHAMP A EXTRAIRE
*  BITPOS  POSITION DU BIT DE DROITE DU CHAMP A EXTRAIRE
*          LE BIT 0 EST LE BIT SIGNE DU PREMIER MOT DU TABLEAU.
*  LG      EST LE NOMBRE DE BITS QU'OCCUPE LE CHAMP. (MAX 32 BITS)
*

*
*
*  MACRO GETBIT(BASE, BITPOS, LG)  OBTENIR UN CHAMP D'UN TABLEAU
*
*  BASE    TABLEAU CONTENANT LE CHAMP A EXTRAIRE
*  BITPOS  POSITION DU BIT DE DROITE DU CHAMP A EXTRAIRE
*          LE BIT 0 EST LE BIT SIGNE DU PREMIER MOT DU TABLEAU.
*  LG      EST LE NOMBRE DE BITS QU'OCCUPE LE CHAMP. (MAX 32 BITS)
*

*
*
*  MACRO INSERT(TABL,KWA,BITPOS,LONG)  INSERER UN CHAMP DANS UN TABLEAU
*
*  TABL    TABLEAU QUI CONTIENDRA LE CHAMP APRES INSERTION
*  KWA     MOT QUI CONTIENT LE CHAMP A INSERER JUSTIFIE A DROITE
*  BITPOS  POSITION DU DERNIER BIT (A DROITE) DU CHAMP A INSERER
*          LE BIT 0 EST LE BIT SIGNE DU PREMIER MOT DU TABLEAU
*  LONG    LONGUEUR EN BIT DU CHAMP A INSERER (PAS PLUS DE 32 BITS)
*

*
*  MACRO PUTBIT(TABL,KWA,BITPOS,LONG)  INSERER UN CHAMP DANS UN TABLEAU
*
*  TABL    TABLEAU QUI CONTIENDRA LE CHAMP APRES INSERTION
*  KWA     MOT QUI CONTIENT LE CHAMP A INSERER JUSTIFIE A DROITE
*  BITPOS  POSITION DU DERNIER BIT (A DROITE) DU CHAMP A INSERER
*          LE BIT 0 EST LE BIT SIGNE DU PREMIER MOT DU TABLEAU
*  LONG    LONGUEUR EN BIT DU CHAMP A INSERER (PAS PLUS DE 32 BITS)
*

*
*  MACRO CLRBIT(TABL,BITPOS,LONG)  METTRE A ZERO UN CHAMP DANS UN TABLEAU
*
*  TABL    TABLEAU
*  BITPOS  POSITION DU DERNIER BIT (A DROITE) DU CHAMP A NETTOYER
*          LE BIT 0 EST LE BIT SIGNE DU PREMIER MOT DU TABLEAU
*  LONG    LONGUEUR EN BIT DU CHAMP A NETTOYER (PAS PLUS DE 32 BITS)
*

*
*  MACRO PUTBITC(TABL,KWA,BITPOS,LONG)  INSERER UN CHAMP DANS UN TABLEAU
*                                       AVEC NETTOYAGE PRELIMINAIRE
*  TABL    TABLEAU QUI CONTIENDRA LE CHAMP APRES INSERTION
*  KWA     MOT QUI CONTIENT LE CHAMP A INSERER JUSTIFIE A DROITE
*  BITPOS  POSITION DU DERNIER BIT (A DROITE) DU CHAMP A INSERER
*          LE BIT 0 EST LE BIT SIGNE DU PREMIER MOT DU TABLEAU
*  LONG    LONGUEUR EN BIT DU CHAMP A INSERER (PAS PLUS DE 32 BITS)
*

*----------------------------------------------------------------------
*

*----------------------------------------------------------------------

*----------------------------------------------------------------------

*
* POSITIONS DES CHAMPS DE L'ENTETE RANDOM
* UTILISABLES AVEC LES MACROS IUNPAK ET INSERT
*

*
*
*  ETIQET  EST L'ETIQUETTE AFFECTEE AU FORMAT
*  LISVAR  EST LA LISTE DE VARIABLES A IMPRIMER AVEC LE MESSAGE D'ERREUR
*  FORMAT  EST LE FORMAT D'IMPRESSION A UTILISER;  CE PARAMETRE CONTIENT
*          LE MESSAGE D'ERREUR ET LES SPECIFICATIONS D'IMPRESSION POUR LISVAR.
*
*         ERREUR MAJEURE FATALE DU SYSTEM
*
*  MACRO FSTSYST(ORIGINE,MESSAGE,CODE,STATUT)
*  CETTE MACRO REALISE LE TRAITEMENT D'ERREUR DU PROGICIEL DES FICHIERS
*  STANDARDS.
*

*
*
*         ERREUR MAJEURE FATALE
*
*  MACRO FSTFATA(ORIGINE,MESSAGE,CODE,STATUT)
*  CETTE MACRO REALISE LE TRAITEMENT D'ERREUR DU PROGICIEL DES FICHIERS
*  STANDARDS.
*

*
*
*       FSTERRO  ERREUR MAJEURE
*
*  MACRO FSTERRO(ORIGINE,MESSAGE,CODE,STATUT)
*  CETTE MACRO REALISE LE TRAITEMENT D'ERREUR DU PROGICIEL DES FICHIERS
*  STANDARDS.
*

*
*
*     WARNING  AVERTISSEMENT MESSAGE
*
*  MACRO FSTWARN(ORIGINE,MESSAGE,CODE,STATUT)
*  CETTE MACRO REALISE LE TRAITEMENT D'ERREUR DU PROGICIEL DES FICHIERS
*  STANDARDS.
*

*
*
*
*      ERREUR MESSAGE TYPE INFORMATIQUE
*
*  MACRO FSTINFO(ORIGINE,MESSAGE,CODE,STATUT)
*  CETTE MACRO REALISE LE TRAITEMENT D'ERREUR DU PROGICIEL DES FICHIERS
*  STANDARDS.
*

*
*    MESSAGE QUI PEUVENT SERVIR AU DEBUGGING
*
*
*  MACRO FSTTRIV(ORIGINE,MESSAGE,CODE,STATUT)
*  CETTE MACRO REALISE LE TRAITEMENT D'ERREUR DU PROGICIEL DES FICHIERS
*  STANDARDS.
*

*
*
*
*
*  MACRO STOVEC(VEC, ADR, DIMVEC) STOCKER UN VECTEUR A UNE ADRESSE
*

*
*  MACRO FETVEC(VEC, ADR, DIMVEC) OBTENIR UN VECTEUR D'UNE ADRESSE (SCM,LCM)
*

*
*  MACRO SEPARAT(INDESC)  UTILISE L'ESPACE DE NRECUP POUR LES SEPARATEURS DE
*                         FIN DE FICHIER
*
*  INDESC    INDICE DU FICHIER SEQUENTIEL STANDARD
*

*
*
*         VALID ARGUMENT %2 > %3  ET  ARGUMENT %2 < %4
*
*  MACRO VALID(NOM,I2,I3,I4)
*
*  NOM = NOM DE LA ROUTINE MAX 7 CARACTERES
*  I2  = ARGUMENT A VERIFIER
*  I3  = LA PLUS PETITE VALEUR DE I2
*  I4  = LA PLUS GRANDE VALEUR DE I2
*
*  NORMALEMENT CETTE MACRO VERIFIER LA VALEUR DE IP1 - IP2   FICHIERS
*  STANDARDS.
*

*
*
*REVISION
*         1.0.1 -  AJOUTER 2 TABLES POUR FICHIER RANDOM/SEQ STANDARD
*                   MODIFICATIONS 22 JANV 89
*
*  COMMON FSTC88   CE BLOC SERT A LA GESTION DES FICHIERS
*  ------ ------   STANDARDS (1988 OU AVANT ) QUI SONT OUVERTS.
*
      COMMON /FSTC88/ IUNTAB(0:41),DIRSIZ(0:40), NBCORR(0:40),NUTIL
     % (0:40), NXTADR(0:40), NUMREC(0:40),PRMECR(0:40), NBECR (0:40)
     %, NBRECR(0:40),NRECUP(0:40), NBEFF (0:40), NBEXT (0:40),MSKPRT
     %( 6,0:40),   SEARCH( 6,0:40),MASQUE( 6,0:40),   CURSLT( 12,0:
     %40),MODIFS(0:40),NXTSEQ(0:40),CURSEQ(0:40),DIAGP, LOGFIL,
     % DEBUG, FASTIO, TOLRNC, MSGLVL,NPAIRE,COPYMOD, CHAINE(0:40),
     % EFFDSC(0:40)
      INTEGER   IUNTAB, DIRSIZ, NUTIL , NXTADR, NUMREC, MSKPRT,
     %MASQUE, SEARCH, NBECR , NBRECR, NRECUP, NBEFF,NBEXT , NBCORR,
     % LOGFIL, DEBUG,  CURSLT, TOLRNC,MSGLVL, NPAIRE, NXTSEQ, CURSEQ
     %, CHAINE, EFFDSC
      LOGICAL   PRMECR, FASTIO, MODIFS, DIAGP, COPYMOD
*
*  IUNTAB  IUNTAB(I) >= 0 EST LE NUMERO D'UNITE DU FICHIER DECRIT
*                         PAR L'ENTREE I DES TABLES DE CE BLOC.
*          IUNTAB(I) = -1 INDIQUE UNE ENTREE LIBRE.
*          IUNTAB(MXSTD1) SERT DE SENTINELLE.
*
*  DIRSIZ  TAILLE DU DIRECTEUR SUR DISQUE, I.E. LE NOMBRE MAXIMUM
*          D'ENTREES QU'IL PEUT CONTENIR. (0 = SEQUENTIEL)
*  NUTIL   NOMBRE D'ENTREES EFFECTIVEMENT UTILISEES DANS L'INDEX
*          (MEMOIRE OU DISQUE).
*  NXTADR  PROCHAINE ADRESSE DISPONIBLE POUR ECRIRE UN ENREGISTREMENT.
*
*  NUMREC  POSITION COURANTE DE RECHERCHE DANS LE DIRECTEUR, UTILISEE
*          POUR PASSER AU 'SUIVANT'.  0 INDIQUE DE PARTIR AU DEBUT,
*          TANDIS QUE -1 INDIQUE QUE L'OPERATION 'SUIVANT' EST ILLEGALE.
*          POUR UN FICHIER SEQUENTIEL, NUMREC = 1 (OU 0).
*          -1 = RIEN DE VALIDE
*           0 = RIEN DE VALIDE / RECHERCHE DU SUIVANT ACCEPTABLE
*           1 = APRES QSTSUI/FSTSUI
*           2 = APRES LECTURE
*           FSTSEL = 0 / FSTPOS = -1 / FSTSKP = -1 / FSTRWD = -1
*
*  MSKPRT  MASQUE PARTIEL COURANT, SPECIFIE PAR LA ROUTINE STDMSQ.
*  MASQUE  MASQUE DE RECHERCHE COURANT.
*  SEARCH  VALEUR DE RECHERCHE COURANTE.
*
*  PRMECR  INDIQUE SI ON PEUT ECRIRE SUR LE FICHIER (.TRUE.) OU NON.
*
*  NBECR   NOMBRE TOTAL D'ECRITURES EFFECTUEES SUR CE FICHIER DEPUIS
*          SA CREATION (INCLUANT RE-ECRITURES ET RECUPERATIONS).
*  NBRECR  NOMBRE TOTAL DE RE-ECRITURES EFFECTUEES.
*  NRECUP  NOMBRE DE FOIS QU'ON A RECUPERE UNE QUELCONQUE ENTREE
*          DETRUITE (Y COMPRIS LORS D'UNE RE-ECRITURE).
*  NBEFF   NOMBRE D'ENREGISTREMENTS EFFACES (Y COMPRIS LORS D'UNE
*          RE-ECRITURE).
*  NBEXT   NOMBRE DE FOIS QU'ON A ETENDU L'INDEX-DISQUE.
*  NBCORR  NOMBRE DE FOIS QU'ON A CORRIGE L'ENTETE DU FICHIER.
*
*  CURSLT  ENTREE DE DIRECTEUR ACTIVE (LA SEULE SI SEQUENTIEL)
*
*  FASTIO  INDIQUE SI ON OPTIMISE L'IO AU PRIX DE LA SECURITE
*  MODIFS  INDIQUE (POUR FASTIO) SI LE DIRECTEUR D'UN FICHIER A ETE MODIFIE
*  DIAGP   INDIQUE S'IL FAUT IMPRIMER CERTAINS MESSAGES.
*  LOGFIL  UNITE FORTRAN QUI RECEVRA LES MESSAGES DE TOUT POIL.
*  DEBUG   OPTION DE TRACE POUR LA ROUTINE STDDBG.
*
*          CE BLOC COMMON CONTIENT AUSSI LE RESULTAT DE L'EXPANSION DE LA
*          DERNIERE "SLOT" ACTIVE (VOIR DECK PCKUNP) DANS LES VARIABLES
*          DATE A UBC. LE TOUT EST AUSSI DISPONIBLE SOUS FORME DE TABLE
*          DANS LE TABLEAU TABPRM (FORMAT HOLLERITH MAX. 4 CARACTERES).
*

      COMMON /FSTC88/ DATE,   DEET,   NPAS,   NI,   NJ,   NK,NBITS,
     %  DATYP,  IP1,   IP2,   IP3,   TYPVAR,NOMVAR, ETIQ14, ETIQ58,
     % GRTYP,  IG1,   IG2,IG3,   IG4,   SWA,   LNG,   DLTF,   UBC,
     %ETIQ56, ETIQ78
      INTEGER   DATE,   DEET,   NPAS,   NI,   NJ,   NK,NBITS,  DATYP
     %,  ETIQ14, NOMVAR, TYPVAR, IP1,IP2,   IP3,   GRTYP,  IG1,
     %   IG2,   IG3,IG4,   SWA,   LNG,   DLTF,   UBC,   ETIQ58,
     %ETIQ56, ETIQ78
      INTEGER   TABPRM(26)
      EQUIVALENCE   (TABPRM(1),DATE)
      INTEGER RDSYNC, RDSYN2, SQSYNC, SQSYN2
      COMMON /FSTSYNC/ RDSYNC, RDSYN2, SQSYNC, SQSYN2



*          L'ENSEMBLE DU PROGICIEL DES  FICHIERS STANDARDS.
*
*

       COMMON /TESTC1/PRNT,DEBREC, FINREC, INCREC,F(15,100)
       COMMON /TESTC2/ WK(10000)
       INTEGER DEBREC, FINREC, INCREC,NNBITS,ITIME,IDBG,IPRM,KEY,IUN
       INTEGER INOMBR,IFRM,IOUV,IVOI,IOPT,IER,INFL,IFOIS,IREC,IOUV1
       INTEGER I,J,K,NB,MIP1,MIP2,MIP3,ILIR,IECR,IMASK,CNBITS,CDATYP
       REAL F, WK, VALEUR
       INTEGER IMEM,FSTSEL,NRECSK,IERR,IWEO
       INTEGER TABNUM(170),IWK(200),IWK2(200),NIOUT,NJOUT,NKOUT,IOUV2
       REAL ERRMAX, WK1(200)
      LOGICAL ERR, REWRIT, LUK ,GET, SET, PRNT
      EXTERNAL FSTINF,FSTINL,FSTNBR,FSTOUV,FSTVOI,FSTLIR,FSTLIS,FSTPRM
      EXTERNAL FSTOPC,FSTFRM,FSTECR,FNOM,FSTOPL,FSTOPR,FSTOPI,QSTDBG
      EXTERNAL FSTMSQ, FSTSUI, FSTCVT , FSTRWD, FSTSKP, FSTLUK
      EXTERNAL CCARD,EXDB,QSTINT,FSTPOS,FSTWEO,FSTAPP
      INTEGER EXDB,EXFIN,INDESC,QSTINT,FSTPOS,FSTWEO,FSTAPP
      INTEGER  FSTMSQ, FSTSUI, FSTCVT, FSTRWD, FSTSKP, FSTLUK, ILUK
      INTEGER FSTINF,FSTINL,FSTNBR,FSTOUV,FSTVOI,FSTLIR,FSTLIS,QSTDBG
      INTEGER FSTOPC,FSTFRM,FSTECR,FNOM,FSTOPL,FSTOPR,FSTOPI,FSTPRM
      INTEGER CNNI,CNNJ,CNNK,LIR,IMSQ,IOK,NEQUIV,JDATE
      INTEGER CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CIP1,CIP2,CIP3,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3

      INTEGER HNOMVAR,HETIKET(2),HTYPVAR,HGRTYP,IVOUT

*

      CHARACTER *8 CETIKET, METIQ

      CHARACTER *2 CNOMVAR
      CHARACTER *1 CTYPVAR,CGRTYP
      DATA HNOMVAR/"YY"/, HTYPVAR/"X"/, HGRTYP/"G"/,

     %     HETIKET(1)/"ETIK"/,HETIKET(2)/"ET89"/

      DATA GET/.TRUE./,SET/.FALSE./,PRNT/.FALSE./
*


*------------------------------------------------------------------
*
      WRITE(6,*)' FNOM POUR TAPE10 - 15 - 20 EXECUTE COMMENCE TEST'
      WRITE(6,*)' ------------------------------------------------'
      WRITE(6,*)' '




      IER = FNOM(10,'TAPE10','STD+RND',0)
      IER = FNOM(15,'TAPE15','STD+SEQ',0)
      IER = FNOM(20,'TAPE20','STD+SEQ+UNF+FTN',0)
      IER = FNOM(6,'TAPE6','SEQ+FMT',0)


*




      JDATE= EXDB('TSTC910','V 0.1.0', 'NON')



*
      WRITE(6,1)
1     FORMAT(1H1,30X,'STDTEST 1.0',//,
     %   ' POUR CONFIRMER REUSSITE DU TEST,',/,
     %   ' VOIR MESSAGE *STDTEST REUSSI* A LA FIN',/,
     %   ' ET VERIFIER VISUELLEMENT TESTS PRECEDES DE ---',//)
*
*
*     INITIALISATION AVANT POUR LES MESSAGES DE MSGLVL
      MSGLVL=2
           NPAIRE = 20
*
*
*
      WRITE(6,*)' NPAIRE INITIALISE A 20 DANS STDMAIN  ',NPAIRE
*
*
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' '
      WRITE(6,*) '---TEST 0.0 (FSTNBR) RND FICHIER INEXISTANT-10-'
*
      INOMBR = FSTNBR(10)
      WRITE(6,*)' RESULTAT DE INOMBR= FSTNBR INOMBR=',INOMBR
      WRITE(6,*)' VERIFIER CODE  SI ERREUR APRES FSTNBR(10)'
      IF( (INOMBR .LT.0))THEN
         CALL ERREUR
         STOP
      ENDIF
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 1.0 (FSTOUV) OUVRIR 3 FICHIERS EN CREATION---'


      IOUV  = FSTOUV(10, 'RND' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(10,RND)  ****'
      ENDIF
      REWIND(20)
      IOUV1 = FSTOUV(20, 'SEQ/FTN')
      IF((IOUV1.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(20,SEQ)  ****'
      ENDIF
      IOUV2 = FSTOUV(15, 'SEQ' )
      IF((IOUV2.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ)  ****'
      ENDIF


*
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 1.1 (FSTOUV) OUVRIR 3 FICHIERS DEJA OUVERT---'
      IOUV  = FSTOUV(10, 'RND' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(10,RND)  DEJA OUVERT  ****'
      ENDIF
      IOUV1 = FSTOUV(15, 'SEQ' )
      INDESC = QSTINT(15)
      WRITE(6,*)' NBCORR APRES FSTOUV(15,SEQ) = ',NBCORR(INDESC)
      IF((IOUV1.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ)  DEJA OUVERT  ****'
      ENDIF
      REWIND(20)
      IOUV2 = FSTOUV(20, 'SEQ/FTN' )
      IF((IOUV2.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(20,SEQ)  DEJA OUVERT  ****'
      ENDIF
      IF(IOUV.LT.0.AND.IOUV1.LT.0.AND.IOUV2.LT.0) THEN
      CALL TESTOK
      ELSE
      CALL ERREUR
      STOP
      ENDIF


      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 1.3 (FSTVOI) LISTER FICHIERS SANS ENTREES---'
*
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 1.3', 1)
      ENDIF
      IDBG =QSTDBG('TEST 1.3', 1)
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 2.0 (FSTFRM) [FSTFERM] FERMER LES  FICHIERS---'


      IFRM= FSTFRM(10)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(10)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(10)'
         CALL ERREUR
         STOP
      ENDIF
      IFRM= FSTFRM(15)
      INDESC = QSTINT(15)
      WRITE(6,*)' NBCORR APRES FSTFRM(15) = ',NBCORR(INDESC)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(15)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(15)'
         CALL ERREUR
         STOP
      ENDIF
      IFRM= FSTFRM(20)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(20)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(20)'
         CALL ERREUR
         STOP
      ENDIF
      CALL TESTOK


      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 3.0 (FSTOPC) NIVEAU - TOLERANCE DES MESSAGES---'
*
*         INITIALISER DANS FSTOUV DEFAUT MSGLVL = 2  TOLRNC = 6
*
      MSGLVL=0
      TOLRNC=8
      WRITE(6,*)' POUR FSTOPC INITIALISATION MSGLVL= ',MSGLVL
      WRITE(6,*)' POUR FSTOPC INITIALISATION TOLRNC= ',TOLRNC
*
      IER = FSTOPC('MSGLVL','INFORM',.FALSE.)
      IER = FSTOPC('TOLRNC','SYSTEM',.FALSE.)
      IER = FSTOPC('MSGLVL','INFORM',.TRUE.)
      IER = FSTOPC('TOLRNC','SYSTEM',.TRUE.)
      WRITE(6,*)' NEXT FSTOPC( MSGLVL , INFURM ) EREUR INFURM'
      IER = FSTOPC('MSGLVL','INFURM',.FALSE.)
      WRITE(6,*)' NEXT FSTOPC( TOLRNC , INFURM ) EREUR INFURM'
      IER = FSTOPC('TOLRNC','INFURM',.FALSE.)
      WRITE(6,*)' NEXT FSTOPC( TOXRNC , INFURM ) EREUR TOXRNC - INFURM'
      IER = FSTOPC('TOXRNC','INFURM',.FALSE.)
      WRITE(6,*)' NEXT FSTOPC( MSGGLV , INFORM ) EREUR MSGGLV'
      IER = FSTOPC('MSGGLV','INFORM',.TRUE.)
      WRITE(6,*)' FSTOPC VALEUR FINAL MSGLVL= ',MSGLVL
      WRITE(6,*)' FSTOPC VALEUR FINAL TOLRNC= ',TOLRNC
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 3.1 (FSTOPL) MODIFIE CHOIX DE FASTIO TRUE - FALSE ---'
      WRITE(6,670) FASTIO
  670 FORMAT(' IMPRIME DANS FSTOPL FASTIO= ',L8)
      FASTIO=.FALSE.
      IER = FSTOPL('FASTIO',.TRUE.,GET)
      WRITE(6,678) FASTIO
  678 FORMAT(' GET DANS FSTOPL FASTIO= ',L8)
      IER = FSTOPL('FASTIO',.TRUE.,SET)
      WRITE(6,679) FASTIO
      IER = FSTOPL('FASTIO',.TRUE.,SET)
      WRITE(6,679) FASTIO
  679 FORMAT(' SET DANS FSTOPL FASTIO= ',L8)
      IER = FSTOPL('FASTIO',.TRUE.,GET)
      WRITE(6,678) FASTIO
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 3.2 (FSTOPR) MODIFIE LOCATION REELLE---'
      WRITE(6,*)' *****  AUCUNE VALEUR DANS LE COMMON FSTC88 ****'
*     IER = FSTOPR(VALEUR,  100.0,GET)
*     WRITE(6,680) VALEUR
* 680 FORMAT(' GET  DANS FSTOPR VALEUR= ',F7.3)
*     IER = FSTOPR(VALEUR,  100.0,SET)
*     WRITE(6,677) VALEUR
* 677 FORMAT(' SET  DANS FSTOPR VALEUR= ',F7.3)
*     CALL TESTOK
*     WRITE(6,*)' ##################################################'
*     WRITE(6,*)' ##################################################'
*     WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 3.3 (FSTOPI) MODIFIE LOCATION ENTIERE---'
      MSGLVL=123
      WRITE(6,*)' AVANT FSTOPI MSGLVL= ', MSGLVL
      IER = FSTOPI('MSGLVL',222,GET)
      WRITE(6,697) MSGLVL
  697 FORMAT('  GET  APRES FSTOPI MSGLVL= ',I5)
      WRITE(6,*)' AVANT FSTOPI MSGLVL= ', MSGLVL
      IER = FSTOPI('MSGLVL',2,SET)
      WRITE(6,676) MSGLVL
  676 FORMAT('  SET  APRES FSTOPI MSGLVL= ',I5)
*
      CALL TESTOK
*
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*) '---TEST 4.0 (FSTOUV) OUVRE FICHIERS SANS ENTREES---'
      IOUV = FSTOUV(15, 'SEQ')
      INDESC = QSTINT(15)
      WRITE(6,*)' NBCORR APRES FSTOUV(15,SEQ) = ',NBCORR(INDESC)
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(15,SEQ)'
      ENDIF
      IOUV1= FSTOUV(10, 'RND')
      IF((IOUV1.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(10,RND)'
      ENDIF
      REWIND 20
      IOUV2= FSTOUV(20, 'SEQ/FTN')
      IF((IOUV2.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(20,SEQ)'
      ENDIF
         CALL TESTOK
*
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 5.0 (FSTVOI) LISTER FICHIERS SANS ENTREES---'
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         IOK=1
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 5.0', 1)
      ENDIF
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         IOK=1
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ/FTN)'
      IDBG =QSTDBG('TEST 5.0', 1)
      ENDIF
      IDBG =QSTDBG('TEST 5.0', 1)
         CALL TESTOK
      IOUV  = FSTOUV(10, 'RND' )
      IF((IOUV.NE.-18))THEN
         IOK=1
         WRITE(6,*)' ****  ERREUR FSTOUV(10,RND) DEJA OUVERT  ****'
      ENDIF
      REWIND(20)
      IOUV  = FSTOUV(20, 'SEQ/FTN' )
      IF((IOUV.NE.-18))THEN
         IOK=1
         WRITE(6,*)' ****  ERREUR FSTOUV(20,SEQ) DEJA OUVERT ****'
      ENDIF
      REWIND(20)
      IOUV  = FSTOUV(15, 'SEQ' )
      IF((IOUV.NE.-18))THEN
         IOK=1
         WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ) DEJA OUVERT ****'
      ENDIF
         CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*) '---TEST 5.7(FSTECR) ECRIRE 10 RECORDS '
      DEBREC = 1
      FINREC = 10
      INCREC = 1
*                   DATYP=0 BINAIRE  DATYP=1 REEL  DATYP=2 ENTIER
*                                                  DATYP=3 CARACTERE
*
      REWRIT = .FALSE.
      CALL TSTECR(20, REWRIT   )
      CALL TSTECR(10, REWRIT   )
      CALL TSTECR(15, REWRIT   )
*
      INDESC = QSTINT(15)
      WRITE(6,*)' NBCORR APRES TSTECR(15,REWRIT) = ',NBCORR(INDESC)
      IFRM= FSTFRM(20)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(20)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(20)'
      CALL ERREUR
      STOP
      ENDIF
*     IFRM= FSTFRM(15)
*     WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(15)'
*     IF((IFRM.LT.0))THEN
*        WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(15)'
*     CALL ERREUR
*     STOP
*     ENDIF
*------------------- ----------------------------------------------
      WRITE(6,*) '---TEST 5.8(FSTAPP)CHERCHER FIN DU FICHIER SEQ(31)'
      WRITE(6,*)'  AJOUTER 10 RECORDS SUR LE FICHIER   '
      IOUV  = FSTOUV(20, 'SEQ/FTN' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(20,SEQ)  ****'
      CALL ERREUR
      STOP
      ENDIF
*
      IER= FSTRWD(20)
      IER= FSTRWD(15)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      IDBG =QSTDBG('TEST 5.8', 1)
      ENDIF
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 5.8', 1)
      ENDIF
      IER= FSTRWD(20)
      IER= FSTRWD(15)
      IER = 0
      IER = FSTAPP(20, 'SEQ/FTN')
      IF(IER.LT.0) WRITE(6,*)' IER < 0 APRES FSTAPP(20) IER= ',IER
      WRITE(6,*)' VALEUR DE IER APRES IER=FSTAPP(20) IER= ',IER
*
      IER = FSTAPP(15, 'SEQ')
      IF(IER.LT.0) WRITE(6,*)' IER < 0 APRES FSTAPP(15) IER= ',IER
      WRITE(6,*)' VALEUR DE IER APRES IER=FSTAPP(15) IER= ',IER
      IER= FSTRWD(15)
      IER = FSTAPP(15, 'SEQ')
      DEBREC = 11
      FINREC = 20
      INCREC = 1
      REWRIT = .FALSE.
      CALL TSTECR(20, REWRIT  )
      CALL TSTECR(10, REWRIT  )
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      IDBG =QSTDBG('TEST 5.8', 1)
      ENDIF
      CALL TSTECR(15, REWRIT  )
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 5.8', 1)
      ENDIF
      IER = FSTAPP(20, 'SEQ/FTN')
      IF(IER.LT.0) WRITE(6,*)' IER < 0 APRES FSTAPP(20) IER= ',IER
      WRITE(6,*)' VALEUR DE IER APRES IER=FSTAPP(20) IER= ',IER
      IER= FSTRWD(15)
      IER = FSTAPP(15, 'SEQ')
      IF(IER.LT.0) WRITE(6,*)' IER < 0 APRES FSTAPP(15) IER= ',IER
      WRITE(6,*)' VALEUR DE IER APRES IER=FSTAPP(15) IER= ',IER
      IFRM= FSTFRM(15)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(15)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(15)'
         CALL ERREUR
         STOP
      ENDIF
      IOUV2 = FSTOUV(15, 'SEQ' )
      IF((IOUV2.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ)  ****'
      ENDIF
      IER = FSTAPP(15, 'SEQ')
      DEBREC = 11
      FINREC = 20
      INCREC = 1
      REWRIT = .FALSE.
      CALL TSTECR(15, REWRIT  )
      IFRM= FSTFRM(15)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(15)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(15)'
         CALL ERREUR
         STOP
      ENDIF
      IFRM= FSTFRM(20)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(20)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(20)'
         CALL ERREUR
         STOP
      ENDIF
      IOUV1 = FSTOUV(20, 'SEQ/FTN')
      IF((IOUV1.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(20,SEQ)  ****'
      ENDIF
      IOUV2 = FSTOUV(15, 'SEQ' )
      IF((IOUV2.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ)  ****'
      ENDIF
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         IOK=1
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ/FTN)'
      IDBG =QSTDBG('TEST 5.8', 1)
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 5.8', 1)
      ENDIF
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 5.9 (FSTECR - FSTFRM - FSTLIR) ECRIRE - LIRE   ---'
      WRITE(6,*)' AJOUTER 10 RECORDS 10 - 15 - 20 ET RELIRE   '
      WRITE(6,*)' IMPRIME PREMIERE ET DERNIERE VALEUR SI PRNT=.TRUE.'
      PRNT = .TRUE.
*
      DEBREC = 21
      FINREC = 31
      INCREC = 1
*                   DATYP=0 BINAIRE  DATYP=1 REEL  DATYP=2 ENTIER
*                                                  DATYP=3 CARACTERE
*
      REWRIT = .TRUE.
      CALL TSTECR(10, REWRIT )
      CALL TSTECR(15, REWRIT  )
      CALL TSTECR(20, REWRIT   )
      IFRM= FSTFRM(10)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(10)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(10)'
      CALL ERREUR
      STOP
      ENDIF
      IFRM= FSTFRM(15)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(15)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(15)'
      CALL ERREUR
      STOP
      ENDIF
      IOUV  = FSTOUV(15, 'SEQ' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ)  IOUV=',IOUV
      CALL ERREUR
      STOP
      ENDIF
      IER = FSTAPP(15, 'SEQ')
      WRITE(6,*)' **************************************'
      WRITE(6,*)' APRES FSTAPP(15) FICHIER FERMER'
      WRITE(6,*)' **************************************'


      IFRM= FSTFRM(20)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(20)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(20)'
      CALL ERREUR
      STOP
      ENDIF
      IFRM= FSTFRM(15)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(15)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(15)'
      CALL ERREUR
      STOP
      ENDIF
*------------------------------------------------------------------
      REWIND(20)
      IOUV  = FSTOUV(20, 'SEQ/FTN' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(20,SEQ)  ****'
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(10)
      IOUV  = FSTOUV(10, 'RND' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(10,RND) IOUV=',IOUV
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(15)
      IOUV  = FSTOUV(15, 'SEQ' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ)  IOUV=',IOUV
      CALL ERREUR
      STOP
      ENDIF
*
      DEBREC = 1
      FINREC = 11
      INCREC = 1
      IER= FSTRWD(10)
      CALL TSTLIR(10 )
      IER= FSTRWD(15)
      CALL TSTLIR(15 )
      IER= FSTRWD(20)
      CALL TSTLIR(20 )
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 6.0 (FSTECR)  ECRIRE 20 RECORD FILE 10 - 15 - 20 ---'
      WRITE(6,*)' ECRIRE PAR DESSUS 10 - 15'
      PRNT = .FALSE.
*
      DEBREC = 31
      FINREC = 40
      INCREC = 1
*                   DATYP=0 BINAIRE  DATYP=1 REEL  DATYP=2 ENTIER
*                                                  DATYP=3 CARACTERE
*
      REWRIT = .TRUE.
      CALL TSTECR(10, REWRIT )
      CALL TSTECR(15, REWRIT  )
      CALL TSTECR(20, REWRIT   )
      WRITE(6,*)' ------------------------------------------------ '
      WRITE(6,*)'  REECRIRE PAR DESSUS AVEC REWRIT=TRUE '
      CALL TSTECR(10, REWRIT )
      CALL TSTECR(15, REWRIT  )
      WRITE(6,*)'  ECRIRE A LA FIN AVEC REWRIT=FALSE RND SEULEMENT '
      DEBREC = 41
      FINREC = 50
      INCREC = 1
      REWRIT = .FALSE.
      CALL TSTECR(10, REWRIT )
      CALL TSTECR(20, REWRIT   )
      CALL STDTEOF(20)
      CALL TSTECR(15, REWRIT  )
      CALL TSTECR(10, REWRIT )
      DEBREC = 51
      FINREC = 60
      INCREC = 1
      CALL TSTECR(10, REWRIT )
      CALL TSTECR(15, REWRIT  )
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*??????????????????????????????????????????????????????????
*          stop
*???????????????????????????????????????????????????????????
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 6.1 (FSTVOI) LISTER FICHIERS AVEC ENTREES---'
      WRITE(6,*) '  FERMER ET REOUVRIR FICHIERS'
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 6.1', 1)
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      IDBG =QSTDBG('TEST 6.1', 1)
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 6.1', 1)
      CALL ERREUR
      STOP
      ENDIF
      IDBG =QSTDBG('TEST 6.1', 1)
      IFRM= FSTFRM(10)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(10)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(10)'
      CALL ERREUR
      STOP
      ENDIF
      IFRM= FSTFRM(15)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(15)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(15)'
      CALL ERREUR
      STOP
      ENDIF
      IFRM= FSTFRM(20)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(20)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(20)'
      CALL ERREUR
      STOP
      ENDIF
*------------------------------------------------------------------
      REWIND(20)
      IOUV  = FSTOUV(20, 'SEQ/FTN' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(20,SEQ)  ****'
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(10)
      IOUV  = FSTOUV(10, 'RND' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(10,RND) IOUV=',IOUV
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(15)
      IOUV  = FSTOUV(15, 'SEQ' )
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ)  IOUV=',IOUV
      CALL ERREUR
      STOP
      ENDIF
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 6.2(FSTVOI)LISTER FICHIERS APRES FSTFRM - FSTOUV---'
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 6.2', 1)
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 6.2', 1)
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      IDBG =QSTDBG('TEST 6.2', 1)
      CALL ERREUR
      STOP
      ENDIF
      IDBG =QSTDBG('TEST 6.2', 1)
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 6.3 (FSTRWD, FSTLIR, FSTINF, FSTSUI) LIRE ET'
      WRITE(6,*) '   VERIFIER  LE CONTENU DES 3 FICHIERS'
      DEBREC = 1
      FINREC = 20
      INCREC = 3
      IER= FSTRWD(10)
      CALL TSTLIR(10 )
      DEBREC = 1
      FINREC = 20
      INCREC = 3
      IER= FSTRWD(15)
      CALL TSTLIR(15 )
      DEBREC = 1
      FINREC = 20
      INCREC = 5
      IER= FSTRWD(20)
      CALL TSTLIR(20 )
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 6.4(FSTFRM,FSTVOI,FSTOUV,FSTLIR,FSTLUK) FERMER LES'
      WRITE(6,*)
     % '   FICHIERS, LES REOUVRIR, LIRE ET REVERIFIER CONTENU---'
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IFRM= FSTFRM(10)
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(10)'
      ENDIF
      IFRM= FSTFRM(15)
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(15)'
      ENDIF
      IFRM= FSTFRM(20)
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(20)'
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, RND)'
      ENDIF
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      ENDIF
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      ENDIF
      IOUV = FSTOUV(10, 'RND')
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(10,RND)'
      ENDIF
      IOUV = FSTOUV(15, 'SEQ')
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(15,SEQ)'
      ENDIF
      IER= FSTRWD(20)
      IOUV = FSTOUV(20, 'SEQ/FTN')
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(20,SEQ)'
      ENDIF
      IER= FSTRWD(20)
      CALL TSTLIR(20 )
      IER= FSTRWD(15)
      CALL TSTLIR(15  )
      IER= FSTRWD(10)
      CALL TSTLIR(10  )
      WRITE(6,*)' LIRE REC IP1=5 ET IP1=10 FSTLUK FILE(10)'
      DEBREC = 5
      FINREC = 10
      INCREC = 5
      IER= FSTRWD(10)
      CALL TSTLIR(10)
      WRITE(6,*)' LIRE REC IP1=5 ET IP1=10 FSTLUK FILE(15)'
      IER= FSTRWD(15)
      CALL TSTLIR(15)
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 6.5(FSTEFF,FSTVOI) EFFACER RECS 3 - 13 FILE 10 - 15'
      WRITE(6,*)'  VOIR FICHIER 10 15 '
      WRITE(6,*)' '
      DEBREC = 3
      FINREC = 30
      INCREC = 10
      IER= FSTRWD(10)
      CALL TSTEFF(10)
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 6.5', 1)
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, RND)'
      IDBG =QSTDBG('TEST 6.5', 1)
      ENDIF
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      IDBG =QSTDBG('TEST 6.5', 1)
      ENDIF
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 6.6 (FSTECR) AJOUTER 30 RECORDS FICHIER 10-15-20'
      DEBREC = 61
      FINREC = 70
      INCREC = 1
      REWRIT = .FALSE.
      CALL TSTECR(10, REWRIT  )
      CALL TSTECR(15, REWRIT  )
      CALL TSTECR(20, REWRIT  )
      IDBG =QSTDBG('TEST 6.6', 1)
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 6.7 (FSTVOI) VOIR LES  FICHIERS   10 - 15'
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 6.7', 1)
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 6.7', 1)
      ENDIF
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ/FTN)'
      IDBG =QSTDBG('TEST 6.7', 1)
      ENDIF
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '


      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*) '---TEST 7.0(FSTEFF,FSTECR,FSTWEO,FSTVOI) EFFACE '
      WRITE(6,*) ' 3 RECORDS ECRIRE AVEC REMPLACEMENTS ET EXTENSIONS '
      WRITE(6,*) ' REWRIT = .TRUE. CALL TSTECR(UNIT, REWRIT) ---'
      WRITE(6,*)' AJOUTER RECORD 61 A 70 '
      DEBREC = 4
      FINREC = 30
      INCREC = 10
      IER= FSTRWD(10)
      CALL TSTEFF(10)
      REWRIT = .TRUE.
      DEBREC = 71
      FINREC = 80
      INCREC = 1
      CALL TSTECR(10, REWRIT )
      CALL TSTECR(15, REWRIT )
      CALL TSTECR(20, REWRIT )
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 7.0', 1)
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 7.0', 1)
      ENDIF
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      IDBG =QSTDBG('TEST 7.0', 1)
      ENDIF
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*) '---TEST 7.1 (FSTLIR) RELIRE ET VERIFIER---'
      DEBREC = 1
      FINREC = 70
      INCREC = 1
      IER= FSTRWD(10)
      CALL TSTLIR(10  )
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IER= FSTRWD(15)
      CALL TSTLIR(15  )
      IER= FSTRWD(20)
      CALL TSTLIR(20  )
      CALL TESTOK


      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 7.2(FSTFRM,FSTOUV,FSTLIR,FSTECR,FSTLIR) FERMER '
      WRITE(6,*) ' OUVRIR, LIRE ECRIRE ETENDRE FICHIER REVERIFIER---'
      IFRM = FSTFRM(10)
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR AVEC FSTFRM UNIT 10'
      CALL ERREUR
      STOP
      ENDIF
      IOUV = FSTOUV(10,'RND')
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(10,RND)'
      CALL ERREUR
      STOP
      ENDIF
      IFRM = FSTFRM(15)
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR AVEC FSTFRM UNIT 10'
      CALL ERREUR
      STOP
      ENDIF
      IOUV = FSTOUV(15,'SEQ')
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(10,RND)'
      CALL ERREUR
      STOP
      ENDIF
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      REWRIT = .FALSE.
      IER= FSTRWD(10)
      CALL TSTLIR(10)
      IER= FSTRWD(15)
      CALL TSTLIR(15)
      IER= FSTRWD(20)
      IER= FSTRWD(15)
      REWRIT = .TRUE.
*
      CALL TSTECR(10, REWRIT )
      CALL TSTECR(15, REWRIT )
      CALL TSTECR(20, REWRIT )
      IFRM = FSTFRM(10)
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR AVEC FSTFRM UNIT 10'
      CALL ERREUR
      STOP
      ENDIF
      IOUV = FSTOUV(10,'RND')
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(10,RND)'
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 7.2', 1)
      CALL ERREUR
      STOP
      ENDIF
      IFRM = FSTFRM(15)
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR AVEC FSTFRM UNIT 10'
      CALL ERREUR
      STOP
      ENDIF
      IOUV = FSTOUV(15,'SEQ')
      IF((IOUV.LT.0))THEN
         WRITE(6,*)' ERREUR FSTOUV(15,SEQ)'
      CALL ERREUR
      STOP
      ENDIF
*
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IER= FSTRWD(10)
      CALL TSTLIR(10)
      IER= FSTRWD(15)
      CALL TSTLIR(15)
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 7.2', 1)
      CALL ERREUR
      STOP
      ENDIF
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 8.0(FSTLIR,FSTINF,FSTVOI) LIRE  AVEC WILDCARD---'
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 8.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IER= FSTRWD(10)
      CALL TSTLIR(10)
      DO 23046 IFOIS = 1, 10
         IREC = FSTLIR(F, 10, NI, NJ, NK, -1, ' ', IFOIS, -1, -1,
     %    ' ',' ')
23046 CONTINUE
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 8.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 8.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IER= FSTRWD(10)
      CALL TSTLIR(10)
      DO 24046 IFOIS = 1, 10
         IREC = FSTLIR(F, 10, NI, NJ, NK, -1, ' ', IFOIS, -1, -1,
     %    ' ',' ')
24046 CONTINUE
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 8.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IER= FSTRWD(15)
      CALL TSTLIR(15)
      IER= FSTRWD(15)
      DO 25046 IFOIS = 1, 10
         IREC = FSTLIR(F, 15, NI, NJ, NK, -1, ' ', IFOIS, -1, -1,
     %    ' ',' ')
25046 CONTINUE
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 8.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      IDBG =QSTDBG('TEST 8.0', 1)
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IER= FSTRWD(15)
      CALL TSTLIR(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 8.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 9.0 (FSTINL) LIRE FICHIER 10 - 15  IMPRIME TABLEAU'
*
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 9.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IER= FSTRWD(10)
      CALL TSTLIR(10)
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 9.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      DEBREC = 1
      FINREC = 10
      INCREC = 1
      IER= FSTRWD(10)
      CALL TSTLIR(10)
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 1.3', 1)
      CALL ERREUR
      STOP
      ENDIF
      INFL =  FSTINL(10, NI, NJ, NK, -1, ' ', -1, -1, -1,
     %         ' ', '  ', TABNUM, NB,170)
*
       WRITE(6,*)' VALEUR DE NB APRES FSTINL(10)TABNUM(NB) NB=',NB
            WRITE(6,21) (TABNUM(K), K=1,NB)
21          FORMAT(' TABNUM=',7(10I10/))
*
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 9.1 (FSTPRM) LIRE ENTETE 10-15-20  AVEC KEY=FSTINF'


      KEY = FSTINF(10,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
       WRITE(6,*)' FICHIER 10 VALEUR DE KEY APRES KEY=FSTINF= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
                WRITE(6,*)' NI DE FSTPRM =',CNI
                WRITE(6,*)' NJ DE FSTPRM =',CNJ
                WRITE(6,*)' NK DE FSTPRM =',CNK
       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
                WRITE(6,*)' NI DE FSTLIR =',CNNI
                WRITE(6,*)' NJ DE FSTLIR =',CNNJ
                WRITE(6,*)' NK DE FSTLIR =',CNNK
*--------------------------------------------------------------------
*     IER= FSTRWD(15)
      IER= FSTPOS(15, 0)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ', 2, 4, 6,' ',' ')
       WRITE(6,*)' FICHIER 15 VALEUR DE KEY APRES KEY=FSTINF= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
                WRITE(6,*)' NI DE FSTPRM =',CNI
                WRITE(6,*)' NJ DE FSTPRM =',CNJ
                WRITE(6,*)' NK DE FSTPRM =',CNK
      IER= FSTPOS(15, 0)
       LIR = FSTLIR(F, 15,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
                WRITE(6,*)' NI DE FSTLIR =',CNNI
                WRITE(6,*)' NJ DE FSTLIR =',CNNJ
                WRITE(6,*)' NK DE FSTLIR =',CNNK
      IER= FSTPOS(15, 0)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ', 1, 1, 3,' ',' ')
       WRITE(6,*)' FICHIER 15 VALEUR DE KEY APRES KEY=FSTINF= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
                WRITE(6,*)' NI DE FSTPRM =',CNI
                WRITE(6,*)' NJ DE FSTPRM =',CNJ
                WRITE(6,*)' NK DE FSTPRM =',CNK
      IER= FSTPOS(15, 0)
       LIR = FSTLIR(F, 15,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
                WRITE(6,*)' NI DE FSTLIR =',CNNI
                WRITE(6,*)' NJ DE FSTLIR =',CNNJ
                WRITE(6,*)' NK DE FSTLIR =',CNNK
*--------------------------------------------------------------------
       REWIND(20)
      KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
       WRITE(6,*)' FICHIER 20 VALEUR DE KEY APRES KEY=FSTINF= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
                WRITE(6,*)' NI DE FSTPRM =',CNI
                WRITE(6,*)' NJ DE FSTPRM =',CNJ
                WRITE(6,*)' NK DE FSTPRM =',CNK
         REWIND(20)
       LIR = FSTLIR(F, 20,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
                WRITE(6,*)' NI DE FSTLIR =',CNNI
                WRITE(6,*)' NJ DE FSTLIR =',CNNJ
                WRITE(6,*)' NK DE FSTLIR =',CNNK
*
*-----------------------------------------------------------------
*
       WRITE(6,*)' *****SPECIAL LIRE 10 PREMIERS RECS FSTINL****'
       WRITE(6,*)' EN UTILISANT FSTPRM DANS LOOP 1 A NB'
      INFL =  FSTINL(10, NI, NJ, NK, -1, ' ', -1, -1, -1,
     %         ' ', '  ', TABNUM, NB,170)
*
       WRITE(6,*)' VALEUR DE NB APRES FSTINL(10)TABNUM(NB) NB=',NB
            WRITE(6,21) (TABNUM(K), K=1,NB)
       DO 54321 I=1,10
*
      IPRM=FSTPRM(TABNUM(I),CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,
     %      CDATYP,CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
54321 CONTINUE
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*) '---TEST 9.2 (FSTMSQ) BATIR MASQUE PARTIEL'
      WRITE(6,*) ' LIRE MASQUE PARTIEL---'
      IER= FSTRWD(10)
      KEY = FSTINF(10,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      DO 54326 IMASK=0,1000,100
      IMSQ=FSTMSQ(10,IMASK, -1, -1, '        ',.FALSE.)
      WRITE(6,*)' POSITION IP1 BIT= ',IMASK
      KEY = FSTSUI(10,NIOUT,NJOUT,NKOUT)
      WRITE(6,*)' KEY APRES FSTSUI DO 54328 KEY= ',KEY
      WRITE(6,*)' KEY APRES FSTSUI DO 54326 KEY= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,-1,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
       WRITE(6,*)' VALEUR DE CIP1 APRES FSTLIR= ',CIP1
54326 CONTINUE
      IER= FSTRWD(10)
      WRITE(6,*)' ******** TEST 10 FSTSUI CONSECUTIFS *********'
      KEY = FSTINF(10,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      DO 54328 IMASK=1,10
      KEY = FSTSUI(10,NIOUT,NJOUT,NKOUT)
      WRITE(6,*)' KEY APRES FSTSUI DO 54328 KEY= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       LIR = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
54328 CONTINUE
      WRITE(6,*)' '
      WRITE(6,*)' INITIALISE BIT DE IP1 1000-2000-3000.......   '
      WRITE(6,*)' '
      WRITE(6,*)' '
      KEY = FSTINF(10,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      DO 54327 IMASK=1000,10000,1000
      IMSQ=FSTMSQ(10,IMASK, -1, -1, '********',.FALSE.)
      WRITE(6,*)' POSITION IP1 BIT= ',IMASK
      KEY = FSTSUI(10,NIOUT,NJOUT,NKOUT)
      WRITE(6,*)' KEY APRES FSTSUI DO 54327 KEY= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,-1,'XXXXXXXX',
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
       WRITE(6,*)' VALEUR DE CIP1 APRES FSTLIR= ',CIP1
54327 CONTINUE
*---------------------------------------------------------------------
      WRITE(6,*)' '
      WRITE(6,*)' INITIALISE BIT DE IP1 0                    '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' FSTMSQ(10, 0, .......'
*
      IMSQ=FSTMSQ(10, 0, -1, -1, '        ',.FALSE.)
*
      INFL =  FSTINL(10, NI, NJ, NK, -1, 'ETIKET89',-1, -1, -1,
     %         ' ', '  ', TABNUM, NB,170)
       WRITE(6,*)' VALEUR DE NB APRES FSTINL(10)TABNUM(NB) NB=',NB
*
       WRITE(6,*)' *****SPECIAL LIRE 10 PREMIERS RECS FSTINL****'
       DO 54320 I=1,10
*
      IPRM=FSTPRM(TABNUM(I),CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,
     %     CDATYP,CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,-1,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
54320 CONTINUE
            WRITE(6,21) (TABNUM(K), K=1,NB)
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' ---------------------------------------------  '
      WRITE(6,*)' INITIALISE BIT DE IP1 1   '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' FSTMSQ(10, 1, .......'
*
      IMSQ=FSTMSQ(10, 1, -1, -1, '        ',.FALSE.)
*
      INFL =  FSTINL(10, NI, NJ, NK, -1, 'ETIKET89',-1, -1, -1,
     %         ' ', '  ', TABNUM, NB,170)
       WRITE(6,*)' VALEUR DE NB APRES FSTINL(10)TABNUM(NB) NB=',NB
*
       WRITE(6,*)' *****SPECIAL LIRE 10 PREMIERS RECS FSTINL****'
       DO 54322 I=1,10
*
      IPRM=FSTPRM(TABNUM(I),CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,
     %     CDATYP,CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,-1,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
54322 CONTINUE
            WRITE(6,21) (TABNUM(K), K=1,NB)
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' ---------------------------------------------  '
      WRITE(6,*)' INITIALISE BIT DE IP1 32767   '
      WRITE(6,*)' '
      WRITE(6,*)' '
*
      WRITE(6,*)' FSTMSQ(10, 32767.....'
      IMSQ=FSTMSQ(10, 32767, -1, -1, '        ',.FALSE.)
*
      INFL =  FSTINL(10, NI, NJ, NK, -1, 'ETIKET89',-1, -1, -1,
     %         ' ', '  ', TABNUM, NB,170)
       WRITE(6,*)' VALEUR DE NB APRES FSTINL(10)TABNUM(NB) NB=',NB
*
       WRITE(6,*)' *****SPECIAL LIRE 10 PREMIERS RECS FSTINL****'
       DO 54323 I=1,10
*
      IPRM=FSTPRM(TABNUM(I),CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,
     %    CDATYP,CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,-1,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
54323 CONTINUE
            WRITE(6,21) (TABNUM(K), K=1,NB)
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' ---------------------------------------------  '
      WRITE(6,*)' INITIALISE BIT DE IP1 -1   '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' FSTMSQ(10,-1, .......'
      IMSQ=FSTMSQ(10, -1, -1, -1, '        ',.FALSE.)
*
      INFL =  FSTINL(10, NI, NJ, NK, -1, 'ETIKET89',-1, -1, -1,
     %         ' ', '  ', TABNUM, NB,170)
*
       WRITE(6,*)' VALEUR DE IP1 DANS FSTLIR IP1=1'
*
       DO 54324 I=1,NB
*
      IPRM=FSTPRM(TABNUM(I),CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,
     %     CDATYP,CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)


       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,-1,CETIKET,
     %     CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
54324 CONTINUE
            WRITE(6,21) (TABNUM(K), K=1,NB)
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' ---------------------------------------------  '
      WRITE(6,*)' INITIALISE BIT DE IP1 0  ETIKET=[ ** *   ]   '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' FSTMSQ(10, 0, .......'
*
       MIP3=4096
      IMSQ=FSTMSQ(10,10,100,MIP3, ' ** *   ',.FALSE.)
      IF(IMSQ.NE.-24) THEN
         CALL ERREUR
         STOP
       ENDIF
*
      INFL =  FSTINL(10, NI, NJ, NK, -1, 'ETIKET89',-1, -1, -1,
     %         ' ', '  ', TABNUM, NB,170)
       WRITE(6,*)' VALEUR DE NB APRES FSTINL(10)TABNUM(NB) NB=',NB
*
       WRITE(6,*)' *****SPECIAL LIRE 10 PREMIERS RECS FSTINL****'
       DO 54329 I=1,10
*
      IPRM=FSTPRM(TABNUM(I),CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,
     %     CDATYP,CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       LIR = FSTLIR(F, 10,CNNI,CNNJ,CNNK,-1,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
54329 CONTINUE
            WRITE(6,21) (TABNUM(K), K=1,NB)
      WRITE(6,*)' ---------------------------------------------  '


      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*) '---TEST 9.3 (FSTMSQ) IMPRIME MASQUE PARTIEL'


      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' ---------------------------------------------  '
      WRITE(6,*)' FSTMSQ(IUN, MIP1,MIP2,MIP3,METIQ, .TRUE.)  '
      WRITE(6,*)' '
      WRITE(6,*)' '
      IER= FSTRWD(15)
      IMSQ=FSTMSQ(15, MIP1, MIP2, MIP3, METIQ, .TRUE.)
*
*
      WRITE(6,63) MIP1, MIP2, MIP3, METIQ
63    FORMAT(' MIP1,MIP2,MIP3=',3I6,' METIQ=',A8)
*
            WRITE(6,21) (TABNUM(K), K=1,NB)


      IER= FSTRWD(15)
      IMSQ=FSTMSQ(15, MIP1, MIP2, MIP3, METIQ, .TRUE.)
*
      WRITE(6,63) MIP1, MIP2, MIP3, METIQ
      WRITE(6,*)' IMSQ DE FSTMSQ(15......  TRUE)',IMSQ
      IF(IMSQ.NE.0) THEN
         CALL ERREUR
         STOP
       ENDIF


*
      WRITE(6,63) MIP1, MIP2, MIP3, METIQ
*
            WRITE(6,21) (TABNUM(K), K=1,NB)
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*) '---TEST 9.4(FSTRWD-FSTPOS-FSTSKP) REWIND FICHIER'
      WRITE(6,*)' FSTSKP SKIP + AVANT  - RECULE NREC # DE RECS  '
      WRITE(6,*)' FSTPOS POSITIONNE FICHIER AVEC KEY  '
      IER= FSTRWD(15)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ', 2,  4, 6,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      IER= FSTRWD(15)
*     IER= FSTPOS(15, KEY)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      IER= FSTRWD(15)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      WRITE(6,*)' CALCUL POSITION DU REC AVEC FSTPOS,KEY= ',KEY
      IER= FSTPOS(15, 0 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      KEY = FSTSUI(15,NIOUT,NJOUT,NKOUT)
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      WRITE(6,*)'***** FAIRE UN FSTLIR APRES UN REWIND *********   '
      IER= FSTRWD(15)
       LIR = FSTLIR(F, 15,CNNI,CNNJ,CNNK,-1   ,' ',
     %         -1,-1,-1,' ',' ')
*    AJOUTER FSTPOS POUR METTRE LE FICHIER AU DEBUT
*
      IER= FSTPOS(15, 0 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',10,100,30,' ',' ')
       WRITE(6,*)' FICHIER 15 VALEUR DE KEY  KEY=FSTINF= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
                WRITE(6,*)' NI DE FSTPRM =',CNI
                WRITE(6,*)' NJ DE FSTPRM =',CNJ
                WRITE(6,*)' NK DE FSTPRM =',CNK
      IER = FSTSKP(15,-1 )
       LIR = FSTLIR(F, 15,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
                WRITE(6,*)' NI DE FSTLIR =',CNNI
                WRITE(6,*)' NJ DE FSTLIR =',CNNJ
                WRITE(6,*)' NK DE FSTLIR =',CNNK
*
         IER= FSTRWD(15)
       LIR = FSTLIR(F, 15,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
         IER= FSTRWD(15)
       LIR = FSTLIR(F, 15,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
         IER= FSTRWD(15)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',11,121,33,' ',' ')
       WRITE(6,*)' FICHIER 15 VALEUR DE KEY  KEY=FSTINF= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
      IER = FSTSKP(15,-1 )
       LIR = FSTLIR(F, 15,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
         IER= FSTRWD(15)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',10,100,30,' ',' ')
       WRITE(6,*)' FICHIER 15 VALEUR DE KEY  KEY=FSTINF= ',KEY
      IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
      IER = FSTSKP(15,-1 )
       LIR = FSTLIR(F, 15,CNNI,CNNJ,CNNK,CDATE,CETIKET,
     %         CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR)
*
      IER = FSTSKP(15,-1 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15,-2 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15,-3 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15,-4 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15,-5 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15,-50)


      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15, 1 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15, 2 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15, 3 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15,-3 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15, 4 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15, 5 )
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
        WRITE(6,*)' TEST 9.4 VALEUR ILUK APRES FSTSKP =',ILUK
      IER = FSTSKP(15, 50)


      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ',' ')
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
*-------------------------------------------------------------------
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      IDBG =QSTDBG('TEST 9.4', 1)
      ENDIF
      IER= FSTRWD(20)
      KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
       WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      IER= FSTRWD(20)
      IER = FSTSKP(20, 5 )
      KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
       WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)


      IER = FSTSKP(20, 5 )
      KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
       WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(20, -5)
      KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
       WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(20, 5 )
      KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
       WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(20,- 5 )
      KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
       WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
*     IER = FSTSKP(20, 5 )
*     KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
*      WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
*     ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
*     IER = FSTSKP(20, 5 )
*     KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
*      WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
*     ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)
*     IER = FSTSKP(20, 40)
*     KEY = FSTINF(20,NIOUT,NJOUT,NKOUT,-1,' ',-1, -1,-1,' ',' ')
*      WRITE(6,*)' FICHIER 20 VALEUR DE KEY  KEY=FSTINF= ',KEY
*     ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)


      WRITE(6,*)' TERMINE EXECUTION DE FSTSKP - FSTLUK X FOIS    '
      WRITE(6,*)' '
      WRITE(6,*)' ---------------------------------------------  '
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*) '---TEST 9.5(FSTVOI-FSTINF-FSTLUK-FSTSEL,FSTSKP)'
      WRITE(6,*)' FSTRWD - FSTPOS - FSTSKP - FSTSEL - FSTLIS  )'
      WRITE(6,*)' FSTPOS POSITIONNE FICHIER AVEC KEY  '
      NRECSK = 9
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 9.5', 1)
      CALL ERREUR
      STOP
      ENDIF
      IER= FSTRWD(15)
      IMEM= FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','UU')
      IPRM=FSTPRM(IMEM,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
     %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
     %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
     %            EXTRA1,EXTRA2,EXTRA3)
       WRITE(6,*)' UU APRES FSTINF VALEUR DE IMEM= ',IMEM
      ILUK = FSTLUK(F,IMEM,NIOUT,NJOUT,NKOUT)
      KEY = FSTINF(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','VV')
       WRITE(6,*)' VV APRES FSTINF VALEUR DE KEY = ',KEY
      ILUK = FSTLUK(F,KEY,NIOUT,NJOUT,NKOUT)


       WRITE(6,*)' AVANT FSTPOS VALEUR DE IMEM= ',IMEM
      IER= FSTPOS(15, IMEM)
      IER = FSTSKP(15, 1)
      IER = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','UU')
       WRITE(6,*)' UU  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(15, NRECSK)
      KEY = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','VV ')
       WRITE(6,*)' VV  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
*------------------------------------------------------------------
       WRITE(6,*)' AVANT FSTPOS VALEUR DE IMEM= ',IMEM
      IER= FSTPOS(15, IMEM)
      IER = FSTSKP(15, 2)
      IER = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','UU')
       WRITE(6,*)' UU  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(15, NRECSK)
      KEY = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','VV ')
       WRITE(6,*)' VV  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)


       WRITE(6,*)' AVANT FSTPOS VALEUR DE IMEM= ',IMEM
      IER= FSTPOS(15, IMEM)
      IER = FSTSKP(15, 3)
      IER = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','UU')
       WRITE(6,*)' UU  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(15, NRECSK)
      KEY = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','VV ')
       WRITE(6,*)' VV  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
       WRITE(6,*)' AVANT FSTPOS VALEUR DE IMEM= ',IMEM
      IER= FSTPOS(15, IMEM)
      IER = FSTSKP(15, 4)
      IER = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','UU')
       WRITE(6,*)' UU  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(15, NRECSK)
      KEY = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','VV ')
       WRITE(6,*)' VV  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
       WRITE(6,*)' AVANT FSTPOS VALEUR DE IMEM= ',IMEM
      IER= FSTPOS(15, IMEM)
      IER = FSTSKP(15, 5)
      IER = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','UU')
       WRITE(6,*)' UU  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(15, NRECSK)
      KEY = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','VV ')
       WRITE(6,*)' VV  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
       WRITE(6,*)' ******  ICI ON FAIT UN FSTRWD(15) ******'
      IER= FSTRWD(15)
      IER = FSTSKP(15, 6)
      IER = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','UU')
       WRITE(6,*)' UU  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
      IER = FSTSKP(15, NRECSK)
      KEY = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','VV ')
       WRITE(6,*)' VV  APRES FSTSEL  '
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)
*------------------------------------------------------------------
      write(6,*)' relire le premier record'
      IER = FSTSKP(15, -(NRECSK+7+1))
      IER = FSTSEL(15,NIOUT,NJOUT,NKOUT,-1,' ',-1,-1,-1,' ','UU')
      ILUK = FSTLIS(F,15,NIOUT,NJOUT,NKOUT)


      WRITE(6,*)' CALCUL POSITION DU REC AVEC FSTPOS '
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 10.0 (FSTVOI) LISTER FICHIERS 20-10-15 ---'
      IER= FSTRWD(20)
      IVOI = FSTVOI(20,'SEQ/FTN')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(20, SEQ)'
      IDBG =QSTDBG('TEST 10.0', 1)
      ENDIF
      IER= FSTRWD(10)
      IVOI = FSTVOI(10,'RND')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(10, RND)'
      IDBG =QSTDBG('TEST 10.0', 1)
      ENDIF
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 10.0', 1)
      CALL ERREUR
      STOP
      ENDIF
      IDBG =QSTDBG('TEST 10.0', 1)
      CALL TESTOK
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)
     % '---TEST 10.1(FSTECR) ECRIRE 10 RECORDS FILE 10-15-20 ---'
      WRITE(6,*)' FSTWEO - FSTVOI SUR FICHIER 15     '
      DEBREC = 1
      FINREC = 10
      INCREC = 1
*                   DATYP=0 BINAIRE  DATYP=1 REEL  DATYP=2 ENTIER
*                                                  DATYP=3 CARACTERE
*
      REWRIT = .FALSE.
      IERR = FSTRWD(10)
      CALL TSTECR(10, REWRIT )
      IERR = FSTRWD(15)
      CALL TSTECR(15, REWRIT  )
      IERR = FSTRWD(20)
      CALL TSTECR(20, REWRIT   )
      IER= FSTRWD(15)
      IVOI = FSTVOI(15,'SEQ')
      IF((IVOI.LT.0 ))THEN
         WRITE(6,*)' ERREUR DANS FSTVOI(15, SEQ)'
      IDBG =QSTDBG('TEST 10.0', 1)
      CALL ERREUR
      STOP
      ENDIF


      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
      WRITE(6,*)' ##################################################'
*
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*)' '
      WRITE(6,*) '--- TEST 11.0   (FSTEOF, FSTWEO)'
      IERR = FSTRWD(20)
      CALL STDTEOF(20)
      IERR = FSTRWD(15)
      CALL STDTEOF(15)
      IERR = FSTFRM(10)
      IERR = FSTFRM(15)
      IERR = FSTFRM(20)
*
*
5555  WRITE(6,*) ' *STDTEST REUSSI*'
      JDATE = EXFIN('STDTEST', 'OK', 'NON')
      STOP
      END
      SUBROUTINE TSTECR(IUN, REWRIT)
      IMPLICIT LOGICAL(A-Z)
**
      COMMON /TESTC1/PRNT, DEBREC, FINREC, INCREC,F(15,100)
      COMMON /TESTC2/ WK(10000)
      INTEGER DEBREC, FINREC, INCREC, ILIR,DATE,NI,NJ,NK,IP1,IP2,IP3
      INTEGER DATEV,I1,I2,I3,I4,I,J,IREC,DEET,NPAS,DATYP,DATNEG
      INTEGER NNI,NNJ,NNK,IECR,IUN,KEY,IEFF,FSTEFF,ICVT
      REAL F, WK
      REAL VAL
      EXTERNAL FSTECR,FSTLIR,INCDAT,QSTINT,FSTLUK,FSTSUI,FSTINF,FSTEFF
      INTEGER FSTECR,FSTLIR,QSTINT,INDESC,FSTLUK,ILUK,FSTSUI,FSTINF
      EXTERNAL QSTSUI,FSTCVT

      INTEGER QSTSUI,HNOM(2),HETIKET(2),HMETIQ,HTYPVAR,HGRTYP,FSTCVT

      LOGICAL ECR, REWRIT, LUK, EFF, LIR, GETSET, PRNT
      CHARACTER *2 CNOMVA1, CNOMVA2, CNOMVAR

      CHARACTER *8 CETIKET, METIQ

      CHARACTER *1 CTYPVAR,CGRTYP
*      SAVE HNOM, HTYPVAR, HGRTYP, HETIKET, GETSET
      SAVE
      DATA HNOM(1)/"UU"/,HNOM(2)/"VV"/,HTYPVAR/"X"/,HGRTYP/"G"/,

     %     HETIKET(1)/"ETIK"/,HETIKET(2)/"ET89"/

      ECR = .TRUE.
100   CONTINUE
      GETSET = .TRUE.

      ICVT = FSTCVT(HNOM(1),HTYPVAR,HETIKET(1),HGRTYP,CNOMVA1,CTYPVAR,
     %              CETIKET,CGRTYP, GETSET)

      CNOMVAR = CNOMVA1
      WRITE(6,*)' CTYPVAR,CNOMVA1,CETIKET,CGRTYP APRES FSTCVT'
      WRITE(6,681) CTYPVAR,CNOMVA1,CETIKET,CGRTYP

      ICVT = FSTCVT(HNOM(2),HTYPVAR,HETIKET(1),HGRTYP,CNOMVA2,CTYPVAR,
     %              CETIKET,CGRTYP, GETSET)

      WRITE(6,*)' CTYPVAR,CNOMVA2,CETIKET,CGRTYP APRES FSTCVT'
      WRITE(6,681) CTYPVAR,CNOMVA2,CETIKET,CGRTYP
  681 FORMAT(' CTYPVAR= ',A1,'  CNOMVAR= ',A2,' CETIKET= ',A8,
     %        ' CGRTYP= ',A1)
      DO 23000 IREC = DEBREC, FINREC, INCREC
         CNOMVAR = CNOMVA1
         IF(IREC .GT.10 .AND. IREC.LT.21 ) CNOMVAR = CNOMVA2
         IF(IREC .GT.30 .AND. IREC.LT.41 ) CNOMVAR = CNOMVA2
         IF(IREC .GT.50 .AND. IREC.LT.61 ) CNOMVAR = CNOMVA2
         IF(IREC .GT.70 .AND. IREC.LT.81 ) CNOMVAR = CNOMVA2
         IF(IREC .GT.90 .AND. IREC.LT.101) CNOMVAR = CNOMVA2
         DATE = 220589003
         DEET = 2
         NPAS = 900
         CALL INCDAT(DATEV, DATE, (DEET*NPAS + 1800)/3600)
         NI = 15
         NJ = IREC
*         NK = 0
         NK = 1
         IP1 = IREC
         IP2 = IREC**2
         IP3 = 3*IREC
         I1 = 1
         I2 = 2
         I3 = 3
         I4 = 4
         IF( (.NOT. ECR))THEN
                IF(EFF) THEN
                   DATNEG = -1
                   KEY = FSTINF(IUN,NNI,NNJ,NNK,DATNEG,CETIKET,
     %                    IP1,IP2,IP3,CTYPVAR,CNOMVAR)
                   IEFF = FSTEFF(KEY)
*                  IF(IEFF .GT. 0) WRITE(6,605)IUN,KEY
                    WRITE(6,605)IUN,KEY,IEFF
  605 FORMAT(' EFFACE - IUN=',I3,',KEY',I6,' EST EFFACE',
     %  ' IEFF DE IEFF=FSTEFF = ',I4)
                   IF((IEFF.LT. 0))THEN
                      CALL ERREUR
                      STOP
                   ENDIF
                ENDIF
                IF(LUK) THEN
                   KEY = FSTINF(IUN,NNI,NNJ,NNK,DATEV,CETIKET,
     %                    IP1,IP2,IP3,CTYPVAR,CNOMVAR)
                   WRITE(6,*)' VALEUR DE KEY SPECIAL= ',KEY
                   ILUK = FSTLUK(F, KEY, NNI, NNJ, NNK)
                   IF(ILUK .LT. 0) THEN
                      CALL ERREUR
                      STOP
                   ENDIF
                ENDIF
             IF(LIR) THEN
                ILIR = FSTLIR(F, IUN, NNI, NNJ, NNK, DATEV,CETIKET,
     %          IP1, IP2, IP3, CTYPVAR,CNOMVAR)
                IF( (NI .NE. NNI .OR. NJ .NE. NNJ))THEN
                   WRITE(6,11) NI, NNI, NJ, NNJ
11             FORMAT(' NI=',I5,' NNI=',I5,' NJ=',I5,' NNJ=',I5)
                ENDIF
             ENDIF
          ENDIF
         DO 23008 J = 1, NJ
            DO 23010 I = 1, NI
               VAL = IUN*1000000 + IREC*10000 + I*100 + J
               IF( (ECR))THEN
                  F(I,J) = VAL
                  IF(PRNT) THEN
                  IF(J.EQ.1.AND.I.EQ.1)THEN
                     WRITE(6,10) F(I,J), VAL
                  ENDIF
                  ENDIF
                  IF(PRNT) THEN
                     IF(J.EQ.NJ.AND.I.EQ.NI)THEN
                        WRITE(6,12) F(NI,NJ), VAL
                     ENDIF
                  ENDIF
               ELSE
                  IF(PRNT) THEN
                     IF(J.EQ.1.AND.I.EQ.1)THEN
                        WRITE(6,10) F(I,J), VAL
                     ENDIF
                     IF(J.EQ.NJ.AND.I.EQ.NI)THEN
                           WRITE(6,12) F(NI,NJ), VAL
                     ENDIF
                  ENDIF
                  IF(EFF) GO TO 12345
                  IF(ILIR.LT.0) GO TO 12345
                  IF(LUK) GO TO 12345
                  IF(PRNT) THEN
                     IF(J.EQ.1.AND.I.EQ.1)THEN
                        IF( (F(I,J) .NE. VAL))THEN
                           WRITE(6,10) F(I,J), VAL
10                         FORMAT(' F(1,1)=',F12.1,' VAL=',F12.1)
                           CALL ERREUR
                           STOP
                        ENDIF
                     ENDIF
                  ENDIF
                  IF(PRNT) THEN
                     IF(J.EQ.NJ.AND.I.EQ.NI)THEN
                        IF( (F(I,J) .NE. VAL))THEN
                           WRITE(6,12) F(NI,NJ), VAL
12                         FORMAT(' F(NI,NJ)=',F12.1,' VAL=',F12.1)
                           CALL ERREUR
                           STOP
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
12345         CONTINUE
23010       CONTINUE
23008    CONTINUE
         IF( (ECR))THEN
          DATYP =1
            IECR = FSTECR(F, WK, 0, IUN, DATE, DEET, NPAS, NI, NJ,
     %       NK,   IP1, IP2, IP3, CTYPVAR,CNOMVAR,CETIKET,CGRTYP,I1,
     %       I2,I3,I4, DATYP, REWRIT )
         ENDIF
23000 CONTINUE
      RETURN
      ENTRY TSTLIR(IUN)
      EFF = .FALSE.
      LUK = .FALSE.
      ECR = .FALSE.
      LIR = .TRUE.
      GO TO 100
      ENTRY TSTLUK(IUN)
      LUK = .TRUE.
      LIR = .FALSE.
      ECR = .FALSE.
      EFF = .FALSE.
      GO TO 100
      ENTRY TSTEFF(IUN)
      LUK = .FALSE.
      LIR = .FALSE.
      ECR = .FALSE.
      EFF = .TRUE.
      GO TO 100
      END
      SUBROUTINE TESTOK
*   TRAITE LE RESULTAT D'UN TEST


      WRITE(6,10)
10    FORMAT(' ---REUSSI',/)
      RETURN
      ENTRY ERREUR
      WRITE(6,*) '$$$$$$$$$$$$$$$$   ERREUR  $$$$$$$$$$$$$$$$$$'
      RETURN
      END


      SUBROUTINE STDTEOF(IUN)
      IMPLICIT LOGICAL(A-Z)
**
      COMMON /TESTC1/PRNT,  DEBREC, FINREC, INCREC,   F(15,100)
      COMMON /TESTC2/ WK(10000)
      INTEGER DEBREC, FINREC, INCREC, DATE,NI,NJ,NK,IP1,IP2,IP3,ISUI
      INTEGER DATEV,I1,I2,I3,I4,I,J,IREC,DEET,NPAS,DATYP,DATNEG
      INTEGER NNI,NNJ,NNK,IECR,IUN,IEOF,IWEO,NIOUT,NJOUT,NKOUT,KEY
      REAL F, WK
      REAL VAL
      EXTERNAL FSTECR,FSTLIR,INCDAT,FSTSUI,FSTINF,FSTWEO,FSTEOF
      EXTERNAL FSTFRM, FSTOUV, FSTVOI,QSTDBG,FSTRWD,FSTLUK,FSTLIS
      INTEGER FSTWEO,FSTEOF,FSTOUV,FSTFRM,IOUV,IFRM,IVOI,IDBG,IERR
      INTEGER FSTECR,FSTLIR,INDESC,ILUK,FSTSUI,FSTINF,NIVEAU(15)
      INTEGER FSTVOI, FSTRWD, LIR, FSTLUK, FSTLIS
*
      SAVE NIVEAU
      DATA NIVEAU/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
      DO 23000 IREC = 1 ,15, 1
         DATE = 220589003
         DEET = 2
         NPAS = 900
         CALL INCDAT(DATEV, DATE, (DEET*NPAS + 1800)/3600)
         NI = 15
         NJ = IREC
         NK = 0
         IP1 = IREC
         IP2 = IREC**2
         IP3 = 3*IREC
         I1 = 1
         I2 = 2
         I3 = 3
         I4 = 4
         DO 23008 J = 1, NJ
            DO 23010 I = 1, NI
               VAL = IUN*1000000 + IREC*10000 + I*100 + J
                  F(I,J) = VAL
23010       CONTINUE
23008    CONTINUE
          DATYP =1
            IECR = FSTECR(F, WK, 0, IUN, DATE, DEET, NPAS, NI, NJ,
     %       NK,   IP1, IP2, IP3, 'X', 'YY','ETIKET89','G',I1,
     %       I2,I3,I4, DATYP, .FALSE.)
*
            IECR = FSTECR(F, WK, 0, IUN, DATE, DEET, NPAS, NI, NJ,
     %       NK,   IP1+100, IP2, IP3, 'X', 'YY','ETIKET89','G',I1,
     %       I2,I3,I4, DATYP, .FALSE.)
*
            IECR = FSTECR(F, WK, 0, IUN, DATE, DEET, NPAS, NI, NJ,
     %       NK,   IP1+101, IP2, IP3, 'X', 'YY','ETIKET89','G',I1,
     %       I2,I3,I4, DATYP, .FALSE.)
*
          IWEO = FSTWEO(IUN,NIVEAU(IREC))
           WRITE(6,*)' NIVEAU= ',NIVEAU(IREC)
           WRITE(6,*)' VALEUR IWEO = FSTWEO...= ',IWEO
*
23000 CONTINUE
*
      IFRM= FSTFRM(IUN)
      WRITE(6,*)' VERIFIER CODE D ERREUR APRES FSTFRM(XX)'
      IF((IFRM.LT.0))THEN
         WRITE(6,*)' ERREUR DANS EXECUTION FSTFRM(XX)'
      ENDIF
      IF(IUN.EQ.20) THEN
         IOUV  = FSTOUV(IUN, 'SEQ/FTN' )
         IF((IOUV.LT.0))THEN
            WRITE(6,*)' ****  ERREUR FSTOUV(20,SEQ)  ****'
         ENDIF
      ENDIF
      IF(IUN.EQ.15) THEN
         IOUV  = FSTOUV(IUN, 'SEQ' )
         IF((IOUV.LT.0))THEN
            WRITE(6,*)' ****  ERREUR FSTOUV(15,SEQ)  ****'
         ENDIF
      ENDIF
*
      IERR=FSTRWD(IUN)
*
*
      DO 23999 IREC =1,60,1
         DATE = 220589003
         DEET = 2
         NPAS = 900
         CALL INCDAT(DATEV, DATE, (DEET*NPAS + 1800)/3600)
         NI = 15
         NJ = IREC
         NK = 0
         IP1 = IREC
         IP2 = IREC**2
         IP3 = 3*IREC
         I1 = 1
         I2 = 2
         I3 = 3
         I4 = 4
*
*        LIR = FSTLIR(F,IUN, NI, NJ, NK, -1, ' ', -1 , -1, -1,
*    %    ' ',' ')
      KEY = FSTINF(IUN,NIOUT,NJOUT,NKOUT,-1,' ',-1 ,-1 ,-1 ,' ',' ')
*     IPRM=FSTPRM(KEY,CDATE,CDEET,CNPAS,CNI,CNJ,CNK,CNBITS,CDATYP,
*    %            CIP1,CIP2,CIP3,CTYPVAR,CNOMVAR,CETIKET,CGRTYP,
*    %            CIG1,CIG2,CIG3,CIG4,CSWA,CLNG,CDLTF,CUBC,
*    %            EXTRA1,EXTRA2,EXTRA3)
*      LIR = FSTLIS(F, IUN,NIOUT,NJOUT,NKOUT)
*
         IEOF = FSTEOF(IUN)
         WRITE(6,*)' VALEUR DE FSTEOF 1 A 14  = ',IEOF
*
23999 CONTINUE
         IERR = FSTRWD(IUN)
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
         IVOI = FSTVOI(IUN,'SEQ/FTN')
*
      RETURN
      END


*
*??????????????????????????????????????????????????????????
*          stop
*???????????????????????????????????????????????????????????
*








