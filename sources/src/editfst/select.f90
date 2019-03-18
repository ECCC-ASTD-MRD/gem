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
!** S/R SELECT TRAITER LES DIRECTIVES DE L'USAGER
      SUBROUTINE SELECT
      use configuration
      IMPLICIT NONE 
  
!AUTEUR       YVON R. BOURASSA JUN 86
!REVISION 001   "  "      "    OCT 90 VERSION QLXINS
!         002   "  "      "    JUL 91 DIRECTIVE ZAP
!         003   "  "      "    DEC 91 QLXINX (NEC)
!         004   "  "      "    MAI 91 QLXINX (partout sauf CRAY)
!         005 Mario Lepine     Juil 2001 ip1 en valeurs reels    
!         006 Mario Lepine     Mars 2003 correction valeur hybrid
!         006 Mario Lepine     Nov  2005 remplacement de fstabt par qqexit
!         007 Michel Valin     Fev  2014 ajout de constantes pour les "kind" ip
!
!LANGUAGE FTN77
!
!MODULES
      EXTERNAL QLXINX
      EXTERNAL SAUTSEQ, STDCOPI, WEOFILE, SETPER, QLXINS, ZAP
      EXTERNAL EXCLURE, SEQCOPI, REWINDS, DESIRE, qqexit
      EXTERNAL CRITSUP, SAUTSQI, SQICOPI, READLX

      INTEGER  MOIN1, SAUTSEQ, STDCOPI, WEOFILE, SETPER, SAUTSQI, KERR
      INTEGER  MOIN2, EXCLURE, SEQCOPI, REWINDS, DESIRE, SQICOPI, ZAP
      INTEGER  MOIN3, MOIN4,   CRITSUP, DUMY

      DATA     MOIN1, MOIN2, MOIN3, MOIN4/ -1, -2, -3, -4/
      integer  M1000, M1001, M1002, M1003, M1004
      integer  M1005, M1006, M1010, M1017, M1021
      data     M1000, M1001, M1002, M1003, M1004 / -1000, -1001, -1002, -1003, -1004/
      data     M1005, M1006, M1010, M1017, M1021 / -1005, -1006, -1010, -1017, -1021/
!
!*

!     OPTIONS (cle = valeur)

      CALL QLXINS(DEBUG  , 'DEBUG'  , DUMY, 1, 1) 
      CALL QLXINS(DIAG   , 'DIAG'   , DUMY, 1, 1) 
      CALL QLXINS(ECR    , 'ECR'    , DUMY, 1, 1) 
      CALL QLXINS(EOF    , 'EOF'    , DUMY, 1, 1) 
      CALL QLXINS(MEOF   , 'FINSEQ' , DUMY, 1, 1) 
      CALL QLXINS(FIXD   , 'FIXDATE', DUMY, 1, 1) 
      CALL QLXINS(LIMITE , 'LIMITE' , DUMY, 1, 1) 
      CALL QLXINS(CEOF   , 'SAUVEOF', DUMY, 1, 1) 
      CALL QLXINS(SAUV   , 'SAUVDES', DUMY, 1, 1) 
      CALL QLXINS(VS     , 'VOIRS'  , DUMY, 1, 1) 
      CALL QLXINS(VD     , 'VOIRD'  , DUMY, 1, 1) 
      CALL QLXINS(VD     , 'VOIR'   , DUMY, 1, 1) 
  
!     FONCTIONS ( cle(.....) )

      CALL QLXINX(DESIRE , 'DESIRE',  NP, 107, 2) 
      CALL QLXINX(CRITSUP, 'CRITSUP', NP, 108, 2) 
      CALL QLXINX(EXCLURE, 'EXCLURE', NP, 107, 2) 
      CALL QLXINX(SEQCOPI, 'FTNCOPI', NP, 107, 2) 
      CALL QLXINX(SETPER,  'PERIODE', NP, 104, 2) 
      CALL QLXINX(REWINDS, 'REWINDS', NP, 102, 2) 
      CALL QLXINX(SAUTSEQ, 'SAUTFTN', NP, 103, 2) 
      CALL QLXINX(SAUTSEQ, 'SAUTSEQ', NP, 103, 2) 
      CALL QLXINX(SAUTSQI, 'SAUTSQI', NP, 103, 2) 
      CALL QLXINX(SEQCOPI, 'SEQCOPI', NP, 107, 2) 
      CALL QLXINX(SQICOPI, 'SQICOPI', NP, 107, 2) 
      CALL QLXINX(STDCOPI, 'STDCOPI', NP, 104, 2) 
      CALL QLXINX(WEOFILE, 'STDWEOF', NP, 103, 2) 
      CALL QLXINX(ZAP,     'ZAP',     NP, 107, 2)
  
!     MOTS CLES (CONSTANTES)

      CALL QLXINS(MOIN1  , 'TOUS'   , DUMY, 1, 0) 
      CALL QLXINS(MOIN2  , '@'      , DUMY, 1, 0) 
      CALL QLXINS(MOIN3  , 'DELTA'  , DUMY, 1, 0) 
      CALL QLXINS(MOIN4  , 'COMMUNE', DUMY, 1, 0) 
      CALL QLXINS(.TRUE. , 'OUI'    , DUMY, 1, 0) 
      CALL QLXINS(.FALSE., 'NON'    , DUMY, 1, 0) 

!     SYMBOLES POUR LES UNITES IP1/2/3 (CONSTANTES)

      CALL QLXINS(M1000  , 'METERS' , DUMY, 1, 0)
      CALL QLXINS(M1000  , 'METRES' , DUMY, 1, 0)
      CALL QLXINS(M1001  , 'SIGMA'  , DUMY, 1, 0) 
      CALL QLXINS(M1002  , 'MBAR'   , DUMY, 1, 0) 
      CALL QLXINS(M1003  , 'OTHER'  , DUMY, 1, 0)
      CALL QLXINS(M1003  , 'AUTRE'  , DUMY, 1, 0)
      CALL QLXINS(M1004  , 'MGND'   , DUMY, 1, 0)
      CALL QLXINS(M1004  , 'MSOL'   , DUMY, 1, 0)
      CALL QLXINS(M1005  , 'HYBRID' , DUMY, 1, 0)
      CALL QLXINS(M1006  , 'THETA'  , DUMY, 1, 0)
      CALL QLXINS(M1010  , 'HEURES' , DUMY, 1, 0)
      CALL QLXINS(M1010  , 'HOURS'  , DUMY, 1, 0)
      CALL QLXINS(M1017  , 'INDX'   , DUMY, 1, 0)
      CALL QLXINS(M1021  , 'MPRES'  , DUMY, 1, 0)
  
      CALL READLX(5, DUMY, KERR)  ! on appelle l'interprete READLX
  
      IF(DUMY .LT. 0) THEN
         WRITE(6,*)'  **************************************'
         WRITE(6,*)' *                                      *'
         WRITE(6,*)'*     ERREUR(S) DANS LES DIRECTIVES      *'
         WRITE(6,*)'*      ERROR(S) FOUND IN DIRECTIVES      *'
         WRITE(6,*)' *                                      *'
         WRITE(6,*)'  **************************************'
         CALL qqexit(66) 
      ENDIF
      RETURN
      END 
