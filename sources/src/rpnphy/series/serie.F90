!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------
!**S/P SERIE - PROCESSEUR DE SERIES TEMPORELLES DES MODELES SEF, MC2 ET GEM
      SUBROUTINE SERIE2( INPUNIT,SERSTD,STATUS,ECHO,COMPRESS,HOUR64, &
                         F_DIAG,F_NIG,F_NK_HYBM,F_NK_HYBT)
      use vGrid_Descriptors, only:vgrid_descriptor,vgd_new,vgd_get,vgd_print,vgd_write,VGD_OK
      implicit none
#include <arch_specific.hf>
      INTEGER INPUNIT,SERSTD,STATUS,F_NIG,F_NK_HYBM, &
              F_NK_HYBT
      LOGICAL ECHO,COMPRESS,HOUR64,F_DIAG

!Arguments
!
!          - Input -
! INPUNIT  unit number attached to sequential input file
! SERSTD   unit number attached to standard RPN output file
!
!          - Output -
! STATUS   =1 if no records read
!          =2 if no surface variables
!          =4 if too much space is required
!
!          - Input -
! ECHO     .TRUE. to write out debug statements
!          .FALSE. to not write debug statements
! COMPRESS .TRUE.    compressed output
!          .FALSE. uncompressed output
! HOUR64   .TRUE. to write a 64-bit real (double) for the forecast hour
!          .FALSE. to write a 32-bit real (float) for the forecast hour
! F_DIAG   .TRUE. to run diagnostics (NOT USED)
!          .FALSE. to shut off diagnostic calculation

#include "series.cdk"

      INTEGER,EXTERNAL :: FSTECR,FSTECR_S,INDSERI
      EXTERNAL :: SERPAT2, SERFIN, SERECRI3

!     AUTOMATIC ARRAYS
      integer, dimension(F_nig) :: ig
      real, dimension(F_nk_hybm) :: hybm
      real, dimension(F_nk_hybt) :: hybt
!
!     VARIABLES ALLOUEES DYNAMIQUEMENT
!
      REAL*8 , DIMENSION(:)     , ALLOCATABLE:: HH8
      REAL , DIMENSION(:)       , ALLOCATABLE:: HH, STOK, DIAG
      REAL , DIMENSION(:,:)     , ALLOCATABLE:: PHYSE
      REAL , DIMENSION(:,:,:)   , ALLOCATABLE:: SERS_T
      REAL , DIMENSION(:,:,:,:) , ALLOCATABLE:: SERP_T

      INTEGER, PARAMETER :: STDOUT=6

      INTEGER TMP,AL,MG,Z0,HS,CND,LAT,LON,GL,MAP,SD,NIG
      INTEGER NK_HYBM,NK_HYBT,NCOEF,ERR
      CHARACTER *4 NOMVAR(NVAR), MODELE*3

      REAL DGRW,RGAS,GRAV
      CHARACTER *12 ETIKET, ETIKMAJ
      INTEGER DATE(14),NK,ITYP,IECR,size_vtbl(3)
      CHARACTER *1 TV
      REAL PHS
      INTEGER I,K,L,M,LP,NT,NREC,NSKIP,IP2,NPAK,NPHYE,DATYP
      LOGICAL SATUES, SATUCO

      REAL DEGRAD
      real*8, dimension(:,:,:), pointer :: vtbl
      real  , dimension(:    ), pointer :: wkpt
      real*8, dimension(:    ), pointer :: wkpt8

      TYPE(VGRID_DESCRIPTOR) :: gd
!
!     ---------------------------------------------------------------
!
      REWIND (INPUNIT)
!
      NREC = 0
!
!     CHAMPS PHYSIQUES INVARIANTS

      read (inpunit,end=2)
      NREC = NREC+1
      read (inpunit,end=2) size_vtbl(1),size_vtbl(2),size_vtbl(3)
      NREC = NREC+1
      allocate (vtbl(size_vtbl(1),size_vtbl(2),size_vtbl(3)))
      read (inpunit,end=2) vtbl
      NREC = NREC+1
      if (vgd_new(gd,vtbl) /= vgd_ok) goto 2

      nk= F_nk_hybm; nk_hybm= F_nk_hybm ; nk_hybt= F_nk_hybt ; NCOEF= 2

      nullify(wkpt8,wkpt)
      err= vgd_get(gd,'VCDM - vertical coordinate (m)',wkpt)
      HYBM(1:F_nk_hybm)= wkpt(1:F_nk_hybm) ; deallocate(wkpt);nullify(wkpt)

      err= vgd_get(gd,'VCDT - vertical coordinate (m)',wkpt)
      HYBT(1:F_nk_hybt)= wkpt(1:F_nk_hybt) ; deallocate(wkpt);nullify(wkpt)

      NIG= 4; NSKIP = 9

!VLEE notes
! NSKIP is 9 because format of binary file is such that:
!ser_write1: nstat_g=         589 nsurf=          12 nprof=           0
!ser_write2:size of vtbl()           3         165           1
!ser_write3:write out vtbl
!ser_write4:header of initial physics static surface variables,prof vars
!ser_write5:data of initial physics static surface variables
!ser_write6: nstat_g=         589 nsurf= user requested surface vars nprof =
!                                        user requested profil vars
!ser_write7: size of vtbl()
!ser_write8: write out vtbl
!ser_write9:header of request physics surface variables
!ser_write10:1st write out of surface timeseries: hour, data
!ser_write11:2nd write out of surface timeseries: hour, data
! etc..

      READ (inpunit,end=2) NSTAT, (NAME(L),IJSTAT(L,1),JSTAT(L),L=1,NSTAT), &
                           NPHYE, (NOMVAR(M),M=1,NPHYE)                   , &
                           NPROF, (PROFILS(M,1),M=1,NPROF)                , &
                           (DATE(K),K=1,14),ETIKET,(IG(K),K=1,NIG)        , &
                           DGRW, RGAS, GRAV, SATUES, SATUCO, TSMOYHR, SRWRI
      NREC = NREC + 1

      WRITE ( STDOUT , * ) 'NSTAT = ',NSTAT
      WRITE ( STDOUT , * ) 'IJSTAT = ',(IJSTAT(L,1),L=1,NSTAT)
      WRITE ( STDOUT , * ) 'JSTAT = ',(JSTAT(L),L=1,NSTAT)
      WRITE ( STDOUT , * ) 'NK = ',NK,' HYBRID(M) = ',(HYBM(K),K=1,NK_HYBM)
      WRITE ( STDOUT , * ) 'DGRW = ',DGRW,' RGAS = ',RGAS,' GRAV = ',GRAV
      WRITE ( STDOUT , '(1X,I3,20(1X,A4))' ) &
                     NPHYE,(NOMVAR(M),M=1,NPHYE)
      WRITE ( STDOUT , '(1X,I3,20(1X,A4))' ) &
                     NPROF,(PROFILS(M,1),M=1,NPROF)
      WRITE ( STDOUT , '(1X,9HETIKET = ,A12)' ) ETIKET
      WRITE ( STDOUT , 6060 ) DATE
 6060 FORMAT(22H  DATE-TIME GROUP     ,2X,6I6,6A4,A1,I12)

      WRITE(6,612)SATUES,SATUCO
612   FORMAT(/,10X,'SATUES,SATUCO=',2(L1,2X),/)

!
!--------------------------------------------------------------------
!
      ALLOCATE ( PHYSE(NSTAT,NPHYE) )
!
      READ ( INPUNIT , END=2 ) HEURE,((PHYSE(L,M),L=1,NSTAT),M=1,NPHYE)
      NREC=NREC+1
!
      WRITE ( STDOUT , * ) HEURE,((PHYSE(L,M),L=1,NSTAT),M=1,NPHYE)
!
!     TROUVER LE NOMBRE DE PAS DE TEMPS SAUVES
!
    1 READ ( INPUNIT , END=2 )
      NREC = NREC + 1
      GO TO 1
    2 NT = NREC - NSKIP

!NT is the number of time-series records without the header
!     write (stdout,*) 'serie: NT=',NT
!
      IF ( NT.LE.0 ) THEN
         STATUS = 1
         RETURN
      ENDIF
!
!     POINTEURS DES CHAMPS PHYSIQUES
!
      IF ( NPHYE .LE. 0 ) THEN
         PRINT *,'NPHYE = ',NPHYE
         STATUS = 2
         RETURN
      ELSE
!
            MAP=indseri('MA',nomvar,nphye)
            LAT=indseri('LA',nomvar,nphye)
            LON=indseri('LO',nomvar,nphye)

!!$             NE=indseri('NE',nomvar,nphye)
             GL=indseri('GL',nomvar,nphye)
             MG=indseri('MG',nomvar,nphye)
             HS=indseri('HS',nomvar,nphye)
             Z0=indseri('ZP',nomvar,nphye)
             AL=indseri('AL',nomvar,nphye)
            TMP=indseri('TM',nomvar,nphye)
            CND=indseri('TP',nomvar,nphye)
             SD=indseri('SD',nomvar,nphye)
!
      ENDIF
!
!     VERIFIER SI CERTAINS CHAMPS MANQUENT
!
!
!     CHARGER L'ENTETE
!
      REWIND (INPUNIT)
      DO 5 I=1,NSKIP-1
 5        READ ( INPUNIT )
!So read until we find the Second header which contains user defined VARS
      READ (inpunit,end=2) NSTAT, (NAME(L),IJSTAT(L,1),JSTAT(L),L=1,NSTAT), &
                           NSURF, (SURFACE(M,1),M=1,NSURF)                  , &
                           NPROF, (PROFILS(M,1),M=1,NPROF)                , &
                           (DATE(K),K=1,14),ETIKET,(IG(K),K=1,NIG)        , &
                           DGRW, RGAS, GRAV, SATUES, SATUCO, TSMOYHR, SRWRI

      WRITE ( STDOUT , * ) 'NSTAT = ',NSTAT
      WRITE ( STDOUT , * ) 'IJSTAT = ',(IJSTAT(L,1),L=1,NSTAT)
      WRITE ( STDOUT , * ) 'JSTAT = ',(JSTAT(L),L=1,NSTAT)
      WRITE ( STDOUT , * ) 'NK = ',NK,' HYBRID(M) = ',(HYBM(K),K=1,NK_HYBM)
      WRITE ( STDOUT , * ) 'TSMOYHR = ',TSMOYHR,' SRWRI = ',SRWRI
      WRITE ( STDOUT , * ) 'DGRW = ',DGRW,' RGAS = ',RGAS,' GRAV = ',GRAV
      WRITE ( STDOUT , '(1X,I3,20(1X,A4))' ) &
                     NSURF,(SURFACE(M,1),M=1,NSURF)
      WRITE ( STDOUT , '(1X,I3,20(1X,A4))' ) &
                     NPROF,(PROFILS(M,1),M=1,NPROF)
      WRITE ( STDOUT , '(1X,9HETIKET = ,A12)' ) ETIKET
!
      WRITE ( STDOUT , 6060 ) DATE
!
!     VERIFICATION DE L'ETIQUETTE
!
!     CONVERSION DE L'ETIQUETTE EN MAJUSCULES
      CALL LOW2UP(ETIKET,ETIKMAJ)
!
      MODELE = 'GEM'
      !
!     ALLOUER LA MEMOIRE
!
      ALLOCATE ( HH8(NT)                   )
      ALLOCATE ( HH(NT)                    )
      ALLOCATE ( STOK(NT*NK)               )
      ALLOCATE ( DIAG(NT*NK)               )
      ALLOCATE ( SERS_T(NT,max(1,NSURF),NSTAT)    )
      ALLOCATE ( SERP_T(NK,NT,max(1,NPROF),NSTAT) )
!
!
!     CHERCHER LES HEURES
!
      REWIND INPUNIT
!Rewind so that now we can read and fill up the hour array.
      DO 6 I=1,NSKIP
    6    READ ( INPUNIT )
      DO 7 I=1,NT
        READ ( INPUNIT) HH8(I)
!
 7    continue
      WRITE ( STDOUT , * ) 'HEURES = ',(HH8(I),I=1,NT)
!
!
!     AJOUTER L'HEURE INITIALE
!
      DO 8 I=1,NT
    8   HH8(I) =  HH8(I) + dble(DATE(5))
      WRITE ( STDOUT , * ) 'HEURES = ',(HH8(I),I=1,NT)
      HH = real(HH8)
!
!     WRITE LIST OF FORECAST HOURS
!
      TV='T'
      NPAK=-32
      DATYP=5
      IP2 =FLOAT(DATE(5))*100
      IF ( HOUR64 ) then
        IECR = FSTECR(HH8,STOK,-64,SERSTD,DATE(14),0,0,NT,1,1, &
                      0,0,0,TV,'HH',ETIKET,'T',0,0,0,0,DATYP,.TRUE.)
      ELSE
        IECR = FSTECR(HH,STOK,NPAK,SERSTD,DATE(14),0,0,NT,1,1, &
                      0,0,0,TV,'HH',ETIKET,'T',0,0,0,0,DATYP,.TRUE.)
      ENDIF
!
!     WRITE LIST OF STATION NAMES (64-bit limit in FSTD)
!
      IECR = FSTECR_S(NAME,STOK,-8,SERSTD,DATE(14),0,0, &
                      STN_STRING_LENGTH,NSTAT,1,0,0,0,TV, &
                      'STNS',ETIKET,'T',0,0,0,0,7,.TRUE.)
!
      status = vgd_write(gd,unit=serstd,format='fst')
      if (status /= VGD_OK) then
         write(STDOUT,*) 'WARNING: error writing grid descriptor: ',status
      endif
!
!     OUTPUT A LIST OF HYBRID COORDINATES FOR PLOTTING
!
!
! N.B. Profils are strictly limited to nk-1 levels independently
! of what has been declared in the model run
!
      IECR = FSTECR(HYBT,STOK,NPAK,SERSTD,DATE(14),0,0,F_nk_hybt,1,1, &
                    0,0,0,TV,'SH',ETIKET,'T',0,0,0,0,DATYP,.TRUE.)
      IECR = FSTECR(HYBM,STOK,NPAK,SERSTD,DATE(14),0,0,F_nk_hybt,1,1, &
                    0,0,0,TV,'SV',ETIKET,'T',0,0,0,0,DATYP,.TRUE.)
!
!-----------------------------------------------------------------
!
      DEGRAD = ACOS ( -1.0 )/180.0
!
!
      REWIND INPUNIT
! Rewind and skip all the headers and finally, read each time-serie record
! and write out to standard file
      DO 9 K=1,NSKIP
    9    READ ( INPUNIT )
!
      DO 10 I=1,NT
!
      READ  ( INPUNIT ) HEURE,((SERS(LP,M),LP=1,NSTAT),M=1,NSURF), &
                       (((SERP(K,LP,M),K=1,NK),LP=1,NSTAT),M=1,NPROF)
!
      IF (ECHO) THEN
         WRITE (STDOUT,*) ' HEURE=',HEURE
         DO L=1,NSTAT
         WRITE (STDOUT,*) ' STATION NO.',L
         WRITE (STDOUT,*) ' SURFACES=',(SERS(L,M),M=1,NSURF)
         DO 20 M=1,NPROF
20       WRITE (STDOUT,*) ' PROFIL NO.=',M,'  ',(SERP(K,L,M),K=1,NK)
         END DO
      ENDIF
!
!     TRANSFERER SERIES DE POINTS A TEMPS
!
      CALL SERPAT2 ( SERS_T , SERP_T , I , NT , NK , &
                     NSURF, NPROF, NSTAT )
!
!
10    CONTINUE
!
!
      DO 100 L=1,NSTAT
!
!     BOUCLE SUR LES POINTS (DO 100)
!     L=POINT, I=TEMPS
!
!     FINALISER LES SERIES EN TEMPS
!
      CALL SERFIN ( SERS_T(1,1,L) , SERP_T(1,1,1,L) , SURFACE , PROFILS , &
                    NT , NK , NSURF , NPROF , &
                    DGRW , PHYSE(L,MAP) , PHYSE(L,LAT), PHYSE(L,LON), &
                    DEGRAD, MODELE, IG)
!
!     ACCOUNT FOR ACCUMULATORS IN CASE THEY WERE
!     NOT RESET TO ZERO EVERY TIME SERIES OUTPUT TIME STEP
!
      IF (TSMOYHR > 0 .and. SRWRI > 0) &
         CALL SERACC(SERS_T(1,1,L),HH8,NSURF,NT,SURFACE,TSMOYHR,SRWRI)
!
!     ECRIRE SERIES POUR CE POINT SUR SERSTD
!
      CALL SERECRI3(SERS_T(1,1,L),SERP_T(1,1,1,L),SERSTD,NSURF, &
                  NPROF,NT,SURFACE,PROFILS,L,PHYSE(L,LAT),PHYSE(L,LON), &
                  STOK,DIAG,DATE(14), ETIKET,FLOAT(DATE(5)), NK, &
                  SATUES, SATUCO, COMPRESS)
!
!
!
100   CONTINUE
!
!
      IECR = FSTECR (PHYSE(1,MG),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'MG',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
!    CARACTERISTIQUES DU SOL DANS TEMPORAIRES POUR STOCKAGE
!
       DO 110 l=1,NSTAT
!
       ITYP=PHYSE(L,MG) + 0.5
         IF (ITYP.LE.0) THEN
           PHS=PHYSE(L,GL)
         ELSE
           PHS=PHYSE(L,HS)
           IF (PHS .LT. 0) PHS = 0
         ENDIF
       PHYSE(L,HS)=PHS
!
!!$       IF(ITYP.LE.0) THEN
!!$         IF(PHYSE(L,GL).LT.0.9) THEN
!!$           NTYPES=0
!!$         ELSE
!!$           NTYPES=2
!!$         ENDIF
!!$       ELSE IF(ITYP.GE.1) THEN
!!$         IF(PHYSE(L,NE).GT.0.) THEN
!!$            NTYPES=1
!!$         ELSE
!!$            NTYPES=-1
!!$         ENDIF
!!$       ENDIF
!!$        PHYSE(L,MG)=NTYPES
        PHYSE(L,TMP)=PHYSE(L,TMP)-273.15
        PHYSE(L,Z0)=MAX(0.,PHYSE(L,Z0))
!
  110  CONTINUE
!
!   ECRITURE DES ENREGISTREMENTS SUR serstd
!
!
      IECR = FSTECR (PHYSE(1,LAT),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'^^',ETIKET,'T',0,0,0,0,DATYP,.FALSE.)
!
      IECR = FSTECR (PHYSE(1,LON),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'>>',ETIKET,'T',0,0,0,0,DATYP,.FALSE.)
!
!!$      IECR = FSTECR (PHYSE(1,MG),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
!!$                 1,1,0,IP2,0,TV,'GS',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
      IECR = FSTECR (PHYSE(1,HS),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'HS',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
      IECR = FSTECR (PHYSE(1,AL),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'AL',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
      IECR = FSTECR (PHYSE(1,TMP),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'TP',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
      IECR = FSTECR (PHYSE(1,Z0),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'Z0',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
      IECR = FSTECR (PHYSE(1,CND),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'PS',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
      IECR = FSTECR (PHYSE(1,SD),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'SD',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
      IECR = FSTECR (PHYSE(1,GL),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
                 1,1,0,IP2,0,TV,'GL',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)
!
!
!  ------------------------------------------------------
!
!     DESALLOUER LA MEMOIRE
!
      DEALLOCATE (HH8,stat=iecr)
      DEALLOCATE (HH,stat=iecr)
      DEALLOCATE (STOK,stat=iecr)
      DEALLOCATE (DIAG,stat=iecr)
      DEALLOCATE (PHYSE,stat=iecr)
      DEALLOCATE (SERS_T,stat=iecr)
      DEALLOCATE (SERP_T,stat=iecr)
!
      RETURN

      END
