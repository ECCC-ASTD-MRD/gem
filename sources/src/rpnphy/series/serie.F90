!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

subroutine serie3(INPUNIT,SERSTD,STATUS,ECHO,COMPRESS,HOUR64, &
     F_NIG,F_NK_HYBM,F_NK_HYBT)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_print,vgd_write,VGD_OK
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   !@object Processeur de series temporelles des modeles SEF, MC2 ET GEM
   !@Arguments
   !          - Input -
   ! INPUNIT  unit number attached to sequential input file
   ! SERSTD   unit number attached to standard RPN output file
   !          - Output -
   ! STATUS   =1 if no records read
   !          =2 if no surface variables
   !          =4 if too much space is required
   !          - Input -
   ! ECHO     .TRUE. to write out debug statements
   !          .FALSE. to not write debug statements
   ! COMPRESS .TRUE.    compressed output
   !          .FALSE. uncompressed output
   ! HOUR64   .TRUE. to write a 64-bit real (double) for the forecast hour
   !          .FALSE. to write a 32-bit real (float) for the forecast hour

   integer INPUNIT,SERSTD,STATUS,F_NIG,F_NK_HYBM, &
        F_NK_HYBT
   logical ECHO,COMPRESS,HOUR64

#include "series.cdk"

   integer, external :: FSTECR_S,INDSERI
   external :: SERPAT2

   !     AUTOMATIC ARRAYS
   integer, dimension(F_nig) :: ig
   real, dimension(F_nk_hybm) :: hybm
   real, dimension(F_nk_hybt) :: hybt

   !     VARIABLES ALLOUEES DYNAMIQUEMENT

   real(REAL64) , dimension(:)     , allocatable:: HH8
   real , dimension(:)       , allocatable:: HH, STOK, DIAG
   real , dimension(:,:)     , allocatable:: PHYSE
   real , dimension(:,:,:)   , allocatable:: SERS_T
   real , dimension(:,:,:,:) , allocatable:: SERP_T

   integer, parameter :: STDOUT=6

   integer TMP,AL,MG,Z0,HS,CND,LAT,LON,GL,MAP,SD,NIG
   integer NK_HYBM,NK_HYBT,NCOEF,ERR
   character(len=4) :: NOMVAR(NVAR), MODELE*3

   real DGRW,RGAS,GRAV
   character(len=12) :: ETIKET, ETIKMAJ
   integer DATE(14),NK,ITYP,IECR,size_vtbl(3)
   character(len=1) :: TV
   real PHS
   integer I,K,L,M,LP,NT,NREC,NSKIP,IP2,NPAK,NPHYE,DATYP
   logical SATUES, SATUCO

   real DEGRAD
   real(REAL64), dimension(:,:,:), pointer :: vtbl
   real  , dimension(:    ), pointer :: wkpt
   real(REAL64), dimension(:    ), pointer :: wkpt8

   type(VGRID_DESCRIPTOR) :: gd

   !     ---------------------------------------------------------------

   rewind (INPUNIT)

   NREC = 0

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

   read (inpunit,end=2) NSTAT, (NAME(L),IJSTAT(L,1),JSTAT(L),L=1,NSTAT), &
        NPHYE, (NOMVAR(M),M=1,NPHYE)                   , &
        NPROF, (PROFILS(M,1),M=1,NPROF)                , &
        (DATE(K),K=1,14),ETIKET,(IG(K),K=1,NIG)        , &
        DGRW, RGAS, GRAV, SATUES, SATUCO, TSMOYHR, SRWRI
   NREC = NREC + 1

   write ( STDOUT , * ) 'NSTAT = ',NSTAT
   write ( STDOUT , * ) 'IJSTAT = ',(IJSTAT(L,1),L=1,NSTAT)
   write ( STDOUT , * ) 'JSTAT = ',(JSTAT(L),L=1,NSTAT)
   write ( STDOUT , * ) 'NK = ',NK,' HYBRID(M) = ',(HYBM(K),K=1,NK_HYBM)
   write ( STDOUT , * ) 'DGRW = ',DGRW,' RGAS = ',RGAS,' GRAV = ',GRAV
   write ( STDOUT , '(1X,I3,20(1X,A4))' ) &
        NPHYE,(NOMVAR(M),M=1,NPHYE)
   write ( STDOUT , '(1X,I3,20(1X,A4))' ) &
        NPROF,(PROFILS(M,1),M=1,NPROF)
   write ( STDOUT , '(1X,A9,A12)' ) 'ETIKET = ', ETIKET
   write ( STDOUT , 6060 ) DATE
6060 format("  DATE-TIME GROUP     ",2X,6I6,6A4,A1,I12)

   write(6,612)SATUES,SATUCO
612 format(/,10X,'SATUES,SATUCO=',2(L1,2X),/)


   !--------------------------------------------------------------------

   allocate ( PHYSE(NSTAT,NPHYE) )

   read ( INPUNIT , end=2 ) HEURE,((PHYSE(L,M),L=1,NSTAT),M=1,NPHYE)
   NREC=NREC+1

   write ( STDOUT , * ) HEURE,((PHYSE(L,M),L=1,NSTAT),M=1,NPHYE)

   !     TROUVER LE NOMBRE DE PAS DE TEMPS SAUVES

1  read ( INPUNIT , end=2 )
   NREC = NREC + 1
   GO TO 1
2  NT = NREC - NSKIP

   !NT is the number of time-series records without the header
   !     write (stdout,*) 'serie: NT=',NT

   if ( NT.le.0 ) then
      STATUS = 1
      return
   endif

   !     POINTEURS DES CHAMPS PHYSIQUES

   if ( NPHYE .le. 0 ) then
      print *,'NPHYE = ',NPHYE
      STATUS = 2
      return
   else

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

   endif

   !     VERIFIER SI CERTAINS CHAMPS MANQUENT


   !     CHARGER L'ENTETE

   rewind(INPUNIT)
   do I=1,NSKIP-1
      read(INPUNIT)
   enddo
   !So read until we find the Second header which contains user defined VARS
   read (inpunit,end=2) NSTAT, (NAME(L),IJSTAT(L,1),JSTAT(L),L=1,NSTAT), &
        NSURF, (SURFACE(M,1),M=1,NSURF)                  , &
        NPROF, (PROFILS(M,1),M=1,NPROF)                , &
        (DATE(K),K=1,14),ETIKET,(IG(K),K=1,NIG)        , &
        DGRW, RGAS, GRAV, SATUES, SATUCO, TSMOYHR, SRWRI

   write ( STDOUT , * ) 'NSTAT = ',NSTAT
   write ( STDOUT , * ) 'IJSTAT = ',(IJSTAT(L,1),L=1,NSTAT)
   write ( STDOUT , * ) 'JSTAT = ',(JSTAT(L),L=1,NSTAT)
   write ( STDOUT , * ) 'NK = ',NK,' HYBRID(M) = ',(HYBM(K),K=1,NK_HYBM)
   write ( STDOUT , * ) 'TSMOYHR = ',TSMOYHR,' SRWRI = ',SRWRI
   write ( STDOUT , * ) 'DGRW = ',DGRW,' RGAS = ',RGAS,' GRAV = ',GRAV
   write ( STDOUT , '(1X,I3,20(1X,A4))' ) &
        NSURF,(SURFACE(M,1),M=1,NSURF)
   write ( STDOUT , '(1X,I3,20(1X,A4))' ) &
        NPROF,(PROFILS(M,1),M=1,NPROF)
   write ( STDOUT , '(1X,A9,A12)' ) 'ETIKET = ', ETIKET

   write ( STDOUT , 6060 ) DATE

   !     VERIFICATION DE L'ETIQUETTE

   !     CONVERSION DE L'ETIQUETTE EN MAJUSCULES
   call LOW2UP(ETIKET,ETIKMAJ)

   MODELE = 'GEM'

   !     ALLOUER LA MEMOIRE

   allocate ( HH8(NT)                   )
   allocate ( HH(NT)                    )
   allocate ( STOK(NT*NK)               )
   allocate ( DIAG(NT*NK)               )
   allocate ( SERS_T(NT,max(1,NSURF),NSTAT)    )
   allocate ( SERP_T(NK,NT,max(1,NPROF),NSTAT) )


   !     CHERCHER LES HEURES

   rewind INPUNIT
   !Rewind so that now we can read and fill up the hour array.
   do I=1,NSKIP
      read ( INPUNIT )
   enddo
   do I=1,NT
      read ( INPUNIT) HH8(I)
   enddo

   write ( STDOUT , * ) 'HEURES = ',(HH8(I),I=1,NT)


   !     AJOUTER L'HEURE INITIALE

   do I=1,NT
      HH8(I) =  HH8(I) + dble(DATE(5))
   enddo
   write ( STDOUT , * ) 'HEURES = ',(HH8(I),I=1,NT)
   HH = real(HH8)

   !     WRITE LIST OF FORECAST HOURS

   TV='T'
   NPAK=-32
   DATYP=5
   IP2 =FLOAT(DATE(5))*100
   if ( HOUR64 ) then
      IECR = FSTECR(HH8,STOK,-64,SERSTD,DATE(14),0,0,NT,1,1, &
           0,0,0,TV,'HH',ETIKET,'T',0,0,0,0,DATYP,.true.)
   else
      IECR = FSTECR(HH,STOK,NPAK,SERSTD,DATE(14),0,0,NT,1,1, &
           0,0,0,TV,'HH',ETIKET,'T',0,0,0,0,DATYP,.true.)
   endif

   !     WRITE LIST OF STATION NAMES (64-bit limit in FSTD)

   IECR = FSTECR_S(NAME,STOK,-8,SERSTD,DATE(14),0,0, &
        STN_STRING_LENGTH,NSTAT,1,0,0,0,TV, &
        'STNS',ETIKET,'T',0,0,0,0,7,.true.)

   status = vgd_write(gd,unit=serstd,format='fst')
   if (status /= VGD_OK) then
      write(STDOUT,*) 'WARNING: error writing grid descriptor: ',status
   endif

   !     OUTPUT A LIST OF HYBRID COORDINATES FOR PLOTTING


   ! N.B. Profils are strictly limited to nk-1 levels independently
   ! of what has been declared in the model run

   IECR = FSTECR(HYBT,STOK,NPAK,SERSTD,DATE(14),0,0,F_nk_hybt,1,1, &
        0,0,0,TV,'SH',ETIKET,'T',0,0,0,0,DATYP,.true.)
   IECR = FSTECR(HYBM,STOK,NPAK,SERSTD,DATE(14),0,0,F_nk_hybt,1,1, &
        0,0,0,TV,'SV',ETIKET,'T',0,0,0,0,DATYP,.true.)

   !-----------------------------------------------------------------

   DEGRAD = acos ( -1.0 )/180.0


   rewind INPUNIT
   ! Rewind and skip all the headers and finally, read each time-serie record
   ! and write out to standard file
   do K=1,NSKIP
      read ( INPUNIT )
   enddo

   do I=1,NT

      read  ( INPUNIT ) HEURE,((SERS(LP,M),LP=1,NSTAT),M=1,NSURF), &
           (((SERP(K,LP,M),K=1,NK),LP=1,NSTAT),M=1,NPROF)

      if (ECHO) then
         write (STDOUT,*) ' HEURE=',HEURE
         do L=1,NSTAT
            write (STDOUT,*) ' STATION NO.',L
            write (STDOUT,*) ' SURFACES=',(SERS(L,M),M=1,NSURF)
            do M=1,NPROF
               write (STDOUT,*) ' PROFIL NO.=',M,'  ',(SERP(K,L,M),K=1,NK)
            enddo
         end do
      endif

      !     TRANSFERER SERIES DE POINTS A TEMPS

      call SERPAT2 ( SERS_T , SERP_T , I , NT , NK , &
           NSURF, NPROF, NSTAT )

   enddo



   do L=1,NSTAT

      !     BOUCLE SUR LES POINTS (DO 100)
      !     L=POINT, I=TEMPS

      !     FINALISER LES SERIES EN TEMPS

      call serfin2(SERS_T(1,1,L) , SERP_T(1,1,1,L) , SURFACE , PROFILS , &
           NT , NK , NSURF , NPROF , &
           DGRW , PHYSE(L,LAT), PHYSE(L,LON), &
           DEGRAD, MODELE, IG)

      !     ACCOUNT FOR ACCUMULATORS IN CASE THEY WERE
      !     NOT RESET TO ZERO EVERY TIME SERIES OUTPUT TIME STEP

      if (TSMOYHR > 0 .and. SRWRI > 0) &
           call SERACC(SERS_T(1,1,L),HH8,NSURF,NT,SURFACE,TSMOYHR,SRWRI)

      !     ECRIRE SERIES POUR CE POINT SUR SERSTD

      call serecri4(SERS_T(1,1,L),SERP_T(1,1,1,L),SERSTD,NSURF, &
           NPROF,NT,SURFACE,PROFILS,L,PHYSE(L,LAT),PHYSE(L,LON), &
           STOK,DATE(14), ETIKET,FLOAT(DATE(5)), NK, &
           COMPRESS)


   enddo



   IECR = FSTECR (PHYSE(1,MG),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'MG',ETIKET,'Y',0,IP2,0,0,DATYP,.false.)

   !    CARACTERISTIQUES DU SOL DANS TEMPORAIRES POUR STOCKAGE

   do l=1,NSTAT

      ITYP=PHYSE(L,MG) + 0.5
      if (ITYP.le.0) then
         PHS=PHYSE(L,GL)
      else
         PHS=PHYSE(L,HS)
         if (PHS .lt. 0) PHS = 0
      endif
      PHYSE(L,HS)=PHS

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
      PHYSE(L,Z0)=max(0.,PHYSE(L,Z0))
   enddo


   !   ECRITURE DES ENREGISTREMENTS SUR serstd


   IECR = FSTECR (PHYSE(1,LAT),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'^^',ETIKET,'T',0,0,0,0,DATYP,.false.)

   IECR = FSTECR (PHYSE(1,LON),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'>>',ETIKET,'T',0,0,0,0,DATYP,.false.)

!!$      IECR = FSTECR (PHYSE(1,MG),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
!!$                 1,1,0,IP2,0,TV,'GS',ETIKET,'Y',0,IP2,0,0,DATYP,.FALSE.)

   IECR = FSTECR (PHYSE(1,HS),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'HS',ETIKET,'Y',0,IP2,0,0,DATYP,.false.)

   IECR = FSTECR (PHYSE(1,AL),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'AL',ETIKET,'Y',0,IP2,0,0,DATYP,.false.)

   IECR = FSTECR (PHYSE(1,TMP),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'TP',ETIKET,'Y',0,IP2,0,0,DATYP,.false.)

   IECR = FSTECR (PHYSE(1,Z0),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'Z0',ETIKET,'Y',0,IP2,0,0,DATYP,.false.)

   IECR = FSTECR (PHYSE(1,CND),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'PS',ETIKET,'Y',0,IP2,0,0,DATYP,.false.)

   IECR = FSTECR (PHYSE(1,SD),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'SD',ETIKET,'Y',0,IP2,0,0,DATYP,.false.)

   IECR = FSTECR (PHYSE(1,GL),STOK,NPAK,serstd,date(14),0,0,NSTAT, &
        1,1,0,IP2,0,TV,'GL',ETIKET,'Y',0,IP2,0,0,DATYP,.false.)


   !  ------------------------------------------------------

   !     DESALLOUER LA MEMOIRE

   deallocate (HH8,stat=iecr)
   deallocate (HH,stat=iecr)
   deallocate (STOK,stat=iecr)
   deallocate (DIAG,stat=iecr)
   deallocate (PHYSE,stat=iecr)
   deallocate (SERS_T,stat=iecr)
   deallocate (SERP_T,stat=iecr)

   return

end subroutine serie3
