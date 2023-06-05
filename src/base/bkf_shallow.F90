!##########################################################################
subroutine BKF_SHALLOW6(KLON, KLEV, PDTCONV, &
     &                  PPABS, PZZ, &
     &                  PT, PRV, PRC, PRI, PU, PV, PW, PDMSEDT, &
     &                  PCLOUD,PURC,PURI, &
     &                  PTTENS, PRVTENS, PRCTENS,PRITENS,PUTENS, PVTENS, &
     &                  UMFS, PCH1, PCH1TEN, &
     &                  PKSHAL, PSP, PWSTAR, PDXDY, PMRK2, PKDEEP)
   !##########################################################################

!KICE      => bkf_kice
!PTADJS    => shal_timeconv_sec
!OUVTRANSS => bkf_lshalm
!OCHTRANS  => bkf_lch1conv
!KCH1      => bkf_kcc

   !!**** Interface routine to the fast Meso-NH convection code developed for ECMWF/ARPEGE
   !!     having a structure typical for operational routines
   !!
   !!     Transformations necessary to call deep+ shallow code
   !!     - skip input vertical arrays/levels : bottom=1, top=KLEV
   !!     - transform specific humidities in mixing ratio
   !!
   !!
   !!    PURPOSE
   !!    -------
   !!      The routine interfaces the MNH convection code as developed for operational
   !!      forecast models like ECMWF/ARPEGE or HIRLAM with the typical Meso-NH array structure
   !!      Calls the deep and/or shallow convection routine
   !!
   !!
   !!**  METHOD
   !!    ------
   !!     Returns one tendency for shallow+deep convection but each part can
   !!     be activated/desactivated separately
   !!     For deep convection one can enable up to 3 additional ensemble members
   !!     - this substantially improves the smoothness of the scheme and reduces
   !!       allows for runs with different cloud radii (entrainment rates) and
   !!       reduces the arbitrariness inherent to convective trigger condition
   !!
   !!
   !!
   !!    EXTERNAL
   !!    --------
   !!    CONVECT_SHALLOW
   !!    SU_CONVPAR, SU_CONVPAR1
   !!    SUCST:   ECMWF/ARPEGE routine
   !!
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !!
   !!
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !!
   !!    MODIFICATIONS
   !!    -------------
   !!      Original    11/12/98
   !!      modified    20/03/2002 by P. Marquet : transformed for ARPEGE/Climat
   !!                             (tsmbkind.h, REAL_B, INTEGER_M, _JPRB, "& &",
   !!                              _ZERO_, _ONE_, _HALF_)
   !!      modified    11/04/O2 allow for ensemble of deep updrafts/downdrafts
   !!
   !!    REFERENCE
   !!    ---------
   !!    Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc., Vol 127, pp 869-886:
   !!           A mass flux convection scheme for regional and global models.
   !!
   !-------------------------------------------------------------------------------
   !Revisions in RPN
   ! 001   -PV-jan2015 - do not aggregate deep and shallow tendencies; add output of shallow tendencies
   ! 002   -PV-apr2015 - output cloud top and base heights for deep only
   ! 003   -JY-Jul2015 - separate shallow convection from convection.ftn90 (which has both deep and shallow convection in bechtold scheme) !
   ! note: kcltop and kclbas are non-flipped, must be flipped for use ; chem transport needs a bit of work to get tendency out

   !* DECLARATIONS

   use phy_status, only: phy_error_L
   use YOMCST
   use YOE_CONVPAR_SHAL
   use YOE_CONVPAREXT
   use cnv_options, only: bkf_tperts, bkf_rads, bkf_kice, bkf_lshalm, bkf_lch1conv, bkf_kch
   use ens_perturb, only: ens_nc2d

   implicit none
!!!#include <arch_specific.hf>

   !* Declarations of dummy arguments

   integer,                    intent(IN)   :: KLON   ! horizontal dimension
   integer,                    intent(IN)   :: KLEV   ! vertical dimension
   real,                       intent(IN)   :: PDTCONV! Interval of time between two
                                                      ! calls of the deep convection

   real, dimension(KLON,KLEV), intent(IN)   :: PT     ! grid scale T at time t  (K)
   real, dimension(KLON,KLEV), intent(IN)   :: PRV    ! grid scale water vapor  (kg/kg)
   real, dimension(KLON,KLEV), intent(IN)   :: PRC    ! grid scale r_c (kg/kg)
   real, dimension(KLON,KLEV), intent(IN)   :: PRI    ! grid scale r_i (kg/kg)
   real, dimension(KLON,KLEV), intent(IN)   :: PU     ! grid scale horiz. wind u (m/s)
   real, dimension(KLON,KLEV), intent(IN)   :: PV     ! grid scale horiz. wind v (m/s)
   real, dimension(KLON,KLEV), intent(IN)   :: PW     ! grid scale vertical velocity (m/s)
   real, dimension(KLON,KLEV), intent(IN)   :: PDMSEDT! tendency of grid-scale moist static energy (m^2/s^3)
   real, dimension(KLON,KLEV), intent(IN)   :: PPABS  ! grid scale pressure (Pa)
   real, dimension(KLON,KLEV), intent(IN)   :: PZZ    ! geopotential (m2/s2)

   real, dimension(KLON,KLEV), intent(OUT)  :: PCLOUD   ! cloud fraction shallow
   real, dimension(KLON,KLEV), intent(OUT)  :: PURC     ! grid-scale liquid cond (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT)  :: PURI     ! grid-scale ice cond (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT)  :: PRCTENS  ! convective r_c tendency (1/s)  shallow
   real, dimension(KLON,KLEV), intent(OUT)  :: PRITENS  ! convective r_i tendency (1/s)  shallow
   real, dimension(KLON,KLEV), intent(INOUT):: PTTENS   ! convective temperat. tendency (K/s) shallow
   real, dimension(KLON,KLEV), intent(INOUT):: PRVTENS  ! convective r_v tendency (1/s)  shallow
   real, dimension(KLON,KLEV), intent(INOUT):: PUTENS   ! convecctive u tendency (m/s^2) shallow
   real, dimension(KLON,KLEV), intent(INOUT):: PVTENS   ! convecctive v tendency (m/s^2) shallow
   real, dimension(KLON,KLEV), intent(INOUT):: UMFS     ! shallow updraft mass flux
   real, dimension(KLON),      intent(IN)   :: PDXDY    ! grid cell area
   real, dimension(KLON,ens_nc2d), intent(IN) :: PMRK2  ! Markov chain for SPP
   real, dimension(KLON),      intent(IN)   :: PWSTAR   ! convective velocity scale (m/s)

   ! transport for chemical tracer
   real, dimension(KLON,KLEV,bkf_kch), intent(IN)   :: PCH1     ! grid scale chemical species
   real, dimension(KLON,KLEV,bkf_kch), intent(INOUT):: PCH1TEN  ! chemical convective tendency
   ! (1/s)

   ! Diagnostic variables

   !integer, DIMENSION(KLON),   INTENT(INOUT) :: KCLTOP ! cloud top level (number of model level)
   !integer, DIMENSION(KLON),   INTENT(INOUT) :: KCLBAS ! cloud base level(number of model level)
   ! they are given a value of
   ! 0 if no convection

   real, dimension(KLON),      intent(OUT)  :: PKSHAL ! shallow convective counter
   real, dimension(KLON),      intent(IN)   :: PKDEEP ! deep convective counter

   real, dimension(KLON),      intent(IN)   :: PSP    ! surface pressure

   !* Declarations of local variables :

   integer  :: JI, JK, JKP, JN, JKE, JKR, JKT  ! loop index

   ! Local arrays (upside/down) necessary for change of ECMWF arrays to convection arrays
   ! increase one extra level for Bechtold scheme
   real, dimension(KLON,KLEV+1) :: ZT     ! grid scale T at time t  (K)
   real, dimension(KLON,KLEV+1) :: ZRV    ! grid scale water vapor  (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZRC    ! grid scale r_c mixing ratio (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZRI    ! grid scale r_i mixing ratio (kg/kg)
   real, dimension(KLON,KLEV+1) :: ZU     ! grid scale horiz. wind u (m/s)
   real, dimension(KLON,KLEV+1) :: ZV     ! grid scale horiz. wind v (m/s)
   real, dimension(KLON,KLEV+1) :: ZW     ! grid scale vertical velocity (m/s)
   real, dimension(KLON,KLEV+1) :: ZPABS  ! grid scale pressure (Pa)
   real, dimension(KLON,KLEV+1) :: ZZZ    ! height of model layer (m)

   real, dimension(KLON,KLEV+1,bkf_kch):: ZCH1     ! grid scale chemical species

   real, dimension(KLON) :: &
        ZCRAD, &   ! radius at LCL
        ZDTPERT    ! temp. perturbation at LCL

   ! special for shallow convection

   real, dimension(KLON,KLEV+1) :: &
        ZCLOUD, ZTTENS, ZUTENS, ZVTENS, ZRVTENS, ZRCTENS, ZRITENS, &
        ZUMFS, ZURVS, ZURC, ZURI, ZDMSEDT, ZUDRS

   real, dimension(KLON,KLEV+1,bkf_kch) :: ZCH1TENS

   integer, dimension(KLON) :: ICLBASS, ICLTOPS

   integer                  :: KLEV1  ! vertical dimension (add sfc as first level)
   integer                  :: ITEST  ! Number of column where deep is not active

   logical, dimension(KLON) :: GTRIG  ! 2D logical mask for trigger test

   !* Declarations for additional shal. conv. calls
   integer                  :: KENS      ! number of additional shal. conv. calls
   integer                  :: KTPERTS   ! number of temperature perturbations
                                         ! for the trigger f-tion
   integer                  :: KRADS     ! number of cloud radii at LCL

   integer, parameter       :: KTPERTSX=10 ! max. allowed number of temp. perturbs.
   integer, parameter       :: KRADSX=10   ! max. allowed number of radii

   real, dimension(:,:,:), pointer :: &
        & ZCLOUDE, ZTTENSE, ZRVTENSE, ZRCTENSE, ZRITENSE, &
        & ZUMFSE, ZURVSE, ZURCE, ZURIE, &
        & ZUTENSE, ZVTENSE
   real, dimension(:,:), pointer     :: ZKSHALE
   integer, dimension(:,:), pointer  :: ICLBASSE, ICLTOPSE
   real, dimension(:,:,:,:), pointer :: ZCH1TENSE
   real, dimension(:,:,:), pointer   :: ZUDRSE

   !* Setup fundamental thermodunamical/physical constants using ECMWF/ARPEGE routine

   call SUCST(54,20020211,0,0)
   call SU_CONVPAR()
   call SU_CONVPAR_SHAL()

   !* Additional ensemble array Allocatinon

   ! Determine number of additional shallow convection calls (KENS)
   !#TODO: KENS = ... should be done once in cnv_nml_post
   KTPERTS = 1
   KRADS   = 1
   if (BKF_TPERTS(3) >= 1.E-4) &
        KTPERTS = max(1, min(KTPERTSX, &
        floor(1.01 + (BKF_TPERTS(2)- BKF_TPERTS(1)) / BKF_TPERTS(3)) &
        ))
   if (BKF_RADS(3) >= 1.E0) &
        KRADS   = max(1, min(KRADSX, &
        floor(1.01 + (BKF_RADS(2)  - BKF_RADS(1))   / BKF_RADS(3)) &
        ))
   KENS = KTPERTS * KRADS - 1

   if (KENS > 0) then
      !#TODO: avoid alloc/dealloc every step every slice in the physics
      allocate( ICLBASSE(KLON,KENS) )
      allocate( ICLTOPSE(KLON,KENS) )
      allocate( ZKSHALE(KLON,KENS) )
      allocate( ZCLOUDE(KLON,KLEV+1,KENS) )
      allocate( ZTTENSE(KLON,KLEV+1,KENS) )
      allocate( ZUTENSE(KLON,KLEV+1,KENS) )
      allocate( ZVTENSE(KLON,KLEV+1,KENS) )
      allocate( ZRVTENSE(KLON,KLEV+1,KENS) )
      allocate( ZRCTENSE(KLON,KLEV+1,KENS) )
      allocate( ZRITENSE(KLON,KLEV+1,KENS) )
      allocate( ZCH1TENSE(KLON,KLEV+1,bkf_kch,KENS) )
      allocate( ZUMFSE(KLON,KLEV+1,KENS) )
      allocate( ZURCE(KLON,KLEV+1,KENS) )
      allocate( ZURIE(KLON,KLEV+1,KENS) )
      allocate( ZURVSE(KLON,KLEV+1,KENS) )
      allocate( ZUDRSE(KLON,KLEV+1,KENS) )

      ICLTOPSE(:,:) = 1 ! set default value when no convection
      ICLBASSE(:,:) = 1
      ZKSHALE = 0.
   endif

   !    KCLTOP(:)  = 1 ! set default value when no convection
   !    KCLBAS(:)  = 1 ! can be changed  depending on user

   PCLOUD = 0.0
   PURC   = 0.0
   PURI   = 0.0

   !* Flip arrays upside-down as  first vertical level in convection is 1

   ! Move whole column one level up, and the first level value is surface value
   do JK = 1, KLEV
      JKP = KLEV - JK + 2
      do JI = 1, KLON
         ZPABS(JI,JKP) = PPABS(JI,JK)
         ZZZ(JI,JKP)   = PZZ(JI,JK)
         ZT(JI,JKP)    = PT(JI,JK)
         ZRV(JI,JKP)   = PRV(JI,JK) / ( 1. - PRV(JI,JK) ) ! transform specific humidity
         ZRC(JI,JKP)   = PRC(JI,JK) / ( 1. - PRC(JI,JK) ) ! in mixing ratio
         ZRI(JI,JKP)   = PRI(JI,JK) / ( 1. - PRI(JI,JK) )
         ZU(JI,JKP)    = PU(JI,JK)
         ZV(JI,JKP)    = PV(JI,JK)
         ZW(JI,JKP)    = PW(JI,JK)
         ZDMSEDT(JI,JKP) = PDMSEDT(JI,JK)
      end do
   end do

   ! Set first level value as surface values (not used)
   do JI = 1, KLON
      ZPABS(JI,1) = PSP(JI)
      ZZZ(JI,1) = 0.
      ZT(JI,1)  = ZT(JI,2)
      ZRV(JI,1) = ZRV(JI,2)
      ZRC(JI,1) = 0.
      ZRI(JI,1) = 0.
      ZU(JI,1)  = 0.
      ZV(JI,1)  = 0.
      ZW(JI,1)  = 0.
      ZDMSEDT(JI,1) = ZDMSEDT(JI,2)
   end do

   if (bkf_lch1conv) then
      do JK = 1, KLEV
         JKP = KLEV - JK + 2
         do JN = 1, bkf_kch
            do JI = 1, KLON
               ZCH1(JI,JKP,JN) = PCH1(JI,JK,JN)
            end do
         end do
      end do
   end if

   !* Call shallow convection routine

   ZCRAD(:)   = BKF_RADS(1)
   ZDTPERT(:) = BKF_TPERTS(1)

   !* Base version

   KLEV1 = KLEV+1
   GTRIG = (PKDEEP(:) <= 0.)
   ITEST = count(GTRIG)
   call CONVECT_SHALLOW6(KLON, KLEV1, ITEST, PDTCONV, &
        &                ZPABS, ZZZ,                                 &
        &                ZT, ZRV, ZRC, ZRI, ZDMSEDT,             &
        &                ZTTENS, ZRVTENS, ZRCTENS, ZRITENS,          &
        &                ICLTOPS, ICLBASS, ZUMFS, ZURVS,             &
        &                ZCLOUD,ZURC,ZURI,                           &
        &                ZCH1, ZCH1TENS,             &
        &                ZUDRS, PWSTAR, ZCRAD, ZDTPERT,              &
        &                PDXDY, PMRK2, PKSHAL, GTRIG,                &
        &                ZU, ZV, ZUTENS, ZVTENS)
   if (phy_error_L) return

   !* Calculations for additional ensemble members (if activated)
   ENSEMBLE_CASE: if (KENS > 0) then

      JKE = 0

      do JKR = 1, KRADS
         ZCRAD(:) = BKF_RADS(1) + (JKR-1)*BKF_RADS(3)

         do JKT = 1, KTPERTS
            ZDTPERT(:) = BKF_TPERTS(1) + (JKT-1)*BKF_TPERTS(3)

            if (JKT==1 .and. JKR==1) cycle

            JKE = JKE + 1
            GTRIG = (PKDEEP(:) <= 0.)
            ITEST = count(GTRIG)
            call CONVECT_SHALLOW6(KLON, KLEV1, ITEST, PDTCONV, &
                 &                ZPABS, ZZZ, &
                 &                ZT, ZRV, ZRC, ZRI, ZDMSEDT, &
                 &                ZTTENSE(:,:,JKE), ZRVTENSE(:,:,JKE), ZRCTENSE(:,:,JKE), ZRITENSE(:,:,JKE), &
                 &                ICLTOPSE(:,JKE), ICLBASSE(:,JKE), ZUMFSE(:,:,JKE), ZURVSE(:,:,JKE),  &
                 &                ZCLOUDE(:,:,JKE),ZURCE(:,:,JKE),ZURIE(:,:,JKE),  &
                 &                ZCH1, ZCH1TENSE(:,:,:,JKE),  &
                 &                ZUDRSE(:,:,JKE), PWSTAR, ZCRAD, ZDTPERT, &
                 &                PDXDY, PMRK2, ZKSHALE(:,JKE), GTRIG,  &
                 &                ZU, ZV, ZUTENSE(:,:,JKE), ZVTENSE(:,:,JKE))
            if (phy_error_L) return
            !#TODO: why not accumulate directly here instead of below... this way we would avoid (array size * KENS)
         enddo
      enddo

      !* Add  - if activated - ensemble average values for shallow
      !         convective tendencies

      DOJKE: do JKE = 1, KENS ! sum over additional clouds
         do JK = 1, KLEV1
            do JI = 1, KLON
               ZTTENS(JI,JK)  = ZTTENS(JI,JK)  + ZTTENSE(JI,JK,JKE)
               ZRVTENS(JI,JK) = ZRVTENS(JI,JK) + ZRVTENSE(JI,JK,JKE)
               ZRCTENS(JI,JK) = ZRCTENS(JI,JK) + ZRCTENSE(JI,JK,JKE)
               ZRITENS(JI,JK) = ZRITENS(JI,JK) + ZRITENSE(JI,JK,JKE)
               ZUMFS(JI,JK)   = ZUMFS(JI,JK)   + ZUMFSE(JI,JK,JKE)
               ZURVS(JI,JK)   = ZURVS(JI,JK)   + ZURVSE(JI,JK,JKE)
               ZURC(JI,JK)    = ZURC(JI,JK)    + ZURCE(JI,JK,JKE)
               ZURI(JI,JK)    = ZURI(JI,JK)    + ZURIE(JI,JK,JKE)
               ZCLOUD(JI,JK)  = ZCLOUD(JI,JK)  + ZCLOUDE(JI,JK,JKE)

               ZUDRS(JI,JK)   = ZUDRS(JI,JK)   + ZUDRSE(JI,JK,JKE) ! for ERA40
            end do
         end do
         do JI = 1, KLON
            ICLTOPS(JI)   = max(ICLTOPS(JI), ICLTOPSE(JI,JKE))
            ICLBASS(JI)   = max(ICLBASS(JI), ICLBASSE(JI,JKE))
         end do

         if (bkf_lshalm) then
            do JK = 1, KLEV1
               do JI = 1, KLON
                  ZUTENS(JI,JK) = ZUTENS(JI,JK) + ZUTENSE(JI,JK,JKE)
                  ZVTENS(JI,JK) = ZVTENS(JI,JK) + ZVTENSE(JI,JK,JKE)
               end do
            end do
         end if

         if (bkf_lch1conv)  then
            do JK = 1, KLEV1
               do JI = 1, KLON
                  ZCH1TENS(JI,JK,:) = ZCH1TENS(JI,JK,:) + ZCH1TENSE(JI,JK,:,JKE)
               end do
            end do
         end if

      end do DOJKE

      do JI = 1, KLON
         PKSHAL(JI)    = PKSHAL(JI) + sum( ZKSHALE(JI,:) )
      end do


      do JK = 1, KLEV1
         do JI = 1, KLON
            ZTTENS(JI,JK)  = ZTTENS(JI,JK)  / (KENS+1)
            ZRVTENS(JI,JK) = ZRVTENS(JI,JK) / (KENS+1)
            ZRCTENS(JI,JK) = ZRCTENS(JI,JK) / (KENS+1)
            ZRITENS(JI,JK) = ZRITENS(JI,JK) / (KENS+1)
            ZUMFS(JI,JK)   = ZUMFS(JI,JK)   / (KENS+1)
            ZURVS(JI,JK)   = ZURVS(JI,JK)   / (KENS+1)
            ZURC(JI,JK)    = ZURC(JI,JK)    / (KENS+1)
            ZURI(JI,JK)    = ZURI(JI,JK)    / (KENS+1)
            ZCLOUD(JI,JK)  = ZCLOUD(JI,JK)  / (KENS+1)
         end do
      end do

      if (bkf_lshalm) then
         do JK = 1, KLEV1
            do JI = 1, KLON
               ZUTENS(JI,JK) = ZUTENS(JI,JK) / (KENS+1)
               ZVTENS(JI,JK) = ZVTENS(JI,JK) / (KENS+1)
            end do
         end do
      end if

      if (bkf_lch1conv) then
         do JK = 1, KLEV1
            do JI = 1, KLON
               ZCH1TENS(JI,JK,:) = ZCH1TENS(JI,JK,:) / (KENS+1)
            end do
         end do
      end if

   end if ENSEMBLE_CASE


   !* Reflip arrays to ECMWF/ARPEGE vertical structure
   !  change mixing ratios to specific humidity

   ! shift the whole column one level down
   do JK = 1, KLEV
      JKP = KLEV - JK + 2
      do JI = 1, KLON
         PTTENS(JI,JK) = ZTTENS(JI,JKP)
         PRVTENS(JI,JK) = ZRVTENS(JI,JKP) / ( 1. + ZRV(JI,JKP) ) ** 2
         PRCTENS(JI,JK) = ZRCTENS(JI,JKP) / ( 1. + ZRC(JI,JKP) ) ** 2
         PRITENS(JI,JK) = ZRITENS(JI,JKP) / ( 1. + ZRI(JI,JKP) ) ** 2
         PUTENS(JI,JK) = ZUTENS(JI,JKP)
         PVTENS(JI,JK) = ZVTENS(JI,JKP)
         UMFS(JI,JK) = ZUMFS(JI,JKP)
         PCLOUD(JI,JK) = ZCLOUD(JI,JKP)
         PURC(JI,JK) = ZURC(JI,JKP)
         PURI(JI,JK) = ZURI(JI,JKP)
      end do
   end do

   if (bkf_lch1conv) then
      do JK = 1, KLEV
         JKP = KLEV - JK + 2
         do JN = 1, bkf_kch
            do JI = 1, KLON
               PCH1TEN(JI,JK,JN) = ZCH1TENS(JI,JKP,JN)
            end do
         end do
      end do
   end if

   !* Deallocate local arrays

   if (KENS > 0) then
      deallocate( ICLBASSE )
      deallocate( ICLTOPSE )
      deallocate( ZUMFSE )
      deallocate( ZURVSE )
      deallocate( ZURCE )
      deallocate( ZURIE )
      deallocate( ZCH1TENSE )
      deallocate( ZRCTENSE )
      deallocate( ZRITENSE )
      deallocate( ZCLOUDE )
      deallocate( ZKSHALE )
      deallocate( ZTTENSE )
      deallocate( ZUTENSE )
      deallocate( ZVTENSE )
      deallocate( ZRVTENSE )
      deallocate( ZUDRSE )
   end if

   return
end subroutine BKF_SHALLOW6
