!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html

!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

subroutine convect_shallow6(KLON, KLEV, ITEST, PDTCONV, &
     &                      PPABST, PZZ, &
     &                      PTT, PRVT, PRCT, PRIT, PDMSEDT, &
     &                      PTTEN, PRVTEN, PRCTEN, PRITEN, &
     &                      KCLTOP, KCLBAS, PUMF, PURV, &
     &                      PCLOUD,PURCOUT,PURIOUT, &
     &                      PCH1, PCH1TEN, &
     &                      PUDR, PWSTAR, PCRAD, PDTPERT, &
     &                      PDXDY, PMRK2, PKSHAL, GTRIG, &
     &                      PUT, PVT, PUTEN, PVTEN)

!KICE      => bkf_kice
!PTADJS    => shal_timeconv_sec
!KCH1      => bkf_kcc
!OCH1CONV  => bkf_lch1conv
!OUVCONV   => bkf_lshalm

   !!**** Monitor routine to compute all convective tendencies by calls
   !!     of several subroutines.
   !!
   !!
   !!    PURPOSE
   !!    -------
   !!      The purpose of this routine is to determine the convective
   !!      tendencies. The routine first prepares all necessary grid-scale
   !!      variables. The final convective tendencies are then computed by
   !!      calls of different subroutines.
   !!
   !!
   !!**  METHOD
   !!    ------
   !!      We start by selecting convective columns in the model domain through
   !!      the call of routine TRIGGER_FUNCT. Then, we gather the grid scale
   !!      variables in convective arrays.
   !!      The updraft and downdraft computations are done level by level starting
   !!      at the  bottom and top of the domain, respectively.
   !!      All computations are done on MNH thermodynamic levels. The depth
   !!      of the current model layer k is defined by DP(k)=P(k-1)-P(k)
   !!
   !!
   !!
   !!    EXTERNAL
   !!    --------
   !!    CONVECT_TRIGGER_SHAL
   !!    CONVECT_SATMIXRATIO
   !!    CONVECT_UPDRAFT_SHAL
   !!        CONVECT_CONDENS
   !!        CONVECT_MIXING_FUNCT
   !!    CONVECT_CLOSURE_SHAL
   !!        CONVECT_CLOSURE_THRVLCL
   !!        CONVECT_CLOSURE_ADJUST_SHAL
   !!
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !!      Module YOMCST
   !!          RG                   ! gravity constant
   !!          RPI                  ! number Pi
   !!          RATM                 ! reference pressure
   !!          RD, RV               ! gaz  constants for dry air and water vapor
   !!          RCPD, RCPV           ! specific heat for dry air and water vapor
   !!          RALPW, RBETW, RGAMW  ! constants for water saturation pressure
   !!          RTT                  ! triple point temperature
   !!          RLVTT, RLSTT         ! vaporization, sublimation heat constant
   !!          RCW, RCS             ! specific heat for liquid water and ice
   !!
   !!      Module YOE_CONVPAREXT
   !!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
   !!
   !!      Module YOE_CONVPAR
   !!          XCRAD                ! cloud radius
   !!
   !!
   !!    REFERENCE
   !!    ---------
   !!
   !!      Bechtold, 1997 : Meso-NH scientific  documentation (31 pp)
   !!      Fritsch and Chappell, 1980, J. Atmos. Sci., Vol. 37, 1722-1761.
   !!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
   !!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
   !!
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !!
   !!    MODIFICATIONS
   !!    -------------
   !!      Original    26/03/96
   !!   Peter Bechtold 15/11/96 replace theta_il by enthalpy
   !!         "        10/12/98 changes for ARPEGE
   !---------------------------------------------------------------------------

   !* 0. DECLARATIONS
   !     ------------

   use phy_status, only: phy_error_L
   use YOMCST
   use YOE_CONVPAREXT
   use YOE_CONVPAR_SHAL
   use cnv_options, only: bkf_kice, bkf_lshalm, bkf_lch1conv, bkf_kch
   use ens_perturb, only: ens_nc2d
   use integrals

   implicit none
!!!#include <arch_specific.hf>

   !* 0.1 Declarations of dummy arguments :

   integer,                    intent(IN) :: KLON     ! horizontal dimension
   integer,                    intent(IN) :: KLEV     ! vertical dimension
   integer,                    intent(IN) :: ITEST    ! Number of column where deep is not active
   real,                       intent(IN) :: PDTCONV  ! Interval of time between two
                                                      ! calls of the deep convection
   real, dimension(KLON),      intent(IN) :: PWSTAR   ! convective velocity scale (m/s)
   real, dimension(KLON),      intent(IN) :: PCRAD    ! cloud radius at the LCL (m)
   real, dimension(KLON),      intent(IN) :: PDTPERT  ! temp. perturbation (K)
   real, dimension(KLON),      intent(IN) :: PDXDY    ! grid cell area
   real, dimension(KLON,ens_nc2d), intent(IN) :: PMRK2! Markov chain for SPP
   real, dimension(KLON,KLEV), intent(IN) :: PTT      ! grid scale temperature at (K)
   real, dimension(KLON,KLEV), intent(IN) :: PRVT     ! grid scale water vapor(kg/kg)
   real, dimension(KLON,KLEV), intent(IN) :: PRCT     ! grid scale r_c  (kg/kg)"
   real, dimension(KLON,KLEV), intent(IN) :: PRIT     ! grid scale r_i (kg/kg)"
   real, dimension(KLON,KLEV), intent(IN) :: PDMSEDT  ! tendency of grid scale moist static energy (m^2/s^3)
   real, dimension(KLON,KLEV), intent(IN) :: PPABST   ! grid scale pressure (Pa)
   real, dimension(KLON,KLEV), intent(IN) :: PZZ      ! height of model layer (m)

   real, dimension(KLON,KLEV), intent(INOUT):: PTTEN  ! convective temperature
                                                      ! tendency (K/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRVTEN ! convective r_v tendency (1/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRCTEN ! convective r_c tendency (1/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRITEN ! convective r_i tendency (1/s)
   integer, dimension(KLON),   intent(INOUT):: KCLTOP ! cloud top level
   integer, dimension(KLON),   intent(INOUT):: KCLBAS ! cloud base level
                                                      ! they are given a value of
                                                      ! 0 if no convection
   real, dimension(KLON,KLEV), intent(OUT)  :: PUMF   ! updraft mass flux (kg/s m2)  !RON -why inout?
   real, dimension(KLON,KLEV), intent(OUT)  :: PURV   ! updraft water vapor (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT)  :: PURCOUT  ! normalized mixing ratio of updraft cloud water (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT)  :: PURIOUT  ! normalized mixing ratio of updraft cloud ice   (kg/kg)

   real, dimension(KLON,KLEV,bkf_kch), intent(IN)  :: PCH1    ! grid scale chemical species
   real, dimension(KLON,KLEV,bkf_kch), intent(INOUT):: PCH1TEN ! species conv. tendency (1/s)

   ! for ERA40
   real, dimension(KLON,KLEV), intent(OUT) :: PUDR   ! updraft detrainment rate (kg/s m3)

   real, dimension(KLON),      intent(OUT) :: PKSHAL ! shallow convective counter
   logical, dimension(KLON),   intent(INOUT):: GTRIG ! 2D logical mask for trigger test

   real, dimension(KLON,KLEV), intent(IN)  :: PUT    ! grid scale horiz. wind u (m/s)
   real, dimension(KLON,KLEV), intent(IN)  :: PVT    ! grid scale horiz. wind v (m/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PUTEN  ! convective u tendency (m/s^2)
   real, dimension(KLON,KLEV), intent(INOUT):: PVTEN  ! convective v tendency (m/s^2)
   real, dimension(KLON,KLEV), intent(INOUT):: PCLOUD ! shallow cloud fraction (%)

   !* 0.2a Declarations of local variables :

   integer :: ICONV           ! number of convective columns
   integer :: JI                 ! horizontal loop index
   integer :: JK            ! vertical loop index
   real    :: ZEPS, ZEPSA, ZEPSB      ! R_d / R_v, R_v / R_d, RCPV / RCPD - ZEPSA
   real    :: ZCPORD, ZRDOCP          ! C_p/R_d,  R_d/C_p

   real,    dimension(KLON,KLEV)  :: ZTHT, ZSTHV, ZSTHES, ZSMSE  ! grid scale theta, theta_v, moist static energy

   integer, dimension(ITEST) :: &
        ISDPL, &  ! index for parcel departure level
        ISPBL, &  ! index for source layer top
        ISLCL     ! index for lifting condensation level

   real, dimension(ITEST) :: &
        ZSTHLCL, &   ! updraft theta at LCL
        ZSTLCL, &    ! updraft temp. at LCL
        ZSRVLCL, &   ! updraft rv at LCL
        ZSWLCL, &    ! updraft w at LCL
        ZSPLCL, &    ! LCL pressure
        ZSZLCL, &    ! LCL height
        ZSTHVELCL, & ! envir. theta_v at LCL
        ZSMSEELCL, & ! envir. MSE at LCL
        ZSDXDY       ! grid area (m^2)

   real, dimension(ITEST,KLEV) :: &
        ZZ, &      ! height of model layer (m)
        ZPRES, &   ! grid scale pressure
        ZTHEST, &  ! grid scale saturated theta_e
        ZTH, &     ! grid scale theta
        ZTHV, &    ! grid scale theta_v
        ZRV, &     ! grid scale water vapor (kg/kg)
        ZMSE       ! grid scale moist static energy (m^2/s^2)

   real, dimension(ITEST,ens_nc2d) :: &
        ZSMRK2     ! Markov chains for SPP
   
   logical, dimension(ITEST) :: GTRIG1  ! logical mask for convection
   real,    dimension(KLON)  :: tmp1d   ! work variables

   !---------------------------------------------------------------------------

   JCVEXB = 0 !#TODO: remove, part of YOE_CONVPAREXT
   JCVEXT = 0 !#TODO: remove, part of YOE_CONVPAREXT

   !* 0.7 Reset convective tendencies to zero if convective
   !      counter becomes negative
   !      -------------------------------------------------

   !PV always set to zero (in other words, refresh is imposed)
   !PV if .not.refresh option is desired, some coding necesssary for pkdeep

   PTTEN(:,:)  = 0.
   PRVTEN(:,:) = 0.
   PRCTEN(:,:) = 0.
   PRITEN(:,:) = 0.
   PUTEN(:,:)  = 0.
   PVTEN(:,:)  = 0.
   PUMF(:,:)   = 0.
   PURV(:,:)   = 0.
   PURCOUT(:,:)  = 0.
   PURIOUT(:,:)  = 0.
   PUDR(:,:)   = 0.
   PCLOUD (:,:) =0.

   KCLTOP(:)  = 1
   KCLBAS(:)  = 1
   PKSHAL(:)  = 0.

   if (ITEST == 0) return

   !PV - init of pch1ten should probably be revised if tracers are activated

   if (bkf_lch1conv) then
      do JK = 1, KLEV
         do JI = 1, bkf_kch
            where(GTRIG(:)) PCH1TEN(:,JK,JI) = 0.
         enddo
      enddo
   endif

   !* 1. Initialize  local variables
   !     ----------------------------

   ZEPS   = RD / RV
   ZEPSA  = RV / RD
   ZEPSB  = RCPV / RCPD - ZEPSA
   ZCPORD = RCPD / RD
   ZRDOCP = RD / RCPD

   !* 1.1  Set up grid scale theta, theta_v, theta_es
   !       ------------------------------------------

   do JK = 1, KLEV
      where (PPABST(:,JK) > 40.E2)
         ZTHT(:,JK)  = PTT(:,JK) * ( RATM / PPABST(:,JK) ) ** ZRDOCP
         ZSTHV(:,JK) = ZTHT(:,JK) * ( 1. + ZEPSA * PRVT(:,JK) ) / &
              & ( 1. + PRVT(:,JK) + PRCT(:,JK) + PRIT(:,JK) )

         ! use conservative Bolton (1980) formula for theta_e
         ! it is used to compute CAPE for undilute parcel ascent
         ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here

         tmp1d(:) = exp( RALPW - RBETW / PTT(:,JK) - RGAMW * log( PTT(:,JK) ) )
         tmp1d(:) = min( 1., ZEPS * tmp1d(:) / ( PPABST(:,JK) - tmp1d(:) ) )
         ZSTHES(:,JK) = PTT(:,JK) * ( ZTHT(:,JK) / PTT(:,JK) ) ** &
              &  ( 1. - 0.28 * tmp1d(:) ) &
              &    * exp( ( 3374.6525 / PTT(:,JK) - 2.5403 ) &
              &   * tmp1d(:) * ( 1. + 0.81 * tmp1d(:) ) )
         ZSMSE(:,JK) = RCPD*PTT(:,JK) + RG*PZZ(:,JK) + RLVTT*PRVT(:,JK)
      elsewhere
         ZTHT(:,JK)   = 300.
         ZSTHV(:,JK)  = 300.
         ZSTHES(:,JK) = 400.
         ZSMSE (:,JK) = 30000.
      end where
   enddo

   !* 2. Test for convective columns and determine properties at the LCL
   !     --------------------------------------------------------------

   !# Active Colomns gathering (where GTRIG)
   do JK = 1, KLEV
      ZPRES(1:ITEST,JK)  = pack(PPABST(:,JK), GTRIG(:))
      ZTH(1:ITEST,JK)    = pack(ZTHT(:,JK),   GTRIG(:))
      ZTHV(1:ITEST,JK)   = pack(ZSTHV(:,JK),  GTRIG(:))
      ZTHEST(1:ITEST,JK) = pack(ZSTHES(:,JK), GTRIG(:))
      ZMSE(1:ITEST,JK)   = pack(ZSMSE(:,JK),  GTRIG(:))
      ZRV(1:ITEST,JK)    = max(0., pack(PRVT(:,JK), GTRIG(:)))
      ZZ(1:ITEST,JK)     = pack(PZZ(:,JK),    GTRIG(:))
   enddo
   do JK = 1, ens_nc2d
      ZSMRK2(1:ITEST,JK) = pack(PMRK2(:,JK), GTRIG(:))
   enddo
   ZSDXDY(1:ITEST) = pack(PDXDY(:), GTRIG(:))
!!$   PDTPERT(1:ITEST) = pack(PDTPERT(:), GTRIG(:)) !#No need to be packed, constant=BKF_TPERTS(1)

   !* 2.3 Test for convective columns and determine properties at the LCL
   !      --------------------------------------------------------------

   ISLCL(:) = 2   ! initialize DPL PBL and LCL
   ISDPL(:) = 1
   ISPBL(:) = 1

   call CONVECT_TRIGGER_SHAL4(ITEST, KLEV, &
        &  ZPRES, ZTH, ZTHV, ZTHEST, ZMSE, &
        &  ZRV, ZZ, PDTPERT, &
        &  ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSPLCL, ZSZLCL, &
        &  ZSTHVELCL, ZSMSEELCL, ISLCL, ISDPL, ISPBL, GTRIG1)
   if (phy_error_L) return
   !#TODO: should be (unpacked) GTRIG1 = GTRIG1 .AND. GTRIG

   !* 3. After the call of TRIGGER_FUNCT we do calculus only in
   !     convective columns. This corresponds to a GATHER operation.
   !     --------------------------------------------------------------

   ICONV = count(GTRIG1(1:ITEST))
   if (ICONV == 0) return ! no convective column has been found, exit CONVECT_SHALLOW

   ! Continue convect_shallow work on only ICONV (active) points
   call convect_shallow_c(KLON, KLEV, ITEST, ICONV, PDTCONV, &
        &                 PPABST, PZZ, &
        &                 PTT, PRVT, PRCT, PRIT, PDMSEDT, &
        &                 PTTEN, PRVTEN, PRCTEN, PRITEN, &
        &                 KCLTOP, KCLBAS, PUMF, PURV, &
        &                 PCLOUD,PURCOUT,PURIOUT, &
        &                 PCH1, PCH1TEN, &
        &                 PUDR, PWSTAR, PCRAD, &
        &                 ZSDXDY, ZSMRK2, PKSHAL, GTRIG, GTRIG1, &
        &                 PUT, PVT, PUTEN, PVTEN, &
        &                 ZTHT, ZSTHV, ZSTHES, ZSMSE, &
        &                 ISDPL, ISPBL, ISLCL, &
        &                 ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSPLCL, &
        &                 ZSZLCL, ZSTHVELCL, ZSMSEELCL)

   return
end subroutine convect_shallow6


!############################################################################
subroutine convect_shallow_c(KLON, KLEV, ITEST, ICONV, PDTCONV, &
     &                      PPABST, PZZ, &
     &                      PTT, PRVT, PRCT, PRIT, PDMSEDT, &
     &                      PTTEN, PRVTEN, PRCTEN, PRITEN, &
     &                      KCLTOP, KCLBAS, PUMF, PURV, &
     &                      PCLOUD,PURCOUT,PURIOUT, &
     &                      PCH1, PCH1TEN, &
     &                      PUDR, PWSTAR, PCRAD, &
     &                      ZSDXDY, ZSMRK2, PKSHAL, GTRIG, GTRIG1, &
     &                      PUT, PVT, PUTEN, PVTEN, &
     &                      ZTHT, ZSTHV, ZSTHES, ZSMSE, &
     &                      ISDPL, ISPBL, ISLCL, &
     &                      ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSPLCL, &
     &                      ZSZLCL, ZSTHVELCL, ZSMSEELCL)
   !############################################################################

   !* 0. DECLARATIONS
   !     ------------

   use phy_status, only: phy_error_L
   use YOMCST
   use YOE_CONVPAREXT
   use YOE_CONVPAR_SHAL
   use cnv_options, only: bkf_kice, shal_timeconv_sec, bkf_lshalm, bkf_lch1conv, bkf_kch
   use ens_perturb, only: ens_nc2d
   use integrals

   implicit none
!!!#include <arch_specific.hf>

   !* 0.1 Declarations of dummy arguments :

   integer,                    intent(IN) :: KLON     ! horizontal dimension
   integer,                    intent(IN) :: KLEV     ! vertical dimension
   integer,                    intent(IN) :: ITEST, ICONV    ! Number of column where deep is not active
   real,                       intent(IN) :: PDTCONV  ! Interval of time between two
                                                      ! calls of the deep convection
   real, dimension(KLON),      intent(IN) :: PWSTAR   ! convective velocity scale (m/s)
   real, dimension(KLON),      intent(IN) :: PCRAD    ! cloud radius at the LCL (m)
   real, dimension(ITEST),     intent(IN) :: ZSDXDY   ! grid cell area
   real, dimension(ITEST,ens_nc2d), intent(IN) :: ZSMRK2 ! Markov chain for SPP
   real, dimension(KLON,KLEV), intent(IN) :: PTT      ! grid scale temperature at (K)
   real, dimension(KLON,KLEV), intent(IN) :: PRVT     ! grid scale water vapor(kg/kg)
   real, dimension(KLON,KLEV), intent(IN) :: PRCT     ! grid scale r_c  (kg/kg)"
   real, dimension(KLON,KLEV), intent(IN) :: PRIT     ! grid scale r_i (kg/kg)"
   real, dimension(KLON,KLEV), intent(IN) :: PDMSEDT  ! tendency of grid scale moist static energy (m^2/s^3)
   real, dimension(KLON,KLEV), intent(IN) :: PPABST   ! grid scale pressure (Pa)
   real, dimension(KLON,KLEV), intent(IN) :: PZZ      ! height of model layer (m)

   real, dimension(KLON,KLEV), intent(INOUT):: PTTEN  ! convective temperature
                                                      ! tendency (K/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRVTEN ! convective r_v tendency (1/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRCTEN ! convective r_c tendency (1/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PRITEN ! convective r_i tendency (1/s)
   integer, dimension(KLON),   intent(INOUT):: KCLTOP ! cloud top level
   integer, dimension(KLON),   intent(INOUT):: KCLBAS ! cloud base level
                                                      ! they are given a value of
                                                      ! 0 if no convection
   real, dimension(KLON,KLEV), intent(OUT)  :: PUMF   ! updraft mass flux (kg/s m2)  !RON -why inout?
   real, dimension(KLON,KLEV), intent(OUT)  :: PURV   ! updraft water vapor (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT)  :: PURCOUT  ! normalized mixing ratio of updraft cloud water (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT)  :: PURIOUT  ! normalized mixing ratio of updraft cloud ice   (kg/kg)

   real, dimension(KLON,KLEV,bkf_kch), intent(IN)  :: PCH1    ! grid scale chemical species
   real, dimension(KLON,KLEV,bkf_kch), intent(INOUT):: PCH1TEN ! species conv. tendency (1/s)

   ! for ERA40
   real, dimension(KLON,KLEV), intent(OUT) :: PUDR   ! updraft detrainment rate (kg/s m3)

   real, dimension(KLON),      intent(OUT) :: PKSHAL ! shallow convective counter
   logical, dimension(KLON),   intent(INOUT):: GTRIG ! 2D logical mask for trigger test
   logical, dimension(ITEST),  intent(INOUT):: GTRIG1 ! 2D logical mask for trigger test

   real, dimension(KLON,KLEV), intent(IN)  :: PUT    ! grid scale horiz. wind u (m/s)
   real, dimension(KLON,KLEV), intent(IN)  :: PVT    ! grid scale horiz. wind v (m/s)
   real, dimension(KLON,KLEV), intent(INOUT):: PUTEN  ! convective u tendency (m/s^2)
   real, dimension(KLON,KLEV), intent(INOUT):: PVTEN  ! convective v tendency (m/s^2)
   real, dimension(KLON,KLEV), intent(INOUT):: PCLOUD ! shallow cloud fraction (%)

   real, dimension(KLON,KLEV), intent(INOUT):: ZTHT, ZSTHV, ZSTHES, ZSMSE  ! grid scale theta, theta_v, moist static energy

   integer, dimension(ITEST), intent(INOUT) :: &
        ISDPL, &  ! index for parcel departure level
        ISPBL, &  ! index for source layer top
        ISLCL     ! index for lifting condensation level

   real, dimension(ITEST), intent(INOUT) :: &
        ZSTHLCL, &   ! updraft theta at LCL
        ZSTLCL, &    ! updraft temp. at LCL
        ZSRVLCL, &   ! updraft rv at LCL
        ZSWLCL, &    ! updraft w at LCL
        ZSPLCL, &    ! LCL pressure
        ZSZLCL, &    ! LCL height
        ZSTHVELCL, & ! envir. theta_v at LCL
        ZSMSEELCL    ! envir. MSE at LCL

   !* 0.2 Declarations of local variables :

   integer :: ICONV1                  ! number of convective columns
   integer :: JI, JL                  ! horizontal loop index
   integer :: JN                      ! number of tracers
   integer :: JK, JKP, JKM            ! vertical loop index
   integer :: IFTSTEPS                ! only used for chemical tracers
   real    :: ZEPS                    ! R_d / R_v
   real    :: ZRDOCP                  ! R_d/C_p

   real, dimension(KLON) :: ZWORK2, ZWORK2B ! work array

   integer, dimension(ICONV) :: &
        IDPL, &   ! index for parcel departure level
        IPBL, &   ! index for source layer top
        ILCL, &   ! index for lifting condensation level
        IETL, &   ! index for zero buoyancy level
        ICTL, &   ! index for cloud top level
        ILFS      ! index for level of free sink

   ! grid scale variables
   real, dimension(ICONV,KLEV) :: &
        ZZ, &      ! height of model layer (m)
        ZPRES, &   ! grid scale pressure
        ZDPRES, &  ! pressure difference between
        ZTT, &     ! temperature
        ZTH, &     ! grid scale theta
        ZTHV, &    ! grid scale theta_v
        ZTHL, &    ! grid scale enthalpy (J/kg)
        ZTHES, &
        ZRW, &     ! grid scale total water (kg/kg)
        ZRV, &     ! grid scale water vapor (kg/kg)
        ZRC, &     ! grid scale cloud water (kg/kg)
        ZRI, &     ! grid scale cloud ice (kg/kg)
        ZDMSEDT, & ! tendency of grid scale moist static energy (m^2/s^3)
        ZMSE       ! grid scale moist static energy (m^2/s^2)

   real, dimension(ICONV,ens_nc2d) :: &
        ZMRK2      ! Markov chains for SPP
   
   real, dimension(ICONV) :: &
        ZDXDY, &   ! grid area (m^2)
        ZWSTAR, &  ! convective velocity scale (m/s)
        ZCRAD      ! cloud radius at the LCL (m)

   ! updraft variables
   real, dimension(ICONV,KLEV) :: &
        ZUMF, &    ! updraft mass flux (kg/s)
        ZUER, &    ! updraft entrainment (kg/s)
        ZUDR, &    ! updraft detrainment (kg/s)
        ZUTHL, &   ! updraft enthalpy (J/kg)
        ZUTHV, &   ! updraft theta_v (K)
        ZWU, &     ! updraft vertical velocity (m/s)
        ZURW, &    ! updraft total water (kg/kg)
        ZURC, &    ! updraft cloud water (kg/kg)
        ZURI       ! updraft cloud ice   (kg/kg)

   real, dimension(ICONV) :: &
        ZMFLCL, &  ! cloud base unit mass flux(kg/s)
        ZCAPE, &   ! available potent. energy
        ZTHLCL, &  ! updraft theta at LCL
        ZTLCL, &   ! updraft temp. at LCL
        ZMSELCL, & ! updraft MSE at LCL
        ZRVLCL, &  ! updraft rv at LCL
        ZWLCL, &   ! updraft w at LCL
        ZPLCL, &   ! LCL pressure
        ZZLCL, &   ! LCL height
        ZTHVELCL, &! envir. theta_v at LCL
        ZMSEELCL   ! envir. MSE at LCL

   ! downdraft variables
   real, dimension(ICONV,KLEV) :: &
        ZDMF, &    ! downdraft mass flux (kg/s)
        ZDER, &    ! downdraft entrainment (kg/s)
        ZDDR       ! downdraft detrainment (kg/s)

   ! closure variables
   real, dimension(ICONV,KLEV) :: ZLMASS  ! mass of model layer (kg)
   real, dimension(ICONV) :: ZTIMEC  ! advective time period

   real, dimension(ICONV,KLEV) :: &
        ZTHC, &    ! conv. adj. grid scale theta
        ZRVC, &    ! conv. adj. grid scale r_w
        ZRCC, &    ! conv. adj. grid scale r_c
        ZRIC, &    ! conv. adj. grid scale r_i
        ZWSUB      ! envir. compensating subsidence (Pa/s)

   integer, dimension(KLON)  :: IINDEX  ! hor.index
   integer, dimension(ICONV) :: IJINDEX ! hor.index

   real :: ZCPH     ! specific heat C_ph
   real :: ZLV, ZLS ! latent heat of vaporis., sublim.
   real :: ZW1      ! work variables

   ! Chemical Tracers:
   real,  dimension(ICONV,KLEV,bkf_kch) :: & 
        ZCH1, &    ! grid scale chemical specy (kg/kg)
        ZCH1C      ! conv. adjust. chemical specy 1
   real,  dimension(ICONV,bkf_kch) :: ZWORK3  ! conv. adjust. chemical specy 1
                                              ! for U, V transport:

   real, dimension(ICONV,KLEV) :: &
        ZU, &      ! grid scale horiz. u component on theta grid
        ZV, &      ! grid scale horiz. v component on theta grid
        ZUC, &     ! horizontal wind u (m/s)
        ZVC        ! horizontal wind v (m/s)

   real,     dimension(ICONV) :: ZTIMC   ! adjust fractional convective time step
   integer,  dimension(ICONV) :: ITSTEP  ! fractional convective time step
   logical,  dimension(ICONV) :: GWORK1  ! flag for convectiton-activated column with MASS CONSERVED

   !* 1. Initialize  local variables
   !     ----------------------------

   ZEPS   = RD / RV
   ZRDOCP = RD / RCPD

   !* 3.1 Gather grid scale and updraft base variables in
   !      arrays using mask GTRIG
   !      ---------------------------------------------------

   GTRIG(1:KLON)    = unpack(GTRIG1(1:ITEST), GTRIG(1:KLON), FIELD=.false.)

   do JK = 1, KLEV
      ZZ(1:ICONV,JK)     = pack(PZZ(1:KLON,JK),    GTRIG(1:KLON))
      ZPRES(1:ICONV,JK)  = pack(PPABST(1:KLON,JK), GTRIG(1:KLON))
      ZTT(1:ICONV,JK)    = pack(PTT(1:KLON,JK),    GTRIG(1:KLON))
      ZTH(1:ICONV,JK)    = pack(ZTHT(1:KLON,JK),   GTRIG(1:KLON))
      ZTHES(1:ICONV,JK)  = pack(ZSTHES(1:KLON,JK), GTRIG(1:KLON))
      ZRV(1:ICONV,JK)    = max(0., pack(PRVT(1:KLON,JK), GTRIG(1:KLON)))
      ZRC(1:ICONV,JK)    = max(0., pack(PRCT(1:KLON,JK), GTRIG(1:KLON)))
      ZRI(1:ICONV,JK)    = max(0., pack(PRIT(1:KLON,JK), GTRIG(1:KLON)))
      ZTHV(1:ICONV,JK)   = pack(ZSTHV(1:KLON,JK),  GTRIG(1:KLON))
      ZU(1:ICONV,JK)     = pack(PUT(1:KLON,JK),    GTRIG(1:KLON))
      ZV(1:ICONV,JK)     = pack(PVT(1:KLON,JK),    GTRIG(1:KLON))
      ZMSE(1:ICONV,JK)   = pack(ZSMSE(1:KLON,JK),  GTRIG(1:KLON))
      ZDMSEDT(1:ICONV,JK)= pack(PDMSEDT(1:KLON,JK),GTRIG(1:KLON))
   enddo

   do JK = 1, ens_nc2d
      ZMRK2(1:ICONV,JK)     = pack(ZSMRK2(1:ITEST,JK), GTRIG1(1:ITEST))
   enddo
   
   IDPL(1:ICONV)      = pack(ISDPL(1:ITEST),   GTRIG1(1:ITEST))
   IPBL(1:ICONV)      = pack(ISPBL(1:ITEST),   GTRIG1(1:ITEST))
   ILCL(1:ICONV)      = pack(ISLCL(1:ITEST),   GTRIG1(1:ITEST))
   ZTHLCL(1:ICONV)    = pack(ZSTHLCL(1:ITEST), GTRIG1(1:ITEST))
   ZTLCL(1:ICONV)     = pack(ZSTLCL(1:ITEST),  GTRIG1(1:ITEST))
   ZRVLCL(1:ICONV)    = pack(ZSRVLCL(1:ITEST), GTRIG1(1:ITEST))
   ZWLCL(1:ICONV)     = pack(ZSWLCL(1:ITEST),  GTRIG1(1:ITEST))
   ZPLCL(1:ICONV)     = pack(ZSPLCL(1:ITEST),  GTRIG1(1:ITEST))
   ZZLCL(1:ICONV)     = pack(ZSZLCL(1:ITEST),  GTRIG1(1:ITEST))
   ZTHVELCL(1:ICONV)  = pack(ZSTHVELCL(1:ITEST), GTRIG1(1:ITEST))
   ZMSEELCL(1:ICONV)  = pack(ZSMSEELCL(1:ITEST), GTRIG1(1:ITEST))
   ZDXDY(1:ICONV)     = pack(ZSDXDY(1:ITEST), GTRIG1(1:ITEST))
   ZWSTAR(1:ICONV)    = pack(PWSTAR(1:ITEST), GTRIG1(1:ITEST))
   ZCRAD(1:ICONV)     = pack(PCRAD(1:ITEST),  GTRIG1(1:ITEST))

   GTRIG1(1:ICONV)    = pack(GTRIG1(1:ITEST), GTRIG1(1:ITEST))

   !* 3.2 Compute pressure difference
   !      ---------------------------------------------------

   ZDPRES(1:ICONV,1) = 0.
   do JK = 1 + 1, KLEV
      ZDPRES(1:ICONV,JK)  = ZPRES(1:ICONV,JK-1) - ZPRES(1:ICONV,JK)
   enddo

   !* 3.3 Compute environm. enthalpy, MSE and total water = r_v + r_i + r_c
   !      ----------------------------------------------------------

   do JK = 1, KLEV, 1
      do JI = 1, ICONV
         ZRW(JI,JK)  = ZRV(JI,JK) + ZRC(JI,JK) + ZRI(JI,JK)
         ZCPH        = RCPD + RCPV * ZRW(JI,JK)
         ZLV         = RLVTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ! compute L_v
         ZLS         = RLSTT + ( RCPV - RCS ) * ( ZTT(JI,JK) - RTT ) ! compute L_i
         ZTHL(JI,JK) = ZCPH * ZTT(JI,JK) + ( 1. + ZRW(JI,JK) ) * RG * ZZ(JI,JK) &
              & - ZLV * ZRC(JI,JK) - ZLS * ZRI(JI,JK)
      enddo
   enddo


   !* 4. Compute updraft properties
   !     ----------------------------

   !* 4.1 Set mass flux at LCL ( here a unit mass flux with w = 1 m/s )
   !      -------------------------------------------------------------

   ! Initial single-cloud updraft estimate
   do JI = 1, ICONV
      JK = ILCL(JI) - 1
      ZMFLCL(JI) = ZPRES(JI,JK) / ( RD * ZTT(JI,JK) * &
           & ( 1. + ZEPS * ZRVLCL(JI) ) ) * RPI * ZCRAD(JI) * ZCRAD(JI)
   enddo

!!$     ! Initial estimate following C. Jones
!!$     ZMFLCL = (0.03*MAX(ZWSTAR,XWSTARMIN))*ZDXDY

   call CONVECT_UPDRAFT_SHAL3(ICONV, KLEV, &
        & bkf_kice, ZPRES, ZDPRES, ZZ, ZTT, ZTHL, ZTHV, ZTHES, ZRV, ZRW, &
        & ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL, &
        & ZMFLCL, ZCRAD, GTRIG1, ILCL, IDPL, IPBL, &
        & ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZWU, &
        & ZURW, ZURC, ZURI, ZCAPE, ICTL, IETL)

   ! Compute derived variables at the LCL
   ZMSELCL(1:ICONV) = RCPD*ZTLCL(1:ICONV) + RG*ZZLCL(1:ICONV) + RLVTT*ZRVLCL(1:ICONV)  !MSE of updraft at LCL

   !* 4.2 In routine UPDRAFT GTRIG1 has been set to false when cloud
   !      thickness is smaller than 3 km
   !      -----------------------------------------------------------

   ICONV1 = count(GTRIG1(1:ICONV))

   IF_ICONV1: if ( ICONV1 > 0 )  then

      ! downdraft variables

      ZDMF(:,:) = 0.
      ZDER(:,:) = 0.
      ZDDR(:,:) = 0.
      ILFS(:)   = 1  !#TODO: only in CONVECT_CHEM_TRANSPORT, make it local to that s/r
      do JK = 1, KLEV
         ZLMASS(1:ICONV,JK)  = ZDXDY(1:ICONV) * ZDPRES(1:ICONV,JK) / RG  ! mass of model layer
      enddo
      ZLMASS(1:ICONV,1) = ZLMASS(1:ICONV,2)

      ! closure variables

      ZTIMEC(:) = shal_timeconv_sec

      !* 7. Determine adjusted environmental values assuming
      !     that all available buoyant energy must be removed
      !     within an advective time step ZTIMEC.
      !     ---------------------------------------------------
      call CONVECT_CLOSURE_SHAL5(ICONV, KLEV, &
           & ZPRES, ZDPRES, ZZ, ZDXDY, ZMRK2, ZCRAD, ZLMASS, &
           & ZTHL, ZTH, ZRW, ZRC, ZRI, ZDMSEDT, GTRIG1, &
           & ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB, &
           & ILCL, IDPL, IPBL, ICTL, &
           & ZPLCL, ZMSELCL, ZMSEELCL, &
           & ZUMF, ZUER, ZUDR, ZUTHL, ZURW, &
           & ZURC, ZURI, ZCAPE, ZTIMEC, PDTCONV, IFTSTEPS,&
           & ZTIMC,ITSTEP,GWORK1)
      if (phy_error_L) return

      !* 8. Determine the final grid-scale (environmental) convective
      !     tendencies and set convective counter
      !     --------------------------------------------------------

      !* 8.1 Grid scale tendencies
      !      ---------------------

      ! in order to save memory, the tendencies are temporarily stored
      ! in the tables for the adjusted grid-scale values

      do JK = 1, KLEV
         do JI = 1, ICONV
            ZTHC(JI,JK) = ( ZTHC(JI,JK) - ZTH(JI,JK) ) / ZTIMEC(JI) &
                 & * ( ZPRES(JI,JK) / RATM ) ** ZRDOCP  ! change theta in temperature
            ZRVC(JI,JK) = ( ZRVC(JI,JK) - ZRW(JI,JK) + ZRC(JI,JK) + ZRI(JI,JK) ) &
                 &                            / ZTIMEC(JI)

            ZRCC(JI,JK) = ( ZRCC(JI,JK) - ZRC(JI,JK) ) / ZTIMEC(JI)
            ZRIC(JI,JK) = ( ZRIC(JI,JK) - ZRI(JI,JK) ) / ZTIMEC(JI)
         enddo
      enddo

      !* 8.2 Apply conservation correction
      !      -----------------------------

      ! adjustment at cloud top to smooth discontinuous profiles at PBL inversions
      ! (+ - - tendencies for moisture )

      do JI = 1, ICONV
         JK = ICTL(JI)
         JKM= max(1,ICTL(JI)-1)
         JKP= max(1,ICTL(JI)-2)
         ZRVC(JI,JKM) = ZRVC(JI,JKM) + .5 * ZRVC(JI,JK)
         ZRCC(JI,JKM) = ZRCC(JI,JKM) + .5 * ZRCC(JI,JK)
         ZRIC(JI,JKM) = ZRIC(JI,JKM) + .5 * ZRIC(JI,JK)
         ZTHC(JI,JKM) = ZTHC(JI,JKM) + .5 * ZTHC(JI,JK)
         ZRVC(JI,JKP) = ZRVC(JI,JKP) + .3 * ZRVC(JI,JK)
         ZRCC(JI,JKP) = ZRCC(JI,JKP) + .3 * ZRCC(JI,JK)
         ZRIC(JI,JKP) = ZRIC(JI,JKP) + .3 * ZRIC(JI,JK)
         ZTHC(JI,JKP) = ZTHC(JI,JKP) + .3 * ZTHC(JI,JK)
         ZRVC(JI,JK)  = .2 * ZRVC(JI,JK)
         ZRCC(JI,JK)  = .2 * ZRCC(JI,JK)
         ZRIC(JI,JK)  = .2 * ZRIC(JI,JK)
         ZTHC(JI,JK)  = .2 * ZTHC(JI,JK)
      end do

      ! Compute moisture integral

      ZWORK2(:) = 0.
      ZWORK2B(:) = 0.
      do JI = 1, ICONV
         do JK = 2, ICTL(JI)
            JKP = JK + 1
            ZW1 = 0.5 *  (ZPRES(JI,max(JK-1,2)) - ZPRES(JI,JKP)) / RG
            ZWORK2(JI) = ZWORK2(JI) + ( ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK) ) *   & ! moisture
                 &                   ZW1
            ZWORK2B(JI) = ZWORK2B(JI) + ZW1
         enddo
      enddo

      ! Compute moisture correction factor (integrates to zero)
      do JI = 1, ICONV
         if ( ICTL(JI) > 2 ) then
            ZWORK2(JI) =  ZWORK2(JI) / ZWORK2B(JI)
         endif
      enddo
      ! Apply moisture correction
      do JI = 1, ICONV
         if ( ICTL(JI) > 2 ) then
            do JK = ICTL(JI), 2, -1
               ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)
            enddo
         endif
      enddo
      ! Compute enthalpy integral
      zwork2(:) = 0.
      zwork2b(:) = 0.
      do ji=1,iconv
         do jk=2,ICTL(JI)
            jkp = jk + 1
            ZW1 = 0.5 *  (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
            ZWORK2(JI) = ZWORK2(JI) + ( &
                 &  ( RCPD + RCPV * ZRW(JI,JK) )* ZTHC(JI,JK) &
                 &  + RCPV * ZTT(JI,JK) * (ZRVC(JI,JK)+ZRCC(JI,JK)+ZRIC(JI,JK)) &
                 &  - ( RLVTT + ( RCPV - RCW ) * ( ZTT(JI,JK) - RTT ) ) * ZRCC(JI,JK) &
                 &  - ( RLSTT + ( RCPV - RCS ) * ( ZTT(JI,JK) - RTT ) ) * ZRIC(JI,JK) &
                 &                        ) * ZW1
            ZWORK2B(JI) = ZWORK2B(JI) + ZW1
         enddo
      enddo
      ! Compute enthalpy correction factor (integrates to zero)
      do JI = 1, ICONV
         if ( ICTL(JI) > 2 ) then
            ZWORK2(JI)= -ZWORK2(JI) / ZWORK2B(JI)
         endif
      enddo
      ! Apply enthalpy correction
      do JI = 1, ICONV
         if ( ICTL(JI) > 2 ) then
            do JK = ICTL(JI), 2, -1
               ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2(JI) / ( RCPD + RCPV * ZRW(JI,JK) )
            enddo
         endif
      enddo

      ! extend tendencies to first model level

      ! DO JI = 1, ICONV
      !    ZWORK2(JI) = ZDPRES(JI,2) + ZDPRES(JI,3)
      !    ZTHC(JI,1)  = ZTHC(JI,2) * ZDPRES(JI,3)/ZWORK2(JI)
      !    ZTHC(JI,2)= ZTHC(JI,2) * ZDPRES(JI,2)/ZWORK2(JI)
      !    ZRVC(JI,1)  = ZRVC(JI,2) * ZDPRES(JI,3)/ZWORK2(JI)
      !    ZRVC(JI,2)= ZRVC(JI,2) * ZDPRES(JI,2)/ZWORK2(JI)
      ! ENDDO

      ! execute a "scatter"/unpack command to store the tendencies in
      ! the final 2D tables

      do JK = 1, KLEV
         PTTEN(1:KLON,JK)  = unpack(ZTHC(1:ICONV,JK), GTRIG(1:KLON), PTTEN(1:KLON,JK))
         PRVTEN(1:KLON,JK) = unpack(ZRVC(1:ICONV,JK), GTRIG(1:KLON), PRVTEN(1:KLON,JK))
         PRCTEN(1:KLON,JK) = unpack(ZRCC(1:ICONV,JK), GTRIG(1:KLON), PRCTEN(1:KLON,JK))
         PRITEN(1:KLON,JK) = unpack(ZRIC(1:ICONV,JK), GTRIG(1:KLON), PRITEN(1:KLON,JK))
         zwork2(1:ICONV)   = ZURW(:,JK) - ZURC(:,JK) - ZURI(:,JK)
         PURV(1:KLON,JK)   = unpack(zwork2(1:ICONV),  GTRIG(1:KLON), PURV(1:KLON,JK))
      enddo

      ! Cloud base and top levels
      ! -------------------------

      ILCL(1:ICONV) = min( ILCL(1:ICONV), ICTL(1:ICONV) )
      zwork2(1:ICONV) = 0.
      where (GWORK1(1:ICONV)) zwork2(1:ICONV) = 1.
      PKSHAL(1:KLON) = unpack(zwork2(1:ICONV), GTRIG(1:KLON), field=0.)
      KCLTOP(1:KLON) = unpack(ICTL(1:ICONV),   GTRIG(1:KLON), KCLTOP(1:KLON))
      KCLBAS(1:KLON) = unpack(ILCL(1:ICONV),   GTRIG(1:KLON), KCLBAS(1:KLON))

      !* 8.7 Compute convective tendencies for Tracers
      !      ------------------------------------------

      IF_BKF_LCH1CONV: if (bkf_lch1conv) then

         do JN = 1, bkf_kch
            do JK = 1, KLEV
               ZCH1(1:ICONV,JK,JN) = pack(PCH1(1:KLON,JK,JN), GTRIG(1:KLON))
            enddo
         enddo

         call CONVECT_CHEM_TRANSPORT1(ICONV, KLEV, bkf_kch, ZCH1, ZCH1C, &
              &  IDPL, IPBL, ILCL, ICTL, ILFS, ILFS, &
              &  ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR, &
              &  ZTIMEC, ZDXDY, ZDMF(:,1), ZLMASS, ZWSUB)

         do JN = 1, bkf_kch
            do JK = 1, KLEV
               ZCH1C(1:ICONV,JK,JN) = ( ZCH1C(1:ICONV,JK,JN)- ZCH1(1:ICONV,JK,JN) ) / ZTIMEC(1:ICONV)
            enddo
         enddo

         !* 8.8 Apply conservation correction
         !      -----------------------------

         ! Compute vertical integrals

         ZWORK3(:,:) = 0.
         do JI = 1, ICONV
            do JK = 1, ICTL(JI)+1
               JKP = JK + 1
               ZWORK3(JI,:) = ZWORK3(JI,:) + ZCH1C(JI,JK,:) * &
                    &   (ZPRES(JI,JK) - ZPRES(JI,JKP)) / RG
            enddo
         enddo

         ! Mass error (integral must be zero)

         do JI = 1, ICONV
            JKP = ICTL(JI) + 1
            if ( ICTL(JI) > 2 ) then
               ZW1 = RG / ( ZPRES(JI,1) - ZPRES(JI,JKP) - &
                    & 0.5*(ZDPRES(JI,2) - ZDPRES(JI,JKP+1)) )
               ZWORK3(JI,:) = ZWORK3(JI,:) * ZW1
            endif
         enddo

         ! Apply uniform correction but assure positive mass at each level

         do JI = 1, ICONV
            if ( ICTL(JI) > 2 ) then
               do JK = ICTL(JI), 1, -1
                  ZCH1C(JI,JK,:) = ZCH1C(JI,JK,:) - ZWORK3(JI,:)
                  ! ZCH1C(JI,JK,:) = MAX( ZCH1C(JI,JK,:), -ZCH1(JI,JK,:)/ZTIMEC(JI) )
               enddo
            endif
         enddo

         ! extend tendencies to first model level

         !  DO JI = 1, ICONV
         !     ZWORK2(JI) = ZDPRES(JI,2) + ZDPRES(JI,3)
         !  ENDDO
         !  DO JN = 1, bkf_kch
         !  DO JI = 1, ICONV
         !    ZCH1(JI,1,JN)  = ZCH1(JI,2,JN) * ZDPRES(JI,3)/ZWORK2(JI)
         !    ZCH1(JI,2,JN)= ZCH1(JI,2,JN) * ZDPRES(JI,2)/ZWORK2(JI)
         !  ENDDO
         !  ENDDO

         do JN = 1, bkf_kch
            do JK = 1, KLEV
               PCH1TEN(1:KLON,JK,JN) = unpack(ZCH1C(1:ICONV,JK,JN), GTRIG(1:KLON), PCH1TEN(1:KLON,JK,JN))
            enddo
         enddo

      endif IF_BKF_LCH1CONV

      !* 8.9  Compute convective tendencies for wind
      !       --------------------------------------

      IF_BKF_LSHALM: if (bkf_lshalm) then

         call CONVECT_UV_TRANSPORT_SHAL1(ICONV, KLEV, ZU, ZV, ZUC, ZVC, &
              &  IDPL, IPBL, ILCL, ICTL,ZUMF, ZUER, ZUDR, &
              &  ZDXDY, ZLMASS, ZWSUB, &
              &  ZTIMC, ITSTEP, GWORK1)

         do JK = 1, KLEV
            ZUC(1:ICONV,JK) = ( ZUC(1:ICONV,JK)- ZU(1:ICONV,JK) ) / ZTIMEC(1:ICONV)
            ZVC(1:ICONV,JK) = ( ZVC(1:ICONV,JK)- ZV(1:ICONV,JK) ) / ZTIMEC(1:ICONV)
         enddo

         !* 8.91 Apply conservation correction
         !       -----------------------------

         ! Compute vertical integrals

         ZWORK2(:) = 0.
         ZWORK2B(:)= 0.
         do JI = 1, ICONV
            do JK = 2, ICTL(JI)
               JKP = JK + 1
               ZW1 = 0.5 *  (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / RG
               ZWORK2(JI) = ZWORK2(JI) + ZUC(JI,JK) * ZW1
               ZWORK2B(JI)= ZWORK2B(JI)+ ZVC(JI,JK) * ZW1
            enddo

         enddo

         !  error (integral must be zero)

         do JI = 1, ICONV
            JKP = ICTL(JI) + 1
            ZW1 = RG / ( ZPRES(JI,1) - ZPRES(JI,JKP) - &
                 & 0.5 * (ZDPRES(JI,2) - ZDPRES(JI,JKP+1)) )
            ZWORK2(JI) = ZWORK2(JI) * ZW1
            ZWORK2B(JI)= ZWORK2B(JI)* ZW1
         enddo

         ! Apply uniform correction

         do JI = 1, ICONV
            if ( ICTL(JI) > 2 ) then
               do JK = ICTL(JI), 1, -1
                  ZUC(JI,JK) = ZUC(JI,JK) - ZWORK2(JI)
                  ZVC(JI,JK) = ZVC(JI,JK) - ZWORK2B(JI)
               enddo
            endif
         enddo
 
         ! extend tendencies to first model level

         ! DO JI = 1, ICONV
         !    ZWORK2(JI) = ZDPRES(JI,2) + ZDPRES(JI,3)
         !    ZUC(JI,1)  = ZUC(JI,2) * ZDPRES(JI,3)/ZWORK2(JI)
         !    ZUC(JI,2)= ZUC(JI,2) * ZDPRES(JI,2)/ZWORK2(JI)
         !    ZVC(JI,1)  = ZVC(JI,2) * ZDPRES(JI,3)/ZWORK2(JI)
         !    ZVC(JI,2)= ZVC(JI,2) * ZDPRES(JI,2)/ZWORK2(JI)
         ! ENDDO


         do JK = 1, KLEV
            PUTEN(1:KLON,JK)  = unpack(ZUC(1:ICONV,JK), GTRIG(1:KLON), PUTEN(1:KLON,JK))
            PVTEN(1:KLON,JK)  = unpack(ZVC(1:ICONV,JK), GTRIG(1:KLON), PVTEN(1:KLON,JK))
         enddo

      endif IF_BKF_LSHALM

      !* 9. Write up- and downdraft mass fluxes and unpack
      !     ----------------------------------------------

      do JK = 1, KLEV
         ZUMF(1:ICONV,JK)  = ZUMF(1:ICONV,JK) / ZDXDY(1:ICONV) ! Mass flux per unit area
      enddo
      do JK = 2, KLEV
         ZUDR(1:ICONV,JK)  = ZUDR(1:ICONV,JK) / ( ZDXDY(1:ICONV) * ZDPRES(1:ICONV,JK) )! detrainment for ERA40
      enddo

      do JI = 1, KLON
         IINDEX(JI) = JI
      enddo
      IJINDEX(1:ICONV) = pack(IINDEX(1:KLON),    GTRIG(1:KLON))

      ZWORK2(:) = 1.
      do JK = 1, KLEV
         do JI = 1, ICONV
            JL = IJINDEX(JI)
            if ( KCLTOP(JL) <= 2 ) ZWORK2(JL) = 0.
            PUMF(JL,JK) = ZUMF(JI,JK) * ZWORK2(JL)

            !  Diagnose cloud coverage
            ZW1 = ZUTHV(JI,JK)*RD /ZPRES(JI,JK)
            if (ZWU(JI,JK) > 1.E-2 .and. JK >= ILCL(JI)) then
               ZW1 = ZUTHV(JI,JK)*RD /ZPRES(JI,JK)
               PCLOUD(JL,JK) = PUMF(JL,JK)/( ZW1* ZWU(JI,JK))
               PCLOUD(JL,JK) = min(1.,max(0.0, PCLOUD(JL,JK)))
            endif

            ! Cloud liquid and ice water mixing ratio in updraft, normalized by convective cloud
            PURCOUT(JL,JK) = ZURC(JI,JK)*PCLOUD(JL,JK)
            PURIOUT(JL,JK) = ZURI(JI,JK)*PCLOUD(JL,JK)

            PUDR(JL,JK) = ZUDR(JI,JK) * ZWORK2(JL)
            PURV(JL,JK) = PURV(JL,JK) * ZWORK2(JL)
         enddo
      enddo

   endif IF_ICONV1
 
   return
end subroutine convect_shallow_c
