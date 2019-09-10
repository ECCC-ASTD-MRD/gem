!######################################################################
subroutine CONVECT_TRIGGER_SHAL4(KLON, KLEV, &
     &  PPRES, PTH, PTHV, PTHES, PMSE, &
     &  PRV, PZ, PDTPERT, &
     &  PTHLCL, PTLCL, PRVLCL, PWLCL, PPLCL, PZLCL, &
     &  PTHVELCL, PMSEELCL, KLCL, KDPL, KPBL, OTRIG)
   !######################################################################

   !!**** Determine convective columns as well as the cloudy values of theta,
   !!     and qv at the lifting condensation level (LCL)
   !!
   !!    PURPOSE
   !!    -------
   !!      The purpose of this routine is to determine convective columns
   !!
   !!
   !!
   !!**  METHOD
   !!    ------
   !!      Computations are done at every model level starting from bottom.
   !!      The use of masks allows to optimise the inner loops (horizontal loops).
   !!      What we look for is the undermost unstable level at each grid point.
   !!
   !!
   !!
   !!    EXTERNAL
   !!    --------
   !!     Routine CONVECT_SATMIXRATIO
   !!
   !!
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !!      Module YOMCST
   !!          RG                 ! gravity constant
   !!          RATM               ! Reference pressure
   !!          RD, RV           ! Gaz  constants for dry air and water vapor
   !!          RCPD               ! Cpd (dry air)
   !!          RTT                ! triple point temperature
   !!          RBETW, RGAMW      ! constants for vapor saturation pressure
   !!
   !!      Module YOE_CONVPAR
   !!          XA25               ! reference grid area
   !!          XZLCL              ! maximum height difference between
   !!                             ! the surface and the DPL
   !!          XZPBL              ! minimum mixed layer depth to sustain convection
   !!          XCDEPTH            ! minimum necessary cloud depth
   !!          XCDEPTH_D          ! maximum allowed cloud depth
   !!          XDTPERT            ! add small Temp peturbation
   !!          XNHGAM             ! coefficient for buoyancy term in w eq.
   !!                             ! accounting for nh-pressure
   !!
   !!      Module YOE_CONVPAREXT
   !!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
   !!
   !!    REFERENCE
   !!    ---------
   !!
   !!      Book2 of documentation ( routine TRIGGER_FUNCT)
   !!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
   !!
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !!
   !!    MODIFICATIONS
   !!    -------------
   !!      Original    07/11/95
   !!   Last modified  20/03/97  Select first departure level
   !!                            that produces a cloud thicker than XCDEPTH
   !-------------------------------------------------------------------------------

   !*       0.    DECLARATIONS
   !              ------------

   use YOMCST
   use YOE_CONVPAR_SHAL
   use YOE_CONVPAREXT

   implicit none
!!!#include <arch_specific.hf>

   !*       0.1   Declarations of dummy arguments :

   integer, intent(IN)                      :: KLON      ! horizontal loop index
   integer, intent(IN)                      :: KLEV      ! vertical loop index
   real,    dimension(KLON,KLEV),intent(IN) :: PTH, PTHV ! theta, theta_v
   real,    dimension(KLON,KLEV),intent(IN) :: PTHES     ! envir. satur. theta_e
   real,    dimension(KLON,KLEV),intent(IN) :: PMSE      ! envir. moist static energy
   real,    dimension(KLON,KLEV),intent(IN) :: PRV       ! vapor mixing ratio
   real,    dimension(KLON,KLEV),intent(IN) :: PPRES     ! pressure
   real,    dimension(KLON,KLEV),intent(IN) :: PZ        ! height of grid point (m)
   real,    dimension(KLON),     intent(IN) :: PDTPERT   ! temp.perturbation (K)

   real,    dimension(KLON),     intent(OUT):: PTHLCL    ! theta at LCL
   real,    dimension(KLON),     intent(OUT):: PTLCL     ! temp. at LCL
   real,    dimension(KLON),     intent(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
   real,    dimension(KLON),     intent(OUT):: PWLCL     ! parcel velocity at  LCL
   real,    dimension(KLON),     intent(OUT):: PPLCL     ! pressure at LCL (Pa)
   real,    dimension(KLON),     intent(OUT):: PZLCL     ! height at LCL (m)
   real,    dimension(KLON),     intent(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
   real,    dimension(KLON),     intent(OUT):: PMSEELCL  ! environm. MSE at LCL (K)
   logical, dimension(KLON),  intent(OUT)  :: OTRIG     ! logical mask for convection
   integer, dimension(KLON),  intent(INOUT):: KLCL    ! contains vert. index of LCL
   integer, dimension(KLON),  intent(INOUT):: KDPL    ! contains vert. index of DPL
   integer, dimension(KLON),  intent(INOUT):: KPBL    ! contains index of source layer top

   !*       0.2   Declarations of local variables :

   integer :: JK, JKP, JKM, JL          ! vertical loop index
   integer :: JI                                  ! horizontal loop index
   real    :: ZEPS, ZEPSA                         ! R_d / R_v, R_v / R_d
   real    :: ZCPORD, ZRDOCP                      ! C_pd / R_d, R_d / C_pd

   integer, dimension(KLON) :: &
        IDPL, IPBL, ILCL, &  ! locals for KDPL, ...
        ITOP                 ! work array to store highest test layer
   real,    dimension(KLON) :: &
        ZTHLCL, ZTLCL, ZRVLCL, &            ! locals for PTHLCL,PTLCL
        ZWLCL, ZZLCL, ZTHVELCL, ZMSEELCL,&  ! PRVLCL, ....
        ZZDPL, &     ! height of DPL
        ZTHVLCL, &   ! theta_v at LCL = mixed layer value
        ZTMIX, &     ! mixed layer temperature
        ZEVMIX , &   ! mixed layer water vapor pressure
        ZDPTHMIX, ZPRESMIX, & ! mixed layer depth and pressure
        ZCAPE, &     ! convective available energy (m^2/s^2/g)
        ZCAP, &      ! pseudo for CAPE
        ZTHEUL, &    ! updraft equiv. pot. temperature (K)
        ZLV, ZCPH, & ! specific heats of vaporisation, dry air
        ZDP, &       ! pressure between LCL and model layer
        ZTOP, &      ! estimated cloud top (m)
        ZWORK1, ZWORK2, ZWORK3    ! work arrays
   logical, dimension(KLON) :: GTRIG2          ! local arrays for OTRIG
   logical, dimension(KLON) :: GWORK1          ! work array

   !-------------------------------------------------------------------------------

   !* Initialize local variables
   !  --------------------------

   ZEPS       = RD  / RV
   ZEPSA      = RV  / RD
   ZCPORD     = RCPD / RD
   ZRDOCP     = RD  / RCPD

   OTRIG(:)   = .false.

   IDPL(:)    = KDPL(:)
   IPBL(:)    = KPBL(:)
   ILCL(:)    = KLCL(:)
   ITOP(:)    = 1 !revive itop

   PWLCL(:)   = 0.
   ZWLCL(:)   = 0.
   PTHLCL(:)  = 1.
   PTHVELCL(:)= 1.
   PMSEELCL(:)= 1.
   PTLCL(:)   = 1.
   PRVLCL(:)  = 0.
   PWLCL(:)   = 0.
   PZLCL(:)   = PZ(:,1)
   ZZDPL(:)   = PZ(:,1)
   GTRIG2(:)  = .true.

   !* 1. Determine highest necessary loop test layer
   !     -------------------------------------------

   do JK = 2, KLEV - 2
      where (PZ(:,JK) - PZ(:,1) <= 5.E3) ITOP(:) = JK
   enddo

   !* 2. Enter loop for convection test
   !     ------------------------------

   JKP = minval(IDPL(:)) + 1  !IDPL was set to 2 in convect_shallow
   if (JKP.ne.2) then
      call physeterror('convect_trigger_shal', 'IDPL is no longer 2; non MPI bit reprod')
      return
   endif

   GWORK1(:) = (ZZDPL(:) - PZ(:,1) < XZLCL)
   ! we exit the trigger test when the center of the mixed layer is more
   ! than 1500 m  above soil level.
   where (GWORK1(:))
      ZDPTHMIX(:) = 0.
      ZPRESMIX(:) = 0.
      ZTHLCL(:)   = 0.
      ZRVLCL(:)   = 0.
      ZZDPL(:)    = PZ(:,JKP)
      IDPL(:)     = JKP
   end where

   !* 3. Construct a mixed layer of at least 50 hPa (XZPBL)
   !     ------------------------------------------

   do JK = JKP, KLEV - 1
      JKM = JK + 1
      where (GWORK1(:) .and. ZDPTHMIX(:) < XZPBL)
         IPBL(:)     = JK
         ZWORK1(:)   = PPRES(:,JK) - PPRES(:,JKM)
         ZDPTHMIX(:) = ZDPTHMIX(:) + ZWORK1(:)
         ZPRESMIX(:) = ZPRESMIX(:) + PPRES(:,JK) * ZWORK1(:)
         ZTHLCL(:)   = ZTHLCL(:)   + PTH(:,JK)   * ZWORK1(:)
         ZRVLCL(:)   = ZRVLCL(:)   + PRV(:,JK)   * ZWORK1(:)
      end where
   enddo


   where ( GWORK1(:) )

      ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
      ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:) + PDTPERT(:) ! add small Temp Perturb.
      ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:)
      ZTHVLCL(:)  = ZTHLCL(:) * ( 1. + ZEPSA * ZRVLCL(:) )                 &
           &           / ( 1. + ZRVLCL(:) )

      !* 4.1 Use an empirical direct solution ( Bolton formula )
      !      to determine temperature and pressure at LCL.
      !      Nota: the adiabatic saturation temperature is not
      !            equal to the dewpoint temperature
      !      ----------------------------------------------------


      ZTMIX(:)  = ZTHLCL(:) * ( ZPRESMIX(:) / RATM ) ** ZRDOCP
      ZEVMIX(:) = ZRVLCL(:) * ZPRESMIX(:) / ( ZRVLCL(:) + ZEPS )
      ZEVMIX(:) = max( 1.E-8, ZEVMIX(:) )
      ZWORK1(:) = log( ZEVMIX(:) / 613.3 )
      ! dewpoint temperature
      ZWORK1(:) = ( 4780.8 - 32.19 * ZWORK1(:) ) / ( 17.502 - ZWORK1(:) )
      ! adiabatic saturation temperature
      ZTLCL(:)  = ZWORK1(:) - ( .212 + 1.571E-3 * ( ZWORK1(:) - RTT )      &
           & - 4.36E-4 * ( ZTMIX(:) - RTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
      ZTLCL(:)  = min( ZTLCL(:), ZTMIX(:) )
      PPLCL(:)  = RATM * ( ZTLCL(:) / ZTHLCL(:) ) ** ZCPORD

   end where


   !* 4.2  Correct ZTLCL in order to be completely consistent
   !       with MNH saturation formula
   !       ---------------------------------------------

   call CONVECT_SATMIXRATIO( KLON, PPLCL, ZTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
   where( GWORK1(:) )
      ZWORK2(:) = ZWORK1(:) / ZTLCL(:) * ( RBETW / ZTLCL(:) - RGAMW ) ! dr_sat/dT
      ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
           &     ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) )
      ZTLCL(:)  = ZTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
   end where


   !* 4.3 If ZRVLCL = PRVMIX is oversaturated set humidity
   !      and temperature to saturation values.
   !      ---------------------------------------------

   call CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
   where( GWORK1(:) .and. ZRVLCL(:) > ZWORK1(:) )
      ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( RBETW / ZTMIX(:) - RGAMW ) ! dr_sat/dT
      ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
           &    ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) )
      ZTLCL(:)  = ZTMIX(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
      !        ZRVLCL(:) = ZRVLCL(:) - ZWORK2(:) ! zwork2 is negative; vapor must be removed not added if supersaturated
      ZRVLCL(:) = ZRVLCL(:) + ZWORK2(:)
      PPLCL(:)  = ZPRESMIX(:)
      ZTHLCL(:) = ZTLCL(:) * ( RATM / PPLCL(:) ) ** ZRDOCP
      ZTHVLCL(:)= ZTHLCL(:) * ( 1. + ZEPSA * ZRVLCL(:) )                   &
           &           / ( 1. + ZRVLCL(:) )
   end where


   !* 5.1 Determine  vertical loop index at the LCL and DPL
   !      --------------------------------------------------

   do JK = JKP, KLEV - 1
      where (PPLCL(:) <= PPRES(:,JK) .and. GWORK1(:)) ILCL(:) = JK + 1
   enddo

   !* 5.2 Estimate height, environm. theta_v and MSE at LCL
   !      --------------------------------------------------

   do JI = 1, KLON
      if (GWORK1(JI)) then
         JK   = ILCL(JI)
         JKM  = JK - 1
         ZDP(JI)    = log( PPLCL(JI) / PPRES(JI,JKM) ) /                     &
              & log( PPRES(JI,JK) / PPRES(JI,JKM) )
         ZWORK1(JI) = PTHV(JI,JKM) + ( PTHV(JI,JK) - PTHV(JI,JKM) ) * ZDP(JI)
         ! we compute the precise value of the LCL
         ! The precise height is between the levels ILCL and ILCL-1.
         ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
         ZWORK3(JI) = PMSE(JI,JKM) + ( PMSE(JI,JK) - PMSE(JI,JKM) ) * ZDP(JI)
      end if
   enddo
   where (GWORK1(:))
      ZTHVELCL(:) = ZWORK1(:)
      ZZLCL(:)    = ZWORK2(:)
      ZMSEELCL(:) = ZWORK3(:)
   end where

   !* 6.  Check to see if cloud is bouyant
   !      --------------------------------

   !* 6.1 Compute grid scale vertical velocity perturbation term ZWORK1
   !      -------------------------------------------------------------

   !* 6.2 Compute parcel vertical velocity at LCL
   !      ---------------------------------------
   ZWLCL(:) = 1.

   !* 6.3 Look for parcel that produces sufficient cloud depth.
   !      The cloud top is estimated as the level where the CAPE
   !      is smaller  than a given value (based on vertical velocity eq.)
   !      --------------------------------------------------------------

   where( GWORK1(:) )
      ZTHEUL(:) = ZTLCL(:) * ( ZTHLCL(:) / ZTLCL(:) ) ** &
           &                              ( 1. - 0.28 * ZRVLCL(:) ) &
           &               * exp( ( 3374.6525 / ZTLCL(:) - 2.5403 )     &
           &                      * ZRVLCL(:) * ( 1. + 0.81 * ZRVLCL(:) ) )
   end where

   ZCAPE(:) = 0.
   ZCAP(:)  = 0.
   ZTOP(:)  = 0.
   ZWORK3(:)= 0.
   do JI = 1, KLON
      if ( GWORK1(JI) ) then
         do JL = ILCL(JI), ITOP(JI)
            !DO JL = ILCL(JI), JT
            JK = JL + 1
            ZWORK1(JI) = ( 2. * ZTHEUL(JI) /                                &
                 & ( PTHES(JI,JK) + PTHES(JI,JL) ) - 1. ) * ( PZ(JI,JK) - PZ(JI,JL) )
            if ( JL < ILCL(JI) ) ZWORK1(JI) = 0.
            ZCAP(JI)   = ZCAP(JI) + ZWORK1(JI)
            ZCAPE(JI)  = ZCAPE(JI) + RG * max( 0., ZWORK1(JI) )
            ZWORK2(JI) = XNHGAM * RG * ZCAP(JI) + 1.05 * ZWLCL(JI) * ZWLCL(JI)
            ! the factor 1.05 takes entrainment into account
            ZWORK2(JI) = sign( 1., ZWORK2(JI) )
            ZWORK3(JI) = ZWORK3(JI) + min(0., ZWORK2(JI) )
            ZWORK3(JI) = max( -1., ZWORK3(JI) )
            ! Nota, the factors ZWORK2 and ZWORK3 are only used to avoid
            ! if and goto statements, the difficulty is to extract only
            ! the level where the criterium is first fullfilled
            ZTOP(JI)   = PZ(JI,JL) * 0.5 * ( 1. + ZWORK2(JI) ) * ( 1. + ZWORK3(JI) ) + &
                 & ZTOP(JI)  * 0.5 * ( 1. - ZWORK2(JI) )
         enddo
      end if
   enddo

   ZWORK2(:) = ZTOP(:) - ZZLCL(:)
   where (ZWORK2(:) >= XCDEPTH  .and. GTRIG2(:) .and. ZCAPE(:) > 10.)
      GTRIG2(:)   = .false.
      OTRIG(:)    = .true.       ! we  select the first departure level
      PTHLCL(:)   = ZTHLCL(:)    ! that gives sufficient cloud depth
      PRVLCL(:)   = ZRVLCL(:)
      PTLCL(:)    = ZTLCL(:)
      PWLCL(:)    = ZWLCL(:)
      PZLCL(:)    = ZZLCL(:)
      PTHVELCL(:) = ZTHVELCL(:)
      PMSEELCL(:) = ZMSEELCL(:)
      KDPL(:)     = IDPL(:)
      KPBL(:)     = IPBL(:)
      KLCL(:)     = ILCL(:)
   end where

   return
end subroutine CONVECT_TRIGGER_SHAL4

