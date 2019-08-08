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

subroutine CONVECT_TRIGGER_FUNCT3(KLON, KLEV, &
     & PPRES, PTH, PTHV, PTHES,               &
     & PRV, PW, PZ,                           &
     & PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL,   &
     & PTHVELCL, KLCL, KDPL, KPBL, OTRIG,     &
     & PCAPE, XLAT, MG, MLAC)

   !!**** Determine convective columns as well as the cloudy values of theta,
   !!     and qv at the lifting condensation level (LCL)
   !
   !!    PURPOSE
   !!    -------
   !!      The purpose of this routine is to determine convective columns
   !
   !
   !
   !!**  METHOD
   !!    ------
   !!      Computations are done at every model level starting from bottom.
   !!      The use of masks allows to optimise the inner loops (horizontal loops).
   !!      What we look for is the undermost unstable level at each grid point.
   !
   !
   !
   !!    EXTERNAL
   !!    --------
   !!     Routine CONVECT_SATMIXRATIO
   !
   !
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !!      Module YOMCST
   !!          RG                 ! gravity constant
   !!          RATM               ! Reference pressure
   !!          RD, RV           ! Gaz  constants for dry air and water vapor
   !!          RCPD               ! Cpd (dry air)
   !!          RTT                ! triple point temperature
   !!          RBETW, RGAMW      ! constants for vapor saturation pressure
   !
   !!      Module YOE_CONVPAR
   !!          XA25               ! reference grid area
   !!          XZLCL              ! maximum height difference between
   !!                             ! the surface and the DPL
   !!          XZPBL              ! minimum mixed layer depth to sustain convection
   !!          XWTRIG             ! constant in vertical velocity trigger
   !!          XCDEPTH            ! minimum necessary cloud depth
   !!          XNHGAM             ! coefficient for buoyancy term in w eq.
   !!                             ! accounting for nh-pressure
   !!          XDTHPBL            ! theta perturbation in PBL
   !!          XDRVPBL            ! moisture perturbation in PBL
   !
   !!      Module YOE_CONVPAREXT
   !!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
   !
   !!    REFERENCE
   !!    ---------
   !
   !!      Book2 of documentation ( routine TRIGGER_FUNCT)
   !!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
   !
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !
   !!    MODIFICATIONS
   !!    -------------
   !!      Original    07/11/95
   !!   Last modified  20/03/97  Select first departure level
   !!                            that produces a cloud thicker than XCDEPTH
   !----------------------------------------------------------------------------


   !*       0.    DECLARATIONS
   !              ------------

   use YOMCST
   use YOE_CONVPAR
   use YOE_CONVPAREXT
   use cnv_options

   implicit none
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

   !*       0.1   Declarations of dummy arguments :

   integer, intent(IN)                   :: KLON      ! horizontal loop index
   integer, intent(IN)                   :: KLEV      ! vertical loop index
   real, dimension(KLON,KLEV),intent(IN) :: PTH, PTHV ! theta, theta_v
   real, dimension(KLON,KLEV),intent(IN) :: PTHES     ! envir. satur. theta_e
   real, dimension(KLON,KLEV),intent(IN) :: PRV       ! vapor mixing ratio
   real, dimension(KLON,KLEV),intent(IN) :: PPRES     ! pressure
   real, dimension(KLON,KLEV),intent(IN) :: PZ        ! height of grid point (m)
   real, dimension(KLON,KLEV),intent(IN) :: PW        ! vertical velocity

   real, dimension(KLON),     intent(OUT):: PTHLCL    ! theta at LCL
   real, dimension(KLON),     intent(OUT):: PTLCL     ! temp. at LCL
   real, dimension(KLON),     intent(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
   real, dimension(KLON),     intent(OUT):: PWLCL     ! parcel velocity at  LCL
   real, dimension(KLON),     intent(OUT):: PZLCL     ! height at LCL (m)
   real, dimension(KLON),     intent(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
   logical,   dimension(KLON),  intent(OUT):: OTRIG     ! logical mask for convection
   integer, dimension(KLON),  intent(INOUT):: KLCL    ! contains vert. index of LCL
   integer, dimension(KLON),  intent(INOUT):: KDPL    ! contains vert. index of DPL
   integer, dimension(KLON),  intent(INOUT):: KPBL    ! contains index of source layer top
   real, dimension(KLON),     intent(OUT):: PCAPE     ! CAPE (J/kg) for diagnostics

   real, dimension(KLON),     intent(IN)  :: XLAT    ! Latitude
   real, dimension(KLON),     intent(IN)  :: MG      ! Land_sea mask
   real, dimension(KLON),     intent(IN)  :: MLAC    ! Lake mask


   !*       0.2   Declarations of local variables :

   integer :: JKK, JK, JKP, JKM, JKDL, JL, JKT, JT! vertical loop index
   integer :: JI                                  ! horizontal loop index
   integer :: IIE, IKB, IKE                       ! horizontal + vertical loop bounds
   real    :: ZEPS, ZEPSA                         ! R_d / R_v, R_v / R_d
   real    :: ZCPORD, ZRDOCP                      ! C_pd / R_d, R_d / C_pd

   real    :: WKLCL, WKLCLD
   real, dimension(KLON) :: WKLCLA


   real, dimension(KLON) :: ZTHLCL, ZTLCL, ZRVLCL, & ! locals for PTHLCL,PTLCL
        &  ZWLCL,  ZZLCL, ZTHVELCL  ! PRVLCL, ....
   integer, dimension(KLON) :: IDPL, IPBL, ILCL      ! locals for KDPL, ...
   real, dimension(KLON) :: ZPLCL    ! pressure at LCL
   real, dimension(KLON) :: ZZDPL    ! height of DPL
   real, dimension(KLON) :: ZTHVLCL  ! theta_v at LCL = mixed layer value
   real, dimension(KLON) :: ZTMIX    ! mixed layer temperature
   real, dimension(KLON) :: ZEVMIX   ! mixed layer water vapor pressure
   real, dimension(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
   real, dimension(KLON) :: ZCAPE    ! convective available energy (m^2/s^2/g)
   real, dimension(KLON) :: ZTHEUL   ! updraft equiv. pot. temperature (K)
   real, dimension(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
   real, dimension(KLON) :: ZDP      ! pressure between LCL and model layer
   real, dimension(KLON) :: ZTOP     ! estimated cloud top (m)
   real, dimension(KLON,KLEV):: ZCAP ! CAPE at every level for diagnostics
   real,  dimension(KLON) :: ZWORK1, ZWORK2, ZWORK3 ! work arrays
   logical, dimension(KLON) :: GTRIG, GTRIG2          ! local arrays for OTRIG
   logical, dimension(KLON) :: GWORK1                 ! work array

   !     "RAMP" FOR WKLCL :
   !     ================

   !     WKLCL WILL INCREASE FROM KFCTRIG4(3) TO KFCTRIG4(4)

   WKLCL = KFCTRIG4(4)

   !      Latitudinal ramp for WKLCL :
   !     ============================

   !     WKLCL will take on different values:
   !     over land and lakes: we kee the value set by the "ramp" above
   !     over sea water:
   !       for |lat| >= TRIGLAT(2) we keep value set by the "ramp" above
   !       for |lat| <= TRIGLAT(1) we use the new value KFCTRIGL
   !       and linear interpolation in between TRIGLAT(1) and TRIGLAT(2)

   WKLCLA(:) = WKLCL

   if (KFCTRIGLAT) then


      do JI=1,KLON


         if (abs(XLAT(JI)) .le. TRIGLAT(1).and. &
              MG(JI) .le. 0.001 .and. &
              MLAC(JI) .le. 0.001) then
            WKLCLA(JI)= KFCTRIGL
         else if (abs(XLAT(JI)).gt.TRIGLAT(1) .and. &
              abs(XLAT(JI)).lt.TRIGLAT(2) .and. &
              MG(JI) .le. 0.001 .and. &
              MLAC(JI) .le. 0.001) then
            WKLCLA(JI)= ( ((abs(XLAT(JI))-TRIGLAT(1))/ &
                 (TRIGLAT(2)-TRIGLAT(1)))* &
                 (WKLCL-KFCTRIGL) ) + KFCTRIGL
         else
            WKLCLA(JI)= WKLCL
         endif

      end do

   endif
   !==============


   !----------------------------------------------------------------------------

   !*       0.3    Compute array bounds
   !               --------------------

   IIE = KLON
   IKB = 1 + JCVEXB
   IKE = KLEV - JCVEXT


   !*       1.     Initialize local variables
   !               --------------------------

   ZEPS       = RD  / RV
   ZEPSA      = RV  / RD
   ZCPORD     = RCPD / RD
   ZRDOCP     = RD  / RCPD

   OTRIG(:)   = .false.

   IDPL(:)    = KDPL(:)
   IPBL(:)    = KPBL(:)
   ILCL(:)    = KLCL(:)

   PWLCL(:)   = _ZERO_
   ZWLCL(:)   = _ZERO_
   PTHLCL(:)  = _ONE_
   PTHVELCL(:)= _ONE_
   PTLCL(:)   = _ONE_
   PRVLCL(:)  = _ZERO_
   PWLCL(:)   = _ZERO_
   PZLCL(:)   = PZ(:,IKB)
   ZZDPL(:)   = PZ(:,IKB)
   GTRIG2(:)  = .true.
   ZCAP(:,:)  = _ZERO_



   !       1.     Determine highest necessary loop test layer
   !              -------------------------------------------

   JT = IKE - 2
   do JK = IKB + 1, IKE - 2
      if ( PZ(1,JK) - PZ(1,IKB) < 12.E3 ) JT = JK
   enddo


   !*       2.     Enter loop for convection test
   !               ------------------------------

   JKP  = minval( IDPL(:) ) + 1
   JKT  = JT

   do JKK = JKP, JKT

      GWORK1(:) = ZZDPL(:) - PZ(:,IKB) < XZLCL .and. GTRIG2(:)
      ! we exit the trigger test when the center of the mixed layer is more
      ! than 3500 m  above soil level.
      where ( GWORK1(:) )
         ZDPTHMIX(:) = _ZERO_
         ZPRESMIX(:) = _ZERO_
         ZTHLCL(:)   = _ZERO_
         ZRVLCL(:)   = _ZERO_
         ZZDPL(:)    = PZ(:,JKK)
         IDPL(:)     = JKK
      end where


      !*       3.     Construct a mixed layer of at least 60 hPa (XZPBL)
      !               ------------------------------------------

      do JK = JKK, IKE - 1
         JKM = JK + 1
         do JI = 1, IIE
            if ( GWORK1(JI) .and. ZDPTHMIX(JI) < XZPBL ) then
               IPBL(JI)     = JK
               ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
               ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
               ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
               ZTHLCL(JI)   = ZTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
               ZRVLCL(JI)   = ZRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
            endif
         enddo
         if ( minval ( ZDPTHMIX(:) ) >= XZPBL ) exit
      enddo


      where ( GWORK1(:) )

         ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
         ! ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:)
         ! ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:)
         ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:) + XDTHPBL
         ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:) + XDRVPBL
         ZTHVLCL(:)  = ZTHLCL(:) * ( _ONE_ + ZEPSA * ZRVLCL(:) )              &
              &           / ( _ONE_ + ZRVLCL(:) )

         !*       4.1    Use an empirical direct solution ( Bolton formula )
         !               to determine temperature and pressure at LCL.
         !               Nota: the adiabatic saturation temperature is not
         !                     equal to the dewpoint temperature
         !               ----------------------------------------------------


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
         ZPLCL(:)  = RATM * ( ZTLCL(:) / ZTHLCL(:) ) ** ZCPORD

      end where


      !*       4.2    Correct ZTLCL in order to be completely consistent
      !               with MNH saturation formula
      !               ---------------------------------------------

      call CONVECT_SATMIXRATIO( KLON, ZPLCL, ZTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
      where( GWORK1(:) )
         ZWORK2(:) = ZWORK1(:) / ZTLCL(:) * ( RBETW / ZTLCL(:) - RGAMW ) ! dr_sat/dT
         ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
              &     ( _ONE_ + ZLV(:) / ZCPH(:) * ZWORK2(:) )
         ZTLCL(:)  = ZTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)

      end where


      !*       4.3    If ZRVLCL = PRVMIX is oversaturated set humidity
      !               and temperature to saturation values.
      !               ---------------------------------------------

      call CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
      where( GWORK1(:) .and. ZRVLCL(:) > ZWORK1(:) )
         ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( RBETW / ZTMIX(:) - RGAMW ) ! dr_sat/dT
         ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
              &    ( _ONE_ + ZLV(:) / ZCPH(:) * ZWORK2(:) )
         ZTLCL(:)  = ZTMIX(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
         ZRVLCL(:) = ZRVLCL(:) + ZWORK2(:)
         ! Jing, bug fixed. should be ZRVLCL+ZWORK2, ZWORK2 here is negative value

         ZPLCL(:)  = ZPRESMIX(:)
         ZTHLCL(:) = ZTLCL(:) * ( RATM / ZPLCL(:) ) ** ZRDOCP
         ZTHVLCL(:)= ZTHLCL(:) * ( _ONE_ + ZEPSA * ZRVLCL(:) )                &
              &           / ( _ONE_ + ZRVLCL(:) )
      end where


      !*        5.1   Determine  vertical loop index at the LCL and DPL
      !               --------------------------------------------------

      do JK = JKK, IKE - 1
         do JI = 1, IIE
            if ( ZPLCL(JI) <= PPRES(JI,JK) .and. GWORK1(JI) ) ILCL(JI) = JK + 1
         enddo
      enddo


      !*        5.2   Estimate height and environm. theta_v at LCL
      !               --------------------------------------------------

      do JI = 1, IIE
         if ( GWORK1(JI) ) then
            JK   = ILCL(JI)
            JKM  = JK - 1
            ZDP(JI)    = log( ZPLCL(JI) / PPRES(JI,JKM) ) /                   &
                 & log( PPRES(JI,JK) / PPRES(JI,JKM) )
            ZWORK1(JI) = PTHV(JI,JKM) + ( PTHV(JI,JK) - PTHV(JI,JKM) ) * ZDP(JI)
            ! we compute the precise value of the LCL
            ! The precise height is between the levels ILCL and ILCL-1.
            ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
         end if
      enddo
      where( GWORK1(:) )
         ZTHVELCL(:) = ZWORK1(:)
         ZZLCL(:)    = ZWORK2(:)
      end where


      !*       6.     Check to see if cloud is bouyant
      !               --------------------------------

      !*      6.1    Compute grid scale vertical velocity perturbation term ZWORK1
      !               ----------------------------------------------------------

      !  normalize w grid scale to a 25 km refer. grid
      do JI = 1, IIE
         if ( GWORK1(JI) )  then
            JK  = ILCL(JI)
            JKM = JK - 1
            JKDL= IDPL(JI)
            !PV following two commented lines is what was used in original bkf code;  estimate of Wlcl (first line)
            ! is incomprehensible to me...
            !Jing        ZWORK1(JI) =  ( PW(JI,JK) + PW(JI,JKDL)*ZZLCL(JI)/MAX(_ONE_,PZ(JI,JKDL)) ) * _HALF_  &
            !Jing                   &       * SQRT( PDXDY(JI) / XA25 )
            !PV commented line is what KF uses
            ! ZWORK1(JI) =  ( PW(JI,JKM)  + ( PW(JI,JK) - PW(JI,JKM) ) * ZDP(JI) )  &

            !Jing, uncomment the above to use same as KF
            !      WKLCLD = min(0.02,WKLCLA(JI))
            WKLCLD = WKLCLA(JI)
            ZWORK1(JI) =   PW(JI,JKM)  + ( PW(JI,JK) - PW(JI,JKM) ) * ZDP(JI) -WKLCLD


            ! compute sign of normalized grid scale w
            ZWORK2(JI) = sign( _ONE_, ZWORK1(JI) )
            !PV eq5 of Bechtoldetal2001: compute perturbation theta to promote/suppress triggering
            ZWORK1(JI) = XWTRIG * ZWORK2(JI) * abs( ZWORK1(JI) ) ** 0.333     &
                 &        * ( RATM / ZPLCL(JI) ) ** ZRDOCP
         end if
      enddo

      !*       6.2    Compute parcel vertical velocity at LCL
      !               ---------------------------------------

      do JI = 1, IIE
         JKDL = IDPL(JI)
         ZWORK3(JI) = RG * ZWORK1(JI) * ( ZZLCL(JI) - PZ(JI,JKDL) )       &
              &   / ( PTHV(JI,JKDL) + ZTHVELCL(JI) )
      enddo

      ! DO JI = 1, IIE
      !    JKDL = IDPL(JI)
      !    JK   = ILCL(JI)
      !    ZWORK4(JI) = RG/RCPD * _HALF_ * ( PHSFLX(JI,JK) + PHSFLX(JI,JKDL) ) &
      !               &   * ( ZZLCL(JI) - PZ(JI,JKDL) ) / ZTHVELCL(JI)
      !    ZWORK4(JI) = 3. * MAX( 1.E-3, ZWORK4(JI) ) ** .3333
      ! ENDDO

      where( GWORK1(:) )
         ZWLCL(:)  = _ONE_ + _HALF_ * ZWORK2(:) * sqrt( abs( ZWORK3(:) ) )
         ! ZWLCL(:)  = ZWORK4(:) + .25 * ZWORK2(:) * SQRT( ABS( ZWORK3(:) ) ) ! UPG PB
         GTRIG(:)  = ZTHVLCL(:) - ZTHVELCL(:) + ZWORK1(:) > _ZERO_ .and.       &
              & ZWLCL(:) > _ZERO_
      end where



      !*       6.3    Look for parcel that produces sufficient cloud depth.
      !               The cloud top is estimated as the level where the CAPE
      !               is smaller  than a given value (based on vertical velocity eq.)
      !               ----------------------------------------------------------

      where( GWORK1(:) )
         ZTHEUL(:) = ZTLCL(:) * ( ZTHLCL(:) / ZTLCL(:) ) **                &
              &            ( _ONE_ - 0.28 * ZRVLCL(:) )                    &
              &          * exp( ( 3374.6525 / ZTLCL(:) - 2.5403 ) *   &
              &                 ZRVLCL(:) * ( _ONE_ + 0.81 * ZRVLCL(:) ) )
      end where

      ZCAPE(:) = _ZERO_
      ZTOP(:)  = _ZERO_
      ZWORK3(:)= _ZERO_
      JKM = minval( ILCL(:) )
      do JL = JKM, JT
         JK = JL + 1
         do JI = 1, IIE
            if ( GWORK1(JI) ) then
               ZWORK1(JI) = ( _TWO_ * ZTHEUL(JI) /                             &
                    & ( PTHES(JI,JK) + PTHES(JI,JL) ) - _ONE_ ) * ( PZ(JI,JK) - PZ(JI,JL) )
               if ( JL < ILCL(JI) ) ZWORK1(JI) = _ZERO_
               ZCAPE(JI)  = ZCAPE(JI) + ZWORK1(JI)
               ZCAP(JI,JKK) = ZCAP(JI,JKK) + RG * max( _ZERO_, ZWORK1(JI) ) ! actual CAPE
               ZWORK2(JI) = XNHGAM * RG * ZCAPE(JI) + 1.05 * ZWLCL(JI) * ZWLCL(JI)
               ! the factor 1.05 takes entrainment into account
               ZWORK2(JI) = sign( _ONE_, ZWORK2(JI) )
               ZWORK3(JI) = ZWORK3(JI) + min(_ZERO_, ZWORK2(JI) )
               ZWORK3(JI) = max( -_ONE_, ZWORK3(JI) )
               ! Nota, the factors ZWORK2 and ZWORK3 are only used to avoid
               ! if and goto statements, the difficulty is to extract only
               ! the level where the criterium is first fullfilled
               ZTOP(JI)   = PZ(JI,JL) * _HALF_ * ( _ONE_ + ZWORK2(JI) ) * ( _ONE_ + ZWORK3(JI) ) + &
                    & ZTOP(JI)  * _HALF_ * ( _ONE_ - ZWORK2(JI) )
            end if
         enddo
      enddo


      where( ZTOP(:) - ZZLCL(:)  >=  XCDEPTH  .and. GTRIG(:) .and. GTRIG2(:) )
         GTRIG2(:)   = .false.
         OTRIG(:)    = GTRIG(:)     ! we  select the first departure level
         PTHLCL(:)   = ZTHLCL(:)    ! that gives sufficient cloud depth
         PRVLCL(:)   = ZRVLCL(:)
         PTLCL(:)    = ZTLCL(:)
         PWLCL(:)    = ZWLCL(:)
         PZLCL(:)    = ZZLCL(:)
         PTHVELCL(:) = ZTHVELCL(:)
         KDPL(:)     = IDPL(:)
         KPBL(:)     = IPBL(:)
         KLCL(:)     = ILCL(:)
      end where

   enddo

   do JI = 1, IIE
      PCAPE(JI) = maxval( ZCAP(JI,:) ) ! maximum CAPE for diagnostics
   enddo

end subroutine CONVECT_TRIGGER_FUNCT3

