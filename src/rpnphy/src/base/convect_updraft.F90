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

subroutine CONVECT_UPDRAFT(KLON, KLEV,                   &
     &  KICE, PPRES, PDPRES, PZ, PTHL, PTHV, PTHES, PRW, &
     &  PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,   &
     &  PMFLCL, OTRIG, KLCL, KDPL, KPBL,                 &
     &  PUMF, PUER, PUDR, PUTHL, PUTHV, PWU,             &
     &  PURW,  PURC, PURI, PURR, PURS, PUPR,             &
     &  PUTPR, PCAPE, KCTL, KETL )

   !!**** Compute updraft properties from DPL to CTL.
   !
   !
   !!    PURPOSE
   !!    -------
   !!      The purpose of this routine is to determine updraft properties
   !!      ( mass flux, thermodynamics, precipitation )
   !
   !
   !!**  METHOD
   !!    ------
   !!      Computations are done at every model level starting from bottom.
   !!      The use of masks allows to optimise the inner loops (horizontal loops).
   !
   !
   !
   !!    EXTERNAL
   !!    --------
   !!     Routine CONVECT_MIXING_FUNCT
   !!     Routine CONVECT_CONDENS
   !
   !
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !!      Module YOMCST
   !!          RG                 ! gravity constant
   !!          RATM               ! reference pressure
   !!          RD, RV           ! gaz  constants for dry air and water vapor
   !!          RCPD, RCPV, RCW    ! Cp of dry air, water vapor and liquid water
   !!          RTT                ! triple point temperature
   !!          RLVTT              ! vaporisation heat at RTT
   !
   !
   !!      Module YOE_CONVPAR
   !!          XA25               ! reference grid area
   !!          XCRAD              ! cloud radius
   !!          XCDEPTH            ! minimum necessary cloud depth
   !!          XENTR              ! entrainment constant
   !!          XRCONV             ! constant in precipitation conversion
   !!          XNHGAM             ! coefficient for buoyancy term in w eq.
   !!                             ! accounting for nh-pressure
   !!          XTFRZ1             ! begin of freezing interval
   !!          XTFRZ2             ! begin of freezing interval
   !
   !!     Module YOE_CONVPAREXT
   !!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
   !
   !!    REFERENCE
   !!    ---------
   !
   !!      Book1,2 of documentation ( routine CONVECT_UPDRAFT)
   !!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol.
   !!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
   !
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !
   !!    MODIFICATIONS
   !!    -------------
   !!      Original    07/11/95
   !!   Last modified  10/12/97
   !----------------------------------------------------------------------------


   !*       0.    DECLARATIONS
   !              ------------

   use YOMCST
   use YOE_CONVPAR
   use YOE_CONVPAREXT

   implicit none
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

   !*       0.1   Declarations of dummy arguments :

   integer, intent(IN)                    :: KLON  ! horizontal dimension
   integer, intent(IN)                    :: KLEV  ! vertical dimension
   integer, intent(IN)                    :: KICE  ! flag for ice ( 1 = yes,
   !                0 = no ice )
   real, dimension(KLON,KLEV), intent(IN) :: PTHL  ! grid scale enthalpy (J/kg)
   real, dimension(KLON,KLEV), intent(IN) :: PTHV  ! grid scale theta_v
   real, dimension(KLON,KLEV), intent(IN) :: PTHES ! grid scale saturated theta_e
   real, dimension(KLON,KLEV), intent(IN) :: PRW   ! grid scale total water
   ! mixing ratio
   real, dimension(KLON,KLEV), intent(IN) :: PPRES ! pressure (P)
   real, dimension(KLON,KLEV), intent(IN) :: PDPRES! pressure difference between
   ! bottom and top of layer (Pa)
   real, dimension(KLON,KLEV), intent(IN) :: PZ    ! height of model layer (m)
   real, dimension(KLON),     intent(IN) :: PTHLCL ! theta at LCL
   real, dimension(KLON),     intent(IN) :: PTLCL  ! temp. at LCL
   real, dimension(KLON),     intent(IN) :: PRVLCL ! vapor mixing ratio at  LCL
   real, dimension(KLON),     intent(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
   real, dimension(KLON),     intent(IN) :: PMFLCL ! cloud  base unit mass flux
   ! (kg/s)
   real, dimension(KLON),     intent(IN) :: PZLCL  ! height at LCL (m)
   real, dimension(KLON),     intent(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
   logical, dimension(KLON),  intent(INOUT):: OTRIG! logical mask for convection
   integer, dimension(KLON),  intent(IN) :: KLCL   ! contains vert. index of LCL
   integer, dimension(KLON),  intent(IN) :: KDPL   ! contains vert. index of DPL
   integer, dimension(KLON),  intent(IN) :: KPBL   !  " vert. index of source layertop


   integer, dimension(KLON),  intent(OUT):: KCTL   ! contains vert. index of CTL
   integer, dimension(KLON),  intent(OUT):: KETL   ! contains vert. index of
   !equilibrium (zero buoyancy) level
   real, dimension(KLON,KLEV), intent(OUT):: PUMF  ! updraft mass flux (kg/s)
   real, dimension(KLON,KLEV), intent(OUT):: PUER  ! updraft entrainment (kg/s)
   real, dimension(KLON,KLEV), intent(OUT):: PUDR  ! updraft detrainment (kg/s)
   real, dimension(KLON,KLEV), intent(OUT):: PUTHL ! updraft enthalpy (J/kg)
   real, dimension(KLON,KLEV), intent(OUT):: PUTHV ! updraft theta_v (K)
   real, dimension(KLON,KLEV), intent(OUT):: PWU   ! updraft vertical velocity (m/s)
   real, dimension(KLON,KLEV), intent(OUT):: PURW  ! updraft total water (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT):: PURC  ! updraft cloud water (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT):: PURI  ! updraft cloud ice   (kg/kg)
   real, dimension(KLON,KLEV), intent(OUT):: PURR  ! liquid precipit. (kg/kg)
   ! produced in  model layer
   real, dimension(KLON,KLEV),   intent(OUT)::PURS ! solid precipit. (kg/kg)
   ! produced in  model layer
   real, dimension(KLON,KLEV),   intent(OUT)::PUPR ! updraft precipitation in
   ! flux units (kg water / s)
   real, dimension(KLON),     intent(OUT):: PUTPR  ! total updraft precipitation
   ! in flux units (kg water / s)
   real, dimension(KLON),     intent(OUT):: PCAPE  ! available potent. energy

   !*       0.2   Declarations of local variables :

   integer :: IIE, IKB, IKE  ! horizontal and vertical loop bounds
   integer :: JI             ! horizontal loop index
   integer :: JK, JKP, JKM, JK1, JK2, JKMIN  ! vertical loop index
   real    :: ZEPSA, ZCVOCD  ! R_v / R_d, C_pv / C_pd
   real    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd

   real, dimension(KLON)    :: ZUT             ! updraft temperature (K)
   real, dimension(KLON)    :: ZUW1, ZUW2      ! square of updraft vert.
   ! velocity at levels k and k+1
   real, dimension(KLON)    :: ZE1,ZE2,ZD1,ZD2 ! fractional entrainm./detrain
   ! rates at levels k and k+1
   real, dimension(KLON)    :: ZMIXF           ! critical mixed fraction
   real, dimension(KLON)    :: ZCPH            ! specific heat C_ph
   real, dimension(KLON)    :: ZLV, ZLS        ! latent heat of vaporis., sublim.
   real, dimension(KLON)    :: ZURV            ! updraft water vapor at level k+1
   real, dimension(KLON)    :: ZPI             ! Pi=(P0/P)**(Rd/Cpd)
   real, dimension(KLON)    :: ZTHEUL          ! theta_e for undilute ascent
   real, dimension(KLON)    :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
   real, dimension(KLON)    :: ZWORK4, ZWORK5, ZWORK6  ! work arrays
   integer, dimension(KLON) :: IWORK           ! wok array
   logical, dimension(KLON)      :: GWORK1, GWORK2, GWORK4, GWORK5 ! work arrays
   logical, dimension(KLON,KLEV) :: GWORK6       ! work array


   !----------------------------------------------------------------------------

   !        0.3   Set loop bounds
   !              ---------------

   IKB = 1 + JCVEXB
   IKE = KLEV - JCVEXT
   IIE = KLON


   !*       1.     Initialize updraft properties and local variables
   !               -------------------------------------------------

   ZEPSA      = RV  / RD
   ZCVOCD     = RCPV / RCPD
   ZCPORD     = RCPD / RD
   ZRDOCP     = RD  / RCPD

   PUMF(:,:)  = _ZERO_
   PUER(:,:)  = _ZERO_
   PUDR(:,:)  = _ZERO_
   PUTHL(:,:) = _ZERO_
   PUTHV(:,:) = _ZERO_
   PWU(:,:)   = _ZERO_
   PURW(:,:)  = _ZERO_
   PURC(:,:)  = _ZERO_
   PURI(:,:)  = _ZERO_
   PUPR(:,:)  = _ZERO_
   PURR(:,:)  = _ZERO_
   PURS(:,:)  = _ZERO_
   PUTPR(:)   = _ZERO_
   ZUW1(:)    = PWLCL(:) * PWLCL(:)
   ZUW2(:)    = _ZERO_
   ZE1(:)     = _ONE_
   ZD1(:)     = _ZERO_
   PCAPE(:)   = _ZERO_
   KCTL(:)    = IKB
   KETL(:)    = KLCL(:)
   GWORK2(:)  = .true.
   GWORK5(:)  = .true.
   ZPI(:)     = _ONE_
   ZWORK3(:)  = _ZERO_
   ZWORK4(:)  = _ZERO_
   ZWORK5(:)  = _ZERO_
   ZWORK6(:)  = _ZERO_
   GWORK1(:)  = .false.
   GWORK4(:)  = .false.


   !*       1.1    Compute undilute updraft theta_e for CAPE computations
   !               Bolton (1980) formula.
   !               Define accurate enthalpy for updraft
   !               -----------------------------------------------------

   ZTHEUL(:) = PTLCL(:) * ( PTHLCL(:) / PTLCL(:) ) ** &
        &            ( 1. - 0.28 * PRVLCL(:) ) &
        &          * exp( ( 3374.6525 / PTLCL(:) - 2.5403 )      &
        &          * PRVLCL(:) * ( _ONE_ + 0.81 * PRVLCL(:) ) )


   ZWORK1(:) = ( RCPD + PRVLCL(:) * RCPV ) *  PTLCL(:)   &
        & + ( _ONE_ + PRVLCL(:) ) * RG * PZLCL(:)


   !*       2.     Set updraft properties between DPL and LCL
   !               ------------------------------------------

   JKP = maxval( KLCL(:) )
   JKM = minval( KDPL(:) )
   do JK = JKM, JKP
      do JI = 1, IIE
         if ( JK >= KDPL(JI) .and. JK < KLCL(JI) ) then
            PUMF(JI,JK)  = PMFLCL(JI)
            PUTHL(JI,JK) = ZWORK1(JI)
            PUTHV(JI,JK) = PTHLCL(JI) * ( _ONE_ + ZEPSA * PRVLCL(JI) ) /        &
                 &            ( _ONE_ + PRVLCL(JI) )
            PURW(JI,JK)  = PRVLCL(JI)
         endif
      enddo
   enddo


   !*       3.     Enter loop for updraft computations
   !               ------------------------------------

   JKMIN = minval( KLCL(:) - 1 )
   do JK = max( IKB + 1, JKMIN ), IKE - 1
      ZWORK6(:) = _ONE_
      JKP = JK + 1

      GWORK4(:) = JK >= KLCL(:) - 1
      GWORK1(:) = GWORK4(:) .and. GWORK2(:) ! this mask is used to confine
      ! updraft computations between the LCL and the CTL

      where( JK == KLCL(:) - 1 ) ZWORK6(:) = _ZERO_ ! factor that is used in buoyancy
      ! computation at first level above LCL


      !*       4.     Estimate condensate, L_v L_i, Cph and theta_v at level k+1
      !               ----------------------------------------------------------

      ZWORK1(:) = PURC(:,JK) + PURR(:,JK)
      ZWORK2(:) = PURI(:,JK) + PURS(:,JK)
      call CONVECT_CONDENS1(KLON, KICE, PPRES(:,JKP), PUTHL(:,JK), PURW(:,JK), &
           & ZWORK1, ZWORK2, PZ(:,JKP), ZUT, ZURV,     &
           & PURC(:,JKP), PURI(:,JKP), ZLV, ZLS, ZCPH)


      ZPI(:) = ( RATM / PPRES(:,JKP) ) ** ZRDOCP
      where ( GWORK1(:) )

         PUTHV(:,JKP) = ZPI(:) * ZUT(:) * ( _ONE_ + ZEPSA * ZURV(:) ) &
              &       / ( _ONE_ + PURW(:,JK) )


         !*       5.     Compute square of vertical velocity using entrainment
         !               at level k
         !               -----------------------------------------------------

         ZWORK3(:) = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -         &
              &      ( _ONE_ - ZWORK6(:) ) * PZLCL(:)        ! level thickness
         ZWORK4(:) = PTHV(:,JK) * ZWORK6(:) +                   &
              &  ( _ONE_ - ZWORK6(:) ) * PTHVELCL(:)
         ZWORK5(:) = _TWO_ * ZUW1(:) * PUER(:,JK) / max( .1, PUMF(:,JK) )
         ZUW2(:)   = ZUW1(:) + ZWORK3(:) * XNHGAM * RG *        &
              &   ( ( PUTHV(:,JK) + PUTHV(:,JKP) ) /       &
              &   ( ZWORK4(:) + PTHV(:,JKP) ) - _ONE_ )    & ! buoyancy term
              & - ZWORK5(:)                                  ! entrainment term


         !*       6.     Update total precipitation: dr_r=(r_c+r_i)*exp(-rate*dz)
         !               -------------------------------------------------------

         !                    compute level mean vertical velocity
         ZWORK2(:)   = _HALF_ *                                        &
              &      ( sqrt( max( 1.E-2, ZUW2(:) ) ) +                 &
              &        sqrt( max( 1.E-2, ZUW1(:) ) ) )
         PURR(:,JKP) = _HALF_ * ( PURC(:,JK) + PURC(:,JKP) + PURI(:,JK) + PURI(:,JKP) )&
              &     * ( _ONE_ - exp( - XRCONV  * ZWORK3(:) / ZWORK2(:) ) )
         PUPR(:,JKP) = PURR(:,JKP) * PUMF(:,JK) ! precipitation rate ( kg water / s)
         PUTPR(:)    = PUTPR(:) + PUPR(:,JKP)   ! total precipitation rate
         ZWORK2(:)   = PURR(:,JKP) / max( 1.E-8, PURC(:,JKP) + PURI(:,JKP) )
         PURR(:,JKP) = ZWORK2(:) * PURC(:,JKP)          ! liquid precipitation
         PURS(:,JKP) = ZWORK2(:) * PURI(:,JKP)          ! solid precipitation


         !*       7.     Update r_c, r_i, enthalpy, r_w  for precipitation
         !               -------------------------------------------------------

         PURW(:,JKP)  = PURW(:,JK) - PURR(:,JKP) - PURS(:,JKP)
         PURC(:,JKP)  = PURC(:,JKP) - PURR(:,JKP)
         PURI(:,JKP)  = PURI(:,JKP) - PURS(:,JKP)
         PUTHL(:,JKP) = ( RCPD + PURW(:,JKP) * RCPV ) * ZUT(:)             &
              & + ( _ONE_ + PURW(:,JKP) ) * RG * PZ(:,JKP)                 &
              & - ZLV(:) * PURC(:,JKP) - ZLS(:) * PURI(:,JKP)

         ZUW1(:)      = ZUW2(:)
         PWU(:,JKP)   = sqrt(max(ZUW1(:),1.E-2))

      end where


      !*       8.     Compute entrainment and detrainment using conservative
      !               variables adjusted for precipitation ( not for entrainment)
      !               ----------------------------------------------------------

      !*       8.1    Compute critical mixed fraction by estimating unknown
      !               T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
      !               We determine the zero crossing of the linear curve
      !               evaluating the derivative using ZMIXF=0.1.
      !               -----------------------------------------------------

      ZMIXF(:)  = 0.1   ! starting value for critical mixed fraction
      !Jing    ZWORK1(:) = ZMIXF(:) * PTHL(:,JKP)                                     &
      !Jing              &      + ( _ONE_ - ZMIXF(:) ) * PUTHL(:,JKP) ! mixed enthalpy
      !Jing    ZWORK2(:) = ZMIXF(:) * PRW(:,JKP)                                      &
      !Jing              &      + ( _ONE_ - ZMIXF(:) ) * PURW(:,JKP)  ! mixed r_w

      !Jing, CONVECT_CONDENS uses  r_w, r_c and r_i at lower level value to give a first temperature estimate, then use iteration to calculate T,r_v(water saturation mixing ratio), r_c(cloud water mixing ratio), r_i(cloud ice mixing ratio), Lv, Ls, CPH at upper levels.
      ZWORK1(:) = ZMIXF(:) * PTHL(:,JK)                                     &
           &      + ( _ONE_ - ZMIXF(:) ) * PUTHL(:,JK) ! mixed enthalpy
      ZWORK2(:) = ZMIXF(:) * PRW(:,JK)                                      &
           &      + ( _ONE_ - ZMIXF(:) ) * PURW(:,JK)  ! mixed r_w

      call CONVECT_CONDENS1(KLON, KICE, PPRES(:,JKP), ZWORK1, ZWORK2,       &
           & PURC(:,JK), PURI(:,JK), PZ(:,JKP), ZUT, &
           & ZWORK3, ZWORK4, ZWORK5, ZLV, ZLS, ZCPH)
      !  ZWORK3:water vapor mixing ratio, ZWORK4/5-cloud water/ICE mixing ratio, ZWORK4
      !        put in enthalpy and r_w and get T r_c, r_i (ZUT, ZWORK4-5)

      ! compute theta_v of mixture
      !Jing    ZWORK3(:) = ZUT(:) * ZPI(:) * ( _ONE_ + ZEPSA * (                         &
      !Jing              & ZWORK2(:) - ZWORK4(:) - ZWORK5(:) ) ) / ( _ONE_ + ZWORK2(:) )
      !Jing, (ZWORK2-ZWORK4-ZWORK5) is (r_w - r_c - r_i = r_v), however, ZWORK2 is JK, ZWORK4&ZWORK5 is JKP, they should be from same lvl JKP
      !Jing, therefore, use ZWORK3 (updraft water vapor at JKP)
      ZWORK3(:) = ZUT(:) * ZPI(:) * ( _ONE_ + ZEPSA * (                        &
           & ZWORK3(:) ) ) / ( _ONE_ + ZWORK2(:) )

      ! compute final value of critical mixed fraction using theta_v
      ! of mixture, grid-scale and updraft
      ZMIXF(:) = max( _ZERO_, PUTHV(:,JKP) - PTHV(:,JKP) ) * ZMIXF(:) /        &
           &                ( PUTHV(:,JKP) - ZWORK3(:) + 1.E-10 )

      ZMIXF(:) = max( _ZERO_, min( _ONE_, ZMIXF(:) ) )


      !*       8.2     Compute final midlevel values for entr. and detrainment
      !                after call of distribution function
      !                -------------------------------------------------------


      call CONVECT_MIXING_FUNCT ( KLON, ZMIXF, 1, ZE2, ZD2 )
      !       Note: routine MIXING_FUNCT returns fractional entrainm/detrainm. rates

      ZWORK1(:) = XENTR * PMFLCL(:) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
      !*MOD
      !Jingorg       ZWORK1(:) = XENTR * RG / XCRAD * PUMF(:,JK) * ( PZ(:,JKP) - PZ(:,JK) )
      ! ZWORK1(:) = XENTR * pumf(:,jk) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
      !*MOD
      ZWORK2(:) = _ZERO_
      where ( GWORK1(:) ) ZWORK2(:) = _ONE_
      !PV -following line is commented to get BKF to look like kfc
      !  ze2=_HALF_; zd2=_HALF_  ! modif entrainment=detrainment, this avoids
      ! too large mass flux values at upper levels
      where ( PUTHV(:,JKP) > PTHV(:,JKP) )
         PUER(:,JKP) = _HALF_ * ZWORK1(:) * ( ZE1(:) + ZE2(:) ) * ZWORK2(:)
         PUDR(:,JKP) = _HALF_ * ZWORK1(:) * ( ZD1(:) + ZD2(:) ) * ZWORK2(:)
      elsewhere
         PUER(:,JKP) = _ZERO_
         PUDR(:,JKP) = ZWORK1(:) * ZWORK2(:)
      end where




      !*       8.3     Determine equilibrium temperature level
      !                --------------------------------------

      where ( PUTHV(:,JKP) > PTHV(:,JKP) .and. JK > KLCL(:) + 1 &
           & .and. GWORK1(:) )
         KETL(:) = JKP            ! equilibrium temperature level
      end where

      !*       8.4     If the calculated detrained mass flux is greater than
      !                the total updraft mass flux, or vertical velocity is
      !                negative, all cloud mass detrains at previous model level,
      !                exit updraft calculations - CTL is attained
      !                -------------------------------------------------------

      where( GWORK1(:) )                                                   &
           & GWORK2(:) = PUMF(:,JK) - PUDR(:,JKP) > 10. .and. ZUW2(:) > _ZERO_
      where ( GWORK2(:) ) KCTL(:) = JKP   ! cloud top level
      GWORK1(:) = GWORK2(:) .and. GWORK4(:)

      if ( count( GWORK2(:) ) == 0 ) exit


      !*       9.   Compute CAPE for undilute ascent using theta_e and
      !             theta_es instead of theta_v. This estimation produces
      !             a significantly larger value for CAPE than the actual one.
      !             ----------------------------------------------------------

      where ( GWORK1(:) )

         ZWORK3(:)   = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -                      &
              & ( _ONE_ - ZWORK6(:) ) *  PZLCL(:)              ! level thickness
         ZWORK2(:)   = PTHES(:,JK) + ( _ONE_ - ZWORK6(:) ) *                   &
              &  ( PTHES(:,JKP) - PTHES(:,JK) ) / ( PZ(:,JKP) - PZ(:,JK) ) *          &
              &  ( PZLCL(:) - PZ(:,JK) ) ! linear interpolation for theta_es at LCL
         ! ( this is only done for model level just above LCL

         ZWORK1(:) = ( _TWO_ * ZTHEUL(:) ) / ( ZWORK2(:) + PTHES(:,JKP) ) - _ONE_
         PCAPE(:)  = PCAPE(:) + RG * ZWORK3(:) * max( _ZERO_, ZWORK1(:) )


         !*       10.   Compute final values of updraft mass flux, enthalpy, r_w
         !              at level k+1
         !              --------------------------------------------------------

         PUMF(:,JKP)  = PUMF(:,JK) - PUDR(:,JKP) + PUER(:,JKP)
         PUMF(:,JKP)  = max( PUMF(:,JKP), 0.1 )
         PUTHL(:,JKP) = ( PUMF(:,JK) * PUTHL(:,JK) +                              &
              &   PUER(:,JKP) * PTHL(:,JK) - PUDR(:,JKP) * PUTHL(:,JK) )  &
              &  / PUMF(:,JKP) + PUTHL(:,JKP) - PUTHL(:,JK)
         PURW(:,JKP)  = ( PUMF(:,JK) * PURW(:,JK) +                               &
              &   PUER(:,JKP) * PRW(:,JK) - PUDR(:,JKP) * PURW(:,JK) )    &
              &  / PUMF(:,JKP) - PURR(:,JKP) - PURS(:,JKP)

         ZE1(:) = ZE2(:) ! update fractional entrainment/detrainment
         ZD1(:) = ZD2(:)

      end where

   enddo

   !*       12.1    Set OTRIG to False if cloud thickness < XCDEPTH
   !                or CAPE < 1
   !                ------------------------------------------------

   do JI = 1, IIE
      JK  = KCTL(JI)
      OTRIG(JI) = PZ(JI,JK) - PZLCL(JI) >= XCDEPTH               &
           & .and. PCAPE(JI) > _ONE_
   enddo
   where( .not. OTRIG(:) )
      KCTL(:) = IKB
   end where
   KETL(:) = max( KETL(:), KLCL(:) + 2 )
   KETL(:) = min( KETL(:), KCTL(:) )


   !*       12.2    If the ETL and CTL are the same detrain updraft mass
   !                flux at this level
   !                -------------------------------------------------------

   ZWORK1(:) = _ZERO_
   where ( KETL(:) == KCTL(:) ) ZWORK1(:) = _ONE_

   do JI = 1, IIE
      JK = KETL(JI)
      PUDR(JI,JK)   = PUDR(JI,JK) +                                    &
           &       ( PUMF(JI,JK) - PUER(JI,JK) )  * ZWORK1(JI)
      PUER(JI,JK)   = PUER(JI,JK) * ( _ONE_ - ZWORK1(JI) )
      PUMF(JI,JK)   = PUMF(JI,JK) * ( _ONE_ - ZWORK1(JI) )
      JKP = KCTL(JI) + 1
      PUER(JI,JKP)  = _ZERO_ ! entrainm/detr rates have been already computed
      PUDR(JI,JKP)  = _ZERO_ ! at level KCTL+1, set them to zero
      PURW(JI,JKP)  = _ZERO_
      PURC(JI,JKP)  = _ZERO_
      PURI(JI,JKP)  = _ZERO_
      PUTHL(JI,JKP) = _ZERO_
      PURC(JI,JKP+1)= _ZERO_
      PURI(JI,JKP+1)= _ZERO_
   enddo

   !*       12.3    Adjust mass flux profiles, detrainment rates, and
   !                precipitation fallout rates to reflect linear decrease
   !                in mass flux between the ETL and CTL
   !                -------------------------------------------------------

   ZWORK1(:) = _ZERO_
   JK1 = minval( KETL(:) )
   JK2 = maxval( KCTL(:) )
   do JK = JK1, JK2
      do JI = 1, IIE
         if( JK > KETL(JI) .and. JK <= KCTL(JI) ) then
            ZWORK1(JI) = ZWORK1(JI) + PDPRES(JI,JK)
         endif
      enddo
   enddo

   do JI = 1, IIE
      JK = KETL(JI)
      ZWORK1(JI) = PUMF(JI,JK) / max( _ONE_, ZWORK1(JI) )
   enddo

   do JK = JK1 + 1, JK2
      JKP = JK - 1
      do JI = 1, IIE
         if ( JK > KETL(JI) .and. JK <= KCTL(JI) ) then
            ! PUTPR(JI)    = PUTPR(JI) - ( PURR(JI,JK) + PURS(JI,JK) ) * PUMF(JI,JKP)
            PUTPR(JI)    = PUTPR(JI) - PUPR(JI,JK)
            PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
            PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
            PUPR(JI,JK)  = PUMF(JI,JKP) * ( PURR(JI,JK) + PURS(JI,JK) )
            PUTPR(JI)    = PUTPR(JI) + PUPR(JI,JK)
         endif
      enddo
   enddo

   !         12.4   Set mass flux and entrainment in the source layer.
   !                Linear increase throughout the source layer.
   !                -------------------------------------------------------

   !IWORK(:) = MIN( KPBL(:), KLCL(:) - 1 )
   IWORK(:) = KPBL(:)
   do JI = 1, IIE
      JK  = KDPL(JI)
      JKP = IWORK(JI)
      !          mixed layer depth
      ZWORK2(JI) = PPRES(JI,JK) - PPRES(JI,JKP) + PDPRES(JI,JK)
   enddo


   JKP = maxval( IWORK(:) )
   do JK = JKM, JKP
      do JI = 1, IIE
         if ( JK >= KDPL(JI)  .and. JK <= IWORK(JI) ) then
            PUER(JI,JK) = PUER(JI,JK) + PMFLCL(JI) * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
            PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
         endif
      enddo
   enddo


   !*       13.   If cloud thickness is smaller than  3 km, no
   !              convection is allowed
   !              Nota: For technical reasons, we stop the convection
   !                    computations in this case and do not go back to
   !                    TRIGGER_FUNCT to look for the next unstable LCL
   !                    which could produce a thicker cloud.
   !              ---------------------------------------------------

   GWORK6(:,:) = spread( OTRIG(:), DIM=2, NCOPIES=KLEV )
   where ( .not. OTRIG(:) ) PUTPR(:) = _ZERO_
   where ( .not. GWORK6(:,:) )
      PUMF(:,:)  = _ZERO_
      PUDR(:,:)  = _ZERO_
      PUER(:,:)  = _ZERO_
      PWU(:,:)   = _ZERO_
      PUTHL(:,:) = PTHL(:,:)
      PURW(:,:)  = PRW(:,:)
      PUPR(:,:)  = _ZERO_
      PURC(:,:)  = _ZERO_
      PURI(:,:)  = _ZERO_
      PURR(:,:)  = _ZERO_
      PURS(:,:)  = _ZERO_
   end where

end subroutine CONVECT_UPDRAFT

