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

subroutine CONVECT_UPDRAFT_SHAL3(KLON, KLEV, &
     & KICE, PPRES, PDPRES, PZ, PTT, PTHL, PTHV, PTHES, PRV, PRW, &
     & PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,  &
     & PMFLCL, PCRAD, OTRIG, KLCL, KDPL, KPBL,         &
     & PUMF, PUER, PUDR, PUTHL, PUTHV, PWU,            &
     & PURW, PURC, PURI, PCAPE, KCTL, KETL )
   !############################################################################

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
   !!      Module YOE_CONVPAR_SHAL
   !!          XA25               ! reference grid area
   !!          XCRAD              ! cloud radius
   !!          XCDEPTH            ! minimum necessary cloud depth
   !!          XCDEPTH_D          ! maximum allowed   cloud depth
   !!          XENTR              ! entrainment constant
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
   !---------------------------------------------------------------------------


   !*       0.    DECLARATIONS
   !              ------------

   use YOMCST
   use YOE_CONVPAR_SHAL
   use YOE_CONVPAREXT
   use cnv_options

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
   real, dimension(KLON,KLEV), intent(IN) :: PRV   ! grid scale water vapour mixing ratio
   real, dimension(KLON,KLEV), intent(IN) :: PRW   ! grid scale total water
   ! mixing ratio
   real, dimension(KLON,KLEV), intent(IN) :: PPRES ! pressure (P)
   real, dimension(KLON,KLEV), intent(IN) :: PDPRES! pressure difference between
   ! bottom and top of layer (Pa)
   real, dimension(KLON,KLEV), intent(IN) :: PZ    ! height of model layer (m)
   real, dimension(KLON,KLEV), intent(IN) :: PTT   ! grid scale temperature (K)
   real, dimension(KLON),     intent(IN) :: PTHLCL ! theta at LCL
   real, dimension(KLON),     intent(IN) :: PTLCL  ! temp. at LCL
   real, dimension(KLON),     intent(IN) :: PRVLCL ! vapor mixing ratio at  LCL
   real, dimension(KLON),     intent(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
   real, dimension(KLON),     intent(IN) :: PMFLCL ! cloud  base unit mass flux
   ! (kg/s)
   real, dimension(KLON),     intent(IN) :: PCRAD  ! cloud base radius (m)
   real, dimension(KLON),     intent(IN) :: PZLCL  ! height at LCL (m)
   real, dimension(KLON),     intent(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
   logical, dimension(KLON),  intent(INOUT):: OTRIG  ! logical mask for convection
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
   real, dimension(KLON),     intent(OUT):: PCAPE  ! available potent. energy

   !*       0.2   Declarations of local variables :

   integer :: IIE, IKB, IKE  ! horizontal and vertical loop bounds
   integer :: JI             ! horizontal loop index
   integer :: JK, JKP   ! vertical loop index
   real    :: ZEPSA, ZCVOCD  ! R_v / R_d, C_pv / C_pd
   real    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
   real    :: ZRH            ! relative humidity

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
   real, dimension(KLON)    :: ZPLCL           ! pressure at LCL
   real, dimension(KLON)    :: ZDFRAC          ! fractional detrainment rate
   real, dimension(KLON)    :: ZEFRAC          ! fractional entrainment rate
   real, dimension(KLON)    :: ZQVSLCL         ! saturation specific humidity at LCL
   real, dimension(KLON,KLEV) :: ZQVS          ! saturation specific humidity of environment

   real, dimension(KLON)    :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
   real, dimension(KLON)    :: ZWORK4, ZWORK5, ZWORK6  ! work arrays
   logical, dimension(KLON) :: GWORK1, GWORK2, GWORK5
   ! work arrays
   logical, dimension(KLON,KLEV) :: GWORK6       ! work array


   !----------------------------------------------------------------------------

   !        0.3   Set loop bounds
   !              ---------------

   IKB = 1 + JCVEXB
   IKE = KLEV - JCVEXT
   IIE = KLON


   !*       1.     Initialize updraft properties and local variables
   !               -------------------------------------------------

   ZEPSA      = RV / RD
   ZCVOCD     = RCPV / RCPD
   ZCPORD     = RCPD / RD
   ZRDOCP     = RD / RCPD

   PUMF(:,:)  = _ZERO_
   PUER(:,:)  = _ZERO_
   PUDR(:,:)  = _ZERO_
   PUTHL(:,:) = _ZERO_
   PUTHV(:,:) = _ZERO_
   PWU(:,:)   = _ZERO_
   PURW(:,:)  = _ZERO_
   PURC(:,:)  = _ZERO_
   PURI(:,:)  = _ZERO_
   ZUW1(:)    = PWLCL(:) * PWLCL(:)
   ZUW2(:)    = _ZERO_
   ZE1(:)     = _ZERO_
   ZE2(:)     = _ZERO_
   ZD1(:)     = _ZERO_
   ZD2(:)     = _ZERO_
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

   !*       1.1    Compute undilute updraft theta_e for CAPE computations
   !               Bolton (1980) formula.
   !               Define accurate enthalpy for updraft
   !               -----------------------------------------------------

   ZTHEUL(:) = PTLCL(:) * ( PTHLCL(:) / PTLCL(:) ) ** ( 1. - 0.28 * PRVLCL(:) )  &
        &          * exp( ( 3374.6525 / PTLCL(:) - 2.5403 )       &
        &          * PRVLCL(:) * ( _ONE_ + 0.81 * PRVLCL(:) ) )


   ZWORK1(:) = ( RCPD + PRVLCL(:) * RCPV ) * PTLCL(:)                            &
        & + ( _ONE_ + PRVLCL(:) ) * RG * PZLCL(:)

   zplcl = RATM * (ptlcl/pthlcl)**(1./ZRDOCP)
   call mfoqst3(zqvs,ptt,ppres,iie,ike,iie)
   call mfoqst3(zqvslcl,ptlcl,zplcl,iie,1,iie)

   !*       2.     Set updraft properties between DPL and LCL
   !               ------------------------------------------

   do JI = 1, IIE
      do JK = KDPL(JI), KLCL(JI)
         !IF ( JK >= KDPL(JI) .AND. JK < KLCL(JI) ) THEN
         if ( JK < KLCL(JI) ) then
            PUMF(JI,JK)  = PMFLCL(JI)
            PUTHL(JI,JK) = ZWORK1(JI)
            PUTHV(JI,JK) = PTHLCL(JI) * ( _ONE_ + ZEPSA * PRVLCL(JI) ) /         &
                 &            ( _ONE_ + PRVLCL(JI) )
            PURW(JI,JK)  = PRVLCL(JI)
            PWU(JI,JK) = PWLCL(JI)
         endif
      enddo
   enddo

   !*       3.     Enter loop for updraft computations
   !               ------------------------------------

   !JKMIN = MINVAL( KLCL(:) - 1 ) !viv definitely not MPI reprod
   !JKMIN = 0

   ILOOP: do JI = 1, IIE
      UPDRAFT_LOOP: do JK = max( IKB + 1, KLCL(JI)-1 ), IKE - 1
         JKP = JK + 1
         ZWORK6(JI) = _ONE_

         GWORK1(JI) = GWORK2(JI) ! this mask is used to confine
         ! updraft computations between the LCL and the CTL

         if ( JK == KLCL(JI) - 1 ) ZWORK6(JI) = _ZERO_
         ! factor that is used in buoyancy
         ! computation at first level above LCL

         !*       4.     Estimate condensate, L_v L_i, Cph and theta_v at level k+1
         !               -------------------------------------------------------

         ZWORK1(JI) = PURC(JI,JK)
         ZWORK2(JI) = PURI(JI,JK)
         call CONVECT_CONDENS1(1, KICE, PPRES(JI,JKP), PUTHL(JI,JK), PURW(JI,JK),&
              & ZWORK1(JI), ZWORK2(JI), PZ(JI,JKP), ZUT(JI), ZURV(JI),   &
              & PURC(JI,JKP), PURI(JI,JKP), ZLV(JI), ZLS(JI), ZCPH(JI))


         ZPI(JI) = ( RATM / PPRES(JI,JKP) ) ** ZRDOCP

         if ( GWORK1(JI) ) then

            PUTHV(JI,JKP) = ZPI(JI) * ZUT(JI) * ( _ONE_ + ZEPSA * ZURV(JI) )  &
                 / ( _ONE_ + PURW(JI,JK) )

            !*       5.     Compute square of vertical velocity using entrainment
            !               at level k
            !               ----------------------------------------------------

            ZWORK3(JI) = PZ(JI,JKP) - PZ(JI,JK) * ZWORK6(JI) -      &
                 ( _ONE_ - ZWORK6(JI) ) * PZLCL(JI)        ! level thickness
            ZWORK4(JI) = PTHV(JI,JK) * ZWORK6(JI) +                 &
                 ( _ONE_ - ZWORK6(JI) ) * PTHVELCL(JI)
            ZWORK5(JI) = _TWO_ * ZUW1(JI) * PUER(JI,JK) / max( .1, PUMF(JI,JK) )
            ZUW2(JI)   = ZUW1(JI) + ZWORK3(JI) * XNHGAM * RG *      &
                 ( ( PUTHV(JI,JK) + PUTHV(JI,JKP) ) /       &
                 ( ZWORK4(JI) + PTHV(JI,JKP) ) - _ONE_ )    & ! buoyancy term
                 - ZWORK5(JI)                               ! entrainment term

            !*       6.     Update total precipitationJI dr_r=(r_c+r_i)*exp(-rate*dz)
            !               ----------------------------------------------------
            !                    compute level mean vertical velocity
            ZWORK2(JI)   = _HALF_ *                                 &
                 ( sqrt( max( 1.E-2, ZUW2(JI) ) ) +  &
                 sqrt( max( 1.E-2, ZUW1(JI) ) ) )

            !*       7.     Update r_c, r_i, enthalpy, r_w  for precipitation
            !               ----------------------------------------------------

            PURW(JI,JKP)  = PURW(JI,JK)
            PUTHL(JI,JKP) = PUTHL(JI,JK)
            ZUW1(JI)      = ZUW2(JI)
            PWU(JI,JKP)   =sqrt(max(ZUW1(JI),1.E-2))
         endif


         !*       8.     Compute entrainment and detrainment using conservative
         !               variables adjusted for precipitation ( not for entrainment)
         !               -------------------------------------------------------

         !*       8.1    Compute critical mixed fraction by estimating unknown
         !               T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
         !               We determine the zero crossing of the linear curve
         !               evaluating the derivative using ZMIXF=0.1.
         !               -----------------------------------------------------

         MIXING_FRACTION: if (any((/bkf_entrains,bkf_detrains/) == 'BECHTOLD01')) then

            ZMIXF(JI)  = 0.1   ! starting value for critical mixed fraction
            ZWORK1(JI) = ZMIXF(JI) * PTHL(JI,JKP)                                &
                 + ( _ONE_ - ZMIXF(JI) ) * PUTHL(JI,JKP) ! mixed enthalpy
            ZWORK2(JI) = ZMIXF(JI) * PRW(JI,JKP)                                 &
                 + ( _ONE_ - ZMIXF(JI) ) * PURW(JI,JKP)  ! mixed r_w

            call CONVECT_CONDENS1(1, KICE, PPRES(JI,JKP), ZWORK1(JI), ZWORK2(JI),&
                 PURC(JI,JKP), PURI(JI,JKP), PZ(JI,JKP), ZUT(JI),&
                 ZWORK3(JI), ZWORK4(JI), ZWORK5(JI), ZLV(JI), ZLS(JI), ZCPH(JI) )
            !        put in enthalpy and r_w and get T r_c, r_i (ZUT, ZWORK4-5)

            ! compute theta_v of mixture
            ZWORK3(JI) = ZUT(JI) * ZPI(JI) * ( _ONE_ + ZEPSA * (  &
                 ZWORK2(JI) - ZWORK4(JI) - ZWORK5(JI) ) ) / ( _ONE_ + ZWORK2(JI) )

            ! compute final value of critical mixed fraction using theta_v
            ! of mixture, grid-scale and updraft
            ZMIXF(JI) = max( _ZERO_, PUTHV(JI,JKP) - PTHV(JI,JKP) ) * ZMIXF(JI) /&
                 ( PUTHV(JI,JKP) - ZWORK3(JI) + 1.E-10 )
            ZMIXF(JI) = max( _ZERO_, min( _ONE_, ZMIXF(JI) ) )

            call CONVECT_MIXING_FUNCT ( 1, ZMIXF(JI), 1, ZE2(JI), ZD2(JI) )
            !  Note: routine MIXING_FUNCT returns fractional entrainm/detrainm. rates

         endif MIXING_FRACTION

         !*       8.2     Compute final midlevel values for entr. and detrainment
         !                ------------------------------------------------------

         ! Compute adjustment mask

         ZWORK2(JI) = _ZERO_
         if ( GWORK1(JI) ) ZWORK2(JI) = _ONE_

         ! Buoyant parcel fractional entrainment/detrainment rates
         ENTRAINMENT: select case (bkf_entrains)
         case ('BECHTOLD01')
            zefrac(JI) = XENTR * RG / PCRAD(JI) * _HALF_ * (ze1(JI)+ze2(JI))
         case ('BECHTOLD08')
            zefrac(JI) = 1.8e-3 * (1.3-prv(JI,jkp)/zqvs(JI,jkp)) * (zqvs(JI,jkp)/zqvslcl(JI))**3
         case ('DEROOY11')
            zefrac(JI) = 1./((pz(JI,jk)+pz(JI,jkp))/2. - pzlcl(JI) + 500)
         case ('SIEBESMA03')
            zefrac(JI) = 1.0/((pz(JI,jk)+pz(JI,jkp))/2.)
         end select ENTRAINMENT
         DETRAINMENT: select case (bkf_detrains)
         case ('BECHTOLD01')
            zdfrac(JI) = XENTR * RG / PCRAD(JI) * _HALF_ * (zd1(JI)+zd2(JI))
         case ('CUIJPERS95')
            zdfrac(JI) = 2.75e-3
         case ('DEROOY10')
            zdfrac(JI) = zefrac(JI)
         end select DETRAINMENT

         ! Correction for negatively buoyant parcels
         if (.not. gwork1(ji)) cycle
         NEG_BUOYANT: if (puthv(ji,jkp) <= pthv(ji,jkp)) then
            zefrac(ji) = _ZERO_
            select case (bkf_detrains)
            case ('BECHTOLD01')
               zdfrac(ji) = XENTR * RG / PCRAD(ji)
            case DEFAULT
               zrh = min(prv(ji,jkp)/zqvs(ji,jkp),1.) ! do not allow environmental supersaturation
               zdfrac(ji) = zdfrac(ji) + (1./(pz(ji,jkp)-pz(ji,jk))) * &
                    max(((1.6-zrh)*(pwu(ji,jk)/pwu(ji,jkp)))**(-1)-1., 0.)
            end select
         endif NEG_BUOYANT

         ! Compute entrainment/detrainment-related updraft mass flux changes
         puer(JI,jkp) = zwork2(JI) * pumf(JI,jk) * (pz(JI,jkp)-pz(JI,jk)) * zefrac(JI)
         pudr(JI,jkp) = zwork2(JI) * pumf(JI,jk) * (pz(JI,jkp)-pz(JI,jk)) * zdfrac(JI)

         !*       8.3     Determine equilibrium temperature level
         !                --------------------------------------

         if (GWORK1(JI) ) then

            ! equilibrium temperature level
            if ( PUTHV(JI,JKP) > PTHV(JI,JKP) .and. JK > KLCL(JI) + 1 ) KETL(JI) = JKP

            !*       8.4     If the calculated detrained mass flux is greater than
            !                the total updraft mass flux, or vertical velocity is
            !                negative, all cloud mass detrains at previous model level,
            !                exit updraft calculations - CTL is attained
            !                -------------------------------------------------------

            GWORK2(JI) = PUMF(JI,JK) - PUDR(JI,JKP) > 10. .and. ZUW2(JI) > _ZERO_
            if ( GWORK2(JI) ) KCTL(JI) = JKP   ! cloud top level
            !else, CTL is attained and GWORK2 is set to False.
         endif

         !GWORK1(JI) = GWORK2(JI)
         if (.not.GWORK2(JI)) exit


         !*       9.   Compute CAPE for undilute ascent using theta_e and
         !             theta_es instead of theta_v. This estimation produces
         !             a significantly larger value for CAPE than the actual one.
         !             ----------------------------------------------------------

         if ( GWORK1(JI) ) then

            ZWORK3(JI)   = PZ(JI,JKP) - PZ(JI,JK) * ZWORK6(JI) -              &
                 ( _ONE_ - ZWORK6(JI) ) *  PZLCL(JI)         ! level thickness
            ! linear interpolation for theta_es at LCL
            ! ( this is only done for model level just above LCL
            ZWORK2(JI)   = PTHES(JI,JK) + ( _ONE_ - ZWORK6(JI) ) *            &
                 ( PTHES(JI,JKP) - PTHES(JI,JK) ) / ( PZ(JI,JKP) - PZ(JI,JK) ) *  &
                 ( PZLCL(JI) - PZ(JI,JK) )
            ZWORK1(JI) = ( _TWO_ * ZTHEUL(JI) ) / ( ZWORK2(JI) + PTHES(JI,JKP) ) - _ONE_
            PCAPE(JI)  = PCAPE(JI) + RG * ZWORK3(JI) * max( _ZERO_, ZWORK1(JI) )

            !*       10.   Compute final values of updraft mass flux, enthalpy, r_w
            !              at level k+1
            !              --------------------------------------------------------

            PUMF(JI,JKP)  = PUMF(JI,JK) - PUDR(JI,JKP) + PUER(JI,JKP)
            PUMF(JI,JKP)  = max( PUMF(JI,JKP), 0.1 )
            PUTHL(JI,JKP) = ( PUMF(JI,JK)  * PUTHL(JI,JK) +                             &
                 PUER(JI,JKP) * PTHL(JI,JK) - PUDR(JI,JKP) * PUTHL(JI,JK) ) &
                 / PUMF(JI,JKP)
            PURW(JI,JKP)  = ( PUMF(JI,JK)  * PURW(JI,JK) +                              &
                 PUER(JI,JKP) * PRW(JI,JK) - PUDR(JI,JKP) * PURW(JI,JK) )   &
                 / PUMF(JI,JKP)

            ZE1(JI) = ZE2(JI) ! update fractional entrainment/detrainment
            ZD1(JI) = ZD2(JI)

         endif

      enddo UPDRAFT_LOOP
   enddo ILOOP

   !*       12.1    Set OTRIG to False if cloud thickness < 0.5km
   !                or > 3km (deep convection) or CAPE < 1
   !                ------------------------------------------------

   do JI = 1, IIE
      JK  = KCTL(JI)
      ZWORK1(JI) = PZ(JI,JK) - PZLCL(JI)
      OTRIG(JI) = ZWORK1(JI) >= XCDEPTH  .and. ZWORK1(JI) < XCDEPTH_D &
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

   do JI = 1, IIE
      do JK = KETL(JI), KCTL(JI)
         if( JK > KETL(JI) ) then
            ZWORK1(JI) = ZWORK1(JI) + PDPRES(JI,JK)
         endif
      enddo
      ZWORK1(JI) = PUMF(JI,KETL(JI)) / max( _ONE_, ZWORK1(JI) )

      do JK = KETL(JI) + 1, KCTL(JI)
         JKP = JK - 1
         PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
         PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
      enddo
   enddo

   !         12.4   Set mass flux and entrainment in the source layer.
   !                Linear increase throughout the source layer.
   !                -------------------------------------------------------

   ! FIXME - the code above seems like a good idea.  If our LCL is very low (within the departure
   !         layer (KLCL < KPBL)) then we're making this post-hoc adjustment to within-cloud values,
   !         where they really only seem justified in the subcloud layer.  Note that once this is
   !         changed, the conditional for the PWU definition below can be removed.

   do JI = 1, IIE
      !  JK  = KDPL(JI) !  JKP = KPBL(JI)
      !  ZWORK2(JI) = PPRES(JI,JK) - PPRES(JI,JKP) + PDPRES(JI,JK)
      !  mixed layer depth
      ZWORK2(JI) = PPRES(JI,KDPL(JI)) - PPRES(JI,KPBL(JI)) + PDPRES(JI,KDPL(JI))
      do JK = KDPL(JI), KPBL(JI)
         PUER(JI,JK) = PUER(JI,JK) + PMFLCL(JI) * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
         PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
         if (JK < KLCL(JI)) PWU(JI,JK) = PWLCL(JI) * PUMF(JI,JK) / PMFLCL(JI)
      enddo
   enddo


   !*       13.   If cloud thickness is smaller than  .5 km or > 3 km
   !              no shallow convection is allowed
   !              Nota: For technical reasons, we stop the convection
   !                    computations in this case and do not go back to
   !                    TRIGGER_FUNCT to look for the next unstable LCL
   !                    which could produce a thicker cloud.
   !              ---------------------------------------------------

   GWORK6(:,:) = spread( OTRIG(:), DIM=2, NCOPIES=KLEV )
   where ( .not. GWORK6(:,:) )
      PUMF(:,:)  = _ZERO_
      PUDR(:,:)  = _ZERO_
      PUER(:,:)  = _ZERO_
      PWU(:,:)   = _ZERO_
      PUTHL(:,:) = PTHL(:,:)
      PURW(:,:)  = PRW(:,:)
      PURC(:,:)  = _ZERO_
      PURI(:,:)  = _ZERO_
   end where

end subroutine CONVECT_UPDRAFT_SHAL3

