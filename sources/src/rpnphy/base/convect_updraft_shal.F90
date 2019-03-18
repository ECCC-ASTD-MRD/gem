!
 SUBROUTINE CONVECT_UPDRAFT_SHAL3( KLON, KLEV,                                &
                           & KICE, PPRES, PDPRES, PZ, PTT, PTHL, PTHV, PTHES, PRV, PRW,&
                           & PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,  &
                           & PMFLCL, PCRAD, OTRIG, KLCL, KDPL, KPBL,         &
                           & PUMF, PUER, PUDR, PUTHL, PUTHV, PWU,           &
                           & PURW, PURC, PURI, PCAPE, KCTL, KETL )
!#############################################################################

!!**** Compute updraft properties from DPL to CTL.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine updraft properties
!!      ( mass flux, thermodynamics, precipitation )
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!
!!
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_MIXING_FUNCT
!!     Routine CONVECT_CONDENS
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RATM               ! reference pressure
!!          RD, RV           ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV, RCW    ! Cp of dry air, water vapor and liquid water
!!          RTT                ! triple point temperature
!!          RLVTT              ! vaporisation heat at RTT
!!
!!
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
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_UPDRAFT)
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  10/12/97
!-------------------------------------------------------------------------------


!*       0.    DECLARATIONS
!              ------------

USE YOMCST
USE YOE_CONVPAR_SHAL
USE YOE_CONVPAREXT
use cnv_options

IMPLICIT NONE
#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

!*       0.1   Declarations of dummy arguments :

integer, INTENT(IN)                    :: KLON  ! horizontal dimension
integer, INTENT(IN)                    :: KLEV  ! vertical dimension
integer, INTENT(IN)                    :: KICE  ! flag for ice ( 1 = yes,
                                                  !                0 = no ice )
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHL  ! grid scale enthalpy (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHV  ! grid scale theta_v
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRV   ! grid scale water vapour mixing ratio
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water
                                                  ! mixing ratio
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (P)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between
                                                  ! bottom and top of layer (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of model layer (m)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT   ! grid scale temperature (K)
real, DIMENSION(KLON),     INTENT(IN) :: PTHLCL ! theta at LCL
real, DIMENSION(KLON),     INTENT(IN) :: PTLCL  ! temp. at LCL
real, DIMENSION(KLON),     INTENT(IN) :: PRVLCL ! vapor mixing ratio at  LCL
real, DIMENSION(KLON),     INTENT(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
real, DIMENSION(KLON),     INTENT(IN) :: PMFLCL ! cloud  base unit mass flux
                                                  ! (kg/s)
real, DIMENSION(KLON),     INTENT(IN) :: PCRAD  ! cloud base radius (m) 
real, DIMENSION(KLON),     INTENT(IN) :: PZLCL  ! height at LCL (m)
real, DIMENSION(KLON),     INTENT(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(INOUT):: OTRIG  ! logical mask for convection
integer, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! contains vert. index of LCL
integer, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! contains vert. index of DPL
integer, DIMENSION(KLON),  INTENT(IN) :: KPBL   !  " vert. index of source layertop


integer, DIMENSION(KLON),  INTENT(OUT):: KCTL   ! contains vert. index of CTL
integer, DIMENSION(KLON),  INTENT(OUT):: KETL   ! contains vert. index of 
                                                  !equilibrium (zero buoyancy) level
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PUMF  ! updraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PUER  ! updraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PUDR  ! updraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHL ! updraft enthalpy (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHV ! updraft theta_v (K)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PWU   ! updraft vertical velocity (m/s)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PURW  ! updraft total water (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PURC  ! updraft cloud water (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(OUT):: PURI  ! updraft cloud ice   (kg/kg)
real, DIMENSION(KLON),     INTENT(OUT):: PCAPE  ! available potent. energy

!*       0.2   Declarations of local variables :

integer :: IIE, IKB, IKE  ! horizontal and vertical loop bounds
integer :: JI             ! horizontal loop index
integer :: JK, JKP, JKMIN   ! vertical loop index
real    :: ZEPSA, ZCVOCD  ! R_v / R_d, C_pv / C_pd
real    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
real    :: ZRH            ! relative humidity

real, DIMENSION(KLON)    :: ZUT             ! updraft temperature (K)
real, DIMENSION(KLON)    :: ZUW1, ZUW2      ! square of updraft vert.
                                              ! velocity at levels k and k+1
real, DIMENSION(KLON)    :: ZE1,ZE2,ZD1,ZD2 ! fractional entrainm./detrain
                                              ! rates at levels k and k+1
real, DIMENSION(KLON)    :: ZMIXF           ! critical mixed fraction
real, DIMENSION(KLON)    :: ZCPH            ! specific heat C_ph
real, DIMENSION(KLON)    :: ZLV, ZLS        ! latent heat of vaporis., sublim.
real, DIMENSION(KLON)    :: ZURV            ! updraft water vapor at level k+1
real, DIMENSION(KLON)    :: ZPI             ! Pi=(P0/P)**(Rd/Cpd)
real, DIMENSION(KLON)    :: ZTHEUL          ! theta_e for undilute ascent
real, DIMENSION(KLON)    :: ZPLCL           ! pressure at LCL
real, DIMENSION(KLON)    :: ZDFRAC          ! fractional detrainment rate
real, DIMENSION(KLON)    :: ZEFRAC          ! fractional entrainment rate
real, DIMENSION(KLON)    :: ZQVSLCL         ! saturation specific humidity at LCL
real, DIMENSION(KLON,KLEV) :: ZQVS          ! saturation specific humidity of environment

real, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
real, DIMENSION(KLON)    :: ZWORK4, ZWORK5, ZWORK6  ! work arrays
LOGICAL, DIMENSION(KLON) :: GWORK1, GWORK2, GWORK5
                                              ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK6       ! work array


!-------------------------------------------------------------------------------

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
GWORK2(:)  = .TRUE.
GWORK5(:)  = .TRUE.
ZPI(:)     = _ONE_
ZWORK3(:)  = _ZERO_
ZWORK4(:)  = _ZERO_
ZWORK5(:)  = _ZERO_
ZWORK6(:)  = _ZERO_
GWORK1(:)  = .FALSE.

!*       1.1    Compute undilute updraft theta_e for CAPE computations
!               Bolton (1980) formula.
!               Define accurate enthalpy for updraft
!               -----------------------------------------------------

ZTHEUL(:) = PTLCL(:) * ( PTHLCL(:) / PTLCL(:) ) ** ( 1. - 0.28 * PRVLCL(:) )  &
          &          * EXP( ( 3374.6525 / PTLCL(:) - 2.5403 )       &
          &          * PRVLCL(:) * ( _ONE_ + 0.81 * PRVLCL(:) ) )


ZWORK1(:) = ( RCPD + PRVLCL(:) * RCPV ) * PTLCL(:)                            &
          & + ( _ONE_ + PRVLCL(:) ) * RG * PZLCL(:)

zplcl = RATM * (ptlcl/pthlcl)**(1./ZRDOCP)
call mfoqst3(zqvs,ptt,ppres,iie,ike,iie)
call mfoqst3(zqvslcl,ptlcl,zplcl,iie,1,iie)

!*       2.     Set updraft properties between DPL and LCL
!               ------------------------------------------

DO JI = 1, IIE
   DO JK = KDPL(JI), KLCL(JI)
   !IF ( JK >= KDPL(JI) .AND. JK < KLCL(JI) ) THEN
   IF ( JK < KLCL(JI) ) THEN
        PUMF(JI,JK)  = PMFLCL(JI)
        PUTHL(JI,JK) = ZWORK1(JI)
        PUTHV(JI,JK) = PTHLCL(JI) * ( _ONE_ + ZEPSA * PRVLCL(JI) ) /         &
                     &            ( _ONE_ + PRVLCL(JI) )
        PURW(JI,JK)  = PRVLCL(JI)
        PWU(JI,JK) = PWLCL(JI)
   ENDIF
   ENDDO
ENDDO

!*       3.     Enter loop for updraft computations
!               ------------------------------------

!JKMIN = MINVAL( KLCL(:) - 1 ) !viv definitely not MPI reprod
!JKMIN = 0

ILOOP: DO JI = 1, IIE
  UPDRAFT_LOOP: DO JK = MAX( IKB + 1, KLCL(JI)-1 ), IKE - 1
    JKP = JK + 1
    ZWORK6(JI) = _ONE_

    GWORK1(JI) = GWORK2(JI) ! this mask is used to confine
                           ! updraft computations between the LCL and the CTL

    IF ( JK == KLCL(JI) - 1 ) ZWORK6(JI) = _ZERO_ 
    ! factor that is used in buoyancy
    ! computation at first level above LCL

!*       4.     Estimate condensate, L_v L_i, Cph and theta_v at level k+1
!               ----------------------------------------------------------

    ZWORK1(JI) = PURC(JI,JK)
    ZWORK2(JI) = PURI(JI,JK)
    CALL CONVECT_CONDENS( 1, KICE, PPRES(JI,JKP), PUTHL(JI,JK), PURW(JI,JK),&
         & ZWORK1(JI), ZWORK2(JI), PZ(JI,JKP), GWORK1(JI), ZUT(JI), ZURV(JI),   &
         & PURC(JI,JKP), PURI(JI,JKP), ZLV(JI), ZLS(JI), ZCPH(JI) )


    ZPI(JI) = ( RATM / PPRES(JI,JKP) ) ** ZRDOCP

    IF ( GWORK1(JI) ) THEN

      PUTHV(JI,JKP) = ZPI(JI) * ZUT(JI) * ( _ONE_ + ZEPSA * ZURV(JI) )  &
                        / ( _ONE_ + PURW(JI,JK) )

!*       5.     Compute square of vertical velocity using entrainment
!               at level k
!               -----------------------------------------------------

      ZWORK3(JI) = PZ(JI,JKP) - PZ(JI,JK) * ZWORK6(JI) -      &
                    ( _ONE_ - ZWORK6(JI) ) * PZLCL(JI)        ! level thickness
      ZWORK4(JI) = PTHV(JI,JK) * ZWORK6(JI) +                 &
                    ( _ONE_ - ZWORK6(JI) ) * PTHVELCL(JI)
      ZWORK5(JI) = _TWO_ * ZUW1(JI) * PUER(JI,JK) / MAX( .1, PUMF(JI,JK) )
      ZUW2(JI)   = ZUW1(JI) + ZWORK3(JI) * XNHGAM * RG *      &
                   ( ( PUTHV(JI,JK) + PUTHV(JI,JKP) ) /       &
                   ( ZWORK4(JI) + PTHV(JI,JKP) ) - _ONE_ )    & ! buoyancy term
                   - ZWORK5(JI)                               ! entrainment term

!*       6.     Update total precipitationJI dr_r=(r_c+r_i)*exp(-rate*dz)
!               --------------------------------------------------------
!                    compute level mean vertical velocity
      ZWORK2(JI)   = _HALF_ *                                 &
                     ( SQRT( MAX( 1.E-2, ZUW2(JI) ) ) +  &
                       SQRT( MAX( 1.E-2, ZUW1(JI) ) ) )

!*       7.     Update r_c, r_i, enthalpy, r_w  for precipitation
!               -------------------------------------------------------

      PURW(JI,JKP)  = PURW(JI,JK)
      PUTHL(JI,JKP) = PUTHL(JI,JK)
      ZUW1(JI)      = ZUW2(JI)
      PWU(JI,JKP)   =SQRT(MAX(ZUW1(JI),1.E-2))
    ENDIF


!*       8.     Compute entrainment and detrainment using conservative
!               variables adjusted for precipitation ( not for entrainment)
!               -----------------------------------------------------------

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
     
     CALL CONVECT_CONDENS( 1, KICE, PPRES(JI,JKP), ZWORK1(JI), ZWORK2(JI),&
           PURC(JI,JKP), PURI(JI,JKP), PZ(JI,JKP), GWORK1(JI), ZUT(JI),&
           ZWORK3(JI), ZWORK4(JI), ZWORK5(JI), ZLV(JI), ZLS(JI), ZCPH(JI) )
     !        put in enthalpy and r_w and get T r_c, r_i (ZUT, ZWORK4-5)

     ! compute theta_v of mixture
     ZWORK3(JI) = ZUT(JI) * ZPI(JI) * ( _ONE_ + ZEPSA * (  &
           ZWORK2(JI) - ZWORK4(JI) - ZWORK5(JI) ) ) / ( _ONE_ + ZWORK2(JI) )

     ! compute final value of critical mixed fraction using theta_v
     ! of mixture, grid-scale and updraft
     ZMIXF(JI) = MAX( _ZERO_, PUTHV(JI,JKP) - PTHV(JI,JKP) ) * ZMIXF(JI) /&
                           ( PUTHV(JI,JKP) - ZWORK3(JI) + 1.E-10 )
     ZMIXF(JI) = MAX( _ZERO_, MIN( _ONE_, ZMIXF(JI) ) )

     CALL CONVECT_MIXING_FUNCT ( 1, ZMIXF(JI), 1, ZE2(JI), ZD2(JI) )
     !  Note: routine MIXING_FUNCT returns fractional entrainm/detrainm. rates
     
  endif MIXING_FRACTION

!*       8.2     Compute final midlevel values for entr. and detrainment
!                -------------------------------------------------------

  ! Compute adjustment mask

    ZWORK2(JI) = _ZERO_
    IF ( GWORK1(JI) ) ZWORK2(JI) = _ONE_

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

  IF (GWORK1(JI) ) THEN 

     ! equilibrium temperature level
     IF ( PUTHV(JI,JKP) > PTHV(JI,JKP) .AND. JK > KLCL(JI) + 1 ) KETL(JI) = JKP 

!*       8.4     If the calculated detrained mass flux is greater than
!                the total updraft mass flux, or vertical velocity is
!                negative, all cloud mass detrains at previous model level,
!                exit updraft calculations - CTL is attained
!                -------------------------------------------------------

      GWORK2(JI) = PUMF(JI,JK) - PUDR(JI,JKP) > 10. .AND. ZUW2(JI) > _ZERO_
      IF ( GWORK2(JI) ) KCTL(JI) = JKP   ! cloud top level
      !else, CTL is attained and GWORK2 is set to False.
  ENDIF

  !GWORK1(JI) = GWORK2(JI)
  IF (.not.GWORK2(JI)) exit


!*       9.   Compute CAPE for undilute ascent using theta_e and
!             theta_es instead of theta_v. This estimation produces
!             a significantly larger value for CAPE than the actual one.
!             ----------------------------------------------------------

  IF ( GWORK1(JI) ) THEN

    ZWORK3(JI)   = PZ(JI,JKP) - PZ(JI,JK) * ZWORK6(JI) -              &
                  ( _ONE_ - ZWORK6(JI) ) *  PZLCL(JI)         ! level thickness
    ! linear interpolation for theta_es at LCL
    ! ( this is only done for model level just above LCL
    ZWORK2(JI)   = PTHES(JI,JK) + ( _ONE_ - ZWORK6(JI) ) *            &
     ( PTHES(JI,JKP) - PTHES(JI,JK) ) / ( PZ(JI,JKP) - PZ(JI,JK) ) *  &
     ( PZLCL(JI) - PZ(JI,JK) ) 
    ZWORK1(JI) = ( _TWO_ * ZTHEUL(JI) ) / ( ZWORK2(JI) + PTHES(JI,JKP) ) - _ONE_
    PCAPE(JI)  = PCAPE(JI) + RG * ZWORK3(JI) * MAX( _ZERO_, ZWORK1(JI) )

!*       10.   Compute final values of updraft mass flux, enthalpy, r_w
!              at level k+1
!              --------------------------------------------------------

    PUMF(JI,JKP)  = PUMF(JI,JK) - PUDR(JI,JKP) + PUER(JI,JKP)
    PUMF(JI,JKP)  = MAX( PUMF(JI,JKP), 0.1 )
    PUTHL(JI,JKP) = ( PUMF(JI,JK)  * PUTHL(JI,JK) +                             &
                     PUER(JI,JKP) * PTHL(JI,JK) - PUDR(JI,JKP) * PUTHL(JI,JK) ) &
                    / PUMF(JI,JKP)
    PURW(JI,JKP)  = ( PUMF(JI,JK)  * PURW(JI,JK) +                              &
                     PUER(JI,JKP) * PRW(JI,JK) - PUDR(JI,JKP) * PURW(JI,JK) )   &
                   / PUMF(JI,JKP)

    ZE1(JI) = ZE2(JI) ! update fractional entrainment/detrainment
    ZD1(JI) = ZD2(JI)

  ENDIF

  ENDDO UPDRAFT_LOOP
ENDDO ILOOP

!*       12.1    Set OTRIG to False if cloud thickness < 0.5km
!                or > 3km (deep convection) or CAPE < 1
!                ------------------------------------------------

    DO JI = 1, IIE
          JK  = KCTL(JI)
          ZWORK1(JI) = PZ(JI,JK) - PZLCL(JI)
          OTRIG(JI) = ZWORK1(JI) >= XCDEPTH  .AND. ZWORK1(JI) < XCDEPTH_D &
                    & .AND. PCAPE(JI) > _ONE_
    ENDDO
    WHERE( .NOT. OTRIG(:) )
          KCTL(:) = IKB
    END WHERE
KETL(:) = MAX( KETL(:), KLCL(:) + 2 )
KETL(:) = MIN( KETL(:), KCTL(:) )


!*       12.2    If the ETL and CTL are the same detrain updraft mass
!                flux at this level
!                -------------------------------------------------------

ZWORK1(:) = _ZERO_
WHERE ( KETL(:) == KCTL(:) ) ZWORK1(:) = _ONE_

DO JI = 1, IIE
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
ENDDO

!*       12.3    Adjust mass flux profiles, detrainment rates, and
!                precipitation fallout rates to reflect linear decrease
!                in mass flux between the ETL and CTL
!                -------------------------------------------------------

ZWORK1(:) = _ZERO_

DO JI = 1, IIE
   DO JK = KETL(JI), KCTL(JI)
     IF( JK > KETL(JI) ) THEN
        ZWORK1(JI) = ZWORK1(JI) + PDPRES(JI,JK)
     ENDIF
   ENDDO
   ZWORK1(JI) = PUMF(JI,KETL(JI)) / MAX( _ONE_, ZWORK1(JI) )

   DO JK = KETL(JI) + 1, KCTL(JI)
      JKP = JK - 1
      PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
      PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
   ENDDO
ENDDO

!         12.4   Set mass flux and entrainment in the source layer.
!                Linear increase throughout the source layer.
!                -------------------------------------------------------

! FIXME - the code above seems like a good idea.  If our LCL is very low (within the departure
!         layer (KLCL < KPBL)) then we're making this post-hoc adjustment to within-cloud values,
!         where they really only seem justified in the subcloud layer.  Note that once this is
!         changed, the conditional for the PWU definition below can be removed.

DO JI = 1, IIE
!  JK  = KDPL(JI) !  JKP = KPBL(JI)
!  ZWORK2(JI) = PPRES(JI,JK) - PPRES(JI,JKP) + PDPRES(JI,JK)
!  mixed layer depth
   ZWORK2(JI) = PPRES(JI,KDPL(JI)) - PPRES(JI,KPBL(JI)) + PDPRES(JI,KDPL(JI))
   DO JK = KDPL(JI), KPBL(JI)
     PUER(JI,JK) = PUER(JI,JK) + PMFLCL(JI) * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
     PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
     IF (JK < KLCL(JI)) PWU(JI,JK) = PWLCL(JI) * PUMF(JI,JK) / PMFLCL(JI)
   ENDDO
ENDDO


!*       13.   If cloud thickness is smaller than  .5 km or > 3 km
!              no shallow convection is allowed
!              Nota: For technical reasons, we stop the convection
!                    computations in this case and do not go back to
!                    TRIGGER_FUNCT to look for the next unstable LCL
!                    which could produce a thicker cloud.
!              ---------------------------------------------------

GWORK6(:,:) = SPREAD( OTRIG(:), DIM=2, NCOPIES=KLEV )
WHERE ( .NOT. GWORK6(:,:) )
    PUMF(:,:)  = _ZERO_
    PUDR(:,:)  = _ZERO_
    PUER(:,:)  = _ZERO_
    PWU(:,:)   = _ZERO_
    PUTHL(:,:) = PTHL(:,:)
    PURW(:,:)  = PRW(:,:)
    PURC(:,:)  = _ZERO_
    PURI(:,:)  = _ZERO_
END WHERE

END SUBROUTINE CONVECT_UPDRAFT_SHAL3

