!####################################################################
 SUBROUTINE CONVECT_CLOSURE2( KLON, KLEV,                            &
                           & PPRES, PDPRES, PZ, PDXDY, PLMASS,      &
                           & PTHL, PTH, PRW, PRC, PRI, OTRIG1,      &
                           & PTHC, PRWC, PRCC, PRIC, PWSUB,         &
                           & KLCL, KDPL, KPBL, KLFS, KCTL, KML,     &
                           & PUMF, PUER, PUDR, PUTHL, PURW,         &
                           & PURC, PURI, PUPR,                      &
                           & PDMF, PDER, PDDR, PDTHL, PDRW,         &
                           & PTPR, PSPR, PDTEVR,                    &
                           & PCAPE, PTIMEC,                         &
                           & KFTSTEPS,                              &
                           & PDTEVRF, PPRLFLX, PPRSFLX,ZTIMEC,ITSTEP,GWORK5 )
!####################################################################

!!**** Uses modified Fritsch-Chappell closure
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!     (over a time step PTIMEC) environmental values of THETA_l, R_w, R_c, R_i
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PTHC-PTH)/PTIMEC
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
!!
!!    CONVECT_CLOSURE_THRVLCL
!!    CONVECT_CLOSURE_ADJUST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RATM               ! reference pressure
!!          RD, RV           ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV         ! specific heat for dry air and water vapor
!!          RCW, RCS           ! specific heat for liquid water and ice
!!          RTT                ! triple point temperature
!!          RLVTT, RLSTT       ! vaporization, sublimation heat constant
!!
!!      Module YOE_CONVPAR
!!          XA25               ! reference grid area
!!          XSTABT             ! stability factor in time integration
!!          XSTABC             ! stability factor in CAPE adjustment
!!          XMELDPTH           ! allow melting over specific pressure depth
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Peter Bechtold 04/10/97 change for enthalpie, r_c + r_i tendencies
!-------------------------------------------------------------------------------

!*       0.    DECLARATIONS
!              ------------

USE YOMCST
USE YOE_CONVPAR
USE YOE_CONVPAREXT

IMPLICIT NONE
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

!*       0.1   Declarations of dummy arguments :

integer,                   INTENT(IN) :: KLON   ! horizontal dimension
integer,                   INTENT(IN) :: KLEV   ! vertical dimension
integer, DIMENSION(KLON),  INTENT(IN) :: KLFS   ! index for level of free sink
integer, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! index lifting condens. level
integer, DIMENSION(KLON),  INTENT(IN) :: KCTL   ! index for cloud top level
integer, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! index for departure level
integer, DIMENSION(KLON),  INTENT(IN) :: KPBL   ! index for top of source layer
integer, DIMENSION(KLON),  INTENT(IN) :: KML    ! index for melting level
real, DIMENSION(KLON),  INTENT(INOUT) :: PTIMEC ! convection time step
real, DIMENSION(KLON),     INTENT(IN) :: PDXDY  ! grid area (m^2)
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHL   ! grid scale enthalpy (J/kg)
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH    ! grid scale theta
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PRW    ! grid scale total water
                                                  ! mixing ratio
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PRC    ! grid scale r_c
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PRI    ! grid scale r_i
LOGICAL, DIMENSION(KLON),  INTENT(IN) :: OTRIG1   ! logical to keep trace of
                                                  ! convective arrays modified in UPDRAFT


real, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (P)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES ! pressure difference between
                                                   ! bottom and top of layer (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PLMASS ! mass of model layer (kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m)
real, DIMENSION(KLON),     INTENT(IN)  :: PCAPE  ! available potent. energy
integer,                INTENT(OUT)   :: KFTSTEPS! maximum of fract time steps
                                                   ! only used for chemical tracers


real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF  ! updraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUER  ! updraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUDR  ! updraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUPR  ! updraft precipitation in
                                                    ! flux units (kg water / s)
real, DIMENSION(KLON,KLEV), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)

real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF  ! downdraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDER  ! downdraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDDR  ! downdraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PDTHL ! downdraft enthalpy (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)   :: PDRW  ! downdraft total water (kg/kg)
real, DIMENSION(KLON),      INTENT(INOUT):: PTPR  ! total surf precipitation (kg/s)
real, DIMENSION(KLON),      INTENT(OUT)  :: PSPR  ! solid surf precipitation (kg/s)
real, DIMENSION(KLON),      INTENT(INOUT):: PDTEVR! donwndraft evapor. (kg/s)

real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)

real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDTEVRF! downdraft evaporation rate
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRLFLX! liquid precip flux
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRSFLX! solid  precip flux

real, DIMENSION(KLON),INTENT(OUT)        :: ZTIMEC       !fractional convective time step
integer, DIMENSION(KLON),INTENT(OUT)     :: ITSTEP       !# of fractional convective timesteps 
LOGICAL, DIMENSION(KLON),INTENT(OUT)       :: GWORK5       !flag for newly activated columns including MASS CONSERVATION criteria from closure
                                                           !this flag is used in convect_uv_transport.F90 

!*       0.2   Declarations of local variables :

integer :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
integer :: IKS            ! vertical dimension
integer :: JK, JKP, JKMAX ! vertical loop index
integer :: JI             ! horizontal loop index
integer :: JITER          ! iteration loop index
integer :: JSTEP          ! fractional time loop index
 real    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!real    :: ZCVOCD, ZEPSA  ! C_pv / C_pd, R_v / R_d

real, DIMENSION(KLON,KLEV) :: ZTHLC       ! convectively adjusted
                                            ! grid scale enthalpy
real, DIMENSION(KLON,KLEV) :: ZOMG        ! conv. environm. subsidence (Pa/s)
real, DIMENSION(KLON,KLEV) :: ZUMF        ! non-adjusted updraft mass flux
real, DIMENSION(KLON,KLEV) :: ZUER        !   "     updraft  entrainm. rate
real, DIMENSION(KLON,KLEV) :: ZUDR        !   "     updraft  detrainm. rate
real, DIMENSION(KLON,KLEV) :: ZDMF        !   "   downdraft mass flux
real, DIMENSION(KLON,KLEV) :: ZDER        !   "   downdraft  entrainm. rate
real, DIMENSION(KLON,KLEV) :: ZDDR        !   "   downdraft  detrainm. rate
real, DIMENSION(KLON)     :: ZTPR         !   "   total precipitation
real, DIMENSION(KLON)     :: ZDTEVR       !   "   total downdraft evapor.
real, DIMENSION(KLON,KLEV):: ZPRLFLX      !   "   liquid precip flux
real, DIMENSION(KLON,KLEV):: ZPRSFLX      !   "   solid  precip flux
real, DIMENSION(KLON)     :: ZPRMELT      ! melting of precipitation
real, DIMENSION(KLON)     :: ZPRMELTO     ! non-adjusted  "
real, DIMENSION(KLON)     :: ZADJ         ! mass adjustment factor
real, DIMENSION(KLON)     :: ZADJMAX      ! limit value for ZADJ
real, DIMENSION(KLON)     :: ZCAPE        ! new CAPE after adjustment
real, DIMENSION(KLON,KLEV):: ZTIMC        ! 2D work array for ZTIMEC

real, DIMENSION(KLON)     :: ZTHLCL       ! new  theta at LCL
real, DIMENSION(KLON)     :: ZRVLCL       ! new  r_v at LCL
real, DIMENSION(KLON)     :: ZZLCL        ! height of LCL
real, DIMENSION(KLON)     :: ZTLCL        ! temperature at LCL
real, DIMENSION(KLON)     :: ZTELCL       ! envir. temper. at LCL
real, DIMENSION(KLON)     :: ZTHEUL       ! theta_e for undilute ascent
real, DIMENSION(KLON)     :: ZTHES1, ZTHES2! saturation environm. theta_e
real, DIMENSION(KLON,KLEV) :: ZTHMFIN, ZTHMFOUT, ZRWMFIN, ZRWMFOUT
real, DIMENSION(KLON,KLEV) :: ZRCMFIN, ZRCMFOUT, ZRIMFIN, ZRIMFOUT
                                    ! work arrays for environm. compensat. mass flux
real, DIMENSION(KLON)     :: ZPI          ! (P/P00)**R_d/C_pd
real, DIMENSION(KLON)     :: ZLV          ! latent heat of vaporisation
real, DIMENSION(KLON)     :: ZLS          ! latent heat of sublimation
real, DIMENSION(KLON)     :: ZLM          ! latent heat of melting
real, DIMENSION(KLON)     :: ZCPH         ! specific heat C_ph
real, DIMENSION(KLON)     :: ZMELDPTH     ! actual depth of melting layer
integer, DIMENSION(KLON)  :: ICOUNT       ! timestep counter
integer, DIMENSION(KLON)  :: ILCL         ! index lifting condens. level
integer, DIMENSION(KLON)  :: IWORK1       ! work array
real,  DIMENSION(KLON)      :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
real,  DIMENSION(KLON)      :: ZWORK4, ZWORK5          ! work arrays
real,  DIMENSION(KLON,KLEV) :: ZWORK6                  ! work array
LOGICAL, DIMENSION(KLON)      :: GWORK1, GWORK3          ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK4                  ! work array


!-------------------------------------------------------------------------------

!*       0.2    Initialize  local variables
!               ----------------------------


PSPR(:)     = _ZERO_

ZTIMC(:,:)  = _ZERO_
ZTHES2(:)   = _ZERO_

ZWORK1(:)   = _ZERO_
ZWORK2(:)   = _ZERO_
ZWORK3(:)   = _ZERO_
ZWORK4(:)   = _ZERO_
ZWORK5(:)   = _ZERO_
ZWORK6(:,:) = _ZERO_

GWORK1(:)   = .FALSE.
GWORK3(:)   = .FALSE.
GWORK4(:,:) = .FALSE.

ILCL(:)     = KLCL(:)

 ZCPORD    = RCPD / RD
 ZRDOCP    = RD  / RCPD
!ZCVOCD    = RCPV / RCPD
!ZEPSA     = RV  / RD

ZADJ(:)   = _ONE_
ZWORK5(:) = _ONE_
WHERE( .NOT. OTRIG1(:) ) ZWORK5(:) = _ZERO_


!*       0.3   Compute loop bounds
!              -------------------

IIE    = KLON
IKB    = 1 + JCVEXB
IKS    = KLEV
IKE    = KLEV - JCVEXT
JKMAX  = MAXVAL( KCTL(:) )


!*       2.     Save initial mass flux values to be used in adjustment procedure
!               ---------------------------------------------------------------

ZUMF(:,:)  = PUMF(:,:)
ZUER(:,:)  = PUER(:,:)
ZUDR(:,:)  = PUDR(:,:)
ZDMF(:,:)  = PDMF(:,:)
ZDER(:,:)  = PDER(:,:)
ZDDR(:,:)  = PDDR(:,:)
ZTPR(:)    = PTPR(:)
ZDTEVR(:)  = PDTEVR(:)
ZOMG(:,:)  = _ZERO_
PWSUB(:,:) = _ZERO_
ZPRMELT(:) = _ZERO_
PPRLFLX(:,:) = _ZERO_
ZPRLFLX(:,:) = _ZERO_
PPRSFLX(:,:) = _ZERO_
ZPRSFLX(:,:) = _ZERO_


!*       2.1    Some preliminar computations for melting of precipitation
!               used later in section 9 and computation of precip fluxes
!               Precipitation fluxes are updated for melting and evaporation
!               ---------------------------------------------------------


ZWORK1(:) = _ZERO_
ZMELDPTH(:) = _ZERO_
ZWORK6(:,:) = _ZERO_

DO JK = JKMAX + 1, IKB + 1, -1
   ! Nota: PUPR is total precipitation flux, but the solid, liquid
   !       precipitation is stored in units kg/kg; therefore we compute here
   !       the solid fraction of the total precipitation flux.
  DO JI = 1, IIE
     ZWORK2(JI)    = PUPR(JI,JK) / ( PURC(JI,JK) + PURI(JI,JK) + 1.E-8 )
     ZPRMELT(JI)   = ZPRMELT(JI) + PURI(JI,JK) * ZWORK2(JI)
     ZWORK1(JI)    = ZWORK1(JI) + PURC(JI,JK) * ZWORK2(JI) - PDTEVRF(JI,JK)
     ZPRLFLX(JI,JK)= MAX( _ZERO_, ZWORK1(JI) )
     ZPRMELT(JI)   = ZPRMELT(JI) + MIN( _ZERO_, ZWORK1(JI) )
     ZPRSFLX(JI,JK)= ZPRMELT(JI)
     IF ( KML(JI) >= JK .AND. ZMELDPTH(JI) <= XMELDPTH ) THEN
          ZPI(JI)    = ( PPRES(JI,JK) / RATM ) ** ZRDOCP
          ZWORK3(JI) = PTH(JI,JK) * ZPI(JI)            ! temperature estimate
          ZLM(JI)    = RLSTT + ( RCPV - RCS ) * ( ZWORK3(JI) - RTT ) -         &
                   & ( RLVTT + ( RCPV - RCW ) * ( ZWORK3(JI) - RTT ) ) ! L_s - L_v
          ZCPH(JI)   = RCPD + RCPV * PRW(JI,JK)
          ZMELDPTH(JI) = ZMELDPTH(JI) + PDPRES(JI,JK)
          ZWORK6(JI,JK)= ZLM(JI) * PTIMEC(JI) / PLMASS(JI,JK) * PDPRES(JI,JK)
          ZOMG(JI,JK)= _ONE_ ! at this place only used as work variable
     ENDIF
  ENDDO

ENDDO

ZWORK2(:) = _ZERO_
DO JK = JKMAX, IKB + 1, -1
    ZWORK1(:) = ZPRMELT(:) * PDPRES(:,JK) / MAX( XMELDPTH, ZMELDPTH(:) )
    ZWORK2(:) = ZWORK2(:) + ZWORK1(:) * ZOMG(:,JK)
    ZPRLFLX(:,JK) = ZPRLFLX(:,JK) + ZWORK2(:)
    ZPRSFLX(:,JK) = ZPRSFLX(:,JK) - ZWORK2(:)
ENDDO
WHERE( ZPRSFLX(:,:) < _ONE_ ) ZPRSFLX(:,:)=_ZERO_
ZPRMELTO(:) = ZPRMELT(:)


!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------

ZADJMAX(:) = 1000.
IWORK1(:)  = MAX( ILCL(:), KLFS(:) )
JKP        = MINVAL( KDPL(:) )

DO JK = JKP, IKE
  DO JI = 1, IIE
    IF( JK > KDPL(JI) .AND. JK <= IWORK1(JI) ) THEN
        ZWORK1(JI)  = PLMASS(JI,JK) /                                      &
                & ( ( PUER(JI,JK) + PDER(JI,JK) + 1.E-5 ) * PTIMEC(JI) )
        ZADJMAX(JI) = MIN( ZADJMAX(JI), ZWORK1(JI) )
    ENDIF
  ENDDO
ENDDO


GWORK5(:) = OTRIG1(:)  !  
GWORK1(:) = OTRIG1(:)  ! logical array to limit adjustment to not definitively
                       ! adjusted columns

DO JK = IKB, IKE
  ZTHLC(:,JK) = PTHL(:,JK) ! initialize adjusted envir. values
  PRWC(:,JK)  = PRW(:,JK)
  PRCC(:,JK)  = PRC(:,JK)
  PRIC(:,JK)  = PRI(:,JK)
  PTHC(:,JK)  = PTH(:,JK)
ENDDO



DO JITER = 1, 6  ! Enter adjustment loop to assure that all CAPE is
                 ! removed within the advective time interval TIMEC

     ZTIMEC(:) = PTIMEC(:)
     GWORK4(:,:)   = SPREAD( GWORK1(:), DIM=2, NCOPIES=IKS )
     WHERE( GWORK4(:,:) ) PWSUB(:,:) = _ZERO_
     ZOMG(:,:)=_ZERO_

     DO JK = IKB + 1, JKMAX
           JKP = MAX( IKB + 1, JK - 1 )
           WHERE ( GWORK1(:) .AND. JK <= KCTL(:) )


!*       4.     Determine vertical velocity at top and bottom of each layer
!               to satisfy mass continuity.
!               ---------------------------------------------------------------
              ! we compute here Domega/Dp = - g rho Dw/Dz = 1/Dt

             ZWORK1(:)   = - ( PUER(:,JKP) + PDER(:,JKP) -                   &
                         & PUDR(:,JKP) - PDDR(:,JKP) ) / PLMASS(:,JKP)

             PWSUB(:,JK) = PWSUB(:,JKP) - PDPRES(:,JK-1) * ZWORK1(:)
              ! we use PDPRES(JK-1) and not JKP in order to have zero subsidence
              ! at the first layer


!*       5.     Compute fractional time step. For stability or
!               mass conservation reasons one must split full time step PTIMEC)
!               ---------------------------------------------------------------

             ZWORK1(:) = XSTABT * PDPRES(:,JKP) / ( ABS( PWSUB(:,JK) ) + 1.E-10 )
              ! the factor XSTABT is used for stability reasons
             ZTIMEC(:) = MIN( ZTIMEC(:), ZWORK1(:) )

              ! transform vertical velocity in mass flux units
             ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / RG
         END WHERE

     ENDDO


     WHERE( GWORK4(:,:) )
           ZTHLC(:,:) = PTHL(:,:) ! reinitialize adjusted envir. values
           PRWC(:,:)  = PRW(:,:)  ! when iteration criterium not attained
           PRCC(:,:)  = PRC(:,:)
           PRIC(:,:)  = PRI(:,:)
           PTHC(:,:)  = PTH(:,:)
     END WHERE


!        6. Check for mass conservation, i.e. ZWORK1 > 1.E-2
!           If mass is not conserved, the convective tendencies
!           automatically become zero.
!           ----------------------------------------------------

    DO JI = 1, IIE
       JK=KCTL(JI)
       ZWORK1(JI) = PUDR(JI,JK) * PDPRES(JI,JK) / ( PLMASS(JI,JK) + .1 )    &
                  &                                         - PWSUB(JI,JK)
    ENDDO
    WHERE( GWORK1(:) .AND. ABS( ZWORK1(:) ) - .01 > _ZERO_ )
        GWORK1(:) = .FALSE.
        GWORK5(:) = .FALSE.    !output this flag for late use in convect_uv_transport.F90
        PTIMEC(:) = 1.E-1
        ZTPR(:)   = _ZERO_
        ZWORK5(:) = _ZERO_
    END WHERE

    DO JK = IKB, IKE
        PWSUB(:,JK) = PWSUB(:,JK) * ZWORK5(:)
        ZPRLFLX(:,JK) = ZPRLFLX(:,JK) * ZWORK5(:)
        ZPRSFLX(:,JK) = ZPRSFLX(:,JK) * ZWORK5(:)
    ENDDO
    GWORK4(:,1:IKB) = .FALSE.
    GWORK4(:,IKE:IKS) = .FALSE.

    ITSTEP(:) = INT( PTIMEC(:) / ZTIMEC(:) ) + 1
    ZTIMEC(:) = PTIMEC(:) / REAL( ITSTEP(:) ) ! adjust  fractional time step
                                              ! to be an integer multiple of PTIMEC
    ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
    ICOUNT(:) = 0



    KFTSTEPS = MAXVAL( ITSTEP(:) )
    DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop here

         ICOUNT(:) = ICOUNT(:) + 1

             GWORK3(:) =  ITSTEP(:) >= ICOUNT(:) .AND. GWORK1(:)


!*       7.     Assign enthalpy and r_w values at the top and bottom of each
!               layer based on the sign of w
!               ------------------------------------------------------------

             ZTHMFIN(:,:)   = _ZERO_
             ZRWMFIN(:,:)   = _ZERO_
             ZRCMFIN(:,:)   = _ZERO_
             ZRIMFIN(:,:)   = _ZERO_
             ZTHMFOUT(:,:)  = _ZERO_
             ZRWMFOUT(:,:)  = _ZERO_
             ZRCMFOUT(:,:)  = _ZERO_
             ZRIMFOUT(:,:)  = _ZERO_

         DO JK = IKB + 1, JKMAX
           DO JI = 1, IIE
              GWORK4(JI,JK) = GWORK3(JI) .AND. JK <= KCTL(JI)
           ENDDO
           JKP = MAX( IKB + 1, JK - 1 )
           DO JI = 1, IIE
           IF ( GWORK3(JI) ) THEN
               ZWORK1(JI)       = SIGN( _ONE_, ZOMG(JI,JK) )
               ZWORK2(JI)       = _HALF_ * ( _ONE_ + ZWORK1(JI) )
               ZWORK1(JI)       = _HALF_ * ( _ONE_ - ZWORK1(JI) )
               ZTHMFIN(JI,JK)   = - ZOMG(JI,JK) * ZTHLC(JI,JKP) * ZWORK1(JI)
               ZTHMFOUT(JI,JK)  =   ZOMG(JI,JK) * ZTHLC(JI,JK)  * ZWORK2(JI)
               ZRWMFIN(JI,JK)   = - ZOMG(JI,JK) * PRWC(JI,JKP) * ZWORK1(JI)
               ZRWMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRWC(JI,JK)  * ZWORK2(JI)
               ZRCMFIN(JI,JK)   = - ZOMG(JI,JK) * PRCC(JI,JKP) * ZWORK1(JI)
               ZRCMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRCC(JI,JK)  * ZWORK2(JI)
               ZRIMFIN(JI,JK)   = - ZOMG(JI,JK) * PRIC(JI,JKP) * ZWORK1(JI)
               ZRIMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRIC(JI,JK)  * ZWORK2(JI)
           ENDIF
           ENDDO
           DO JI = 1, IIE
           IF ( GWORK3(JI) ) THEN
               ZTHMFIN(JI,JKP)  = ZTHMFIN(JI,JKP)  + ZTHMFOUT(JI,JK) * ZWORK2(JI)
               ZTHMFOUT(JI,JKP) = ZTHMFOUT(JI,JKP) + ZTHMFIN(JI,JK)  * ZWORK1(JI)
               ZRWMFIN(JI,JKP)  = ZRWMFIN(JI,JKP)  + ZRWMFOUT(JI,JK) * ZWORK2(JI)
               ZRWMFOUT(JI,JKP) = ZRWMFOUT(JI,JKP) + ZRWMFIN(JI,JK)  * ZWORK1(JI)
               ZRCMFIN(JI,JKP)  = ZRCMFIN(JI,JKP)  + ZRCMFOUT(JI,JK) * ZWORK2(JI)
               ZRCMFOUT(JI,JKP) = ZRCMFOUT(JI,JKP) + ZRCMFIN(JI,JK)  * ZWORK1(JI)
               ZRIMFIN(JI,JKP)  = ZRIMFIN(JI,JKP)  + ZRIMFOUT(JI,JK) * ZWORK2(JI)
               ZRIMFOUT(JI,JKP) = ZRIMFOUT(JI,JKP) + ZRIMFIN(JI,JK)  * ZWORK1(JI)
           ENDIF
           ENDDO
         ENDDO

         WHERE ( GWORK4(:,:) )

!******************************************************************************

!*       8.     Update the environmental values of enthalpy and r_w at each level
!               NOTA: These are the MAIN EQUATIONS of the scheme
!               -----------------------------------------------------------------


           ZTHLC(:,:) = ZTHLC(:,:) + ZTIMC(:,:) / PLMASS(:,:) * (      &
                      &   ZTHMFIN(:,:) + PUDR(:,:) * PUTHL(:,:)  +     &
                      &   PDDR(:,:) * PDTHL(:,:) - ZTHMFOUT(:,:) -     &
                      & ( PUER(:,:) + PDER(:,:) ) * PTHL(:,:)   )

           PRWC(:,:)  = PRWC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
                      &   ZRWMFIN(:,:) + PUDR(:,:) * PURW(:,:)  +      &
                      &   PDDR(:,:) * PDRW(:,:) - ZRWMFOUT(:,:) -      &
                      & ( PUER(:,:) + PDER(:,:) ) * PRW(:,:)    )

           PRCC(:,:)  = PRCC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
                      &   ZRCMFIN(:,:) + PUDR(:,:) * PURC(:,:)         &
                      &                         - ZRCMFOUT(:,:) -      &
                      & ( PUER(:,:) + PDER(:,:) ) * PRC(:,:)    )

           PRIC(:,:)  = PRIC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
                      &   ZRIMFIN(:,:) + PUDR(:,:) * PURI(:,:)         &
                      &                         - ZRIMFOUT(:,:) -      &
                      & ( PUER(:,:) + PDER(:,:) ) * PRI(:,:)    )


!******************************************************************************

         END WHERE

    ENDDO ! Exit the fractional time step loop


!*           9.    Allow frozen precipitation to melt over a 200 mb deep layer
!                  -----------------------------------------------------------

      DO JK = JKMAX, IKB + 1, -1
            ZTHLC(:,JK) = ZTHLC(:,JK) -                                &
            &  ZPRMELT(:) * ZWORK6(:,JK) / MAX( XMELDPTH, ZMELDPTH(:) )
      ENDDO


!*          10.    Compute final linearized value of theta envir.
!                  ----------------------------------------------

      DO JK = IKB + 1, JKMAX
         DO JI = 1, IIE
         IF( GWORK1(JI) .AND. JK <= KCTL(JI) ) THEN
           ZPI(JI)    = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
           ZCPH(JI)   = RCPD + PRWC(JI,JK) * RCPV
           ZWORK2(JI) = PTH(JI,JK) / ZPI(JI)  ! first temperature estimate
           ZLV(JI)    = RLVTT + ( RCPV - RCW ) * ( ZWORK2(JI) - RTT )
           ZLS(JI)    = RLVTT + ( RCPV - RCS ) * ( ZWORK2(JI) - RTT )
             ! final linearized temperature
           ZWORK2(JI) = ( ZTHLC(JI,JK) + ZLV(JI) * PRCC(JI,JK) + ZLS(JI) * PRIC(JI,JK) &
                      & - (_ONE_ + PRWC(JI,JK) ) * RG * PZ(JI,JK) ) / ZCPH(JI)
           ZWORK2(JI) = MAX( 180., MIN( 340., ZWORK2(JI) ) )
           PTHC(JI,JK)= ZWORK2(JI) * ZPI(JI) ! final adjusted envir. theta
         ENDIF
         ENDDO
      ENDDO


!*         11.     Compute new cloud ( properties at new LCL )
!                     NOTA: The computations are very close to
!                           that in routine TRIGGER_FUNCT
!                  ---------------------------------------------

      CALL CONVECT_CLOSURE_THRVLCL(  KLON, KLEV,                           &
                                  &  PPRES, PTHC, PRWC, PZ, GWORK1,        &
                                  &  ZTHLCL, ZRVLCL, ZZLCL, ZTLCL, ZTELCL, &
                                  &  ILCL, KDPL, KPBL )


       ZTLCL(:)  = MAX( 230., MIN( 335., ZTLCL(:)  ) )  ! set some overflow bounds
       ZTELCL(:) = MAX( 230., MIN( 335., ZTELCL(:) ) )
       ZTHLCL(:) = MAX( 230., MIN( 345., ZTHLCL(:) ) )
       ZRVLCL(:) = MAX(   _ZERO_,  MIN(   _ONE_,   ZRVLCL(:) ) )


!*         12.    Compute adjusted CAPE
!                 ---------------------

       ZCAPE(:)  = _ZERO_
       ZPI(:)    = ZTHLCL(:) / ZTLCL(:)
       ZPI(:)    = MAX( 0.95, MIN( 1.5, ZPI(:) ) )
       ZWORK1(:) = RATM / ZPI(:) ** ZCPORD ! pressure at LCL

       CALL CONVECT_SATMIXRATIO( KLON, ZWORK1, ZTELCL, ZWORK3, ZLV, ZLS, ZCPH )
       ZWORK3(:) = MIN(   .1, MAX(   _ZERO_, ZWORK3(:) ) )

                ! compute theta_e updraft undilute
       ZTHEUL(:) = ZTLCL(:) * ZPI(:) ** ( _ONE_ - 0.28 * ZRVLCL(:) )             &
                 &                * EXP( ( 3374.6525 / ZTLCL(:) - 2.5403 )  &
                 &                * ZRVLCL(:) * ( _ONE_ + 0.81 * ZRVLCL(:) ) )

                ! compute theta_e saturated environment at LCL
       ZTHES1(:) = ZTELCL(:) * ZPI(:) ** ( _ONE_ - 0.28 * ZWORK3(:) )            &
                 &                * EXP( ( 3374.6525 / ZTELCL(:) - 2.5403 ) &
                 &                * ZWORK3(:) * ( _ONE_ + 0.81 * ZWORK3(:) ) )

      DO JK = MINVAL( ILCL(:) ), JKMAX
        JKP = JK - 1
        DO JI = 1, IIE
          ZWORK4(JI) = _ONE_
          IF ( JK == ILCL(JI) ) ZWORK4(JI) = _ZERO_

           ! compute theta_e saturated environment and adjusted values
           ! of theta

          GWORK3(JI)  = JK >= ILCL(JI) .AND. JK <= KCTL(JI) .AND. GWORK1(JI)

          ZPI(JI)     = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
          ZWORK2(JI)  = PTHC(JI,JK) / ZPI(JI)
        ENDDO

        CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZWORK2, ZWORK3, ZLV, ZLS, ZCPH )


        DO JI = 1, IIE
          IF ( GWORK3(JI) ) THEN
              ZTHES2(JI)  = ZWORK2(JI) * ZPI(JI) ** ( _ONE_ - 0.28 * ZWORK3(JI) )  &
                          &        * EXP( ( 3374.6525 / ZWORK2(JI) - 2.5403 ) &
                          &        * ZWORK3(JI) * ( _ONE_ + 0.81 * ZWORK3(JI) ) )

              ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) * ZWORK4(JI) -                       &
                         & ( _ONE_ - ZWORK4(JI) ) * ZZLCL(JI)    ! level thickness
              ZWORK1(JI)  = ( _TWO_ * ZTHEUL(JI) ) / ( ZTHES1(JI) + ZTHES2(JI) ) - _ONE_
              ZCAPE(JI)   = ZCAPE(JI) + RG * ZWORK3(JI) * MAX( _ZERO_, ZWORK1(JI) )
              ZTHES1(JI)  = ZTHES2(JI)
          ENDIF
        ENDDO
      ENDDO


!*         13.     Determine mass adjustment factor knowing how much
!                  CAPE has been removed.
!                  -------------------------------------------------

       WHERE ( GWORK1(:) )
           ZWORK1(:) = MAX( PCAPE(:) - ZCAPE(:), 0.1 * PCAPE(:) )
           ZWORK2(:) = ZCAPE(:) / ( PCAPE(:) + 1.E-8 )

           GWORK1(:) = ZWORK2(:) > 0.1 .OR. ZCAPE(:) == _ZERO_ ! mask for adjustment
       END WHERE

       WHERE ( ZCAPE(:) == _ZERO_ .AND. GWORK1(:) )  ZADJ(:) = ZADJ(:) * _HALF_
       WHERE ( ZCAPE(:) /= _ZERO_ .AND. GWORK1(:) )                                    &
             & ZADJ(:) = ZADJ(:) * XSTABC * PCAPE(:) / ( ZWORK1(:) + 1.E-8 )
       ZADJ(:) = MIN( ZADJ(:), ZADJMAX(:) )


!*         14.     Adjust mass flux by the factor ZADJ to converge to
!                  specified degree of stabilization
!                 ----------------------------------------------------

       CALL CONVECT_CLOSURE_ADJUST( KLON, KLEV, ZADJ,                     &
                                  & PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR,   &
                                  & PDMF, ZDMF, PDER, ZDER, PDDR, ZDDR,   &
                                  & ZPRMELT, ZPRMELTO, PDTEVR, ZDTEVR,    &
                                  & PTPR, ZTPR,                           &
                                  & PPRLFLX, ZPRLFLX, PPRSFLX, ZPRSFLX    )

      IF ( COUNT( GWORK1(:) ) == 0 ) EXIT ! exit big adjustment iteration loop
                                          ! when all columns have reached
                                          ! desired degree of stabilization.

ENDDO  ! end of big adjustment iteration loop

!*         15.     Update compensating environmental motion based on revised
!                  entrainment/detrainment rates so that a consistent value
!                  is available for subsequent transport calculations.
!                  ----------------------------------------------------

     PWSUB(:,:) = _ZERO_
     DO JK = IKB + 1, JKMAX
           JKP = MAX( IKB + 1, JK - 1 )
           WHERE ( GWORK5(:) .AND. JK <= KCTL(:) )

             ZWORK1(:)   = - ( PUER(:,JKP) + PDER(:,JKP) -                   &
                         & PUDR(:,JKP) - PDDR(:,JKP) ) / PLMASS(:,JKP)

             PWSUB(:,JK) = PWSUB(:,JKP) - PDPRES(:,JK-1) * ZWORK1(:)
           END WHERE
     ENDDO

        ! skip adj. total water array  to water vapor
DO JK = IKB, IKE
  PRWC(:,JK) = MAX( _ZERO_, PRWC(:,JK) - PRCC(:,JK) - PRIC(:,JK) )
ENDDO

        ! compute surface solid (ice) precipitation
PSPR(:) = ZPRMELT(:) * ( _ONE_ - ZMELDPTH(:) / XMELDPTH )
PSPR(:) = MAX( _ZERO_ , PSPR(:) )


END SUBROUTINE CONVECT_CLOSURE2

