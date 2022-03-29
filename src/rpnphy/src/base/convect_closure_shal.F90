!#######################################################################
 SUBROUTINE CONVECT_CLOSURE_SHAL5( KLON, KLEV,                         &
                       &   PPRES, PDPRES, PZ, PDXDY, PMRK2, PCRAD, PLMASS,&
                       &   PTHL, PTH, PRW, PRC, PRI, PDMSEDT, OTRIG1,  &
                       &   PTHC, PRWC, PRCC, PRIC, PWSUB,              &
                       &   KLCL, KDPL, KPBL, KCTL,                     &
                       &   PPLCL, PMSELCL, PMSEELCL,                   &
                       &   PUMF, PUER, PUDR, PUTHL, PURW,              &
                       &   PURC, PURI, PCAPE, PTIMEC, PDTCONV, KFTSTEPS,&
                       &   ZTIMEC, ITSTEP, GWORK5 )
!#######################################################################
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
!!    CONVECT_CLOSURE_ADJUST_SHAL
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
!!      Module YOE_CONVPAR_SHAL
!!          XA25               ! reference grid area
!!          XSTABT             ! stability factor in time integration
!!          XSTABC             ! stability factor in CAPE adjustment
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
!!   Peter Bechtold 15/11/96 change for enthalpie, r_c + r_i tendencies
!!      Tony Dore   14/10/96 Initialise local variables
!-------------------------------------------------------------------------------


!*       0.    DECLARATIONS
!              ------------

USE YOMCST
USE YOE_CONVPAR_SHAL
USE YOE_CONVPAREXT
use integrals, only: int_profile,INT_OK
use cnv_options
use ens_perturb, only: ens_nc2d, ens_spp_get

IMPLICIT NONE
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

!*       0.1   Declarations of dummy arguments :

integer,                   INTENT(IN) :: KLON   ! horizontal dimension
integer,                   INTENT(IN) :: KLEV   ! vertical dimension
integer, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! index lifting condens. level
integer, DIMENSION(KLON),  INTENT(IN) :: KCTL   ! index for cloud top level
integer, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! index for departure level
integer, DIMENSION(KLON),  INTENT(IN) :: KPBL   ! index for top of source layer
real, DIMENSION(KLON),  INTENT(INOUT) :: PTIMEC ! convection time step
real, DIMENSION(KLON),     INTENT(IN) :: PDXDY  ! grid area (m^2)
real, DIMENSION(KLON,ens_nc2d), INTENT(IN) :: PMRK2 ! Markov chains for SPP
real, DIMENSION(KLON),     INTENT(IN) :: PCRAD  ! cloud radius (m)
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHL   ! grid scale enthalpy (J/kg)
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH    ! grid scale theta
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PRW    ! grid scale total water
                                                  ! mixing ratio
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PRC    ! grid scale r_c
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PRI    ! grid scale r_i
real, DIMENSION(KLON,KLEV),INTENT(IN) :: PDMSEDT! tendency of grid scale moist static energy (m^2/s^3)
LOGICAL, DIMENSION(KLON),  INTENT(IN)   :: OTRIG1 ! logical to keep trace of
                                                  ! convective arrays modified in UPDRAFT
real, DIMENSION(KLON),     INTENT(IN) :: PPLCL  ! pressure at the lcl (Pa)
real, DIMENSION(KLON),     INTENT(IN) :: PMSELCL! updraft MSE at lcl (m^2/s^2)
real, DIMENSION(KLON),     INTENT(IN) :: PMSEELCL!envir. MSE at the lcl (m^2/s^2)

real, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (P)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES ! pressure difference between
                                                   ! bottom and top of layer (Pa)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PLMASS ! mass of model layer (kg)
real, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m)
real, DIMENSION(KLON),     INTENT(IN)  :: PCAPE  ! available potent. energy
real,                       INTENT(IN) :: PDTCONV! interval between calls to convection (s)
integer,                INTENT(OUT)   :: KFTSTEPS! maximum of fract time steps
                                                   ! only used for chemical tracers


real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF  ! updraft mass flux (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUER  ! updraft entrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUDR  ! updraft detrainment (kg/s)
real, DIMENSION(KLON,KLEV), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
real, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)

real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i
real, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)

real, DIMENSION(KLON),INTENT(OUT)        :: ZTIMEC! adjust fractional convective time step,an integer multiple of PTIMEC 
integer, DIMENSION(KLON),INTENT(OUT)     :: ITSTEP! fractional convective time step
LOGICAL, DIMENSION(KLON),INTENT(OUT)       :: GWORK5! flag for convectiton-activated column with MASS CONSERVED !! 
                                                    ! this flag is used in convect_uv_transport.F90 

!*       0.2   Declarations of local variables :

integer :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
integer :: IKS            ! vertical dimension
integer :: JK, JKP        ! vertical loop index
integer :: JI             ! horizontal loop index
integer :: JITER,NITER    ! iteration loop index and total number
integer :: JSTEP          ! fractional time loop index

 real    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
 !real    :: ZCVOCD, ZEPSA  ! C_pv / C_pd, R_v / R_d
 real    :: detr_cond     ! fraction of detrained condensate

real, DIMENSION(KLON,KLEV) :: ZTHLC       ! convectively adjusted
                                            ! grid scale enthalpy
real, DIMENSION(KLON,KLEV) :: ZOMG        ! conv. environm. subsidence (Pa/s)
real, DIMENSION(KLON,KLEV) :: ZUMF        ! non-adjusted updraft mass flux
real, DIMENSION(KLON,KLEV) :: ZUER        !   "     updraft  entrainm. rate
real, DIMENSION(KLON,KLEV) :: ZUDR        !   "     updraft  detrainm. rate
real, DIMENSION(KLON,KLEV) :: ZURC        !   "     updraft liquid mixing ratio
real, DIMENSION(KLON,KLEV) :: ZURI        !   "     updraft ice mixing ratio
real, DIMENSION(KLON)      :: ZDMSEDT_PBL ! integrated subcloud MSE tendency (kg m/s^5)
real, DIMENSION(KLON)     :: ZADJ         ! mass adjustment factor
real, DIMENSION(KLON)     :: ZADJMAX      ! limit value for ZADJ
real, DIMENSION(KLON)     :: ZCAPE        ! new CAPE after adjustment
real, DIMENSION(KLON,KLEV):: ZTIMC        ! 2D work array for ZTIMEC

real, DIMENSION(KLON)     :: ZTHLCL       ! new  theta at LCL
real, DIMENSION(KLON)     :: ZRVLCL       ! new  r_v at LCL
real, DIMENSION(KLON)     :: ZZLCL        ! height of LCL
real, DIMENSION(KLON)     :: ZMFLCL       ! mass flux at the LCL
real, DIMENSION(KLON)     :: ZTLCL        ! temperature at LCL
real, DIMENSION(KLON)     :: ZTELCL       ! envir. temper. at LCL
!!$real, DIMENSION(KLON)     :: ZMSELCL      ! updraft MSE at the LCL (m^2/s^2)
!!$real, DIMENSION(KLON)     :: ZMSEELCL     ! envir. MSE at the LCL (m^2/s^2)
real, DIMENSION(KLON)     :: ZTHEUL       ! theta_e for undilute ascent
real, DIMENSION(KLON)     :: ZTHES1, ZTHES2! saturation environm. theta_e
real, DIMENSION(KLON,KLEV) :: ZTHMFIN, ZTHMFOUT, ZRWMFIN, ZRWMFOUT
real, DIMENSION(KLON,KLEV) :: ZRCMFIN, ZRCMFOUT, ZRIMFIN, ZRIMFOUT
                                    ! work arrays for environm. compensat. mass flux
 real, DIMENSION(KLON)     :: ZPI          ! (P/P00)**R_d/C_pd
 real, DIMENSION(KLON)     :: ZLV          ! latent heat of vaporisation
 real, DIMENSION(KLON)     :: ZLS          ! latent heat of sublimation
!real, DIMENSION(KLON)     :: ZLM          ! latent heat of melting
 real, DIMENSION(KLON)     :: ZCPH         ! specific heat C_ph
integer, DIMENSION(KLON)  :: ICOUNT       ! timestep counter
integer, DIMENSION(KLON)  :: ILCL         ! index lifting condens. level
integer, DIMENSION(KLON)  :: IWORK1       ! work array
real,  DIMENSION(KLON)      :: ZWORK1, ZWORK2, ZWORK3 ! work arrays
real,  DIMENSION(KLON)      :: ZWORK4, ZWORK5         ! work arrays
LOGICAL, DIMENSION(KLON)      :: GWORK1                 ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK4                 ! work array
real, DIMENSION(KLON)       :: evaps      ! fraction of detrained condensate

!-------------------------------------------------------------------------------

!*       0.2    Initialize  local variables
!               ----------------------------


ZMFLCL(:)   = _ZERO_
ZTIMC(:,:)  = _ZERO_
ZTHES2(:)   = _ZERO_
ZWORK1(:)   = _ZERO_
ZWORK2(:)   = _ZERO_
ZWORK3(:)   = _ZERO_
ZWORK4(:)   = _ZERO_
ZWORK5(:)   = _ZERO_
GWORK1(:)   = .FALSE.
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
!JKMAX  = MAXVAL( KCTL(:) ) !VIV non-MPI bit reprod
GWORK1(:) = OTRIG1(:)  ! logical array to limit adjustment to not definitively
                       ! adjusted columns

!*       0.4   Determine attributes for different closure types
!              ------------------------------------------------
select case (bkf_closures)
case ('CAPE')
   niter = 4
   zadjmax(:) = 10000.
case ('EQUILIBRIUM')
   niter = 1
   zadjmax(:) = 0.1 * pdxdy(:) / (RPI*PCRAD**2)  !limit updraft area to 10% of grid cell
   ! Integrated subcloud MSE tendencies from all other sources
   if (int_profile(zdmsedt_pbl,pdmsedt,ppres,pplcl,ppres(:,ikb)) /= INT_OK) then
      call physeterror('convect_closure_shal', 'vertical interpolation error')
      return
   endif
   ! Mass flux at the LCL
   where (pmselcl(:)-pmseelcl(:) > 10.) &
   zmflcl(:) = max(zdmsedt_pbl(:) / (RG*(pmselcl(:)-pmseelcl(:))) * pdxdy(:),0.)
   ! Adjustment factor for closure
   do ji=1,iie
      if (pumf(ji,ilcl(ji)) > 0.) then
         zadj(ji) = zmflcl(ji) / pumf(ji,ilcl(ji))
      else
         zadj(ji) = 0.
      endif
   enddo
   zadj(:) = min(zadj(:),zadjmax(:))
   ! Adjustment occurs over a model timestep
   ptimec(:) = pdtconv
   where (zadj(:) <= 0.) gwork1(:) = .false.
case DEFAULT
   call physeterror('convect_closure_shal', 'unknown closure type '//trim(bkf_closures))
   return
end select

!*       2.     Save initial mass flux values to be used in adjustment procedure
!               ---------------------------------------------------------------

ZUMF(:,:)  = PUMF(:,:)
ZUER(:,:)  = PUER(:,:)
ZUDR(:,:)  = PUDR(:,:)
ZOMG(:,:)  = _ZERO_
PWSUB(:,:) = _ZERO_


!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------

IWORK1(:) = ILCL(:)
!JKP = MINVAL( KDPL(:) ) !VIV non-MPI bit-reprod
DO JI = 1, IIE
DO JK = KDPL(JI), IKE
    IF( JK > KDPL(JI) .AND. JK <= IWORK1(JI) ) THEN
        ZWORK1(JI)  = PLMASS(JI,JK) / ( ( PUER(JI,JK) + 1.E-5 ) * PTIMEC(JI) )
        ZADJMAX(JI) = MIN( ZADJMAX(JI), ZWORK1(JI) )
    ENDIF
  ENDDO
ENDDO

GWORK5(:) = GWORK1(:) 

DO JK = IKB, IKE
  ZTHLC(:,JK) = PTHL(:,JK) ! initialize adjusted envir. values
  PRWC(:,JK)  = PRW(:,JK)
  PRCC(:,JK)  = PRC(:,JK)
  PRIC(:,JK)  = PRI(:,JK)
  PTHC(:,JK)  = PTH(:,JK)
ENDDO

CLOSURE_ITERATIONS: DO JITER = 1, NITER  ! Enter adjustment loop to assure that all CAPE is
                                         ! removed within the advective time interval TIMEC (for
                                         ! bkf_closures == 'cape'), otherwise no iterations (i.e.
                                         ! niter=1).

     ZTIMEC(:) = PTIMEC(:)
     GWORK4(:,:)   = SPREAD( GWORK1(:), DIM=2, NCOPIES=IKS )
     WHERE( GWORK4(:,:) ) PWSUB(:,:) = _ZERO_
     ZOMG(:,:)=_ZERO_


     ! If this is a single-iteration closure, then use the adjustment factor based
     ! on the computed cloud base mass flux directly.
     if (niter == 1) then
        CALL CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, ZADJ,                     &
                                        & PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR    )
     endif

     DO JI = 1,IIE
        DO JK = IKB + 1, KCTL(JI)
           JKP = MAX( IKB + 1, JK - 1 )
           IF( GWORK1(JI) ) THEN


!*       4.     Determine vertical velocity at top and bottom of each layer
!               to satisfy mass continuity.
!               ---------------------------------------------------------------
              ! we compute here Domega/Dp = - g rho Dw/Dz = 1/Dt

             ZWORK1(JI)   = - ( PUER(JI,JKP) - PUDR(JI,JKP) ) / PLMASS(JI,JKP)

             PWSUB(JI,JK) = PWSUB(JI,JKP) - PDPRES(JI,JK-1) * ZWORK1(JI)
              ! we use PDPRES(JK-1) and not JKP in order to have zero subsidence
              ! at the first layer


!*       5.     Compute fractional time step. For stability or
!               mass conservation reasons one must split full time step PTIMEC)
!               ---------------------------------------------------------------

             ZWORK1(JI) = XSTABT * PDPRES(JI,JKP) / ( ABS( PWSUB(Ji,JK) ) + 1.E-10 )
              ! the factor XSTABT is used for stability reasons

             ZTIMEC(JI) = MIN( ZTIMEC(JI), ZWORK1(JI) )

              ! transform vertical velocity in mass flux units
             ZOMG(JI,JK) = PWSUB(JI,JK) * PDXDY(JI) / RG
             ENDIF
        ENDDO
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
       ZWORK1(JI) = PUDR(JI,JK) * PDPRES(JI,JK) / ( PLMASS(JI,JK) + .1 ) &
                  &                                         - PWSUB(JI,JK)
    ENDDO
    WHERE( GWORK1(:) .AND. ABS( ZWORK1(:) ) - .01 > _ZERO_ )
        GWORK1(:) = .FALSE.
        GWORK5(:) = .FALSE. 
        PTIMEC(:) = 1.E-1
        ZWORK5(:) = _ZERO_
    END WHERE
    DO JK = IKB, IKE
        PWSUB(:,JK) = PWSUB(:,JK) * ZWORK5(:)
    ENDDO
    GWORK4(:,1:IKB) = .FALSE.
    GWORK4(:,IKE:IKS) = .FALSE.

    ITSTEP(:) = INT( PTIMEC(:) / ZTIMEC(:) ) + 1
    ZTIMEC(:) = PTIMEC(:) / REAL( ITSTEP(:) ) ! adjust  fractional time step
                                              ! to be an integer multiple of PTIMEC
    ICOUNT(:) = 0

    !The value KFTSTEPS is also used for CONVECT_CHEM_TRANSPORT
    !when it comes out of this routine -VIV, inversed loop to use itstep(JI)

    KFTSTEPS = MAXVAL( ITSTEP(:) ) !VIV non-MPI bit reprod?
                                   !replace KFTSTEP with ITSTEP(JI) in loop
                       !KFTSTEPS is still used for convect_chem_transport

! Immediate evaporation of detrained condensate is implemented as a shut-off
! of condensate detrainment (total water and enthalpy are detrained as usual,
! which means all detrained total water appears following evaporation as vapour)
! Initialize these vectors in case they are needed before fractional substeps

    detr_cond = 1.
    if (bkf_evaps) detr_cond = 0.
    evaps(:) = ens_spp_get('bkf_evaps', pmrk2, default=detr_cond)    
    do JK = IKB, IKE
       ZURC(:,JK) = evaps(:) * PURC(:,JK)
       ZURI(:,JK) = evaps(:) * PURI(:,JK)
    enddo

    ILOOP: DO JI = 1, IIE

        IF ( GWORK1(JI) ) THEN

          FRACTIONAL_SUBSTEPS: DO JSTEP = 1, ITSTEP(JI) !was KFTSTEPS

!*       7.     Assign enthalpy and r_w values at the top and bottom of each
!               layer based on the sign of w
!               ------------------------------------------------------------

             ZTHMFIN(JI,:)   = _ZERO_
             ZRWMFIN(JI,:)   = _ZERO_
             ZRCMFIN(JI,:)   = _ZERO_
             ZRIMFIN(JI,:)   = _ZERO_
             ZTHMFOUT(JI,:)  = _ZERO_
             ZRWMFOUT(JI,:)  = _ZERO_
             ZRCMFOUT(JI,:)  = _ZERO_
             ZRIMFOUT(JI,:)  = _ZERO_

           DO JK = IKB + 1, KCTL(JI) !Must do this loop first
              JKP = MAX( IKB + 1, JK - 1 )
               ZWORK1(JI)      = SIGN( _ONE_, ZOMG(JI,JK) )
               ZWORK2(JI)      = _HALF_ * ( _ONE_ + ZWORK1(JI) )
               ZWORK1(JI)      = _HALF_ * ( _ONE_ - ZWORK1(JI) )
               ZTHMFIN(JI,JK)  = - ZOMG(JI,JK) * ZTHLC(JI,JKP) * ZWORK1(JI)
               ZTHMFOUT(JI,JK) =   ZOMG(JI,JK) * ZTHLC(JI,JK)  * ZWORK2(JI)
               ZRWMFIN(JI,JK)  = - ZOMG(JI,JK) * PRWC(JI,JKP) * ZWORK1(JI)
               ZRWMFOUT(JI,JK) =   ZOMG(JI,JK) * PRWC(JI,JK)  * ZWORK2(JI)
               ZRCMFIN(JI,JK)  = - ZOMG(JI,JK) * PRCC(JI,JKP) * ZWORK1(JI)
               ZRCMFOUT(JI,JK) =   ZOMG(JI,JK) * PRCC(JI,JK)  * ZWORK2(JI)
               ZRIMFIN(JI,JK)  = - ZOMG(JI,JK) * PRIC(JI,JKP) * ZWORK1(JI)
               ZRIMFOUT(JI,JK) =   ZOMG(JI,JK) * PRIC(JI,JK)  * ZWORK2(JI)
               ZTHMFIN(JI,JKP) = ZTHMFIN(JI,JKP)  + ZTHMFOUT(JI,JK) * ZWORK2(JI)
               ZTHMFOUT(JI,JKP)= ZTHMFOUT(JI,JKP) + ZTHMFIN(JI,JK)  * ZWORK1(JI)
               ZRWMFIN(JI,JKP) = ZRWMFIN(JI,JKP)  + ZRWMFOUT(JI,JK) * ZWORK2(JI)
               ZRWMFOUT(JI,JKP)= ZRWMFOUT(JI,JKP) + ZRWMFIN(JI,JK)  * ZWORK1(JI)
               ZRCMFIN(JI,JKP) = ZRCMFIN(JI,JKP)  + ZRCMFOUT(JI,JK) * ZWORK2(JI)
               ZRCMFOUT(JI,JKP)= ZRCMFOUT(JI,JKP) + ZRCMFIN(JI,JK)  * ZWORK1(JI)
               ZRIMFIN(JI,JKP) = ZRIMFIN(JI,JKP)  + ZRIMFOUT(JI,JK) * ZWORK2(JI)
               ZRIMFOUT(JI,JKP)= ZRIMFOUT(JI,JKP) + ZRIMFIN(JI,JK)  * ZWORK1(JI)
           ENDDO
           DO JK = IKB + 1, KCTL(JI)

!******************************************************************************
!*       8.   Update the environmental values of enthalpy and r_w at each level
!             NOTA: These are the MAIN EQUATIONS of the scheme
!             -----------------------------------------------------------------

           ZTHLC(JI,JK) = ZTHLC(JI,JK) + ZTIMEC(JI) / PLMASS(JI,JK) * (      &
                      &    ZTHMFIN(JI,JK) + PUDR(JI,JK) * PUTHL(JI,JK)       &
                      & - ZTHMFOUT(JI,JK) - PUER(JI,JK) * PTHL(JI,JK)   )
           PRWC(JI,JK)  = PRWC(JI,JK) + ZTIMEC(JI) / PLMASS(JI,JK) *  (      &
                      &    ZRWMFIN(JI,JK) + PUDR(JI,JK) * PURW(JI,JK)        &
                      & - ZRWMFOUT(JI,JK) - PUER(JI,JK) * PRW(JI,JK)    )
           PRCC(JI,JK)  = PRCC(JI,JK) + ZTIMEC(JI) / PLMASS(JI,JK) *  (      &
                      &    ZRCMFIN(JI,JK) + PUDR(JI,JK) * ZURC(JI,JK)        &
                      & - ZRCMFOUT(JI,JK) - PUER(JI,JK) * PRC(JI,JK)    )
           PRIC(JI,JK)  = PRIC(JI,JK) + ZTIMEC(JI) / PLMASS(JI,JK) *  (      &
                      &    ZRIMFIN(JI,JK) + PUDR(JI,JK) * ZURI(JI,JK)        &
                      & - ZRIMFOUT(JI,JK) - PUER(JI,JK) * PRI(JI,JK)    )

!******************************************************************************

           ENDDO

          ENDDO FRACTIONAL_SUBSTEPS
        ENDIF
    ENDDO ILOOP


!*          10.    Compute final linearized value of theta envir.
!                  ----------------------------------------------

      DO JI = 1, IIE
         IF( GWORK1(JI) ) THEN
         DO JK = IKB + 1, KCTL(JI)
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
         ENDDO
         ENDIF
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

       CAPE_CLOSURE: if (bkf_closures == 'CAPE') then

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

          !DO JK = MINVAL( ILCL(:) ), JKMAX !VIV non-MPI bit reprod
          !DO JK = ILCL(JI),KCTL(JI)
          DO JI = 1, IIE
             DO JK = ILCL(JI),KCTL(JI)
             JKP = JK - 1
                ZWORK4(JI) = _ONE_
                IF ( JK == ILCL(JI) ) ZWORK4(JI) = _ZERO_

                ! compute theta_e saturated environment and adjusted values
                ! of theta

                IF ( GWORK1(JI) ) THEN
                   ZPI(JI)     = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
                   ZWORK2(JI)  = PTHC(JI,JK) / ZPI(JI)

                   CALL CONVECT_SATMIXRATIO( 1, PPRES(JI,JK), ZWORK2(JI), &
                                       ZWORK3(JI), ZLV(JI), ZLS(JI), ZCPH(JI) )

                   ZTHES2(JI)  = ZWORK2(JI) * ZPI(JI) ** ( _ONE_ - 0.28 * ZWORK3(JI) )   &
                        &        * EXP( ( 3374.6525 / ZWORK2(JI) - 2.5403 )  &
                        &        * ZWORK3(JI) * ( _ONE_ + 0.81 * ZWORK3(JI) ) )

                   ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) * ZWORK4(JI) -                        &
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
          WHERE ( ZCAPE(:) /= _ZERO_ .AND. GWORK1(:) )                              &
               & ZADJ(:) = ZADJ(:) * XSTABC * PCAPE(:) / ( ZWORK1(:) + 1.E-8 )
          ZADJ(:) = MIN( ZADJ(:), ZADJMAX(:) )

       endif CAPE_CLOSURE

       ! Compute adjustments for multiple iterations
       ITERATING: if (niter > 1) then

!*         14.     Adjust mass flux by the factor ZADJ to converge to
!                  specified degree of stabilization if iterations are
!                  required for this type of closure.
!                 ----------------------------------------------------

          CALL CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, ZADJ,                     &
                                          & PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR    )

!*         15.     Update compensating environmental motion based on revised
!                  entrainment/detrainment rates so that a consistent value
!                  is available for subsequent transport calculations.
!                  ----------------------------------------------------

          DO JI = 1, IIE
          IF ( GWORK5(JI)) THEN
             DO JK = IKB + 1, KCTL(JI)
             JKP = MAX( IKB + 1, JK - 1 )
               ZWORK1(JI)   = - ( PUER(JI,JKP) - PUDR(JI,JKP) ) / PLMASS(JI,JKP)
               PWSUB(JI,JK) = PWSUB(JI,JKP) - PDPRES(JI,JK-1) * ZWORK1(JI)
             ENDDO
          ENDIF
          ENDDO

       endif ITERATING


      IF ( COUNT( GWORK1(:) ) == 0 ) EXIT ! exit big adjustment iteration loop
                                          ! when all columns have reached
                                          ! desired degree of stabilization.

ENDDO CLOSURE_ITERATIONS

        ! skip adj. total water array  to water vapor
DO JK = IKB, IKE
   PRWC(:,JK) = MAX( _ZERO_, PRWC(:,JK) - PRCC(:,JK) - PRIC(:,JK) )
ENDDO

END SUBROUTINE CONVECT_CLOSURE_SHAL5

