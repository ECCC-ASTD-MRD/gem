!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!   ##########################################################################
    SUBROUTINE SNOW_COVER_1LAYER(PTSTEP, PANSMIN, PANSMAX, PTODRY,         &
                                 PRHOSMIN, PRHOSMAX, PRHOFOLD, OALL_MELT,  &
                                 PDRAIN_TIME, PWCRN, PZ0SN, PZ0HSN,        &
                                 PTSNOW, PASNOW, PRSNOW, PWSNOW, PTS_SNOW, &
                                 PESNOW,                                   &
                                 PTG,PABS_SW, PLW1, PLW2,                  &
                                 PTA, PQA, PVMOD, PPS, PRHOA, PSR,         &
                                 PZREF, PUREF,                             &
                                 PRNSNOW, PHSNOW, PLESNOW, PGSNOW, PMELT   )
!   ##########################################################################
!
!!****  *SNOW_COVER_1LAYER*  
!!
!!    PURPOSE
!!    -------
!
!     One layer snow mantel scheme
!         
!     
!!**  METHOD
!     ------
!
!
! The temperature equation is written as:
!
!              b T+ = y
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/09/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,       ONLY : XTT, XCI, XRHOLI, XRHOLW, XCPD, XLSTT, XLMTT, XDAY, XCONDI
USE MODD_SNOW_PAR,   ONLY : XEMISSN, SWE_CRIT
!
USE MODE_THERMOS
!
USE MODI_SURFACE_RI
USE MODI_SURFACE_AERO_COND
!
implicit none
!!!#include <arch_specific.hf>
!
REAL  XUNDEF
PARAMETER (XUNDEF=999.)
!*      0.1    declarations of arguments
!
!
REAL,                 INTENT(IN)    :: PTSTEP   ! time step
REAL,                 INTENT(IN)    :: PANSMIN  ! minimum snow albedo
REAL,                 INTENT(IN)    :: PANSMAX  ! maximum snow albedo
REAL,                 INTENT(IN)    :: PTODRY   ! snow albedo decreasing constant
REAL,                 INTENT(IN)    :: PRHOSMIN ! minimum snow density
REAL,                 INTENT(IN)    :: PRHOSMAX ! maximum snow density
REAL,                 INTENT(IN)    :: PRHOFOLD ! snow density increasing constant
LOGICAL,              INTENT(IN)    :: OALL_MELT! T --> all snow runs off if
                                                ! lower surf. temperature is
                                                ! positive
REAL,                 INTENT(IN)    :: PDRAIN_TIME ! drainage folding time (days)
REAL,                 INTENT(IN)    :: PWCRN    ! critical snow amount necessary
                                                ! to cover the considered surface
REAL,                 INTENT(IN)    :: PZ0SN    ! snow roughness length for momentum
REAL,                 INTENT(IN)    :: PZ0HSN   ! snow roughness length for heat
REAL, DIMENSION(:), INTENT(INOUT) :: PWSNOW   ! snow reservoir (kg/m2)
REAL, DIMENSION(:), INTENT(INOUT) :: PTSNOW   ! snow temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PASNOW   ! snow albedo
REAL, DIMENSION(:), INTENT(INOUT) :: PRSNOW   ! snow density
REAL, DIMENSION(:), INTENT(INOUT) :: PTS_SNOW ! snow surface temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PESNOW   ! snow emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! underlying ground temperature
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW  ! absorbed SW energy (Wm-2)
REAL, DIMENSION(:), INTENT(IN)    :: PLW1     ! LW coef independant of TSNOW
                                              ! (Wm-2)     usually equal to:
                                              !      emis_snow * LW_down
                                              !
REAL, DIMENSION(:), INTENT(IN)    :: PLW2     ! LW coef dependant   of TSNOW
                                              ! (Wm-2 K-4) usually equal to:
                                              ! -1 * emis_snow * stefan_constant
                                              !
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
                                                ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PPS      ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA    ! air density
                                                ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PSR      ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the first
                                              ! atmospheric level (wind)
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT   ! snow melting rate (kg/m2/s)
!
!
!*      0.2    declarations of local variables
!
REAL :: ZEXPL = 0.
REAL :: ZIMPL = 1.
!
REAL, DIMENSION(SIZE(PWSNOW)) :: ZEXNS, ZEXNA, ZDIRCOSZW
REAL, DIMENSION(SIZE(PWSNOW)) :: ZZ0      ! roughness length for momentum
REAL, DIMENSION(SIZE(PWSNOW)) :: ZZ0H     ! roughness length forheat
!
REAL, DIMENSION(SIZE(PWSNOW)) :: ZRI      ! Richardson number
REAL, DIMENSION(SIZE(PWSNOW)) :: ZAC      ! aerodynamical conductance
REAL, DIMENSION(SIZE(PWSNOW)) :: ZB, ZY   ! coefficients in Ts eq.
REAL, DIMENSION(SIZE(PWSNOW)) :: ZWSNOW   ! snow before evolution
REAL, DIMENSION(SIZE(PWSNOW)) :: ZSNOW_HC ! snow heat capacity
REAL, DIMENSION(SIZE(PWSNOW)) :: ZSNOW_TC ! snow thermal conductivity
REAL, DIMENSION(SIZE(PWSNOW)) :: ZSNOW_D  ! snow depth
REAL, DIMENSION(SIZE(PWSNOW)) :: ZMELT    ! snow melting rate (kg/m3/s)
REAL, DIMENSION(SIZE(PWSNOW)) :: ZTS_SNOW ! snow surface temperature
                                          ! at previous time-step
REAL, DIMENSION(SIZE(PWSNOW)) :: ZQSAT    ! specific humidity
!                                         ! for ice
REAL, DIMENSION(SIZE(PWSNOW)) :: ZDQSAT   ! d(specific humidity)/dT
!                                         ! for ice
!
REAL, DIMENSION(SIZE(PWSNOW)) :: ZSR1, ZSR2   ! norm. snow precip.
!
LOGICAL, DIMENSION(SIZE(PWSNOW)) :: GSNOWMASK ! where snow is
!                                             ! at previuos time-step
LOGICAL, DIMENSION(SIZE(PWSNOW)) :: GFLUXMASK ! where fluxes can
!                                             ! be computed at
!                                             ! new time-step
!                                             ! i.e. snow occurence
!                                             ! at previous time-step
!                                             ! OR snow fall
!
REAL :: ZWSNOW_MIN = 0.1 ! minimum value of snow content (kg/m2) for prognostic
!                        ! computations
REAL :: MIN_TDIFF = 0.01 ! minimum temperature diff. for flux calculations
!
!-------------------------------------------------------------------------------
!
ZB=0.
ZY=0.
ZMELT  (:) = 0.
PMELT  (:) = 0.
PRNSNOW(:) = 0.
PHSNOW (:) = 0.
PLESNOW(:) = 0.
PGSNOW (:) = 0.
!
!* snow reservoir before evolution
!
ZWSNOW(:) = PWSNOW(:)
ZTS_SNOW(:) = MIN(XTT,PTG(:))
!
ZSNOW_D (:) = 0.
ZSNOW_TC(:) = 0.
ZSNOW_HC(:) = 0.
!
!-------------------------------------------------------------------------------
!
!*      1.1    most useful masks
!              -----------------
!
!* snow occurence at previous time-step
!
GSNOWMASK(:) = ZWSNOW(:)>0.
!
!* snow occurence during the time-step for fluxes computation
!
GFLUXMASK(:) = GSNOWMASK(:) .OR. PSR(:)>0.
!
!* surface temperature
!
WHERE (GSNOWMASK(:))
  ZTS_SNOW(:) = PTS_SNOW(:)
ENDWHERE

!-------------------------------------------------------------------------------
!
!*      1.2    lower limit of snow cover for prognostic computations:
!              -----------------------------------------------------
!
!              0.1 kg/m2 of snow water content
!
!
WHERE (ZWSNOW(:)<ZWSNOW_MIN)
  PTSNOW  (:) = MIN( PTG(:) , XTT )
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      1.3    drag
!              ----
!
!*      1.3.1  defaults
!              --------
!
!* variation of temperature with altitude is neglected
!
ZEXNS(:) = 1.
ZEXNA(:) = 1.
!
!* slope is neglected in drag computation
!
ZDIRCOSZW(:) = 1.
!
!* roughness length are imposed:
!
ZZ0   (:) = PZ0SN
ZZ0H  (:) = PZ0HSN
!
!
!*      1.3.2   qsat (Tsnow)
!              ------------
!
ZQSAT(:) = QSATI(ZTS_SNOW(:), PPS(:) )
!
!*      1.3.3  Richardson number
!              -----------------
!
!* snow is present on all the considered surface.
!* computation occurs where snow is and/or falls.
!
CALL SURFACE_RI(ZTS_SNOW, ZQSAT, ZEXNS, ZEXNA, PTA, PQA, &
                PZREF, PUREF, ZDIRCOSZW, PVMOD, ZRI      )
!
!*      1.3.4  Aerodynamical conductance
!              -------------------------
!
CALL SURFACE_AERO_COND(ZRI, PZREF, PUREF, PVMOD, ZZ0, ZZ0H, ZAC)
!
!-------------------------------------------------------------------------------

!
!*      2.     snow thermal characteristics
!              ----------------------------
!
!*      2.1    snow heat capacity
!              ------------------
!
WHERE (GSNOWMASK(:))
  ZSNOW_HC(:) = PRSNOW(:) * XCI * XRHOLI / XRHOLW
ELSEWHERE
  ZSNOW_HC(:) = PRHOSMIN * XCI * XRHOLI / XRHOLW
ENDWHERE
!
!*      2.2    snow depth
!              ----------
!
WHERE (GSNOWMASK(:))
  ZSNOW_D(:) = ZWSNOW(:) / PRSNOW(:)
ELSEWHERE
  ZSNOW_D(:) = PTSTEP * PSR(:) / PRHOSMIN
END WHERE
!
!*      2.3    snow thermal conductivity
!              -------------------------
!
WHERE (GSNOWMASK(:))
  ZSNOW_TC(:) = XCONDI * (PRSNOW(:)/XRHOLW)**1.885
ELSEWHERE
  ZSNOW_TC(:) = XCONDI * (PRHOSMIN /XRHOLW)**1.885
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      3.     Snow temperature evolution
!              --------------------------
!
!*      3.1    coefficients from Temperature tendency
!              --------------------------------------
!
WHERE (GSNOWMASK(:) .AND. ZWSNOW(:)>=ZWSNOW_MIN)
!
  ZB(:) = ZB(:) + ZSNOW_D(:) * ZSNOW_HC(:) / PTSTEP
!
  ZY(:) = ZY(:) + ZSNOW_D(:) * ZSNOW_HC(:) / PTSTEP * PTSNOW(:)
!
!*      3.2    coefficients from solar radiation
!              ---------------------------------
!
  ZY(:) = ZY(:) + PABS_SW(:)
!
!*      3.3    coefficients from infra-red radiation
!              -------------------------------------
!
  ZB(:) = ZB(:) - PLW2 * 4. * ZIMPL * PTSNOW(:)**3
!
  ZY(:) = ZY(:) + PLW1 + PLW2 * (ZEXPL-3.*ZIMPL) * PTSNOW(:)**4
!
!
!*      3.4    coefficients from sensible heat flux
!              ------------------------------------
!
  ZB(:) = ZB(:) + XCPD * PRHOA * ZAC *   ZIMPL
!
  ZY(:) = ZY(:) - XCPD * PRHOA * ZAC * ( ZEXPL * PTSNOW(:) - PTA(:) )
!
!*      3.5    dqsat/ dT (Tsnow)
!              -----------------
!
  ZDQSAT(:) = DQSATI(ZTS_SNOW(:), PPS(:), ZQSAT(:) )
!
!*      3.6    coefficients from latent heat flux
!              ----------------------------------
!
  ZB(:) = ZB(:) + XLSTT * PRHOA * ZAC *    ZIMPL * ZDQSAT(:)
!
  ZY(:) = ZY(:) - XLSTT * PRHOA * ZAC * (  ZQSAT(:) - PQA(:)           &
                                         - ZIMPL * ZDQSAT(:)*PTSNOW(:) )
!
!*      3.7    coefficients from conduction flux at snow base
!              ----------------------------------------------
!
  ZB(:) = ZB(:) + ZSNOW_TC/(0.5*ZSNOW_D) *  ZIMPL
!
  ZY(:) = ZY(:) - ZSNOW_TC/(0.5*ZSNOW_D) * (ZEXPL * PTSNOW(:) - PTG(:))
!
!*      3.8    guess of snow temperature before accumulation and melting
!              ---------------------------------------------------------
!
  PTSNOW(:) = ZY(:) / ZB(:)
!
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      4.     Snow melt
!              ---------
!
!*      4.1    melting
!              -------
!
WHERE (GSNOWMASK(:))
!
  ZMELT(:)  = MAX( PTSNOW(:) - XTT , 0. ) * ZSNOW_HC(:) /  XLMTT / PTSTEP
!
  ZMELT(:)  = MIN( ZMELT(:) , ZWSNOW(:) / ZSNOW_D(:) / PTSTEP )
!
  PTSNOW(:) = MIN( PTSNOW(:) , XTT )
!
END WHERE
!
!*      4.2    run-off of all snow if lower surface temperature is positive
!              ------------------------------------------------------------
!
!* this option is used when snow is located on sloping roofs for example.
!
IF (OALL_MELT) THEN
  WHERE ( GSNOWMASK(:) .AND. PTG(:)>XTT .AND. ZWSNOW(:)>=ZWSNOW_MIN )
    PMELT(:) = PMELT(:) + ZWSNOW(:) / PTSTEP
  END WHERE
END IF
!
!*      4.3    output melting in kg/m2/s
!              -------------------------
!
PMELT(:) = ZMELT(:) * ZSNOW_D(:)
!
!-------------------------------------------------------------------------------
!
!*      5.     fluxes
!              ------
!
!*      5.1    net radiation (with Ts lin. extrapolation)
!              -------------
!
WHERE (GFLUXMASK(:))
!
  PRNSNOW(:) = PABS_SW(:) + PLW1(:) + PLW2(:) * PTSNOW(:)**4
!
!
!*      5.2    sensible heat flux
!              ------------------
!
! Only calc. sensible heat flux if diff. is bigger than min_tdiff
  WHERE( ABS(PTSNOW(:)-PTA(:)).GE. MIN_TDIFF )
     PHSNOW(:) = XCPD * PRHOA * ZAC * ( PTSNOW(:) - PTA(:) )
  END WHERE
   
!
!*      5.3    qsat (Tsnow)
!              ------------
!
  ZQSAT(:) = QSATI(PTSNOW(:) , PPS(:) )
!
!*      5.4    latent heat flux
!              ----------------
!
  PLESNOW(:) = XLSTT * PRHOA * ZAC * ( ZQSAT(:) - PQA(:) )
!
!
!*      5.5    Conduction heat flux
!              --------------------
!
! Only calc. conduction heat flux if diff. is bigger than min_tdiff
  WHERE( ABS(PTSNOW(:)-PTG(:)).GE. MIN_TDIFF )
     PGSNOW(:) = ZSNOW_TC(:)/(0.5*ZSNOW_D(:)) * ( PTSNOW(:) - PTG(:) )
  END WHERE
!
END WHERE
!
!*      5.6    If ground T>0°C, Melting is estimated from conduction heat flux
!              ---------------------------------------------------------------
!
WHERE (GFLUXMASK(:) .AND. PTG(:)>XTT)
  PMELT(:) = MAX (PMELT(:), -PGSNOW(:)/XLMTT)
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      6.     reservoir evolution
!              -------------------
!
!*      6.1    snow fall
!              ---------
!
PWSNOW(:) = PWSNOW(:) + PTSTEP * PSR(:)
!
!
!*      6.2    sublimation
!              -----------
!
PLESNOW(:) = MIN( PLESNOW(:), XLSTT*PWSNOW/PTSTEP )
!
PWSNOW(:)  = MAX( PWSNOW(:) - PTSTEP * PLESNOW(:)/XLSTT , 0.)
!
WHERE ( PWSNOW(:)<1.E-8 * PTSTEP ) PWSNOW(:) = 0.
!
!*      6.3    melting
!              -------
!
PMELT(:) = MIN( PMELT(:), PWSNOW/PTSTEP )
!
PWSNOW(:)= MAX( PWSNOW(:) - PTSTEP * PMELT(:) , 0.)
!
WHERE ( PWSNOW(:)<1.E-8 * PTSTEP ) PWSNOW(:) = 0.
!
WHERE(PWSNOW(:)==0.) PGSNOW(:) = MAX ( PGSNOW(:), - PMELT(:)*XLMTT )
!
!
!*      6.4    time dependent drainage
!              -----------------------
!
IF (PDRAIN_TIME>0.) THEN
  WHERE ( PWSNOW(:)>0.)
    PWSNOW(:) = PWSNOW(:) * EXP(-PTSTEP/PDRAIN_TIME/XDAY)
  END WHERE
END IF
!
!*      6.5    melting of last 1mm of snow depth
!              ---------------------------------
!
WHERE ( PWSNOW(:)<ZWSNOW_MIN .AND. PMELT(:)>0. .AND. PSR(:)==0. )
  PMELT(:) = PMELT(:) + PWSNOW(:) / PTSTEP
  PWSNOW(:)=0.
END WHERE
!
WHERE ( PWSNOW(:)<1.E-8 * PTSTEP ) PWSNOW(:) = 0.
!
!-------------------------------------------------------------------------------
!
!*      7.     albedo evolution
!              ----------------
!
!*      7.1    If melting occurs
!              -----------------
!
WHERE ( GSNOWMASK(:) .AND. PMELT>0. )
!
  PASNOW(:) = (PASNOW(:)-PANSMIN)*EXP(-0.01*PTSTEP/3600.) + PANSMIN   &
              + PSR(:)*PTSTEP/PWCRN*PANSMAX
!
END WHERE
!
!*      7.2    If no melting occurs
!              --------------------
!
WHERE ( GSNOWMASK(:) .AND. PMELT==0. )
  PASNOW(:) = PASNOW(:) - PTODRY*PTSTEP/XDAY                          &
              + PSR(:)*PTSTEP/PWCRN*PANSMAX
END WHERE
!
!*      7.3    Limits
!              ------
!
WHERE (  PWSNOW(:)>0. )
  PASNOW(:) = MAX(PASNOW(:),PANSMIN)
  PASNOW(:) = MIN(PASNOW(:),PANSMAX)
END WHERE
!
!*      7.4    fresh snow
!              ----------
!
WHERE (  ZWSNOW(:)==0. .AND. PWSNOW(:)>0. )
  PASNOW(:) = PANSMAX
  PESNOW(:) = XEMISSN
END WHERE
!-------------------------------------------------------------------------------
!
!*      8.     density evolution
!              -----------------
!
!*      8.1    old snow
!              --------
!
WHERE ( GSNOWMASK(:) .AND. PWSNOW(:)>0. )
  ZSR1(:) = MAX( PWSNOW(:) , PSR(:) * PTSTEP )
!
  PRSNOW(:) = (PRSNOW(:)-PRHOSMAX)*EXP(-PRHOFOLD*PTSTEP/3600.) + PRHOSMAX
  PRSNOW(:) = ( (ZSR1(:)-PSR(:)*PTSTEP) * PRSNOW(:)    &
              + (PSR(:)*PTSTEP) * PRHOSMIN ) / ZSR1(:)

! IMPOSE MIN/MAX ON SNOW DENSITY
  PRSNOW(:) = MIN(   MAX(PRSNOW(:),PRHOSMIN)  , PRHOSMAX)
END WHERE
!
!*      8.2    fresh snow
!              ----------
!
WHERE (  ZWSNOW(:)==0. .AND. PWSNOW(:)>0. )
  PRSNOW(:) = PRHOSMIN
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      9.     fresh snow accumulation (if more than 1mm of snow depth)
!              -----------------------
!
WHERE (  ZWSNOW(:)>=ZWSNOW_MIN .AND. PSR(:)>0. .AND. PWSNOW(:)>0. )
  ZSR2(:) = MIN( PWSNOW(:) , PSR(:) * PTSTEP )
!
  PTSNOW(:) =( ( PWSNOW(:) - ZSR2(:) ) *      PTSNOW(:)        &
            +                ZSR2(:)   * MIN( PTA   (:) ,XTT ))&
              /(   PWSNOW(:) )
END WHERE
!
!-------------------------------------------------------------------------------
!
!*     10.     Surface temperature
!              -------------------
!
!* note that if the relation between snow pack temperature and its
!  surface temperature is modified, think to modify it also in
!  subroutine init_snow_lw.f90
!
WHERE (GSNOWMASK(:) )
  PTS_SNOW(:) = PTSNOW(:)
END WHERE
!
!-------------------------------------------------------------------------------
!
!*     11.     bogus values
!              ------------
!
!*     11.1    snow characteristics where snow IS present at current time-step
!              ---------------------------------------------------------------
!
WHERE ( PWSNOW(:)< SWE_CRIT )
  PWSNOW  (:) = 0.0
  PTSNOW  (:) = XUNDEF
  PRSNOW  (:) = XUNDEF
  PASNOW  (:) = XUNDEF
  PTS_SNOW(:) = XUNDEF
  PESNOW  (:) = XUNDEF
END WHERE
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW_COVER_1LAYER
