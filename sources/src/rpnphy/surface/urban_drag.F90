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
    SUBROUTINE URBAN_DRAG2(PTSTEP, PT_CANYON, PQ_CANYON,                    &
                          PTS_ROOF, PTS_ROAD, PTS_WALL, PDELT_SNOW_ROOF,    &
                          PEXNS, PEXNA, PTA, PQA, PPS, PRHOA,               &
                          PZREF, PUREF, PVMOD, PVDIR, PLAT,                 &
                          PZ0_TOWN, PZ0_ROOF, PZ0_ROAD,                     &
                          PBLD, PBLD_HEIGHT, PCAN_HW_RATIO,                 &
                          PWALL_O_ROAD,                                     &
                          PWS_ROOF, PWS_ROAD,                               &
                          PWS_ROOF_MAX, PWS_ROAD_MAX,                       &
                          PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD,   &
                          PCD, PAC_ROOF, PAC_ROOF_WAT,                      &
                          PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,        &
                          PU_CAN, PRI,                                      &
                          ZTRDZT,ZTRFZT,ZUZT,ZVZT                           )
!   ##########################################################################
       use sfclayer_mod, only: sl_sfclayer,SL_OK
       use sfc_options    
!
!!****  *URBAN_DRAG*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the surface drag over artificial surfaces as towns, 
!     taking into account the canyon like geometry of urbanized areas.
!         
!     
!!**  METHOD
!!    ------
!
!
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
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
!!      Original    20/01/98 
!!          01/00 (V. Masson)  separation of skimming, wake and isolated flows
!!          09/00 (V. Masson)  use of Z0 for roads
!!          12/02 (A. Lemonsu) convective speed w* in canyon
!            /04 (A. Lemonsu) z0h=z0m for resistance canyon-atmosphere
!          03/08 (S. Leroyer) debug PU_CAN (1. * H/3)
!          12/08 (S. Leroyer) call flxsurf for new options of z0h applied on roof, road and town
!          04/12 (S. Leroyer) z0h(s) determined in compz0 called in flxsurf4 
!          11/12 (M. Abrahamowicz) Impose min of 0.01 on PU_CAN 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS, ONLY : XLVTT, XPI, XCPD, XG
USE MODD_TOWN, ONLY : XQ_TOWN
!
USE MODE_THERMOS
!
USE MODI_SURFACE_RI
USE MODI_SURFACE_CD
USE MODI_SURFACE_AERO_COND
!
implicit none
#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
!
!
REAL,               INTENT(IN)    :: PTSTEP         ! time-step
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON      ! canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_CANYON      ! canyon air specific humidity.
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROOF       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_SNOW_ROOF! fraction of snow
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS          ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! specific humidity
                                                    ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD          ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PVDIR          ! direction of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PLAT           ! latitude
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA          ! exner function
                                                    ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PZREF          ! reference height of the first
                                                    ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF          ! reference height of the first
                                                    ! atmospheric level (wind)
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_TOWN       ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_ROOF       ! roughness length for momentum for roof
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_ROAD       ! roughness length for momentum for road
REAL, DIMENSION(:), INTENT(IN)    :: PBLD           ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PBLD_HEIGHT    ! h
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO  ! h/W
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_ROAD   ! wall surf. / road surf.
!
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROOF       ! roof water content (kg/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROAD       ! road water content (kg/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROOF_MAX   ! maximum deepness of roof
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROAD_MAX   ! and water reservoirs (kg/m2)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROOF     ! qsat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROAD     ! qsat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROOF     ! water fraction on
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROAD     ! snow-free surfaces
REAL, DIMENSION(:), INTENT(OUT)   :: PCD            ! drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WALL       ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! walls 
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD       ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! roads
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD_WAT   ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! road (for water)
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP        ! aerodynamical conductance
!                                                   ! between canyon top and atm.
REAL, DIMENSION(:), INTENT(OUT)   :: PU_CAN         ! hor. wind in canyon
REAL, DIMENSION(:), INTENT(OUT)   :: PRI            ! Town Richardson number
!
REAL, DIMENSION(:), INTENT(OUT)   :: ztrfzt ! output for diasurf T,q at z=zt over roof
REAL, DIMENSION(:), INTENT(OUT)   :: ztrdzt ! output for diasurf T,q at z=zt over road
REAL, DIMENSION(:), INTENT(OUT)   :: zuzt  ! output for diasurf u,v at z=zt over roof
REAL, DIMENSION(:), INTENT(OUT)   :: zvzt  ! output for diasurf u,v at z=zt over roof
!REAL, DIMENSION(:), INTENT(OUT)   :: zurdzt ! output for diasurf T,q at z=zt over road
!REAL, DIMENSION(:), INTENT(OUT)   :: zvrdzt ! output for diasurf T,q at z=zt over road
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZTS_TOWN     ! town averaged temp.
REAL, DIMENSION(SIZE(PTA)) :: ZQ_TOWN      ! town averaged hum.
REAL, DIMENSION(SIZE(PTA)) :: ZAVDELT_ROOF ! averaged water frac.
REAL, DIMENSION(SIZE(PTA)) :: ZQ_ROOF      ! roof spec. hum.
REAL, DIMENSION(SIZE(PTA)) :: ZDIRCOSZW    ! orography slope cosine (=1 in TEB)
REAL, DIMENSION(SIZE(PTA)) :: ZWAKE        ! reduction of average wind speed
                                           ! in canyon due to direction average.
REAL, DIMENSION(SIZE(PTA)) :: ZW_CAN       ! ver. wind in canyon
REAL, DIMENSION(SIZE(PTA)) :: ZRI          ! Richardson number
REAL, DIMENSION(SIZE(PTA)) :: ZCDN         ! neutral drag coefficient
REAL, DIMENSION(SIZE(PTA)) :: ZLE_MAX      ! maximum latent heat flux available
REAL, DIMENSION(SIZE(PTA)) :: ZLE          ! actual latent heat flux
!
REAL, DIMENSION(SIZE(PTA)) :: ZU_STAR, ZW_STAR !! 
REAL, DIMENSION(SIZE(PTA)) :: ZQ0              !! 
INTEGER                   ::  JLOOP, STAT      !! 
!
REAL,DIMENSION(SIZE(PTA)) :: cmu, ctu, ue, fcor                  ! temp var for sfc layer
REAL,DIMENSION(SIZE(PTA)) :: z0h_roof,z0h_town,z0h_road,z0h_can  ! local thermal roughness
REAL,DIMENSION(SIZE(PTA)) :: ztemp  ! temporary array for -not used- diagnostic var (to avoid print in the listing)
real,    parameter :: XZT = 2.5  ! height for diag calc. different from zt

INTEGER N
   logical, parameter :: TOWN_TDIAGLIM = .false.
!
!-------------------------------------------------------------------------------
!
!
ZDIRCOSZW=1. ! no orography slope effect taken into account in TEB
!
fcor(:)=1.0372462E-04  !really ... we can't do better than this?
N=SIZE(PTA)
!-------------------------------------------------------------------------------
!
!*      1.     roof and road saturation specific humidity
!              ------------------------------------------
!
PQSAT_ROOF(:) =  QSAT(PTS_ROOF(:),PPS(:))
!
PQSAT_ROAD(:) =  QSAT(PTS_ROAD(:),PPS(:))
!
!-------------------------------------------------------------------------------
!
!*      2.     fraction of water on roofs
!              --------------------------
!
PDELT_ROOF=1.
!
!*      2.1    general case
!              ------------
!
WHERE (PQSAT_ROOF(:) >= PQA(:) )
  PDELT_ROOF(:) = (PWS_ROOF(:)/PWS_ROOF_MAX)**(2./3.)
END WHERE
!
!*      2.2    dew deposition on roofs (PDELT_ROOF=1)
!              -----------------------
!
!-------------------------------------------------------------------------------
!
!*      3.     fraction of water on roads
!              --------------------------
!
PDELT_ROAD=1.
!
!*      3.1    general case
!              ------------
!
WHERE (PQSAT_ROAD(:) >= PQ_CANYON(:) )
  PDELT_ROAD(:) = (PWS_ROAD(:)/PWS_ROAD_MAX)**(2./3.)
END WHERE
!
!
!*      3.2    dew deposition on roads (PDELT_ROAD=1)
!              -----------------------
!
!-------------------------------------------------------------------------------
!
!*      4.     Drag coefficient for momentum between roof level and atmosphere
!              ---------------------------------------------------------------
!
!
!*      4.1    Averaged temperature at roof level
!              ----------------------------------
!
ZTS_TOWN(:) = PBLD(:) * PTS_ROOF(:) + (1.-PBLD(:)) * PT_CANYON(:)
!
!*      4.2    Averaged water fraction on roofs
!              -------------------------------
!
ZAVDELT_ROOF(:) = PDELT_ROOF(:) * PDELT_SNOW_ROOF(:)
!
!*      4.3    Roof specific humidity
!              ----------------------
!
ZQ_ROOF(:) = QSAT(PTS_ROOF(:),PPS(:)) * ZAVDELT_ROOF(:)
!
!*      4.4    Averaged Saturation specific humidity
!              -------------------------------------
!
ZQ_TOWN(:) =       PBLD(:)  * ZQ_ROOF(:) &
             + (1.-PBLD(:)) * PQ_CANYON  (:)
!
XQ_TOWN(:) = ZQ_TOWN(:)
!
!
!*      4.5    first guess of z0h
! note : it will be better to use the one computed at the last time step (new argument or bus)
!              -------------------------------------
! First guess of u*
ZU_STAR(:) = 0.4 * PVMOD(:) / LOG( (PUREF+ PBLD_HEIGHT(:)/3.)/PZ0_TOWN(:) )
!
Z0H_TOWN(:)= MAX( PZ0_TOWN(:) * 7.4 * EXP( - 1.29 *( PZ0_TOWN(:)*ZU_STAR(:)/1.461e-5)**0.25), 1.e-20 )
Z0H_ROOF(:)= MAX( PZ0_ROOF(:) * 7.4 * EXP( - 1.29 *( PZ0_ROOF(:)*ZU_STAR(:)/1.461e-5)**0.25), 1.e-20 )
Z0H_ROAD(:)= MAX( PZ0_ROAD(:) * 7.4 * EXP( - 1.29 *( PZ0_ROAD(:)*ZU_STAR(:)/1.461e-5)**0.25), 1.e-20 )
Z0H_CAN(:) = MAX( PZ0_TOWN(:),1.e-20 )
!
!-------------------------------------------------------------------------------
!
!*      5.     Momentum drag coefficient
!              -------------------------
!
! z0h computed with Kanda (2007) formulation - option 8 in compz0
!
stat = sl_sfclayer(PTA/PEXNA,PQA,PVMOD,PVDIR,PUREF+PBLD_HEIGHT/3.,PZREF+PBLD_HEIGHT/3., &
     ZTS_TOWN/PEXNS,ZQ_TOWN,PZ0_TOWN,Z0H_TOWN,PLAT,fcor,optz0=8,coefm=cmu,coeft=ctu,ue=ue,&
     tdiaglim=TOWN_TDIAGLIM)
!
   if (stat /= SL_OK) then
      call physeterror('urban_drag mom.', 'error returned by sl_sfclayer()')
      return
   endif

		PCD(:) = (cmu(:)/ue(:))**2  
!
!
PRI = ZRI
!
!-------------------------------------------------------------------------------
!
!*      6.     Drag coefficient for heat fluxes between roofs and atmosphere
!              -------------------------------------------------------------
!
! z0h computed with Kanda (2007) formulation - option 8 in compz0
!
stat = sl_sfclayer(PTA/PEXNA,PQA,PVMOD,PVDIR,PUREF,PZREF, &
     PTS_ROOF/PEXNS,ZQ_ROOF,PZ0_ROOF,Z0H_ROOF,PLAT,fcor,optz0=8,coefm=cmu,coeft=ctu,ue=ue,&
     tdiaglim=TOWN_TDIAGLIM,   &
     hghtm_diag=zt,hghtt_diag=zt,t_diag=ztrfzt,u_diag=zuzt,v_diag=zvzt)
!
	          PAC_ROOF(:) = (cmu(:)*ctu(:)/ue(:)**2)  * PVMOD(:) 
!
!
!
!
ZLE_MAX(:)     = PWS_ROOF(:) / PTSTEP * XLVTT
ZLE    (:)     =(PQSAT_ROOF(:) - PQA(:))                     &
                 * PAC_ROOF(:) * PDELT_ROOF(:) * XLVTT * PRHOA(:)
!
PAC_ROOF_WAT(:) = PAC_ROOF(:)
!
WHERE (PDELT_ROOF(:)==0.) &
PAC_ROOF_WAT(:)=0.
!
WHERE (ZLE(:)>0.) &
PAC_ROOF_WAT(:) = PAC_ROOF(:) * MIN ( 1. , ZLE_MAX(:)/ZLE(:) )
!-------------------------------------------------------------------------------
!
!*      7.     Drag coefficient for heat fluxes between canyon and atmosphere
!              --------------------------------------------------------------
!
! z0h = z0m - option 1 in compz0
!
stat = sl_sfclayer(PTA/PEXNA,PQA,PVMOD,PVDIR,PUREF+PBLD_HEIGHT/2.,PZREF+PBLD_HEIGHT/2., &
     PT_CANYON/PEXNS,PQ_CANYON,PZ0_TOWN,Z0H_can,PLAT,fcor,optz0=1,coefm=cmu,coeft=ctu,ue=ue,&
     tdiaglim=TOWN_TDIAGLIM)
!
      PAC_TOP(:) = (cmu(:)*ctu(:)/ue(:)**2)  * PVMOD(:) 
!
!-------------------------------------------------------------------------------
!
!*      8.     Drag coefficient for heat fluxes between walls, road and canyon
!              ---------------------------------------------------------------
!
!*      8.1    Horizontal wind speed in canyon
!              -------------------------------
!
!* skimming flow for h/w>1 (maximum effect of direction on wind in the canyon);
!* isolated flow for h/w<0.5 (wind is the same in large streets for all dir.)
!* wake flow between.
!
ZWAKE(:)= 1. + (2./XPI-1.) * 2. * (PCAN_HW_RATIO(:)-0.5)
ZWAKE(:)= MAX(MIN(ZWAKE(:),1.),2./XPI)
!
!* Estimation of canyon wind speed from wind just above roof level
!  (at 1.33h). Wind at 1.33h is estimated using the log law.
!
PU_CAN(:) =  ZWAKE(:) * EXP(-PCAN_HW_RATIO(:)/4.) * PVMOD(:)              &
           * LOG( (             1.* PBLD_HEIGHT(:)/3.) / PZ0_TOWN(:))     &
           / LOG( (PUREF(:)   + 1.* PBLD_HEIGHT(:)/3.) / PZ0_TOWN(:))
!
!* IMPOSE MINIMUM OF 0.01 for PU_CAN
PU_CAN(:) = MAX(PU_CAN(:),0.01)
!
!
!*      8.2    Vertical wind speed in canyon = ustar
!              -----------------------------
!
ZU_STAR(:) = SQRT (PCD(:)) * PVMOD(:)
!
!*      8.3    aerodynamical conductance for roads or walls
!              --------------------------------------------
!
ZW_STAR(:) = 0.
ZQ0(:)     = 0.
!
DO JLOOP=1,3
 !
  ZW_CAN(:)   = ZW_STAR(:) + ZU_STAR(:)
  PAC_WALL(:) = ( 11.8 + 4.2 * SQRT(PU_CAN(:)**2 + ZW_CAN(:)**2) ) &
                  / XCPD / PRHOA(:)
!
!
! z0h computed with Kanda (2007) formulation - option 8 in compz0
stat = sl_sfclayer(PT_CANYON/PEXNS,PQ_CANYON,PU_CAN+ZW_CAN,PVDIR,PBLD_HEIGHT/2.,PBLD_HEIGHT/2., &
     PTS_ROAD/PEXNS,PQ_CANYON,PZ0_ROAD,Z0H_ROAD,PLAT,fcor,optz0=8,coefm=cmu,coeft=ctu,ue=ue, &
     tdiaglim=TOWN_TDIAGLIM,hghtt_diag=zt,t_diag=ztrdzt)
!
		PAC_ROAD(:) = (cmu(:)*ctu(:)/ue(:)**2)  * (PU_CAN(:)+ZW_CAN(:)) 
!
!
  ZQ0(:)     = (PTS_ROAD(:) - PT_CANYON(:)) * PAC_ROAD(:) &
              +(PTS_WALL(:) - PT_CANYON(:)) * PAC_WALL(:) * PWALL_O_ROAD(:)
!
  IF (ALL(ZQ0(:)<0.)) EXIT
  WHERE (ZQ0(:) >= 0.)
    ZW_STAR(:) = ( (XG * PEXNA(:) / PTA(:)) * ZQ0(:) * PBLD_HEIGHT(:)) ** (1/3.)
  ELSEWHERE
    ZW_STAR(:) = 0.
  ENDWHERE
  !
END DO
!
!
!*      8.4    aerodynamical conductance for water limited by available water
!              --------------------------------------------------------------
!
ZLE_MAX(:)     = PWS_ROAD(:) / PTSTEP * XLVTT
ZLE    (:)     = ( PQSAT_ROAD(:) - PQ_CANYON(:) )                   &
                 *   PAC_ROAD(:) * PDELT_ROAD(:) * XLVTT * PRHOA(:)
!
PAC_ROAD_WAT(:) = PAC_ROAD(:)
!
WHERE (PDELT_ROAD(:)==0.) &
PAC_ROAD_WAT(:) = 0.
!
WHERE (ZLE(:)>0.) &
PAC_ROAD_WAT(:) = PAC_ROAD(:) * MIN ( 1. , ZLE_MAX(:)/ZLE(:) )
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_DRAG2
