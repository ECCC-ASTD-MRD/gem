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
!     ###############
      MODULE MODD_TEB
!     ###############
!
!!****  *MODD_TEB - declaration of surface parameters for urban surface
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!	               09/2009   S. Leroyer - optimization for GEM
!                      05/2013   M. Abrahamowicz - add new variables
!
!*       0.   DECLARATIONS
!             ------------
!
!USE MODD_TYPE_SNOW
!
implicit none
!
common/tttd/XZS,XCOVER,XBLD,XCAN_HW_RATIO,XBLD_HEIGHT,XWALL_O_HOR,  &
XZ0_TOWN,XZ0_ROOF,XZ0_ROAD, &
XSVF_ROAD,XSVF_WALL,XALB_ROOF,XEMIS_ROOF,XHC_ROOF,XTC_ROOF,&
XD_ROOF,XALB_ROAD,XEMIS_ROAD,XHC_ROAD,XTC_ROAD,XD_ROAD,XALB_WALL,   &
XEMIS_WALL,XHC_WALL,XTC_WALL,XD_WALL,XH_TRAFFIC,XLE_TRAFFIC,        &
XH_INDUSTRY,XLE_INDUSTRY,XTI_BLD,XTI_ROAD,XWS_ROOF,XWS_ROAD,XT_ROOF,&
XT_ROAD,XT_WALL,XT_CANYON,XQ_CANYON,&
XWSNOW_ROOF,  XTSNOW_ROOF, XRSNOW_ROOF,&
XASNOW_ROOF, XESNOW_ROOF, XTSSNOW_ROOF,&
XWSNOW_ROAD,  XTSNOW_ROAD, XRSNOW_ROAD,&
XASNOW_ROAD, XESNOW_ROAD, XTSSNOW_ROAD,&
XTSTEP,XOUT_TSTEP,              &
LCOVER,NROOF_LAYER,NROAD_LAYER,NWALL_LAYER
!$OMP THREADPRIVATE(/tttd/)
#define ALLOCATABLE POINTER
!
! General surface: 
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XZS           ! orography                        (m)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XCOVER        ! fraction of each ecosystem       (-)
LOGICAL, ALLOCATABLE, DIMENSION(:):: LCOVER        ! GCOVER(i)=T --> ith cover field is not 0.
!
! Number of layers
!
INTEGER                           :: NROOF_LAYER   ! number of layers in roofs
INTEGER                           :: NROAD_LAYER   ! number of layers in roads
INTEGER                           :: NWALL_LAYER   ! number of layers in walls
!
! Geometric Parameters:
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XBLD          ! fraction of buildings            (-)
REAL, ALLOCATABLE, DIMENSION(:)   :: XCAN_HW_RATIO ! canyon    h/W                    (-)
REAL, ALLOCATABLE, DIMENSION(:)   :: XBLD_HEIGHT   ! buildings height 'h'             (m)
REAL, ALLOCATABLE, DIMENSION(:)   :: XWALL_O_HOR   ! wall surf. / hor. surf.          (-)
REAL, ALLOCATABLE, DIMENSION(:)   :: XZ0_TOWN      ! roughness length for momentum    (m)
REAL, ALLOCATABLE, DIMENSION(:)   :: XZ0_ROOF      ! Roof roughness length for momentum  (m)
REAL, ALLOCATABLE, DIMENSION(:)   :: XZ0_ROAD      ! Road roughness length for momentum  (m)
REAL, ALLOCATABLE, DIMENSION(:)   :: XSVF_ROAD     ! road sky view factor             (-)
REAL, ALLOCATABLE, DIMENSION(:)   :: XSVF_WALL     ! wall sky view factor             (-)
!
! Roof parameters
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XALB_ROOF     ! roof albedo                      (-)
REAL, ALLOCATABLE, DIMENSION(:)   :: XEMIS_ROOF    ! roof emissivity                  (-)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XHC_ROOF      ! roof layers heat capacity        (J/K/m3)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XTC_ROOF      ! roof layers thermal conductivity (W/K/m)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XD_ROOF       ! depth of roof layers             (m)
!
!
! Road parameters
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XALB_ROAD     ! road albedo                      (-)
REAL, ALLOCATABLE, DIMENSION(:)   :: XEMIS_ROAD    ! road emissivity                  (-)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XHC_ROAD      ! road layers heat capacity        (J/K/m3)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XTC_ROAD      ! road layers thermal conductivity (W/K/m)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XD_ROAD       ! depth of road layers             (m)
!
! Wall parameters
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XALB_WALL     ! wall albedo                      (-)
REAL, ALLOCATABLE, DIMENSION(:)   :: XEMIS_WALL    ! wall emissivity                  (-)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XHC_WALL      ! wall layers heat capacity        (J/K/m3)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XTC_WALL      ! wall layers thermal conductivity (W/K/m)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XD_WALL       ! depth of wall layers             (m)
!
! anthropogenic fluxes
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic       (W/m2)
REAL, ALLOCATABLE, DIMENSION(:)   :: XLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic       (W/m2)
REAL, ALLOCATABLE, DIMENSION(:)   :: XH_INDUSTRY   ! anthropogenic sensible                   
!                                                  ! heat fluxes due to factories     (W/m2)
REAL, ALLOCATABLE, DIMENSION(:)   :: XLE_INDUSTRY  ! anthropogenic latent
!                                                  ! heat fluxes due to factories     (W/m2)
!
! temperatures for boundary conditions
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XTI_BLD       ! building interior temperature    (K)
REAL, ALLOCATABLE, DIMENSION(:)   :: XTI_ROAD      ! road interior temperature        (K)
!
! Prognostic variables:
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XWS_ROOF      ! roof water reservoir             (kg/m2)
REAL, ALLOCATABLE, DIMENSION(:)   :: XWS_ROAD      ! road water reservoir             (kg/m2)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XT_ROOF       ! roof layer temperatures          (K)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XT_ROAD       ! road layer temperatures          (K)
REAL, ALLOCATABLE, DIMENSION(:,:) :: XT_WALL       ! wall layer temperatures          (K)
!
! Semi-prognostic variables:
!
REAL, ALLOCATABLE, DIMENSION(:)   :: XT_CANYON     ! canyon air temperature           (K)
REAL, ALLOCATABLE, DIMENSION(:)   :: XQ_CANYON     ! canyon air specific humidity     (kg/kg)
!
! Prognostic snow:
REAL, ALLOCATABLE, DIMENSION(:)   :: XWSNOW_ROOF ! snow layers reservoir              (kg/m2)
REAL, ALLOCATABLE, DIMENSION(:)   :: XTSNOW_ROOF ! snow layers temperature            (K)
REAL, ALLOCATABLE, DIMENSION(:)   :: XRSNOW_ROOF ! snow layers density
REAL, ALLOCATABLE, DIMENSION(:)   :: XASNOW_ROOF ! snow albedo
REAL, ALLOCATABLE, DIMENSION(:)   :: XESNOW_ROOF ! snow emissivity
REAL, ALLOCATABLE, DIMENSION(:)   :: XTSSNOW_ROOF! snow surface temperature           (K) 
REAL, ALLOCATABLE, DIMENSION(:)   :: XWSNOW_ROAD ! snow layers reservoir              (kg/m2)
REAL, ALLOCATABLE, DIMENSION(:)   :: XTSNOW_ROAD ! snow layers temperature            (K)
REAL, ALLOCATABLE, DIMENSION(:)   :: XRSNOW_ROAD ! snow layers density
REAL, ALLOCATABLE, DIMENSION(:)   :: XASNOW_ROAD ! snow albedo
REAL, ALLOCATABLE, DIMENSION(:)   :: XESNOW_ROAD ! snow emissivity
REAL, ALLOCATABLE, DIMENSION(:)   :: XTSSNOW_ROAD! snow surface temperature           (K)
!
! Time-step:
!
REAL                              :: XTSTEP        ! time step for TEB
!
REAL                              :: XOUT_TSTEP    ! TEB output writing time step
!
END MODULE MODD_TEB
