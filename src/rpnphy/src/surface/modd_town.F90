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
!     ################
      MODULE MODD_TOWN
!     ################
!
!!****  Initialisation des variables necessaires pour TEB      
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
implicit none 
!
common/ttte/XMASK,XINI_LAT,XINI_LON,XINI_ZS,XINI_TSOIL,XINI_TS, &
XINI_TA,XINI_QA,XTOWN,XQFVH,XQFRE,XQFRH,XQFCE,XQFCH,XQFIE,XQFIH,&
XQ_TOWN,XU_CANYON,XRN_ROOF,XH_ROOF,                             &
XLE_ROOF,XLES_ROOF,XGFLUX_ROOF,XRUNOFF_ROOF,XRN_ROAD,XH_ROAD,   &
XLE_ROAD,XLES_ROAD,XGFLUX_ROAD,XRUNOFF_ROAD,XRN_WALL,XH_WALL,   &
XLE_WALL,XGFLUX_WALL,XRNSNOW_ROOF,XHSNOW_ROOF,XLESNOW_ROOF,     &
XGSNOW_ROOF,XMELT_ROOF,XRNSNOW_ROAD,XHSNOW_ROAD,XLESNOW_ROAD,   &
XGSNOW_ROAD,XMELT_ROAD,XRN,XH,XLE,XGFLUX,XEVAP,XRUNOFF,XCH,XRI,XUSTAR,  &
XTRAD_IN,XTRAD_SUN,XTRAD_SHADE,        &
XTRAD_RFSUN,XTRAD_RFSHADE,             &
XTGLOBE_SUN,XTGLOBE_SHADE,             &
XTGLOBE_RFSUN,XTGLOBE_RFSHADE,         &
XTWETB,XTWETB_ROOF,                    &
XUTCI_IN,XUTCI_OUTSUN,XUTCI_OUTSHADE,  &
XUTCI_RFSUN,XUTCI_RFSHADE,             &
XWBGT_OUTSUN,XWBGT_OUTSHADE,           &
XWBGT_RFSUN,XWBGT_RFSHADE,             &
XUTCIC_IN,XUTCIC_OUTSUN,               &
XUTCIC_OUTSHADE,XUTCIC_RFSUN,          &
XUTCIC_RFSHADE,                        &
XTRFZT,XTRDZT,XURDZU,                  &
XQ1,XQ2,XQ3,XQ4,XQ5,XQ6,XQ7,           &
XQ8,XQ9,XQ10,XQ11,XQ12,XQ13

!$OMP THREADPRIVATE(/ttte/)
!
INTEGER                           :: NNI             ! Number of grid points
REAL, DIMENSION(:)  , POINTER :: XMASK           ! Land/sea maks          
REAL, DIMENSION(:)  , POINTER :: XINI_LAT        ! Latitude              
REAL, DIMENSION(:)  , POINTER :: XINI_LON        ! Longitude             
REAL, DIMENSION(:)  , POINTER :: XINI_ZS         ! Topography            
REAL, DIMENSION(:)  , POINTER :: XINI_TSOIL      ! Soil temperature     
REAL, DIMENSION(:)  , POINTER :: XINI_TS         ! Surface temperature   
REAL, DIMENSION(:)  , POINTER :: XINI_TA         ! Air temperature at first level
REAL, DIMENSION(:)  , POINTER :: XINI_QA         ! Air specific humidity at first level
REAL, DIMENSION(:)  , POINTER :: XTOWN           ! Total fraction of urban cover
!
REAL, DIMENSION(:)  , POINTER :: XQFVH           ! Energy concumption for traffic
REAL, DIMENSION(:)  , POINTER :: XQFRE           ! Elec concumption for residential areas
REAL, DIMENSION(:)  , POINTER :: XQFRH           ! Fuel concumption for residential areas
REAL, DIMENSION(:)  , POINTER :: XQFCE           ! Elec concumption for commercial areas
REAL, DIMENSION(:)  , POINTER :: XQFCH           ! Fuel concumption for commercial areas
REAL, DIMENSION(:)  , POINTER :: XQFIE           ! Elec concumption for industrial areas
REAL, DIMENSION(:)  , POINTER :: XQFIH           ! Fuel concumption for industrial areas
!
REAL, DIMENSION(:)  , POINTER :: XQ_TOWN         ! Town averaged Specific humidity
REAL, DIMENSION(:)  , POINTER :: XU_CANYON       ! Wind in canyon         
REAL, DIMENSION(:)  , POINTER :: XRN_ROOF        ! Net radiation on roof 
REAL, DIMENSION(:)  , POINTER :: XH_ROOF         ! Sensible heat flux on roof
REAL, DIMENSION(:)  , POINTER :: XLE_ROOF        ! Latent heat flux on roof
REAL, DIMENSION(:)  , POINTER :: XLES_ROOF       ! Sublimation heat flux on roof
REAL, DIMENSION(:)  , POINTER :: XGFLUX_ROOF     ! Storage heat flux on roof
REAL, DIMENSION(:)  , POINTER :: XRUNOFF_ROOF    ! Water runoff from roof        
REAL, DIMENSION(:)  , POINTER :: XRN_ROAD        ! Net radiation on road                
REAL, DIMENSION(:)  , POINTER :: XH_ROAD         ! Sensible heat flux on road     
REAL, DIMENSION(:)  , POINTER :: XLE_ROAD        ! Latent heat flux on road            
REAL, DIMENSION(:)  , POINTER :: XLES_ROAD       ! Sublimation heat flux on road
REAL, DIMENSION(:)  , POINTER :: XGFLUX_ROAD     ! Storage heat flux on road       
REAL, DIMENSION(:)  , POINTER :: XRUNOFF_ROAD    ! Water runoff from road 
REAL, DIMENSION(:)  , POINTER :: XRN_WALL        ! Net radiation on wall                
REAL, DIMENSION(:)  , POINTER :: XH_WALL         ! Sensible heat flux on wall     
REAL, DIMENSION(:)  , POINTER :: XLE_WALL        ! Latent heat flux on wall            
REAL, DIMENSION(:)  , POINTER :: XGFLUX_WALL     ! Storage heat flux on wall            
REAL, DIMENSION(:)  , POINTER :: XRNSNOW_ROOF    ! Net radiation over snow            
REAL, DIMENSION(:)  , POINTER :: XHSNOW_ROOF     ! Sensible heat flux over snow  
REAL, DIMENSION(:)  , POINTER :: XLESNOW_ROOF    ! Latent heat flux over snow   
REAL, DIMENSION(:)  , POINTER :: XGSNOW_ROOF     ! Flux under snow             
REAL, DIMENSION(:)  , POINTER :: XMELT_ROOF      ! Snow melt             
REAL, DIMENSION(:)  , POINTER :: XRNSNOW_ROAD    ! Net radiation over snow            
REAL, DIMENSION(:)  , POINTER :: XHSNOW_ROAD     ! Sensible heat flux over snow  
REAL, DIMENSION(:)  , POINTER :: XLESNOW_ROAD    ! Latent heat flux over snow     
REAL, DIMENSION(:)  , POINTER :: XGSNOW_ROAD     ! Flux under snow        
REAL, DIMENSION(:)  , POINTER :: XMELT_ROAD      ! Snow melt 
REAL, DIMENSION(:)  , POINTER :: XRN             ! Net radiation over town 
REAL, DIMENSION(:)  , POINTER :: XH              ! Sensible heat flux over town
REAL, DIMENSION(:)  , POINTER :: XLE             ! Latent heat flux over town     
REAL, DIMENSION(:)  , POINTER :: XGFLUX          ! Storage heat flux over town          
REAL, DIMENSION(:)  , POINTER :: XEVAP           ! Evaporation                    
REAL, DIMENSION(:)  , POINTER :: XRUNOFF         ! Runoff over ground   
REAL, DIMENSION(:)  , POINTER :: XCH             ! Heat drag             
REAL, DIMENSION(:)  , POINTER :: XRI             ! Richardson number               
REAL, DIMENSION(:)  , POINTER :: XUSTAR          ! Friction velocity              
!
REAL, DIMENSION(:) , POINTER    :: XTRAD_IN       ! body MRT inside building (K)
REAL, DIMENSION(:) , POINTER    :: XTRAD_SUN      ! body MRT in the exposed street (K)
REAL, DIMENSION(:) , POINTER    :: XTRAD_SHADE    ! body MRT in the shaded street (K)
REAL, DIMENSION(:) , POINTER    :: XTRAD_RFSUN      ! body MRT on the exposed roof (K)
REAL, DIMENSION(:) , POINTER    :: XTRAD_RFSHADE    ! body MRT on the shaded roof (K)
REAL, DIMENSION(:) , POINTER    :: XUTCI_IN       ! UTCI inside building
REAL, DIMENSION(:) , POINTER    :: XUTCI_OUTSUN   ! UTCI in the exposed street
REAL, DIMENSION(:) , POINTER    :: XUTCI_OUTSHADE ! UTCI in the shaded street
REAL, DIMENSION(:) , POINTER    :: XUTCI_RFSUN    ! UTCI on the exposed roof
REAL, DIMENSION(:) , POINTER    :: XUTCI_RFSHADE  ! UTCI on the shaded roof
REAL, DIMENSION(:) , POINTER    :: XWBGT_OUTSUN    ! WBGT  wet bulb globe temperature in the street (C)
REAL, DIMENSION(:) , POINTER    :: XWBGT_OUTSHADE  ! WBGT  wet bulb globe temperature in the street (C)
REAL, DIMENSION(:) , POINTER    :: XWBGT_RFSUN     ! WBGT  wet bulb globe temperature on the roof (C)
REAL, DIMENSION(:) , POINTER    :: XWBGT_RFSHADE   ! WBGT  wet bulb globe temperature on the roof (C)
REAL, DIMENSION(:) , POINTER    :: XUTCIC_IN       ! UTCI inside building
REAL, DIMENSION(:) , POINTER    :: XUTCIC_OUTSUN   ! UTCI in the exposed street
REAL, DIMENSION(:) , POINTER    :: XUTCIC_OUTSHADE ! UTCI in the shaded street
REAL, DIMENSION(:) , POINTER    :: XUTCIC_RFSUN   ! UTCI in the exposed street
REAL, DIMENSION(:) , POINTER    :: XUTCIC_RFSHADE ! UTCI in the shaded street
!
REAL, DIMENSION(:) , POINTER    :: XTGLOBE_SUN      ! Globe Temperature in the exposed street (K)
REAL, DIMENSION(:) , POINTER    :: XTGLOBE_SHADE    ! Globe Temperature in the shaded street (K)
REAL, DIMENSION(:) , POINTER    :: XTGLOBE_RFSUN    ! Globe Temperature  on the exposed roof (K)
REAL, DIMENSION(:) , POINTER    :: XTGLOBE_RFSHADE  ! Globe Temperature  on the shaded roof  (K)
REAL, DIMENSION(:) , POINTER    :: XTWETB           ! wet-bulb temperature on the ground (K)
REAL, DIMENSION(:) , POINTER    :: XTWETB_ROOF      ! wet-bulb tempertaure over the roof  (K)
!
REAL, DIMENSION(:) , POINTER    :: XTRFZT  ! 
REAL, DIMENSION(:) , POINTER    :: XTRDZT  ! 
REAL, DIMENSION(:) , POINTER    :: XURDZU  ! 
!
REAL, DIMENSION(:) , POINTER    :: XQ1  ! 
REAL, DIMENSION(:) , POINTER    :: XQ2  ! 
REAL, DIMENSION(:) , POINTER    :: XQ3  ! 
REAL, DIMENSION(:) , POINTER    :: XQ4  ! 
REAL, DIMENSION(:) , POINTER    :: XQ5  ! 
REAL, DIMENSION(:) , POINTER    :: XQ6  ! 
REAL, DIMENSION(:) , POINTER    :: XQ7  ! 
REAL, DIMENSION(:) , POINTER    :: XQ8  ! 
REAL, DIMENSION(:) , POINTER    :: XQ9  ! 
REAL, DIMENSION(:) , POINTER    :: XQ10  ! 
REAL, DIMENSION(:) , POINTER    :: XQ11  ! 
REAL, DIMENSION(:) , POINTER    :: XQ12  ! 
REAL, DIMENSION(:) , POINTER    :: XQ13  ! 

!
END MODULE MODD_TOWN
