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
    MODULE MODI_URBAN_HYDRO
!
!
!
INTERFACE
!
!
    SUBROUTINE URBAN_HYDRO(PWS_ROOF_MAX,PWS_ROAD_MAX, PWS_ROOF, PWS_ROAD,  &
                           PRR, PTSTEP, PBLD, PLE_ROOF, PLE_ROAD,          &
                           PRUNOFF_ROOF,                                   &
                           PRUNOFF_ROAD,                                   &
                           PRUNOFF_TOWN                                    )
!
!*      0.1    declarations of arguments 
!
!
REAL, DIMENSION(:), INTENT(IN) :: PWS_ROOF_MAX    ! maximum deepness of roof water reservoir
REAL, DIMENSION(:), INTENT(IN) :: PWS_ROAD_MAX    ! maximum deepness of road water reservoir

REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROOF     ! roof water reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROAD     ! road water reservoir
REAL, DIMENSION(:), INTENT(IN)    :: PRR          ! rain rate
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PBLD         ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROAD     ! latent heat flux over road
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROOF ! runoff (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROAD ! runoff (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_TOWN ! runoff (kg/m2/s)
!
END SUBROUTINE URBAN_HYDRO
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_URBAN_HYDRO
