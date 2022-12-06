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
    MODULE MODI_URBAN_LW_COEF
!
INTERFACE
!
    SUBROUTINE URBAN_LW_COEF(PEMIS_ROAD, PSVF_ROAD, PEMIS_WALL, PSVF_WALL,   &
                             PDN, PDF, PESNOW_ROAD,                          &
                             PLW_W_TO_W, PLW_R_TO_W, PLW_W_TO_R, PLW_R_TO_R, &
                             PLW_S_TO_W, PLW_S_TO_R, PLW_N_TO_W, PLW_N_TO_R  )
!
!*      0.1    declarations of arguments 
!
!
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_ROAD  ! road emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PSVF_ROAD   ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_WALL  ! wall emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PSVF_WALL   ! wall sky view factor
REAL, DIMENSION(:), INTENT(IN)  :: PDN         ! snow-covered surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PDF         ! snow-free surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PESNOW_ROAD ! road snow emissivity
REAL, DIMENSION(:), INTENT(OUT) :: PLW_W_TO_W  ! L.W. interactions
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_W  ! from first Temp.
REAL, DIMENSION(:), INTENT(OUT) :: PLW_W_TO_R  ! on second Temp.
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_R  !
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_W  ! idem. but from
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_R  ! sky rad.
REAL, DIMENSION(:), INTENT(OUT) :: PLW_N_TO_W  ! idem. but from
REAL, DIMENSION(:), INTENT(OUT) :: PLW_N_TO_R  ! snow rad.
!
END SUBROUTINE URBAN_LW_COEF
! 
END INTERFACE
!
END MODULE MODI_URBAN_LW_COEF
