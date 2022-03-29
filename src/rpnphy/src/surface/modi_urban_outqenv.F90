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

!#TODO: transform the s/r into a module instead of making an interface module to avoid non matching interface

    MODULE MODI_URBAN_OUTQENV
!
!
!
INTERFACE
!
!
SUBROUTINE URBAN_OUTQENV(PSCA_SW, PREF_SW_FAC, PREF_SW_GRND, PREF_SW_ROOF,  &
                       PEMIT_LW_FAC, PEMIT_LW_GRND, PEMIT_LW_ROOF, PLW_RAD, &
                       PBLD, PBLD_HEIGHT, PWALL_O_HOR, PDIR_SW, PZENITH,    &
                       PQ1,PQ2,PQ3,PQ4,PQ5,PQ6,PQ7,                         &
                       PQ8,PQ9,PQ10,PQ11,PQ12,PQ13,OPT,OPT_BODY,            &
                       ZEB,ZAB,ZHB                                         )
!
!
!*      0.1    declarations of arguments
INTEGER, INTENT(IN)  :: OPT
INTEGER, INTENT(IN)  :: OPT_BODY
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_GRND !Solar radiation reflected by ground [road + garden] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_FAC !Solar radiation reflected by facade [wall + glazing] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PSCA_SW !Diffuse solar radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PDIR_SW !Direct solar radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH !solar zenithal angle (rad from vert.)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_FAC !Longwave radiation emitted by the facade [wall + glazing] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_GRND !Longwave radiation emitted by the ground [road + garden] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD !Atmospheric longwave radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_ROOF !Solar radiation reflected by the roof (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_ROOF !Longwave radiation emitted by the roof (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PBLD !Building surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PBLD_HEIGHT !Building surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PWALL_O_HOR !Building surface fraction
REAL, DIMENSION(:), INTENT(OUT) :: PQ1
REAL, DIMENSION(:), INTENT(OUT) :: PQ2
REAL, DIMENSION(:), INTENT(OUT) :: PQ3
REAL, DIMENSION(:), INTENT(OUT) :: PQ4
REAL, DIMENSION(:), INTENT(OUT) :: PQ5
REAL, DIMENSION(:), INTENT(OUT) :: PQ6
REAL, DIMENSION(:), INTENT(OUT) :: PQ7
REAL, DIMENSION(:), INTENT(OUT) :: PQ8
REAL, DIMENSION(:), INTENT(OUT) :: PQ9
REAL, DIMENSION(:), INTENT(OUT) :: PQ10
REAL, DIMENSION(:), INTENT(OUT) :: PQ11
REAL, DIMENSION(:), INTENT(OUT) :: PQ12
REAL, DIMENSION(:), INTENT(OUT) :: PQ13
REAL :: ZHB  ! average height of the body
REAL :: ZAB  !absorption coef of solar radiation by the  body       
REAL :: ZEB  !emissivity of  body                              
!
END SUBROUTINE URBAN_OUTQENV
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_URBAN_OUTQENV
