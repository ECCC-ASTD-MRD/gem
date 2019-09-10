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
    MODULE MODI_ROOF_LAYER_E_BUDGET
!
!
!
INTERFACE
!
!
    SUBROUTINE ROOF_LAYER_E_BUDGET(PT_ROOF, PQSAT_ROOF,                        &
                             PTA, PQA, PPS,                                    &
                             PLW, PTSTEP,                                      &
                             PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,          &
                             PTI_BLD, PAC_BLD, PDELT_ROOF,                     &
                             PDELTSNOW_ROOF, PGSNOW_ROOF,                      &
                             PRHOA, PAC_ROOF, PAC_ROOF_WAT,                    &
                             PABS_SW_ROOF, PABS_LW_ROOF                        )
!
!
!*      0.1    declarations of arguments 
!
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF      ! roof layers temperatures
REAL, DIMENSION(:), INTENT(INOUT) :: PQSAT_ROOF     ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! specific humidity
                                                    ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PLW            ! atmospheric infrared radiation
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF     ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF       ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF       ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF        ! depth of roof layers
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD        ! inside building temp.
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD        ! aerodynamical resistance
                                                    ! inside building itself
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF     ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PDELTSNOW_ROOF ! roof snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF    ! roof snow conduction
!                                                   ! heat fluxes at mantel
!                                                   ! base
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF   ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROOF   ! absorbed infra-red rad.
!
END SUBROUTINE ROOF_LAYER_E_BUDGET
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_ROOF_LAYER_E_BUDGET
