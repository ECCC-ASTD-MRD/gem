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
    MODULE MODI_BLD_E_BUDGET
!
!
!
INTERFACE
!
      SUBROUTINE BLD_E_BUDGET(OTI_EVOL, PTSTEP, PBLD, PWALL_O_HOR,           &
                              PRHOA, PT_ROOF, PT_WALL, PTI_BLD, PAC_BLD      )
!
!
LOGICAL,              INTENT(IN)   :: OTI_EVOL      ! true --> internal temp. of
!                                                   !      of buildings evolves
!                                                   ! false--> it is fixed
REAL,                 INTENT(IN)   :: PTSTEP        ! time step
REAL, DIMENSION(:),   INTENT(IN)   :: PBLD          ! building fraction
REAL, DIMENSION(:),   INTENT(IN)   :: PWALL_O_HOR   ! wall surf. / hor. surf.
REAL, DIMENSION(:),   INTENT(IN)   :: PRHOA         ! air density
                                                    ! at the lowest level
REAL, DIMENSION(:,:), INTENT(IN)   :: PT_ROOF       ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(IN)   :: PT_WALL       ! wall layers temperatures
REAL, DIMENSION(:),   INTENT(INOUT):: PTI_BLD       ! building air temperature
REAL, DIMENSION(:),   INTENT(OUT)  :: PAC_BLD       ! aerodynamical conductance
                                                    ! inside building itself
!
END SUBROUTINE BLD_E_BUDGET
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_BLD_E_BUDGET
