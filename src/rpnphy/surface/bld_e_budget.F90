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
    SUBROUTINE BLD_E_BUDGET(OTI_EVOL, PTSTEP, PBLD, PWALL_O_HOR,           &
                            PRHOA, PT_ROOF, PT_WALL, PTI_BLD, PAC_BLD      )
!   ##########################################################################
!
!!****  *BLD_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of the temperature of inside building air
!         
!     
!!**  METHOD
!     ------
!
!     The resistance term between the surfaces and the room is given
!     by a standard value, which mimics both the convection
!     and the radiative interactions in the room.
!     This explains the very low resistance. It is used to compute
!     the evolution of the surfaces only.
!     This resistance value is 0.123 Km/W  (typical for inside surfaces).
!     (ENVIRONMENTAL SCIENCE IN BUILDING, 3rd Edition, Randall McMullan,
!      THE MACMILLAN PRESS Limited).
!
!
!
!     On the contrary, the evolution of the air temperature is mainly
!     governed by the convection (considering the low radiative absorption
!     of the air itself).
!     In order to have a simple formulation, a diurnal cycle is assumed,
!     with a force restore formulation.
!
!     The floor temperature is fixed
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
!!      Original    24/08/00 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XTT, XCPD, XDAY
!
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
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
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTI_BLD)) :: ZT_FLOOR     ! floor temperature
!
REAL                           :: ZTAU         ! temporal filter period
!
INTEGER                        :: IROOF        ! number of roof layers
INTEGER                        :: IWALL        ! number of wall layers
!-------------------------------------------------------------------------------
!
!*      1.   initializations
!            ---------------
!
IROOF = SIZE(PT_ROOF,2)
IWALL = SIZE(PT_WALL,2)
!
ZT_FLOOR (:) = 19. + XTT
!
!
!*      2.   inside conductance FOR SURFACES
!            ------------------
!
!* (normalized by rho Cp for convenience)
!
PAC_BLD(:) = 1. / 0.123 / (XCPD * PRHOA(:))
!
!*      3.   no evolution of interior temperature if OTI_EVOL=.FALSE.
!            --------------------------------------------------------
!
IF (.NOT. OTI_EVOL) RETURN
!
!*      4.   evolution of the internal temperature
!            -------------------------------------
!
ZTAU=XDAY
!
PTI_BLD(:) = PTI_BLD(:) * (ZTAU-PTSTEP)/ZTAU                   &
            + (   PT_ROOF(:,IROOF) * PBLD       (:)            &
                + PT_WALL(:,IWALL) * PWALL_O_HOR(:)            &
                + ZT_FLOOR(:)      * PBLD       (:)         )  &
             /(  2. * PBLD(:)  +  PWALL_O_HOR(:) ) * PTSTEP / ZTAU
!
!*      5.   internal temperature set to a minimum value (heating)
!            -----------------------------------------------------
! 
PTI_BLD(:) = MAX( PTI_BLD(:) , 19. + XTT )
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE BLD_E_BUDGET
