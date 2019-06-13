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
!   ######################################################################

    SUBROUTINE SURFACE_RI(PTG, PQS, PEXNS, PEXNA, PTA, PQA,   &
                             PZREF, PUREF, PDIRCOSZW, PVMOD, PRI )
!   ######################################################################
!
!!****  *SURFACE_RI*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the richardson number near the ground
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
!!    MODD_CST
!!    MODD_GROUND_PAR
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
!!      Original    22/09/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XRV, XRD, XG
!
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQS      ! surface specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS    ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
                                              ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA    ! exner function
                                              ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
!
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
!                                             ! NOTE this is different from ZZREF
!                                             ! ONLY in stand-alone/forced mode,
!                                             ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:), INTENT(IN)    :: PDIRCOSZW! Cosine of the angle between
!                                             ! the normal to the surface and
!                                             ! the vertical
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRI      ! Richardson number
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PTG)) :: ZTHVA, ZTHVS
!-------------------------------------------------------------------------------
!
!       1.     Richardson number
!              -----------------
!                
!                                                 virtual potential        
!                                                 temperature at the 
!                                                 first atmospheric level and
!                                                 at the surface
!
ZTHVA(:)=PTA(:)/PEXNA(:)*( 1.+(XRV/XRD-1.)*PQA(:) )   
ZTHVS(:)=PTG(:)/PEXNS(:)*( 1.+(XRV/XRD-1.)*PQS(:) )
!                                                 
!
                                                ! Richardson's number
PRI(:) = XG * PDIRCOSZW(:) * PUREF(:) * PUREF(:)              &
        * (ZTHVA(:)-ZTHVS(:)) / (0.5 * (ZTHVA(:)+ZTHVS(:)) )  &
        / (PVMOD(:)*PVMOD(:)) /PZREF(:)
!-------------------------------------------------------------------------------
!
END SUBROUTINE SURFACE_RI
