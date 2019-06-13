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
    SUBROUTINE SURFACE_AERO_COND(PRI, PZREF, PUREF, PVMOD, PZ0,&
                                    PZ0H, PAC                     )
!   ######################################################################
!
!!****  *SURFACE_AERO_COND*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the drag coefficients for heat and momentum near the ground
!         
!     
!!**  METHOD
!!    ------
!
!
!
!    1 and 2 : computation of relative humidity near the ground
!
!    3 : richardson number
!
!    4 : the aerodynamical resistance for heat transfers is deduced
!
!    5 : the drag coefficient for momentum ZCD is computed
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
!!                  02/04/01 (P Jabouille) limitation of Z0 with 0.5 PUREF
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XKARMAN
!
USE MODE_THERMOS
!
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN)    :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the wind
                                              ! NOTE this is different from ZZREF
                                              ! ONLY in stand-alone/forced mode,
                                              ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:), INTENT(IN)    :: PZ0      ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H     ! roughness length for heat
!
REAL, DIMENSION(:), INTENT(OUT)   :: PAC      ! aerodynamical conductance
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PRI)) :: ZZ0, ZZ0H, ZMU,          &
                              ZFH, ZCHSTAR, ZPH, ZCDN, &
                              ZSTA, ZDI
!-------------------------------------------------------------------------------
!
!*       4.     Surface aerodynamic resistance for heat transfers
!               -------------------------------------------------
!
ZZ0(:) = MIN(PZ0(:),PUREF(:)*0.5)
ZZ0H(:) = MIN(ZZ0(:),PZ0H(:))
!
ZMU(:) = MAX( LOG( ZZ0(:)/ZZ0H(:) ), 0.0 )
ZFH(:) = LOG( PUREF(:)/ZZ0(:) )/ LOG( PZREF(:)/ ZZ0H(:) )
!
ZCHSTAR(:) = CHSTAR(ZMU(:))
ZPH(:)     = PH(ZMU(:))
!
! 
ZCDN(:) = (XKARMAN/LOG(PUREF(:)/ZZ0(:)))**2.
ZSTA(:) = PRI(:)*PVMOD(:)*PVMOD(:)
!
!
WHERE ( PRI < 0.0 )
  ZDI(:)= 1. / ( PVMOD(:)                                           &
                  +ZCHSTAR(:)*ZCDN(:)*15.                           &
                               *(PZREF(:)/ZZ0H(:))**ZPH(:)        &
                               *ZFH(:) * SQRT(-ZSTA(:))             &
                 )
  PAC(:) = ZCDN(:)*(PVMOD(:)-15.* ZSTA(:)*ZDI(:))*ZFH(:)

ELSEWHERE
  ZDI(:) = SQRT(PVMOD(:) * PVMOD(:) + 5. * ZSTA(:) )
  PAC(:) = ZCDN(:)*PVMOD(:)/(1.+15.*ZSTA(:)*ZDI(:)  &
           / PVMOD(:) /PVMOD(:) /PVMOD(:) )*ZFH(:)
END WHERE
!
!-------------------------------------------------------------------------------
!
!
!
CONTAINS
!
!                                              functions used in the calculation
!                                              the terms Ch
  FUNCTION CHSTAR(X)
  REAL, DIMENSION(:), INTENT(IN)     :: X
  REAL, DIMENSION(SIZE(X))           :: CHSTAR
  !
  CHSTAR = 3.2165 + 4.3431*X + 0.5360*X*X - 0.0781*X*X*X
  !
  END FUNCTION CHSTAR
  !
  !
  FUNCTION PH(X)
  REAL, DIMENSION(:), INTENT(IN)     :: X
  REAL, DIMENSION(SIZE(X))           :: PH
  !
  PH = 0.5802 - 0.1571*X + 0.0327*X*X - 0.0026*X*X*X
  !
  END FUNCTION PH                          
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE SURFACE_AERO_COND
