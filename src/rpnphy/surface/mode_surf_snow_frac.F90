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
!     ##########################
      MODULE MODE_SURF_SNOW_FRAC
!     ##########################
!
!!****  *MODE_SURF_SNOW_FRAC* -  module for routines to compute snow fraction
!!                               for surface schemes
!!
!!    PURPOSE
!!    -------
!    
!      The purpose of this routine is to store here all routines to compute
!     snow fractions for the TEB scheme. This allows to insure a coherent
!     way in retrieving snow fraction or snow contents.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE          
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/03/99
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!-------------------------------------------------------------------------------
!

CONTAINS
!-------------------------------------------------------------------------------
!
!     ###############################################
      FUNCTION SNOW_FRAC_GROUND(PWSNOW) RESULT(PPSNG)
!     ###############################################
!
USE MODD_SNOW_PAR, ONLY : XWCRN
implicit none
!!!#include <arch_specific.hf>
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL, DIMENSION(SIZE(PWSNOW))   :: PPSNG  ! snow fraction over bare ground
!
PPSNG(:) = PWSNOW(:) / (PWSNOW(:)+XWCRN)       ! fraction of ground covered
!
END FUNCTION SNOW_FRAC_GROUND
!
!-------------------------------------------------------------------------------
!
!     ##########################################################
      FUNCTION WSNOW_FROM_SNOW_FRAC_GROUND(PPSNG) RESULT(PWSNOW)
!     ##########################################################
!
USE MODD_SNOW_PAR, ONLY : XWCRN
implicit none
!!!#include <arch_specific.hf>
!
REAL, DIMENSION(:), INTENT(IN)  :: PPSNG  ! snow fraction over bare ground
REAL, DIMENSION(SIZE(PPSNG))    :: PWSNOW ! snow amount over natural areas (kg/m2)
!
PWSNOW(:) = XWCRN * PPSNG(:) / (1. - PPSNG(:))
!
END FUNCTION WSNOW_FROM_SNOW_FRAC_GROUND
!-------------------------------------------------------------------------------
!
!     #########################################################
      FUNCTION SNOW_FRAC_VEG(PWSNOW,PZ0VEG,PRHOS) RESULT(PPSNV)
!     #########################################################
!
implicit none
!!!#include <arch_specific.hf>
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PZ0VEG ! vegetation roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)  :: PRHOS  ! snow density (kg/m3)
REAL, DIMENSION(SIZE(PWSNOW))   :: PPSNV  ! snow fraction over vegetation
!
PPSNV(:) = PWSNOW(:) / (PWSNOW(:)+PRHOS(:)*5.*PZ0VEG(:))
!
END FUNCTION SNOW_FRAC_VEG
!
!-------------------------------------------------------------------------------
!
!     ############################################################
      FUNCTION SNOW_FRAC_NAT(PWSNOW,PPSNG,PPSNV,PVEG) RESULT(PPSN)
!     ############################################################
!
implicit none
!!!#include <arch_specific.hf>
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW ! snow amount over natural areas (kg/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PPSNG  ! snow fraction over bare ground
REAL, DIMENSION(:), INTENT(IN)  :: PPSNV  ! snow fraction over vegetation
REAL, DIMENSION(:), INTENT(IN)  :: PVEG   ! vegetation fraction
REAL, DIMENSION(SIZE(PWSNOW))   :: PPSN   ! snow fraction over natural areas
!
PPSN(:) = ( 1-PVEG(:) )*PPSNG(:) + PVEG(:)*PPSNV(:)
!
END FUNCTION SNOW_FRAC_NAT
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!     ##############################################################
      SUBROUTINE SNOW_FRAC_ROAD(PWSNOW_ROAD,OSNOW,PDN_ROAD,PDF_ROAD)
!     ##############################################################
!
USE MODD_SNOW_PAR, ONLY : XWCRN
!!!#include <arch_specific.hf>
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW_ROAD ! snow amount over roads (kg/m2) 
LOGICAL, DIMENSION(:), INTENT(IN)  :: OSNOW    ! T: snow-fall is occuring
REAL, DIMENSION(:), INTENT(OUT) :: PDN_ROAD    ! snow fraction over roads
REAL, DIMENSION(:), INTENT(OUT) :: PDF_ROAD    ! snow-free fraction over roads
!
PDF_ROAD(:)     = 1.
PDN_ROAD(:)     = 0.
!
! due to the flatness of horizontal surfaces (compared to landscape and
! vegetation), the amount of snow necessary to cover the entire surface XWCRN
! is reduced (equal to 1kg/m2 instead of 10).
!
WHERE (PWSNOW_ROAD(:)>0. .OR. OSNOW)
  PDN_ROAD(:)     = MAX(MIN(PWSNOW_ROAD(:)/(PWSNOW_ROAD(:) + XWCRN*0.1) , 0.7), 0.01)
  PDF_ROAD(:)     = 1.-PDN_ROAD(:)
END WHERE
!
END SUBROUTINE SNOW_FRAC_ROAD
!
!-------------------------------------------------------------------------------
!
!     ##############################################################
      SUBROUTINE SNOW_FRAC_ROOF(PWSNOW_ROOF,OSNOW,PDN_ROOF,PDF_ROOF)
!     ##############################################################
!
USE MODD_SNOW_PAR, ONLY : XWCRN
!!!#include <arch_specific.hf>
!
REAL, DIMENSION(:), INTENT(IN)  :: PWSNOW_ROOF ! snow amount over roofs (kg/m2) 
LOGICAL, DIMENSION(:), INTENT(IN)  :: OSNOW    ! T: snow-fall is occuring
REAL, DIMENSION(:), INTENT(OUT) :: PDN_ROOF    ! snow fraction over roofs
REAL, DIMENSION(:), INTENT(OUT) :: PDF_ROOF    ! snow-free fraction over roofs
!
PDF_ROOF(:)     = 1.
PDN_ROOF(:)     = 0.
!
! due to the flatness of horizontal surfaces (compared to landscape and
! vegetation), the amount of snow necessary to cover the entire surface XWCRN
! is reduced (equal to 1kg/m2 instead of 10).
!
WHERE (PWSNOW_ROOF(:)>0. .OR. OSNOW)
  PDN_ROOF(:)     = MAX(PWSNOW_ROOF(:)/(PWSNOW_ROOF(:) + XWCRN*0.1),0.01)
  PDF_ROOF(:)     = 1.-PDN_ROOF(:)
END WHERE
!
END SUBROUTINE SNOW_FRAC_ROOF
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
END MODULE MODE_SURF_SNOW_FRAC
