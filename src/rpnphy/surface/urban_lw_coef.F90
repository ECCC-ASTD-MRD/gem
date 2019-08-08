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
    SUBROUTINE URBAN_LW_COEF(PEMIS_ROAD, PSVF_ROAD, PEMIS_WALL, PSVF_WALL,   &
                             PDN, PDF, PESNOW_ROAD,                          &
                             PLW_W_TO_W, PLW_R_TO_W, PLW_W_TO_R, PLW_R_TO_R, &
                             PLW_S_TO_W, PLW_S_TO_R, PLW_N_TO_W, PLW_N_TO_R  )
!   ##########################################################################
!
!!****  *URBAN_LW_COEF*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the coefficients before each of the temperatures in the
!     radiative budgets
!         
!     
!!**  METHOD
!     ------
!
! without snow, the radiative budgets read:
!
!   Rn_w = abs_Rg_w 
!  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
!  +         emis_w                       *      SVF_w                * LWR
!  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
!  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
!  +         emis_w            (1-emis_r) *      SVF_r  *      SVF_w  * LWR
!  +         emis_w            (1-emis_w) *      SVF_w  * (1-2*SVF_w) * LWR
!  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
!
!   Rn_r = abs_Rg_r
!  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
!  +         emis_r                       *    SVF_r                  * LWR
!  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
!  +         emis_r            (1-emis_w) * (1-SVF_r)   *      SVF_w  * LWR
!  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
!
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
!!      Original    08/09/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS, ONLY : XSTEFAN
!
implicit none
!!!#include <arch_specific.hf>
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
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PDN)) :: ZAVEMIS_R   ! averaged road emissivity
!-------------------------------------------------------------------------------
!
!
ZAVEMIS_R(:) =     PDF(:)  * PEMIS_ROAD (:) &
              +    PDN(:)  * PESNOW_ROAD(:)
!
!
PLW_W_TO_W(:) =     XSTEFAN * PEMIS_WALL(:)                      &
               * ( -   1.                                        &
                   +       PEMIS_WALL(:)  * (1.-2.*PSVF_WALL(:)) &
                   +       PEMIS_WALL(:)  * (1.-   PSVF_ROAD(:)) &
                     * (1.-ZAVEMIS_R (:)) *        PSVF_WALL(:)  &
                   +       PEMIS_WALL(:)  * (1.-2.*PSVF_WALL(:)) &
                     * (1.-PEMIS_WALL(:)) * (1.-2.*PSVF_WALL(:)) &
                 )
!
PLW_R_TO_W(:) =   XSTEFAN * PEMIS_ROAD(:)                          &
                          * PEMIS_WALL(:)                          &
                          * PSVF_WALL(:)                           &
              * ( 1. + (1.-PEMIS_WALL(:)) * (1.-2.*PSVF_WALL(:)))  &
                          * PDF(:)
!
PLW_S_TO_W(:) =  PEMIS_WALL(:) * PSVF_WALL(:)                    &
               * (    1.                                         &
                   + (1.-ZAVEMIS_R (:)) *        PSVF_ROAD(:)    &
                   + (1.-PEMIS_WALL(:)) * (1.-2.*PSVF_WALL(:))   &
                 )
!
PLW_N_TO_W(:) =   XSTEFAN * PESNOW_ROAD(:)                         &
                          * PEMIS_WALL(:)                          &
                          * PSVF_WALL(:)                           &
              * ( 1. + (1.-PEMIS_WALL(:)) * (1.-2.*PSVF_WALL(:)))  &
                          * PDN(:)
!
PLW_R_TO_R(:) =            XSTEFAN * PEMIS_ROAD(:)                  &
               * (  -   1.                                          &
                    + PDF(:) * PEMIS_ROAD(:) * (1.-PEMIS_WALL(:))   &
                               * (1.-PSVF_ROAD(:)) * PSVF_WALL(:)   &
                 )
!
PLW_W_TO_R(:) =   XSTEFAN * PEMIS_ROAD(:)                           &
                          * PEMIS_WALL(:)                           &
                          * (1. -  PSVF_ROAD(:))                    &
              * ( 1. + (1.-PEMIS_WALL(:)) * (1.-2.*PSVF_WALL(:)) )
!
PLW_S_TO_R(:) =  PEMIS_ROAD(:)* PSVF_ROAD(:)                      &
               + PEMIS_ROAD(:)*(1.-PEMIS_WALL(:))                 &
                                *(1.-PSVF_ROAD(:))*PSVF_WALL(:)
!
PLW_N_TO_R(:) = XSTEFAN * PEMIS_ROAD(:)     * (1.-PEMIS_WALL(:)) &
                        * (1.-PSVF_ROAD(:)) * PSVF_WALL(:)       &
                        * PESNOW_ROAD(:)    * PDN(:)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_LW_COEF
