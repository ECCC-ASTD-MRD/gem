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
!     ######################
      MODULE MODD_SNOW_PAR
!     ######################
!
!!****  *MODD_SNOW_PAR* - declaration of parameters related
!!                          to the snow parameterization
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     parameters related to the surface parameterization of snow.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
implicit none
!
!--------------------------------------------------------------------------------
! Snow on the ground:
!--------------------------------------------------------------------------------
!
! Snow emissivity:
!
REAL, PARAMETER       :: XEMISSN = 1.0  ! (-)
!
! Critical value of the equivalent water content
! of the snow reservoir for snow fractional coverage and albedo computations
!
REAL, PARAMETER       :: XWCRN = 10.0   ! (kg m-2)
!
! Roughness length of pure snow surface 
!
REAL, PARAMETER       :: XZ0SN = 0.001  ! (m)
!
! Roughness length for heat of pure snow surface 
!
REAL, PARAMETER       :: XZ0HSN = 0.0001! (m)
!
! Minimum and maximum values of the albedo of snow:
!
REAL, PARAMETER       :: XANSMIN = 0.50 ! (-)
REAL, PARAMETER       :: XANSMAX = 0.85 ! (-)
!
! Snow aging coefficients (albedo and Force-Restore density):
!
REAL, PARAMETER       :: XANS_TODRY    = 0.008     ! (-) 
REAL, PARAMETER       :: XANS_T        = 0.240     ! (-)
!
! Minimum and maximum values of the density of snow 
! for Force-Restore snow option
!
REAL, PARAMETER       :: XRHOSMIN = 100.  ! (kg m-3)
REAL, PARAMETER       :: XRHOSMAX = 750.  ! (kg m-3) ! WAS 300 kg/m3
!
! Minimum snow water equivalent water content in TEB for all calc.
!
REAL, PARAMETER       :: SWE_CRIT=0.0001

!
! Minimum and maximum values of the density of snow 
! for ISBA-ES snow option
!
REAL, PARAMETER       :: XRHOSMIN_ES =  50.  ! (kg m-3)
REAL, PARAMETER       :: XRHOSMAX_ES = 750.  ! (kg m-3)
!
! ISBA-ES Critical snow depth at which snow grid thicknesses constant
!
REAL, PARAMETER                      :: XSNOWCRITD = 0.03  ! (m)
!                                       
! ISBA-ES Minimum total snow depth for model 
!
REAL, PARAMETER                      :: XSNOWDMIN = 0.000001  ! (m)
!                                       
! Maximum Richardson number limit for very stable conditions using the ISBA-ES 'RIL' option
!
REAL, PARAMETER                      :: X_RI_MAX = 0.20
!                                       
! ISBA-ES Maximum snow liquid water holding capacity (fraction by mass) parameters:
!
REAL, PARAMETER                      :: XWSNOWHOLDMAX2   = 0.10  ! (-) 
REAL, PARAMETER                      :: XWSNOWHOLDMAX1   = 0.03  ! (-)
REAL, PARAMETER                      :: XSNOWRHOHOLD     = 200.0 ! (kg/m3)
!
END MODULE MODD_SNOW_PAR












