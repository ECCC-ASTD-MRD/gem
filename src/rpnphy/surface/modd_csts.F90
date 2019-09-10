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
!     ###############
      MODULE MODD_CSTS      
!     ###############
!
!!****  *MODD_CSTS* - declaration of Physic constants 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the 
!     Physics constants.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!      V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    16/05/94  
!!      J. Stein    02/01/95  add xrholw                    
!!      J.-P. Pinty 13/12/95  add XALPI,XBETAI,XGAMI
!!      J. Stein    25/07/97  add XTH00                    
!!      V. Masson   05/10/98  add XRHOLI
!!      C. Mari     31/10/00  add NDAYSEC
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
implicit none
REAL,SAVE :: XPI                ! Pi
!
REAL,SAVE :: XDAY,XSIYEA,XSIDAY ! day duration, sideral year duration,
                                ! sideral day duration
!
REAL,SAVE :: XKARMAN            ! von karman constant
REAL,SAVE :: XLIGHTSPEED        ! light speed
REAL,SAVE :: XPLANCK            ! Planck constant
REAL,SAVE :: XBOLTZ             ! Boltzman constant 
REAL,SAVE :: XAVOGADRO          ! Avogadro number
!
REAL,SAVE :: XRADIUS,XOMEGA     ! Earth radius, earth rotation
REAL,SAVE :: XG                 ! Gravity constant
!
REAL,SAVE :: XP00               ! Reference pressure
!
REAL,SAVE :: XSTEFAN,XI0        ! Stefan-Boltzman constant, solar constant
!
REAL,SAVE :: XMD,XMV            ! Molar mass of dry air and molar mass of vapor
REAL,SAVE :: XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
REAL,SAVE :: XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
REAL,SAVE :: XRHOLW             ! Volumic mass of liquid water
REAL,SAVE :: XCL,XCI            ! Cl (liquid), Ci (ice)
REAL,SAVE :: XTT                ! Triple point temperature
REAL,SAVE :: XLVTT              ! Vaporization heat constant
REAL,SAVE :: XLSTT              ! Sublimation heat constant
REAL,SAVE :: XLMTT              ! Melting heat constant
REAL,SAVE :: XESTT              ! Saturation vapor pressure  at triple point
                                ! temperature  
REAL,SAVE :: XALPW,XBETAW,XGAMW ! Constants for saturation vapor 
                                !  pressure  function 
REAL,SAVE :: XALPI,XBETAI,XGAMI ! Constants for saturation vapor
                                !  pressure  function over solid ice
REAL, SAVE        :: XTH00      ! reference value  for the potential
                                ! temperature
REAL,SAVE :: XRHOLI             ! Volumic mass of liquid water
REAL,SAVE :: XCONDI             ! thermal conductivity of ice (W m-1 K-1)
!
INTEGER, SAVE :: NDAYSEC        ! Number of seconds in a day
!
END MODULE MODD_CSTS

