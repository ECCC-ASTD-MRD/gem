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
!     ######spl
      MODULE MODE_THERMOS
!     ####################
!
!!****  *MODE_THERMO* -
!!
!!    PURPOSE
!!    -------
!      
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
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/08/94 
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!-------------------------------------------------------------------------------
!
INTERFACE QSAT
  MODULE PROCEDURE QSATW_1D
END INTERFACE
INTERFACE DQSAT
  MODULE PROCEDURE DQSATW_O_DT_1D
END INTERFACE
INTERFACE QSATI
  MODULE PROCEDURE QSATI_1D
END INTERFACE
INTERFACE DQSATI
  MODULE PROCEDURE DQSATI_O_DT_1D
END INTERFACE

CONTAINS
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATW_1D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS
!
implicit none
!!!#include <arch_specific.hf>
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature
                                                        ! (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PP     ! Pressure
                                                        ! (Pa)
REAL, DIMENSION(SIZE(PT,1))                   :: PQSAT  ! saturation vapor 
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1))                   :: ZFOES  ! saturation vapor 
                                                        ! pressure
                                                        ! (Pascal) 
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
ZFOES(:) = EXP( XALPW - XBETAW/PT(:) - XGAMW*LOG(PT(:))  )
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT(:) = XRD/XRV*ZFOES(:)/PP(:)   &
                   / (1.+(XRD/XRV-1.)*ZFOES(:)/PP(:))
!-------------------------------------------------------------------------------
!
END FUNCTION QSATW_1D
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!     ##############################################################
      FUNCTION DQSATW_O_DT_1D(PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS       
!
implicit none
!!!#include <arch_specific.hf>
!
!*       0.1   Declarations of arguments and results
!
!
REAL,    DIMENSION(:), INTENT(IN)             :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:), INTENT(IN)               :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:), INTENT(IN)               :: PQSAT  ! saturation vapor 
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT))                    :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor 
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))                       :: ZFOES  ! saturation vapor 
                                                          ! pressure
                                                          ! (Pascal) 
!
!-------------------------------------------------------------------------------
!
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
ZFOES(:) = PP(:) / (1.+XRD/XRV*(1./PQSAT(:)-1.))
!
!*       2.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
PDQSAT(:) = PQSAT(:) / (1.+(XRD/XRV-1.)*ZFOES(:)/PP(:) ) &
                   * (XBETAW/PT(:)**2 - XGAMW/PT(:))
!
!-------------------------------------------------------------------------------
!
END FUNCTION DQSATW_O_DT_1D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################################
      FUNCTION DQSATI_O_DT_1D(PT,PP,PQSAT) RESULT(PDQSAT)
!     ##############################################################
!
!!****  *QSATW * - function to compute saturation vapor humidity from
!!                 temperature (with respect to ice)
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPW) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAW) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMW) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!
!!      Finally, dqsat / dT  (T) is computed.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPW   : Constant for saturation vapor pressure function
!!        XBETAW  : Constant for saturation vapor pressure function
!!        XGAMW   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS       
!
implicit none
!!!#include <arch_specific.hf>
!
!*       0.1   Declarations of arguments and results
!
!
REAL,    DIMENSION(:), INTENT(IN)               :: PT     ! Temperature
                                                          ! (Kelvin)
REAL,    DIMENSION(:), INTENT(IN)               :: PP     ! Pressure
                                                          ! (Pa)
REAL,    DIMENSION(:), INTENT(IN)               :: PQSAT  ! saturation vapor 
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
REAL,    DIMENSION(SIZE(PT))                    :: PDQSAT ! derivative according
                                                          ! to temperature of
                                                          ! saturation vapor 
                                                          ! specific humidity
                                                          ! with respect to
                                                          ! water (kg/kg))
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT))                       :: ZFOES  ! saturation vapor 
                                                          ! pressure
                                                          ! (Pascal) 
!
!-------------------------------------------------------------------------------
!
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
ZFOES(:) = PP(:) / (1.+XRD/XRV*(1./PQSAT(:)-1.))
!
!*       3.    DERIVATION ACCORDING TO TEMPERATURE
!              -----------------------------------
!
PDQSAT(:) = PQSAT(:) / (1.+(XRD/XRV-1.)*ZFOES(:)/PP(:) ) &
                   * (XBETAI/PT(:)**2 - XGAMI/PT(:))
!
!-------------------------------------------------------------------------------
!
END FUNCTION DQSATI_O_DT_1D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ######################################
      FUNCTION QSATI_1D(PT,PP) RESULT(PQSAT)
!     ######################################
!
!!****  *QSATI * - function to compute saturation vapor humidity from
!!                 temperature
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the saturation vapor 
!     pressure from temperature 
!      
!
!!**  METHOD
!!    ------
!!       Given temperature T (PT), the saturation vapor pressure es(T)
!!    (FOES(PT)) is computed by integration of the Clapeyron equation
!!    from the triple point temperature Tt (XTT) and the saturation vapor 
!!    pressure of the triple point es(Tt) (XESTT), i.e  
!!     
!!         es(T)= EXP( alphaw - betaw /T - gammaw Log(T) )
!!  
!!     with :
!!       alphaw (XALPI) = LOG(es(Tt))+ betaw/Tt + gammaw Log(Tt) 
!!       betaw (XBETAI) = Lv(Tt)/Rv + gammaw Tt
!!       gammaw (XGAMI) = (Cl -Cpv) /Rv
!!
!!      Then, the specific humidity at saturation is deduced.
!!  
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : comtains physical constants
!!        XALPI   : Constant for saturation vapor pressure function
!!        XBETAI  : Constant for saturation vapor pressure function
!!        XGAMI   : Constant for saturation vapor pressure function  
!!      
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/09/98 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS       
!
implicit none
!!!#include <arch_specific.hf>
!
!*       0.1   Declarations of arguments and results
!
!
REAL, DIMENSION(:), INTENT(IN)                :: PT     ! Temperature
                                                        ! (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PP     ! Pressure
                                                        ! (Pa)
REAL, DIMENSION(SIZE(PT,1))                   :: PQSAT  ! saturation vapor 
                                                        ! specific humidity
                                                        ! with respect to
                                                        ! water (kg/kg)
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PT,1))                   :: ZFOES  ! saturation vapor 
                                                        ! pressure
                                                        ! (Pascal) 
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE SATURATION VAPOR PRESSURE
!              ---------------------------------
!
ZFOES(:) = EXP( XALPI - XBETAI/PT(:) - XGAMI*LOG(PT(:))  )
!
!*       2.    COMPUTE SATURATION HUMIDITY
!              ---------------------------
!
PQSAT(:) = XRD/XRV*ZFOES(:)/PP(:)   &
                   / (1.+(XRD/XRV-1.)*ZFOES(:)/PP(:))
!-------------------------------------------------------------------------------
!
END FUNCTION QSATI_1D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END MODULE MODE_THERMOS
