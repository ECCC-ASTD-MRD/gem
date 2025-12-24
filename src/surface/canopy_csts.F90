!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------
MODULE CANOPY_CSTS 


!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the physical constants
!     related to canopy energy and mass balance, intercepted snow
!      
implicit none
!
!       ------------------------------------------------------
!               Used in canopy_met_svs2 and drag_svs2
!       ------------------------------------------------------
!
! Choice of the method to calculate the wind in the forest
! 'VDENS_WCAN' uses the high canopy density (VDENS) in the calculation of WCAN
! 'WEIGHT_AVG' does a weighted average of open and closed forest wind speeds based on VDENS (Mazzotti et al. 2021)
character(10) :: lwind_forest = 'VDENS_WCAN' ! 'VDENS_WCAN' or 'WEIGHT_AVG' 
REAL, PARAMETER :: KEXT = 0.5    ! Vegetation light extinction coefficient

REAL, PARAMETER :: RCHD = 0.67    ! Ratio of displacement height to canopy height
REAL, PARAMETER :: HSUBCANO = 2.0 ! Sub canopy reference height for wind, tair and hu
REAL, PARAMETER :: ZBETA = 0.9 ! Constant used in the canopy wind decay coefficient WCAN ( Marke et al., 2016; Liston and Elder, 2006)

!       ------------------------------------------------------
!                Used in vegi_svs2
!             From Gouttevin et al. (2015)
!!       ------------------------------------------------------


REAL, PARAMETER :: E_LEAF = 0.001  ! Leaf thickness (m)
REAL, PARAMETER :: RHO_BIO = 900. ! Biomass density (kg m-3)
REAL, PARAMETER :: CP_BIO = 2800.   ! Biomass specific heat mass (J kg-1 K-1) 
REAL, PARAMETER :: ZB = 40. / 10000.  ! Basal tree area (m2 m-2)
REAL, PARAMETER :: ZCVAI = 3.6e4  ! Vegetation heat capacity per unit VAI (J/K/m^2)
REAL, PARAMETER ::  CVAI = 4.4 ! Intercepted snow capacity per unit lai (kg/m^2) from Essery et a. (2003)
          ! Note that HP98 are using a temperature-dependant param for falling snow density that could be considered later
          ! The maximum canopy snow interception load in HP98 also depends on the type of trees.
REAL, PARAMETER ::  M_SCAP = 5. ! Intercepted snow coefficient per unit lai (kg/m^2) from Andreadis et a. (2009)
CHARACTER(LEN=4) :: HCANO_HM = 'E03'  ! Select the approach used to compute the mass loss due sublimation of intercepted snow
              ! 'E03': using Essery et al. (2003) approach. Also in FSM2
              ! 'G15': formulation proposed by Gouttevin et al. (2015)

!       ------------------------------------------------------
!                   Used in snow_interception_svs2
!       ------------------------------------------------------

REAL, PARAMETER ::     ALPHA = 5   ! Shape parameter used when computing the averahe sublimation rate
REAL, PARAMETER ::     CLUMPING = 0.5 ! Clumping parameter (switch from LAI to effective LAI) 
                                      ! to account for spatial distribution of patter of leaves. (Chen et al. 1997, Essery et al. 2008)
REAL, PARAMETER ::    ZVENT = 0.75 ! Ratio between ventilation wind speed height and tree height [-]
REAL, PARAMETER :: RADIUS_ICESPH = 5e-4 ! Radius of single 'ideal' ice shpere [m]
REAL, PARAMETER :: ALBEDO_ICESPH = 0.8  ! Albedo of single 'ideal' ice shpere [-]
REAL, PARAMETER :: KS = 0.0114          ! Snow shape coefficient for jack pine. Taken from Pomeroy et al. (1998)
REAL, PARAMETER :: FRACT = 0.37 ! Fractal dimension of intercepted snow [-]
REAL, SAVE :: MWAT = 0.0180153 ! Molecular weight of water [kg mol-1]
REAL, SAVE :: RGAZ = 8.314  !  Universal gas constant [J mol-1 K-1]
REAL, SAVE :: CICE = 2.102e3 !  Heat capacity of ice [J kg-1 K-1]
REAL, PARAMETER  ::  TCNC = 240*3600. ! Canopy unloading time scale for cold snow (s)
REAL, PARAMETER ::  TCNM = 48*3600. ! Canopy unloading time scale for melting snow (s)

!       ------------------------------------------------------
!                             Used in drag_svs2
!       ------------------------------------------------------
LOGICAL :: lres_snca = .false. ! Use or not the resistance from the intercepted snow (from Essery et al. 2003)

!       ------------------------------------------------------
!                               Used in ebudget
!       ------------------------------------------------------
REAL, PARAMETER ::  EMSNV = 1. ! Emissivity of snow below high vegetation                    
REAL, PARAMETER ::  ABARK  = 0.15  !   Albedo of Bark (S. Wang, Ecological Modelling, 2005)
REAL, PARAMETER ::  ALSNV = 0.3 ! Snow-covered canopy albedo (Goottevin et al. 2015)

!       ------------------------------------------------------
!                  Function used in drag_svs2 and canopy_met_svs2
!       ------------------------------------------------------
CONTAINS 

REAL FUNCTION WIND_CANO_COEF(LAIVH,VGH_DENS)
    !    
    ! Compute wind speed attenuation coefficient in the canopy 
    ! From Marke et al., (2016); Liston and Elder (2006)
    !          - Input -
    !   LAIVH: LAI of trees in areas of high vegeation (-)
    !   VGH_DENS: density of trees in areas of high vegeation (-)
    !          - Output -
    !   WIND_CANO_COEF: wind speed attenuation coefficient 

    REAL, INTENT(IN) :: LAIVH,VGH_DENS

    ! Use a minimum value of 1 to ensure a proper reduction of wind speed in the canopy
    ! Set a maximum value of 4 according to Mahat et al (2013), Bonan(1991)
    WIND_CANO_COEF =  MIN(4.,MAX(1.0,ZBETA * LAIVH * VGH_DENS))

    RETURN
END FUNCTION WIND_CANO_COEF

END MODULE CANOPY_CSTS 
