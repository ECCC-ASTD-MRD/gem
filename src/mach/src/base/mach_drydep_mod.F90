!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet / Project : GEM-MACH
! Fichier / File   : mach_drydep_mod.ftn90
! Creation         : A. Robichaud - Dec 2002
! Description      : Declare, as parameters, the tabled values used for the
!                    dry deposition scheme "mach_gas_drydep_solver"
!
! Extra Info       : Modified by A. Kallaur, M. Moran, P.A. Beaulieu for GEM-MACH, Jan 2008
!                    Modified by M. Moran, Aug. 2014
!                    Modified by Verica Savic-Jovcic per instructions of Paul Makar, Nov 2017
!                                (unoxygenated VOCs, HO2, RO2, MCO3 are not allowed to deposit
!                                 and H* for NO corrected to Wesely's value)
!                   : Jack Chen, Feb 2017
!                      - Deposition parameters (alpha/beta/hstar/fzero/ats) definitions now moved 
!                        to the specific gas-phase mechanism package
!                      - This module now saves the generic 'depo' derived type (to be set runtime
!                        by the specified gas-phase mechanism).
!
!                   : Deji Akingunola, July 2019
!                      - Consolidate all the dry deposition parameters used in GEM-MACH
!                         into this module; specifically move 'pllp, 'pgamma' and 'aest'
!                         for aerosol dry deposition from mach_cam_utils_mod here.
!                      -  And rename the module "mach_drydep_mod".
!
!                    Possible future modifications:
!                       - Short term:  Some meteorological input (fb, LAI) taken directly from GEM
!                       - Medium term: LAI from objective analysis (S. Belair)
!                                      Ground water content influence on stomatal resistance (CALDAS) (S. Belair)
! 
!
!------------------ vegetation for dry deposition scheme  -------------
! land use categories  (luc)
! 1  evergreen needleleaf forest (west coast/candian shield)
! 2  evergreen broadleaf forest (only in wrn mexico)
! 3  deciduous needleleaf forest (none in na)
! 4  deciduous broadleaf forest (nw bc and ern na)
! 5  mixed forest
! 6  grassland (prairies)
! 7  crops, mixed farming (prairies)
! 8  desert (wrn us and canada)
! 9  tundra
! 10 dwarf trees, shrubs with ground cover 
! 11 wet land with plants (mainly in nrn canada)
! 12 ice caps and glaciers (very little in na)
! 13 inland water (over canadian shield)
! 14 ocean
! 15 urban
!
!------------------ seasonal categories  (sc)  ------------------------
! 1  midsummer with lush vegetation
! 2  autumn with cropland before harvest
! 3  later autumn after frost, no snow
! 4  winter, snow on ground and subfreezing
! 5  transitional spring with partially green short annuals
!
!------------------ Species list --------------------------
! The list of gas species undergoing deposition are defined in the gas-phase mechanism
! package
!
module mach_drydep_mod
   implicit none
   save

!  Declaration of various constantes and switch settings
!
!  Name       Description
!==============================================================================
!  prandtl    Prandtl number for air
!  b4         Constant for Jarvis scheme
!  ao         Constant for Jarvis scheme
!  bo         Constant for Jarvis scheme
!  co         Constant for Jarvis scheme
!  dzero      Constant for Jarvis scheme
!  en         Constant for Jarvis scheme
!  isimple    Use (1) or not (0) simple scheme for solar radiation effect on
!             stomata aperture. Default is 0 for more complex scheme
!  inew       (=1)  We are extending last paper of ZHANG et al.
!                   to all gas species for non-stomatal resistance
!             (=0)  We are take non-stomatal resistance from ZHANG for Ozone only.
!  insz       (=0)  Bypass non-stomatal resistance (ZHANG et al. 2002) for all
!                   species
!                 

   real(kind=4), parameter :: prandtl   = 0.72
   real(kind=4), parameter :: b4        = 0.6154
   real(kind=4), parameter :: ao        = 7.352e-4
   real(kind=4), parameter :: bo        = 8.748e-4
   real(kind=4), parameter :: co        = 0.205935
   real(kind=4), parameter :: dzero     = 0.6052
   real(kind=4), parameter :: en        = 2.302585
   integer(kind=4)         :: isimple   = 0
   integer(kind=4)         :: inew      = 1
   integer(kind=4)         :: inew2     = 0
   integer(kind=4)         :: insz      = 0

! depo. species pointer defined in gas mechanism package
   type depo
     real(kind=4)     :: alpha
     real(kind=4)     :: beta
     real(kind=4)     :: hstar
     real(kind=4)     :: fzero
     real(kind=4)     :: hstar2
     real(kind=4)     :: fzero2
     real(kind=4)     :: ats
     integer(kind=4)  :: sp_id
     character(len=3) :: sname
     character(len=2) :: descr
   end type depo

   type(depo), pointer :: gas_depo(:) => null()
   integer(kind=4)     :: nb_gas_depo

   integer(kind=4), parameter :: lucprm  = 15 ! 15 land use categories used in Robichaud Dev. vel. scheme
   integer(kind=4), parameter :: nsn = 5      ! 5 Seasonal periods in Robichaud Dev. vel. scheme
   integer(kind=4), parameter :: nmth = 12    ! 12 months of LAI data

   logical(kind=4) :: chm_lu15_out_l  ! Is landuse fraction LU15 requested for output
   
!  RCUT : leaf cuticle resistance 
!  RCONV: buoyant convection in canopy (called r_dc by Wesely, 1989)
!  REXP : exposed surfaces resistance (leaves, twigs, barks, etc.; called r_lu and r_cl by Wesely, 1989)
!  RCAN : in canopy aerodynamic resistance (same as r_ac in Wesely, 1989)

!  Cuticle resistance for SO2 and ozone under dry conditions (from Zhang et al., 2002)

   real(kind=4), dimension(lucprm, nsn) :: rcutd = reshape ((/                                              &
   1000., 1000., 2000., 1200., 1000., 1500., 1500.,9999., 9999., 4000., 6000., 9999., 9999., 9999., 6000.,  &
   1500., 1500., 3000., 2000., 1500., 2000., 2000.,9999., 9999., 1500., 6000., 9999., 9999., 9999., 6000.,  &
   2000., 2000., 8000., 9000., 2000., 3000., 3000.,9999., 9999., 2000., 9000., 9999., 9999., 9999., 9000.,  &
   2000., 2000., 9999., 9999., 2000., 6000., 9999.,9999., 9999., 3000., 9999., 9999., 9999., 9999., 9999.,  &
   1000., 1000., 4000., 2000., 1000., 1500., 1500.,9999., 9999., 4000., 6000., 9999., 9999., 9999., 6000.   &
                                                      /),(/lucprm, nsn/))

!  Rgds  --- dry ground resistance for SO2 (WESELY, 1989; note different interpretation vs. Zhang et al., 2002)

   real(kind=4), dimension(lucprm, nsn) :: rgdso2 = reshape ((/                              &
   500., 500., 500., 500., 100., 350., 150., 1000., 400., 400., 10., 100., 10., 10., 400.,   &
   500., 500., 500., 500., 100., 350., 200., 1000., 400., 400., 10., 100., 10., 10., 400.,   &
   500., 500., 500., 500., 200., 350., 150., 1000., 400., 400., 10., 100., 10., 10., 400.,   &
   100., 100., 100., 100., 100., 100., 100., 1000.,  50.,  50.,100., 200., 10., 10., 100.,   &
   500., 500., 500., 500., 200., 350., 150., 1000., 400., 400., 10., 100., 10., 10., 500.    &
                                                      /),(/lucprm, nsn/))

!  Rgdo  --- dry ground resistance for O3  (WESELY, 1989; note different interpretation vs. Zhang et al., 2002)

   real(kind=4), dimension(lucprm, nsn) ::  rgdo3 =  reshape ((/                                              &
   200. ,  200.,  200.,  200.,  300.,  200.,  150.,  400.,  200.,  200., 1000., 2000., 2000., 2000.,  300.,   &
   200. ,  200.,  200.,  200.,  300.,  200.,  150.,  400.,  200.,  200.,  800., 2000., 2000., 2000.,  300.,   &
   200. ,  200.,  200.,  200.,  300.,  200.,  150.,  400.,  200.,  200., 1000., 2000., 2000., 2000.,  300.,   &
   3500., 3500., 3500., 3500., 3500., 3500., 3500.,  400., 3500., 3500., 3500., 2000., 2000., 2000.,  600.,   &
   200. ,  200.,  200.,  200.,  300.,  200.,  150.,  400.,  200.,  200., 1000., 2000., 2000., 2000.,  300.    &
                                                      /),(/lucprm, nsn/))

!  resistance for other exposed surfaces (ozone) (referred to as r_lu (for upper canopy) and r_clO (for lower
!      canopy) in WESELY, 1989); N.B. values for tundra (LUC9) and dwarf trees and shrubs (LUC10) were obtained
!      using Wesely's class 11 (rocky open spaces with low-growing shrubs) and class 4 (deciduous forest) as
!      proxies

   real(kind=4), dimension(lucprm, nsn) :: rexpo3 = reshape ((/                                               &
   1000., 1000., 1000., 1000., 1000., 1000., 1000., 9999., 1000., 1000., 1000., 9999., 9999., 9999., 9999.,   &
   1000., 1000.,  400.,  400.,  600.,  400.,  400., 9999.,  400.,  400.,  600., 9999., 9999., 9999., 9999.,   &
   1000., 1000.,  400.,  400., 1000.,  400., 1000., 9999.,  600.,  400.,  600., 9999., 9999., 9999., 9999.,   &
   1500., 1500.,  400.,  400.,  600., 1000., 1000., 9999.,  800.,  400.,  600., 9999., 9999., 9999., 9999.,   &
   1500., 1500.,  500.,  500.,  700.,  500., 1000., 9999.,  800.,  500.,  600., 9999., 9999., 9999., 9999.    &
                                                      /),(/lucprm, nsn/))

!  resistance for other exposed surfaces (sulfur dioxide) (referred to as r_lu (for upper canopy) and r_clS 
!      (for lower canopy) in WESELY, 1989); N.B. values for tundra (LUC9) and dwarf trees and shrubs (LUC10)  
!      were obtained using Wesely's class 11 (rocky open spaces with low-growing shrubs) and class 4 
!      (deciduous forest) as proxies

   real(kind=4), dimension(lucprm, nsn) :: rexpso2 = reshape ((/                                              &
   2000., 2000., 2000., 2000., 2000., 2000., 2000., 9999., 4000., 2000., 2500., 9999., 9999., 9999., 9999.,   &
   2000., 2000., 9000., 9000., 4000., 9000., 9000., 9999., 9000., 9000., 9000., 9999., 9999., 9999., 9999.,   &
   3000., 3000., 9000., 9000., 6000., 9000., 9000., 9999., 9000., 9000.,  400., 9999., 9999., 9999., 9999.,   &
    200.,  200., 9000., 9000.,  400., 9999., 9999., 9999., 9000., 9000.,  400., 9999., 9999., 9999., 9999.,   &
   2000., 2000., 4000., 4000., 3000., 4000., 4000., 9999., 8000., 4000., 3000., 9999., 9999., 9999., 9999.    &
                                                       /),(/lucprm, nsn/))

!  canopy resistance (referred to as r_ac, in-canopy aerodynamic resistance, by WESELY, 1989); N.B. values
!     for dwarf trees and shrubs (LUC10) were obtained from Wesely's class 11 (rocky open spaces with 
!     low-growing shrubs)
!
!   real(kind=4), dimension(lucprm, nsn) :: rcanp = reshape ((/                            &
!   2000., 2000., 2000., 2000., 2000., 100., 200., 0., 0., 200., 300., 0., 0., 0., 100.,   &
!   2000., 2000., 1500., 1500., 1700., 100., 150., 0., 0., 140., 200., 0., 0., 0., 100.,   &
!   2000., 2000., 1000., 1000., 1500., 100.,  10., 0., 0., 120., 100., 0., 0., 0., 100.,   &
!   2000., 2000., 1000., 1000., 1500.,  10.,  10., 0., 0.,  50.,  50., 0., 0., 0., 100.,   &
!   2000., 2000., 1200., 1200., 1500.,  80.,  50., 0., 0., 120., 200., 0., 0., 0., 100.    &
!                                                     /),(/lucprm, nsn/))

!  new values for canopy resistance (referred to as r_ac, in-canopy aerodynamic resistance, by Zhang et al., 
!     2003) as a test.  values have been extracted from Table 1 of Zhang et al. (2003) and seasonal variation
!     is based on the range of values in this table plus the seasonal variation of leaf area index given 
!     below in the LAINDEX table

   real(kind=4), dimension(lucprm, nsn) :: rcanp = reshape ((/                  &
   100., 250., 100., 250., 190., 40., 40., 0., 0., 60., 20., 0., 0., 0., 40.,   &
   100., 250.,  85., 190., 150., 40., 40., 0., 0., 60., 20., 0., 0., 0., 40.,   &
   100., 250.,  70., 115., 110., 40., 10., 0., 0., 30., 20., 0., 0., 0., 40.,   &
   100., 250.,  60., 100., 100., 10., 10., 0., 0., 20., 20., 0., 0., 0., 40.,   &
   100., 250.,  60., 100., 100., 30., 20., 0., 0., 30., 20., 0., 0., 0., 40.    &
                                                     /),(/lucprm, nsn/))

!  Stomatal-resistance-related parameters (rsmin from WESELY, 1989: based on r_i); note that 
!     values for tundra (LUC9) and dwarf trees and shrubs (LUC10) were obtained using Wesely's
!     class 11 (rocky open spaces with low-growing shrubs) and class 4 (deciduous forest) as proxies

   real(kind=4), dimension(lucprm, nsn) :: rsmin = reshape ((/                                            &
   130., 130.,   70.,   70., 100.,  120.,   60., 9999.,  150.,   70.,   80., 9999., 9999., 9999., 9999.,  &
   250., 250., 9999., 9999., 800., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.,  &
   250., 250., 9999., 9999., 800., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.,  &
   400., 400., 9999., 9999., 800., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999., 9999.,  &
   250., 250.,  140.,  140., 190.,  240.,  120., 9999.,  300.,  140.,  160., 9999., 9999., 9999., 9999.   &
                                                      /),(/lucprm, nsn/))

!  from Zhang, Brook, and Vet (2002) for 15 land-use categories: 

   real(kind=4), dimension(lucprm) :: tmin = (/ -5.0,    0.0,   -5.0,    0.0, -3.0, &
                                         5.0,    5.0, 9999.0, 9999.0,  0.0, &
                                         5.0, 9999.0, 9999.0, 9999.0,  0.0  /)

   real(kind=4), dimension(lucprm) :: tmax = (/ 40.0,   45.0,   40.0,   45.0, 42.0, &
                                        45.0,   45.0, 9999.0, 9999.0, 43.0, &
                                        45.0, 9999.0, 9999.0, 9999.0, 45.0  /)

   real(kind=4), dimension(lucprm) :: topt = (/ 15.0,   30.0,   15.0,   27.0, 21.0, &
                                        27.0,   25.0, 9999.0, 9999.0, 21.5, &
                                        25.0, 9999.0, 9999.0, 9999.0, 22.0  /)

!  Leaf Area Index from Zhang et al. (2002)

   real(kind=4), dimension(lucprm, nsn) :: laindex = reshape ((/                                             &
   5.30,  4.50,  1.10,  3.40,  4.50,  2.00,  2.00,  0.00, 0.00,  2.50,  0.20,  0.00,  0.00,  0.00,  0.30,    &
   5.30,  4.50,  0.80,  1.90,  3.50,  1.50,  1.50,  0.00, 0.00,  2.50,  0.20,  0.00,  0.00,  0.00,  0.20,    &
   4.70,  4.50,  0.30,  0.10,  2.30,  1.00,  1.00,  0.00, 0.00,  1.50,  0.10,  0.00,  0.00,  0.00,  0.10,    &
   5.50,  4.50,  0.00,  0.00,  2.30,  0.50,  0.00,  0.00, 0.00,  1.20,  0.00,  0.00,  0.00,  0.00,  0.00,    &
   5.50,  4.50,  0.00,  0.00,  2.30,  0.50,  0.00,  0.00, 0.00,  0.50,  0.10,  0.00,  0.00,  0.00,  0.20     &
                                                        /),(/lucprm, nsn/))

!  Revised LAI from Shailesh Kharol(2020), by month starting in January, for the forest categories:  
!  Evergreen Needleleaf Forest,  Evergreen Broadleaf Forest,  Deciduous Needleleaf Forest,
!   Deciduous Broadleaf Forest,  Mixed Forest,Grassland,  "Crops, mixed farming",
!   Desert,Tundra,  "Dwarf trees, shrubs with ground cover",  Wetland with plants,
!   Ice caps and glaciers,  Inland Water,  Ocean,  Urban
   real(kind=4), dimension(lucprm, nmth) :: laindex_sat = reshape ((/ &
   1.40,  3.10,  0.10,  1.10,  1.50,  0.50,  0.50,  0.00,  0.30,  0.30,  0.80,  0.00,  0.00,  0.00,  0.50,   & 
   1.60,  3.10,  0.10,  1.10,  1.50,  0.50,  0.50,  0.00,  0.40,  0.30,  0.90,  0.00,  0.00,  0.00,  0.50,   & 
   1.80,  3.60,  0.10,  1.10,  1.90,  0.50,  0.50,  0.00,  0.50,  0.40,  0.90,  0.00,  0.00,  0.00,  0.60,   & 
   1.90,  4.60,  0.30,  1.10,  2.00,  0.50,  0.60,  0.00,  0.50,  0.50,  1.70,  0.00,  0.00,  0.00,  0.70,   & 
   2.00,  5.80,  1.00,  2.30,  2.10,  0.50,  0.90,  0.00,  0.50,  0.80,  1.70,  0.00,  0.00,  0.00,  1.10,   & 
   2.40,  5.80,  1.50,  5.50,  4.10,  0.50,  1.60,  0.00,  0.60,  0.90,  1.80,  0.00,  0.00,  0.00,  1.20,   & 
   2.90,  5.80,  2.30,  5.60,  5.20,  1.00,  2.70,  0.00,  0.90,  0.90,  1.80,  0.00,  0.00,  0.00,  1.20,   & 
   2.90,  5.60,  1.90,  5.70,  5.10,  2.00,  2.40,  0.00,  1.00,  0.61,  1.50,  0.00,  0.00,  0.00,  1.20,   & 
   2.10,  5.40,  1.30,  5.40,  4.10,  2.00,  1.70,  0.00,  0.50,  0.51,  1.10,  0.00,  0.00,  0.00,  1.30,   & 
   1.80,  4.40,  0.80,  2.90,  1.90,  1.50,  0.80,  0.00,  0.50,  0.41,  0.90,  0.00,  0.00,  0.00,  0.80,   & 
   1.60,  4.40,  0.70,  1.10,  1.90,  1.00,  0.60,  0.00,  0.30,  0.40,  0.80,  0.00,  0.00,  0.00,  0.60,   & 
   1.40,  3.10,  0.60,  1.10,  1.30,  1.00,  0.50,  0.00,  0.20,  0.30,  0.70,  0.00,  0.00,  0.00,  0.60    & 
                                    /), (/lucprm, nmth/))
!
!  Surface Roughness Length [m] from Zhang et al. (2002)
!  Note that urban area roughness length was increased from 1.0 to 2.0

   real(kind=4), dimension(lucprm, nsn) :: zz0 = reshape ((/                                                     &
   0.80,  2.65,  0.85,  1.05,  1.15,  0.10,  0.10,  0.04, 0.03,  0.10,  0.03,  0.01,  0.0001,  0.0001 ,  2.00,   &
   0.90,  2.65,  0.90,  1.05,  1.15,  0.10,  0.10,  0.04, 0.03,  0.10,  0.03,  0.01,  0.0001,  0.0001,   2.00,   &
   0.90,  2.65,  0.80,  0.95,  1.15,  0.05,  0.02,  0.04, 0.03,  0.10,  0.02,  0.01,  0.0001,  0.0001,   2.00,   &
   0.90,  2.65,  0.55,  0.55,  1.15,  0.02,  0.02,  0.04, 0.03,  0.10,  0.02,  0.01,  0.0001,  0.0001,   2.00,   &
   0.90,  2.65,  0.60,  0.75,  1.15,  0.05,  0.05,  0.04, 0.03,  0.10,  0.03,  0.01,  0.0001,  0.0001,   2.00    &
                                                    /),(/lucprm, nsn/))

!  parameters for dry deposition velocity of PM
   real(kind=4), dimension(lucprm, nsn) :: pllp = reshape ((/                                          &
   2.00, 5.00, 2.00,  5.00, 5.00, 2.00,  2.00, -1.00, -1.00, 10.00, 10.00, -1.00, -1.00, -1.00, 10.00, &
   2.00, 5.00, 2.00,  5.00, 5.00, 2.00,  2.00, -1.00, -1.00, 10.00, 10.00, -1.00, -1.00, -1.00, 10.00, &
   2.00, 5.00, 5.00, 10.00, 5.00, 5.00, 10.00, -1.00, -1.00, 10.00, 10.00, -1.00, -1.00, -1.00, 10.00, &
   2.00, 5.00, 5.00, 10.00, 5.00, 5.00, 10.00, -1.00, -1.00, 10.00, 10.00, -1.00, -1.00, -1.00, 10.00, &
   2.00, 5.00, 2.00,  5.00, 5.00, 2.00,  2.00, -1.00, -1.00, 10.00, 10.00, -1.00, -1.00, -1.00, 10.00  &
                                                            /), (/lucprm, nsn/))
   real(kind=4), dimension(lucprm) :: aest = (/ 1.00, 0.60, 1.10, 0.80, 0.80, 1.20, 1.20, 50.00, 50.00, &
                                                1.30, 2.00, 50.00, 100.00, 100.00, 1.50 /)
   real(kind=4), dimension(lucprm) :: pgamma = (/ -.56, -.58, -.56, -.56, -.56, -.54, -.54, -.54, -.54, &
                                                  -.54, -.54, -.54, -.50, -.50, -.56 /)

end module mach_drydep_mod
