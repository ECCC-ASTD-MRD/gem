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
! Projet/Project : GEM-MACH
! Fichier/File   : mach_soa_odum_mod.ftn90
! Creation       : Jack C. and Deji A. 2020
! Description    : contains parameters and routines needed for calls to
!                  mach_soa_odum
!                  handles different versions of saprc mechanisms
!
! snippets of this may be better defined/combined with mechanism specific
!   modules i.e mach_pkg_saprc07_mod.ftn90
!
! 2021-06 - Craig updated cmatx and alpha coefficients
!============================================================================

module mach_soa_odum_mod
   use chm_utils_mod,    only: dp, sp, i4
   use mach_pkg_gas_mod, only: nvsoa

   implicit none

   public  :: map_soa_dvoc_svoc_saprc

!! parameters specific to SOA method but independent of mechanism
   integer(i4), parameter :: nvbin = 5  ! number of volitility bins
   integer(i4), parameter :: ncond = 6  ! number of oxidation conditions
   integer(i4), parameter :: nspc  = 14 ! total SOA precursor species tracked
   integer(i4), parameter :: lo_OH = 1, hi_OH = 2, lo_O3 = 3, hi_O3 = 4, &
                             lo_NO3 = 5, hi_NO3 = 6
   integer(i4), parameter :: mstep = 10 ! number of increment MOi steps

!! expandable SOA precursor species considered here, in alphabatical order
!!! this is to be mapped from input precursor gas_voc
   integer(i4), parameter :: soa_alk4 = 1,  soa_alk5 = 2,  soa_ole1 = 3,  &
                             soa_ole2 = 4,  soa_aro1 = 5,  soa_aro2 = 6,  &
                             soa_isop = 7,  soa_apin = 8,  soa_bpin = 9,  &
                             soa_limo = 10, soa_benz = 11, soa_sesq = 12, &
                             soa_ivoc = 13, soa_naph = 14

   integer(i4), parameter :: i_kOH = 1, i_kO3 = 2, i_kNO3 = 3

   real(sp), parameter :: h_vap = 37000.d0
   real(dp), parameter, dimension(nvbin, ncond, nspc) :: cmatx = reshape([ &
!! parameters for each SOA precursors (alphabatical order)
   ! C_xxx volatility bins define the effective saturation concentrations
   ! Col-nvbin: C1, C2, C3, C4, C5
   ! Row-ncond: hi_OH, lo_O3, hi_O3, lo_NO3, hi_NO3
   !soa_alk4 based on Tsimpidi 2010
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_alk5 based on Tsimpidi 2010
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_ole1
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_ole2
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_aro1
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_aro2
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_isop
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_apin
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_bpin
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_limo
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_benz
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_sesq
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_ivoc
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
   !soa_naph
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0, &
                   0.1d0, 1.0d0, 10.0d0, 100.0d0, 1000.0d0], &
                   [nvbin, ncond, nspc])
!
   real(dp), parameter, dimension(nvbin, ncond, nspc) :: alpha = reshape([ &
   ! A_xxx mass stoichiometric coefficients for products of a precursor reaction
   ! Col-nvbin: C1, C2, C3, C4, C5
   ! Row-ncond: lowNOx_OH, hiNOx_OH, lo_O3, hi_O3, lo_NO3, hi_NO3
   !soa_alk4
                  0.0d0, 0.0d0, 0.075d0, 0.0d0, 0.0d0, & !lowNOx OH, Tsi2010, SOAyield=0.0068
                  0.0d0, 0.0d0, 0.038d0, 0.0d0, 0.0d0, & !high NOx OH, Tsi2010, SOAyield=0.0035
                  0.0d0, 0.0d0,   0.0d0, 0.0d0, 0.0d0, & !no reaction
                  0.0d0, 0.0d0,   0.0d0, 0.0d0, 0.0d0, & !no reaction
                  0.0d0, 0.0d0,   0.0d0, 0.0d0, 0.0d0, & !no reaction
                  0.0d0, 0.0d0,   0.0d0, 0.0d0, 0.0d0, & !no reaction
   !soa_alk5 based on Tsimpidi 2010 den=1.5
                  0.0d0, 0.0d0, 0.3d0,  0.0d0, 0.0d0, & !lowNOx OH, Tsi2010, SOAyield=0.027
                  0.0d0, 0.0d0, 0.15d0, 0.0d0, 0.0d0, & !highNOx OH, Ma, SOAyield=0.078
                  0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, &
                  0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, &
                  0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, &
                  0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, &
   !soa_ole1
               0.0d0, 0.014d0, 0.0d0, 0.098d0, 0.088d0, & ! Ma et al, SOAyield=0.0081 @ 1ug/m3
               0.0d0, 0.014d0, 0.0d0, 0.098d0, 0.088d0, &
               0.0d0, 0.014d0, 0.0d0, 0.098d0, 0.088d0, &
               0.0d0, 0.014d0, 0.0d0, 0.098d0, 0.088d0, &
               0.0d0, 0.014d0, 0.0d0, 0.098d0, 0.088d0, &
               0.0d0, 0.014d0, 0.0d0, 0.098d0, 0.088d0, &
   !soa_ole2
               0.0d0, 0.052d0, 0.000d0, 0.183d0, 0.157d0, & ! Ma et al, SOAyield=0.028 @ 1ug/m3
               0.0d0, 0.052d0, 0.000d0, 0.183d0, 0.157d0, &
               0.0d0, 0.052d0, 0.000d0, 0.183d0, 0.157d0, &
               0.0d0, 0.052d0, 0.000d0, 0.183d0, 0.157d0, &
               0.0d0, 0.052d0, 0.000d0, 0.183d0, 0.157d0, &
               0.0d0, 0.052d0, 0.000d0, 0.183d0, 0.157d0, &
   !soa_aro1
               0.0d0, 0.457d0, 0.000d0, 0.502d0, 0.251d0, &  !lowNOx, OH, Ma, SOAyield=0.23 @1ug/m3
               0.0d0, 0.256d0, 0.002d0, 0.431d0, 0.202d0, & !highNOx, OH, Ma, SOAyield=0.13 @1ug/m3
               0.0d0, 0.457d0, 0.000d0, 0.502d0, 0.251d0, &
               0.0d0, 0.256d0, 0.002d0, 0.431d0, 0.202d0, &
               0.0d0, 0.457d0, 0.000d0, 0.502d0, 0.251d0, &
               0.0d0, 0.256d0, 0.002d0, 0.431d0, 0.202d0, &
   !soa_aro2
               0.0d0, 0.539d0, 0.000d0, 0.491d0,0.255d0, & !lowNOx, OH, Ma, SOAyield=0.27 @1ug/m3
               0.0d0, 0.310d0, 0.000d0, 0.420d0,0.209d0, &  !highNOx, OH, Ma, SOAyield=0.16
               0.0d0, 0.539d0,   0.0d0, 0.491d0,0.255d0, &
               0.0d0, 0.310d0,   0.0d0, 0.420d0,0.209d0, &
               0.0d0, 0.539d0,   0.0d0, 0.491d0,0.255d0, &
               0.0d0, 0.310d0,   0.0d0, 0.420d0,0.209d0, &
   !soa_isop
               0.0d0,    0.0510d0,  0.0d0,   0.003d0,  0.0000d0, & !lowNOx, OH, , SOAyield=0.025 @1ug/m3
               0.0d0,    0.034d0,   0.000d0, 0.0050d0, 0.0000d0, & !highNOx, OH, , SOAyield = 0.017
               0.0000d0, 0.0510d0,  0.00d0,  0.0030d0, 0.0000d0, & !lowNOx, O3
               0.0000d0, 0.034d0,   0.000d0, 0.0050d0, 0.0000d0, & !high NOx, O3
               0.000d0,  0.291d0,   0.000d0, 0.007d0,  0.018d0,  & !lowNOx, NO3, Ma,  SOAyield=0.15
               0.373d0,  0.033d0,   0.000d0, 0.941d0,  0.0d0,    & !highNOx, NO3, Boyd2015, SOAyield=0.36
   !soa_apin
               0.000d0,  0.337d0,  0.000d0,  0.539d0,  0.295d0, & ! lowNOx, OH, Ma, SOAyield=0.174
               0.000d0,  0.210d0,  0.000d0,  0.348d0,  0.297d0, &  !highNOx, OH, Ma, SOA yield=0.11
               0.00d0,   0.337d0,  0.00d0,   0.539d0,  0.295d0, & !lowNOx, O3
               0.0000d0, 0.210d0,  0.000d0,  0.348d0,  0.297d0, & !highNOx, O3
               0.000d0,  0.399d0,  0.783d0,  0.235d0,  0.000d0, & !lowNOx, NO3, Ma , SOAyield=0.27
               0.000d0,  0.399d0,  0.783d0,  0.235d0,  0.000d0, & !
   !soa_bpin
               0.000d0,  0.337d0,  0.000d0,  0.539d0,  0.295d0, & ! lowNOx, OH, Ma, SOAyield=0.174
               0.000d0,  0.210d0,  0.000d0,  0.348d0,  0.297d0, &  !highNOx, OH, Ma, SOA yield=0.11
               0.00d0,   0.337d0,  0.00d0,   0.539d0,  0.295d0, & !lowNOx, O3
               0.0000d0, 0.210d0,  0.000d0,  0.348d0,  0.297d0, & !highNOx, O3
               0.000d0,  0.399d0,  0.783d0,  0.235d0,  0.000d0, & !lowNOx, NO3, Ma , SOAyield=0.27
               0.000d0,  0.399d0,  0.783d0,  0.235d0,  0.000d0, & !
   !soa_limo
               0.000d0,  0.337d0,  0.000d0,  0.539d0,  0.295d0, & ! lowNOx, OH, Ma, SOAyield=0.174
               0.000d0,  0.210d0,  0.000d0,  0.348d0,  0.297d0, &  !highNOx, OH, Ma, SOA yield=0.11
               0.00d0,   0.337d0,  0.000d0,  0.539d0,  0.295d0, & !lowNOx, O3
               0.0000d0, 0.210d0,  0.000d0,  0.348d0,  0.297d0, & !highNOx, O3
               0.000d0,  0.399d0,  0.783d0,  0.235d0,  0.000d0, & !lowNOx, NO3, Ma , SOAyield=0.27
               0.000d0,  0.399d0,  0.783d0,  0.235d0,  0.000d0, & !
   !soa_benz low NOx Ng 2007, high NOx Pye/Ng 2007
               0.37d0,  0.0d0, 0.0d0,   0.0d0, 0.0d0, &  !SOAyield=0.34
               0.0d0, 0.078d0, 0.0d0, 0.793d0, 0.0d0, &  !SOAyield=0.047
               0.0d0,   0.0d0, 0.0d0,   0.0d0, 0.0d0, &
               0.0d0,   0.0d0, 0.0d0,   0.0d0, 0.0d0, &
               0.0d0,   0.0d0, 0.0d0,   0.0d0, 0.0d0, &
               0.0d0,   0.0d0, 0.0d0,   0.0d0, 0.0d0, &
   !soa_sesq
               0.0d0, 0.463d0, 0.072d0, 0.953d0, 0.401d0, & !lowNOx, OH, Tsi, SOAyield=0.25
               0.0d0, 0.075d0, 0.150d0, 0.750d0, 0.900d0, & !highNOx, OH, Hayes report, Ma, SOAyield=0.11
               0.0d0, 0.463d0, 0.072d0, 0.953d0, 0.401d0, & !lowNOx, OH, Tsi, SOAyield=0.25
               0.0d0, 0.075d0, 0.150d0, 0.750d0, 0.900d0, & !highNOx, OH, Hayes report, Ma, SOAyield=0.11
               0.0d0, 0.399d0, 0.783d0, 0.235d0, 0.000d0, & !lowNOx, NO3, Ma, SOAyield=0.27
               0.0d0, 0.399d0, 0.783d0, 0.235d0, 0.000d0, &
   ! soa_ivoc
               0.0d0, 0.088d0,  0.1420d0, 0.82d0, 0.60d0, & !lowNOx, OH, 2x Presto , SOAyield=0.066
               0.0d0, 0.044d0,  0.071d0,  0.41d0, 0.30d0, & !highNOx, OH, Presto2010 C13, SOAyield=0.033
               0.0d0, 0.088d0,  0.1420d0, 0.82d0, 0.60d0, &
               0.0d0, 0.044d0,  0.071d0,  0.41d0, 0.30d0, &
               0.0d0, 0.0880d0, 0.1420d0, 0.82d0, 0.60d0, &
               0.0d0, 0.044d0,  0.071d0,  0.41d0, 0.30d0, &
   ! soa_naph
               0.73d0,0.0d0,   0.0d0,  0.0d0, 0.0d0, &  !lowNOx, OH, Chan2009, SOAyield=0.70
               0.0d0, 0.039d0, 0.296d0, 0.235d0, 0.0d0, & !highNOx, OH, Chan2009 SOAyield=0.049
               0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
               0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
               0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
               0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], &
               [nvbin, ncond, nspc])
   save

!===============================================================================
 CONTAINS

   subroutine map_soa_dvoc_svoc_saprc(dvoc_ppm, dvoc_soa, tt, gni, gnk)
!! Routine that maps input dvoc_gas with precursor dvoc_soa considered for
!! yield calculation.
!! The order of input gas_xxx species parameter is needed to match input array.
!! convert input dvoc_gas from ppmv to ug/m3 and output to dvoc_soa
      use chm_nml_mod,      only: chm_pkg_gas_s, nk_start
      use mach_pkg_gas_mod, only: nsp_soa_gases, soa_gas_idx, gas_species
      use chm_consphychm_mod,   only: avno
      use chm_species_info_mod, only: sm

      implicit none
      integer(kind=4), intent (in) :: gni, gnk
      real(dp),        intent (in) :: dvoc_ppm(nsp_soa_gases, gni, gnk)
      real(dp),        intent(out) :: dvoc_soa(nvsoa, gni, gnk)
      real(sp),        intent (in) :: tt(gni, gnk) ! temp in Kevin

! local variables
 !! precursor species index number order as in soa_gas_idx defined in
 !   mach_pkg_saprc_mod.ftn90
      integer(kind=4) :: gas_alk4, gas_alk5, gas_ole1, gas_ole2, gas_aro1, &
                         gas_aro2, gas_isop, gas_cres, gas_terp, gas_benz, &
                         gas_sesq, gas_ivoc, i, k, jj, isp
      real(dp)        :: dvoc_gas(nsp_soa_gases, gni, gnk) ! in ug/m3

      ! convert from ppmv to ug/m3
      do jj = 1, nsp_soa_gases
         do k = nk_start, gnk
            do i = 1, gni
               isp = soa_gas_idx(jj)
               dvoc_gas(isp, i, k) = dvoc_ppm(isp, i, k) * 1.d-12 / avno * &
                                     sm(gas_species(isp))%mol_wt
            end do
         end do
      end do


      if (chm_pkg_gas_s(1:7) == 'SAPRC07') then
         gas_alk4 = 1; gas_ole1 = 2; gas_ole2 = 3; gas_aro1 = 4
         gas_aro2 = 5; gas_isop = 6; gas_cres = 7; gas_terp = 8
      else if (chm_pkg_gas_s(1:7) == 'SAPRC11') then
         gas_alk4 = 1; gas_alk5 = 2; gas_ole1 = 3; gas_ole2 = 4; gas_aro1 = 5
         gas_aro2 = 6; gas_isop = 7; gas_terp = 8; gas_benz = 9; gas_sesq = 10
         gas_ivoc = 11
      end if
      ! The above can be done dynamically as follows;
      ! gas_alk4 = findloc(soa_gas_idx, sp_ALK4, dim=1)
      ! gas_alk5 = findloc(soa_gas_idx, sp_ALK5, dim=1)
      ! ...
      ! gas_sesq = findloc(soa_gas_idx, sp_SESQ, dim=1)

      dvoc_soa = 0.d0
      do k = nk_start, gnk
         do i = 1, gni
            dvoc_soa(soa_ole1, i, k) = dvoc_gas(gas_ole1, i, k)
            dvoc_soa(soa_ole2, i, k) = dvoc_gas(gas_ole2, i, k)
            dvoc_soa(soa_aro1, i, k) = dvoc_gas(gas_aro1, i, k)
            dvoc_soa(soa_aro2, i, k) = dvoc_gas(gas_aro2, i, k)
!           dvoc_soa(soa_naph, i, k) = 0.d0 !! not tracked in SAPRC
            dvoc_soa(soa_isop, i, k) = dvoc_gas(gas_isop, i, k)
      !! TERP split equally to apin, bpin, limo
            dvoc_soa(soa_apin, i, k) = dvoc_gas(gas_terp, i, k) * 0.333d0
            dvoc_soa(soa_bpin, i, k) = dvoc_soa(soa_apin, i, k)
            dvoc_soa(soa_limo, i, k) = dvoc_soa(soa_apin, i, k)
         end do
      end do
!
      if (chm_pkg_gas_s(1:7) == 'SAPRC07') then
         do k = nk_start, gnk
            do i = 1, gni
      !! ALK5 split equally to ALK4, ALK5 for SAPRC07C(Carter, 2010)
               dvoc_soa(soa_alk4, i, k) = dvoc_gas(gas_alk4, i, k) * 0.5d0
               dvoc_soa(soa_alk5, i, k) = dvoc_soa(soa_alk4, i, k)
      !! BENZ not in SAPRC07CS split from cres and aro1 (Carter, 2010)
               dvoc_soa(soa_benz, i, k) = 0.45d0 * dvoc_gas(gas_aro1, i, k) + &
                                          0.145d0 * dvoc_gas(gas_cres, i, k)
      !! sesq from terp with temp correction based on (Helmig et al., ES&T 2007)
      ! Eos,Eot are basal emission rate at Ts=30oC, Bs, Bt are temp response
      ! factor in decC Eos=0.302, Eot=2.085, Bt=0.11, Bs=0.17
      ! sesq= terp * Eos/Ett*exp((T-30)*(Bs-Bt))
               dvoc_soa(soa_sesq, i, k) = dvoc_soa(gas_terp, i, k) * 0.14484d0 * &
                                          dble(exp((tt(i, k) - 303.15) * 0.06))
            end do
         end do
      else if (chm_pkg_gas_s(1:7) == 'SAPRC11') then
         do k = nk_start, gnk
            do i = 1, gni
               dvoc_soa(soa_alk4, i, k) = dvoc_gas(gas_alk4, i, k)
               dvoc_soa(soa_alk5, i, k) = dvoc_gas(soa_alk5, i, k)
               dvoc_soa(soa_benz, i, k) = dvoc_gas(gas_benz, i, k)
               dvoc_soa(soa_sesq, i, k) = dvoc_gas(gas_sesq, i, k)
               dvoc_soa(soa_ivoc, i, k) = dvoc_gas(gas_ivoc, i, k)
            end do
         end do
      end if

      return
   end subroutine map_soa_dvoc_svoc_saprc

end module mach_soa_odum_mod
