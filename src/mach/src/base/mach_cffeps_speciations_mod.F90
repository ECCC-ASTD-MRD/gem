!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2018 - Air Quality Research Division &
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
! Fichier / File   : mach_cffeps_emiss_speciations.ftn90
! Creation         : A. Akingunola, J. Chen, P. Makar, and Kerry Anderson
!                     - Fall 2018
! Description      : Emissions' speciations for CFFEPS
!
!============================================================================
!
module mach_cffeps_speciations_mod

  use mach_cffeps_mod, only: cffeps_nspec
  implicit none
!
  real(kind=4), dimension(3, cffeps_nspec), parameter :: emission_factors = &
!                              F,       S,       R
!PM10
                reshape((/  16.05,   27.38,   40.75,   &
!PM2.5
                            13.60,   23.20,   34.53,   &
!!PM
!                            23.0,    34.0,    34.0,    &
!CO
                            83.0,    135.0,   248.0,   &
!CO2
                            1662.33, 1600.0,  1383.0,  &
!CH4
                            3.23,    7.32,    9.94,    &
!NOX
                            1.83,    2.0,     0.45,    &
!NH3
                            0.99,    1.50,    1.94,    &
!SO2
                            0.93,    1.06,    1.76,    &
!NMHC
                            19.85,   33.87,   56.08/), (/3, cffeps_nspec/))
!
  integer, parameter :: i_pm10 = 1, i_pm25 = 2, i_co = 3, i_co2 = 4, i_ch4 = 5, &
                         i_nox = 6, i_nh3 = 7, i_so2 = 8, i_nmhc = 9
!
! GEM-MACH speciated species
  integer(kind=4), parameter :: ngas  = 7
  integer(kind=4), parameter :: npm   = 6
  integer(kind=4), parameter :: nvocs = 14
  character(len=4), dimension(ngas), parameter :: me_species_gas =    &
                (/ 'ECO ', 'ECO2', 'ENO ', 'ENO2', 'ENH3', 'ESO2', 'ESO4'/)
  character(len=4), dimension(npm),  parameter :: me_species_pm =   &
                (/'EAM ', 'ECM ', 'EEC ', 'ENT ', 'EPC ', 'ESU '/)
!#     Methane (CH4) - EC1, not usually used in ADOM2
!#     Ethane (C2H6) - EC2, not usually used in ADOM2
  character(len=4), dimension(nvocs), parameter :: me_species_voc = &
                (/'EA2 ', 'EA3 ', 'EALD', 'EARO', 'EC38', 'ECRE', 'EETH', &
                  'EHCH', 'EISO', 'EMEK', 'ETOL', 'EC1 ', 'EC2 ', 'EOTH'/)
  character(len=4), dimension(:), allocatable, save :: fire_me_species
  integer(kind=4), save :: nb_fire_me_species

!# Speciation profiles VOCs
!# VOCs - flaming (from profile #95425 in EPA SPECIATEv4.5)
   real(kind=4), parameter :: voc2togflm = 1.116755
   real(kind=4), parameter, dimension(nvocs) :: vocflm = &
               (/0.082981979, 0.061976018, 0.029332253, 0.013065163, &
                 0.055802877, 0.013031624, 0.03140907,  0.047365044, &
                 0.00177427,  0.0118993,   0.038453197, 0.08083481,  &
                 0.01459366,  0.516154304 /)
!# VOC speciation - smoldering (from profile #95428 in EPA SPECIATEv4.5)
   real(kind=4), parameter :: voc2togsml = 1.0664
   real(kind=4), parameter, dimension(nvocs) :: vocsml = &
               (/0.048898572, 0.071567121, 0.018099495, 0.006872343, &
                 0.024771925, 0.006245644, 0.006480751, 0.008376069, &
                 0.000410902, 0.008586111, 0.01630913,  0.043683492, &
                 0.011210157, 0.627109175/)

!# Speciation profiles for NOx and SO2
   real(kind=4), parameter :: sno2 = 0.1, sno = 0.586956522,     &
                              sso4 = 0.015467172
!# Speciation profiles for PM - 2bin only
   real(kind=4), parameter, dimension(npm) :: spm = &
                                    (/0.00879149, 0.0973774,  0.09488849, &
                                      0.001323,   0.78500909, 0.0126105/)
   !sam = 0.00879149, scm = 0.0973774, sec = 0.09488849, &
   !snt = 0.001323,   spc = 0.78500909, ssu = 0.0126105
!# 2bin to 12bin
   real(kind=4), parameter :: bin_fac(12) = &
                (/0.000005, 0.001054, 0.041466, 0.328735, 0.490694,    &
                  0.127949, 0.009841, 0.000253, 0.000002, 0.000000001, &
                  0.0, 0.0/)
  contains
!============================================================================
!

  subroutine cffeps_emiss_speciations(emiss, mj_emis, isize)

   use mach_cffeps_mod, only: emissions_type
   implicit none

   type(emissions_type),                        intent(in)  :: emiss
   integer(kind=4),                             intent(in)  :: isize
   real(kind=4), dimension(nb_fire_me_species), intent(out) :: mj_emis

! convert metric ton/hr to grams/sec
   real(kind=4), parameter :: conv = 1000000.0 / 3600.0

   integer(kind=4)  :: nsp, sn, nbin, emis_nb_bins
   real(kind=4) :: pmf, pmc, pm, co, co2, ch4, nox, nh3, so2, nmhc, nmhc_flm, &
                   nmhc_sml
!
! Sum Flaming, Smoldering, Residual emissions
! later may want to break up flaming and smoldering+residual for separate input to GEM
   pmf = (emiss%f(1) * emission_factors(1, i_pm25) + &
          emiss%s(1) * emission_factors(2, i_pm25) + &
          emiss%r(1) * emission_factors(3, i_pm25)) * conv
   pmc = (emiss%f(1) * emission_factors(1, i_pm10) + &
          emiss%s(1) * emission_factors(2, i_pm10) + &
          emiss%r(1) * emission_factors(3, i_pm10)) * conv - pmf
   pm  = pmf + pmc
   co  = (emiss%f(1) * emission_factors(1, i_co) + &
          emiss%s(1) * emission_factors(2, i_co) + &
          emiss%r(1) * emission_factors(3, i_co)) * conv
   co2 = (emiss%f(1) * emission_factors(1, i_co2) + &
          emiss%s(1) * emission_factors(2, i_co2) + &
          emiss%r(1) * emission_factors(3, i_co2)) * conv
   ch4 = (emiss%f(1) * emission_factors(1, i_ch4) + &
          emiss%s(1) * emission_factors(2, i_ch4) + &
          emiss%r(1) * emission_factors(3, i_ch4)) * conv
   nox = (emiss%f(1) * emission_factors(1, i_nox) + &
          emiss%s(1) * emission_factors(2, i_nox) + &
          emiss%r(1) * emission_factors(3, i_nox)) * conv
   nh3 = (emiss%f(1) * emission_factors(1, i_nh3) + &
          emiss%s(1) * emission_factors(2, i_nh3) + &
          emiss%r(1) * emission_factors(3, i_nh3)) * conv
   so2 = (emiss%f(1) * emission_factors(1, i_so2) + &
          emiss%s(1) * emission_factors(2, i_so2) + &
          emiss%r(1) * emission_factors(3, i_so2)) * conv
   nmhc_flm = emiss%f(1) * emission_factors(1, i_nmhc) * conv
   nmhc_sml = (emiss%s(1) * emission_factors(2, i_nmhc) + &
               emiss%r(1) * emission_factors(3, i_nmhc)) * conv
   nmhc = nmhc_flm + nmhc_sml

!# Speciate the emissions

   mj_emis(1)  = co
   mj_emis(2)  = co2
   mj_emis(3)  = sno  * nox
   mj_emis(4)  = sno2 * nox
   mj_emis(5)  = nh3
   mj_emis(6)  = so2
   mj_emis(7)  = sso4 * so2

   do sn = 1, nvocs
      mj_emis(ngas + sn) = vocflm(sn) * voc2togflm * nmhc_flm + &
                           vocsml(sn) * voc2togsml * nmhc_sml
   end do
! Override with specific emission factor for CH4
   mj_emis(ngas + 12)  = ch4

   nsp = ngas + nvocs
   emis_nb_bins = max(2, (isize - 2)) ! no emissions for bins B and C
   do sn = 1, npm
      if (isize == 2) then
         mj_emis(nsp + 1) = spm(sn) * pmf
         mj_emis(nsp + 2) = spm(sn) * pmc
         nsp = nsp + 2
      else if (isize == 12) then
         do nbin = 1, emis_nb_bins
            nsp = nsp + 1
            mj_emis(nsp) = spm(sn) * bin_fac(nbin) * pm
         end do
      end if
   end do

   return
  end subroutine cffeps_emiss_speciations

end module mach_cffeps_speciations_mod

