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
! Fichier/File   : mach_gas_lookup_table_rates_mod.ftn90
! Creation       : P. Makar  , GEM-MACH, Feb 2007
!                  Pudykiewicz / A. Kallaur, CHRONOS , 1997-2000
! Description    : Gas phase chemical mechanism - version ii
!                  (including isoprene chemistry), version 4.11.88
!
! Extra info     : For CTM see : J. Pudykiewicz et. al. "Semi-Lagrangian Modelling of
!                  Tropospheric Ozone" Tellus 49B, pp 231-248 1997
!
! Modifications  : Convert into a module and modify split out the evaluation
!                  of JNO2 and JO3 into a separate subroutine, which can be
!                  called elsewhere.
!                  D. Akingunola, Oct. 2020
!
!                : Move mach_adom2_chemi (again) here from mach_pkg_adom2_mod,
!                  since the initializations in the routine are only needed for
!                  the Young & Boris solver. Now initialzed (called) from the
!                  chm_init function
!
!============================================================================
!
module mach_adom2_rates_mod
   use chm_nml_mod,        only: nk_start, chm_messy_jval_l
   use mach_pkg_adom2_mod, only: nreac, nprcf

   implicit none

   integer(kind=4), parameter     :: nzenth = 12
   integer(kind=4), parameter     :: nzhite = 4
! Values of zenith angles (degrees)
   real(kind=4), parameter, dimension(nzenth) :: zenith = (/ 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 78.0, 86.0, 90.0, 180.0 /)

! Elevations for which type 2 data are given (m)
   real(kind=4), parameter, dimension(nzhite) :: zhite = (/ 0.0, 1000.0, 4200.0, 14000.0 /)

! Rate constant data for the version II mechanism
! a)Type 0 rate constant: constant value (in molecules-cc-sec units)
   integer(kind=4), parameter :: ntype0 = 29
   integer(kind=4), parameter, dimension(ntype0) :: indx0 = (/ 11, 17, 18, 21, 46, 47, 60, 65, 74, 75, 78, 79,   &
                                                    82, 83, 92, 93, 94, 95, 96, 97, 98, 99, 102, 104, &
                                                    105, 106, 108, 111, 112/)
   integer(kind=4), parameter, dimension(ntype0) :: jord0 = (/ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
                                                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
                                                    2, 2, 1, 2, 2/)
   real(kind=4), parameter, dimension(ntype0) :: arat0 = (/ 1.00e-21, 2.20e-10, 2.90e-11, 1.00e-24, 1.00e-15, 3.00e-12, &
                                         5.30e-12, 1.70e-11, 1.00e-15, 3.00e-12, 1.00e-15, 3.00e-12,         &
                                         1.00e-15, 3.00e-12, 1.00e-13, 1.00e-13, 2.30e-17, 2.30e-17,         &
                                         2.50e-14, 2.50e-14, 2.50e-14, 2.50e-14, 3.00e-11, 4.00e-11,         &
                                         2.20e-11, 1.50e-11, 1.00e-03, 1.80e-11, 3.20e-13/)

! b)Type 1 rate constant: function of temperature
   integer(kind=4), parameter :: ntype1 = 48
   real(kind=8), parameter, dimension(ntype1) :: arat1 = (/ 6.50d-12, 1.80d-12, 1.20d-13, 1.70d-11, 3.30d-39, 2.50d-14,  &
                                                 9.40d-15, 1.60d-12, 3.70d-12, 1.30d-12, 1.10d-14, 2.20d-13,  &
                                                 1.90d-33, 3.10d-34, 3.30d-12, 2.27d-13, 1.90d-33, 3.10d-34,  &
                                                 4.20d-12, 1.75d-13, 1.60d-11, 2.80d-12, 5.60d-12, 1.40d-12,  &
                                                 4.20d-12, 1.75d-13, 2.00d+16, 1.20d-11, 3.00d-13, 2.40d-12,  &
                                                 1.70d-11, 2.19d-11, 4.20d-12, 1.75d-13, 4.20d-12, 1.75d-13,  &
                                                 4.20d-12, 1.75d-13, 2.15d-12, 1.20d-14, 1.04d-11, 3.70d-12,  &
                                                 2.10d-12, 1.75d-13, 1.50d-11, 7.00d-15, 4.60d-11, 4.00d-12 /)
   real(kind=4), parameter, dimension(ntype1) :: brat1= (/ 120.0, -1370.0,  -2450.0,   150.0,   529.0, -1229.0,          &
                                                778.0,  -942.0,    240.0,   380.0,  -502.0,   619.0,          &
                                                982.0,  2818.0,   -200.0,   771.0,   982.0,  2818.0,          &
                                                180.0,  1000.0,   -110.0, -2518.0,   311.0, -1867.0,          &
                                                180.0,  1000.0, -13542.0,  -745.0, -1427.0, -1710.0,          &
                                              -1232.0,  -709.0,    180.0,  1000.0,   180.0,  1000.0,          &
                                                180.0,  1000.0,    411.0, -2634.0,  -792.0, -2925.0,          &
                                                322.0,  1000.0,    500.0, -1900.0,   230.0,   180.0/)
   integer(kind=4), parameter, dimension(ntype1) :: indx1 = (/  3,  5,  6,  7,  8, 12,  24,  26,  27,  31,  32,  33,     &
                                                     34, 35, 38, 39, 40, 41,  44,  45,  51,  52,  54,  56,    &
                                                     57, 59, 61, 63, 66, 67,  68,  71,  72,  73,  76,  77,    &
                                                     80, 81, 84, 85, 86, 87, 100, 107, 109, 110, 113, 114/)
   integer(kind=4), parameter, dimension(ntype1) :: jord1 = (/  2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2,        &
                                                     3, 3, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2,        &
                                                     2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,        &
                                                     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2/)

! c)Type 2 rate constant: photolysis of no2 and o3 ---> o*sd
! ---   function of zenith angle and altitude
   integer(kind=4), parameter         :: ntype2 = 2
   integer(kind=4), dimension(ntype2) :: indx2 = (/1, 16/)
   integer(kind=4), dimension(ntype2) :: jord2 = (/1, 1/)
! Specify type 2 photolytic rates for k1 & k16 at nzhite elevations
   real(kind=4), parameter, dimension(nzenth, nzhite, ntype2) :: ratz2 = reshape( &
                                       (/ .898e-02, .891e-02, .869e-02, .830e-02, .769e-02, .679e-02,   &
                                          .548e-02, .360e-02, .178e-02, .381e-03, .100e-09, .100e-09,   &
                                          .109e-01, .108e-01, .107e-01, .103e-01, .978e-02, .893e-02,   &
                                          .761e-02, .548e-02, .301e-02, .191e-02, .100e-09, .100e-09,   &
                                          .127e-01, .127e-01, .126e-01, .123e-01, .118e-01, .111e-01,   &
                                          .980e-02, .759e-02, .473e-02, .343e-02, .100e-09, .100e-09,   &
                                          .152e-01, .152e-01, .151e-01, .147e-01, .141e-01, .133e-01,   &
                                          .117e-01, .912e-02, .568e-02, .419e-02, .100e-09, .100e-09,   &
                                          .378e-04, .366e-04, .330e-04, .272e-04, .203e-04, .129e-04,   &
                                          .636e-05, .199e-05, .459e-06, .541e-07, .100e-12, .100e-12,   &
                                          .561e-04, .539e-04, .497e-04, .416e-04, .316e-04, .207e-04,   &
                                          .104e-04, .335e-05, .775e-06, .812e-07, .100e-12, .100e-12,   &
                                          .786e-04, .759e-04, .705e-04, .602e-04, .468e-04, .316e-04,   &
                                          .168e-04, .584e-05, .125e-05, .108e-06, .100e-12, .100e-12,   &
                                          .102e-03, .984e-04, .916e-04, .783e-04, .609e-04, .412e-04,   &
                                          .218e-04, .758e-05, .161e-05, .135e-06, .100e-12, .100e-12/), &
                                          (/nzenth, nzhite, ntype2/) )

! d)Type 4 rate constant: other photolysis reactions; proportional
! ---   to no2 photolysis rates as a function of zenith angle
   integer(kind=4), parameter                          :: ntype4 = 14
   integer(kind=4), parameter, dimension(ntype4)       :: indx4  = (/ 13, 14, 15, 20, 23, 30, 37, 48, 49, 50, 55, 62, 64, 103/)
   real(kind=4), parameter,  dimension(nzenth, ntype4) :: ratz4 = reshape( &
! ---  1- reaction 13: no3 + hv ---> no
               (/ 2.08e+00, 2.09e+00, 2.11e+00, 2.14e+00, 2.20e+00, 2.36e+00, &
                  2.61e+00, 3.08e+00, 3.75e+00, 2.57e+00, 0.00e-01, 0.00e-01, &
! ---  2- reaction 14: no3 + hv ---> no2
                  1.89e+01, 1.89e+01, 1.91e+01, 1.95e+01, 2.01e+01, 2.13e+01, &
                  2.34e+01, 2.74e+01, 3.28e+01, 2.37e+01, 0.00e-01, 0.00e-01, &
! ---  3- reaction 15: o3 + hv ---> o
                  5.19e-02, 5.20e-02, 5.23e-02, 5.26e-02, 5.36e-02, 5.60e-02, &
                  6.01e-02, 6.90e-02, 8.24e-02, 6.10e-02, 0.00e-01, 0.00e-01, &
! ---  4- reaction 20: hono photolysis
                  1.81e-01, 1.81e-01, 1.81e-01, 1.80e-01, 1.79e-01, 1.78e-01, &
                  1.75e-01, 1.70e-01, 1.67e-01, 1.73e-01, 0.00e-01, 0.00e-01, &
! ---  5- reaction 23: hno3 photolysis
                  6.49e-05, 6.41e-05, 6.08e-05, 5.51e-05, 4.77e-05, 3.85e-05, &
                  2.82e-05, 1.80e-05, 1.15e-05, 7.30e-06, 0.00e-01, 0.00e-01, &
! ---  6- reaction 30: hno4 photolysis
                  9.03e-04, 8.95e-04, 8.62e-04, 8.05e-04, 7.27e-04, 6.25e-04, &
                  5.01e-04, 3.64e-04, 2.68e-04, 1.90e-04, 0.00e-01, 0.00e-01, &
! ---  7- reaction 37: h2o2 photolysis
                  8.38e-04, 8.33e-04, 8.14e-04, 7.79e-04, 7.30e-04, 6.62e-04, &
                  5.76e-04, 4.72e-04, 3.97e-04, 3.52e-04, 0.00e-01, 0.00e-01, &
! ---  8- reaction 48: rooh photolysis
                  8.23e-04, 8.18e-04, 8.00e-04, 7.67e-04, 7.21e-04, 6.57e-04, &
                  5.75e-04, 4.76e-04, 4.04e-04, 3.64e-04, 0.00e-01, 0.00e-01, &
! ---  9- reaction 49: hcho + hv ---> radicals
                  3.35e-03, 3.33e-03, 3.23e-03, 3.07e-03, 2.83e-03, 2.51e-03, &
                  2.10e-03, 1.61e-03, 1.24e-03, 9.25e-04, 0.00e-01, 0.00e-01, &
! --- 10- reaction 50: hcho + hv ---> stable
                  4.79e-03, 4.77e-03, 4.71e-03, 4.60e-03, 4.42e-03, 4.17e-03, &
                  3.81e-03, 3.33e-03, 2.96e-03, 2.79e-03, 0.00e-01, 0.00e-01, &
! --- 11- reaction 55: ald2 photolysis
                  6.61e-04, 6.54e-04, 6.22e-04, 5.68e-04, 4.96e-04, 4.05e-04, &
                  2.99e-04, 1.90e-04, 1.20e-04, 7.13e-05, 0.00e-01, 0.00e-01, &
! --- 12- reaction 62: mek photolysis
                  1.76e-04, 1.74e-04, 1.68e-04, 1.56e-04, 1.39e-04, 1.18e-04, &
                  9.29e-05, 6.55e-05, 4.68e-05, 3.28e-05, 0.00e-01, 0.00e-01, &
! --- 13- reaction 64: mgly photolysis
                  1.67e-02, 1.67e-02, 1.68e-02, 1.69e-02, 1.73e-02, 1.78e-02, &
                  1.87e-02, 2.03e-02, 2.21e-02, 1.99e-02, 0.00e-01, 0.00e-01, &
! --- 14- reaction 103: dial photolysis
                  5.89e-02, 5.87e-02, 5.80e-02, 5.68e-02, 5.49e-02, 5.22e-02, &
                  4.84e-02, 4.33e-02, 3.96e-02, 3.93e-02, 0.00e-01, 0.00e-01/), &
                 (/nzenth, ntype4/) )

! e)Type 5 rate constant: temperature and pressure dependent
   integer(kind=4), parameter                          :: ntype5 = 7
   integer(kind=4), parameter, dimension(ntype5)       :: indx5  = (/4, 9, 19, 22, 28, 43, 58/)
   integer(kind=4), parameter, dimension(ntype5)       :: jord5  = (/2, 2, 2, 2, 2, 2, 2/)
   real(kind=4),    parameter, dimension(ntype5)       :: a5     = (/0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.19/)
   real(kind=4),    parameter, dimension(ntype5)       :: b5     = (/8.10e-27, 9.86e-20, 1.93e-24, 2.20e-22, 1.52e-23, 4.48e-23, 6.29e-19/)
   real(kind=4),    parameter, dimension(ntype5)       :: c5     = (/-2.0, -4.3, -2.6, -3.2, -3.2, -3.3, -4.1/)
   real(kind=4),    parameter, dimension(ntype5)       :: d5     = (/2.20e-11, 2.60e-11, 2.60e-10, 4.00e-08, 1.38e-08, 1.50e-12, 4.92e-03/)
   real(kind=4),    parameter, dimension(ntype5)       :: e5     = (/0.0, -0.5, -0.5, -1.3, -1.4, 0.0, -3.6/)
   integer(kind=4), parameter                          :: npm1   = 4
   integer(kind=4), parameter                          :: ntm1   = 5
! unit conversion factor: 1/1013.25 is from mb to atm, 0.01/1013.25 is from Pa to atm
   real(kind=4),    parameter                          :: inv_norm_press = (0.01 / 1013.25)

!  for alkanes rates
   integer(kind=4), parameter :: npres  = 5      ! 5 ranges of pressure
   integer(kind=4), parameter :: ntemp  = 6      ! 6 ranges of temperature
   integer(kind=4), parameter :: ncoeff = 7      ! 7 types of coefficents
   real(kind=4)    :: xx, y, z, xc, yc, zc
   real(kind=4)    :: balka(ncoeff, ntemp, npres)

   save


  contains

!============================================================================!
! Description    : Calculates the clear sky correction factor of Chang et al.,
!                  JGR 92, pp 14, 681-14, 700, 1987
! Description    : Calculates the clear sky correction factor of Chang et al.,
!                  JGR 92, pp 14, 681-14, 700, 1987
!
! Extra info     : corr(x, y, z) = ( 1 + a(x, y) * ( Fcld (x, y, z) - 1) )
!                  a(x, y) is the fraction of the sky which is cloudy
!                  For solar zenith angles (chi) less than 60 degrees,
!                  Fcld = 1 + (1 - tr) * cos (chi) above cloud layer
!                       = 1.4 cos(chi)             in cloud layer
!                       = 1.6 tr * cos(chi)        below cloud layer
!                  For solar zenith angles greater than 60 degrees, the 60 degree
!                  value is used.
!                  "tr" is the energy transmission coefficient for normally incident light
!                  tr= (5 - exp(-tau) ) / (4 + 3*tau*(1 - f) )
!                  f = scattering phase function asymmetry factor = 0.86
!                  tau = cloud optical depth
!                      = sum(i)  3 L(i) dz(i) / (2 rho_h2o r)
!                  That is, the total optical depth is the sum of the optical depths
!                  of each cloudy layer.
!
!                  Important Note: Chang et al's correction is for a single cloud
!                  layer, while the information available from the current
!                  meteorological model allows for multiple cloud layers.  This version
!                  of the correction factor assumes an "average cloud"; summing the
!                  optical depths over height to get the net value, then assuming the
!                  existance of a single cloud over the given cloudy region.  No
!                  attempt in this version of the correction factor to account for
!                  multiple cloud layer effects.
! Arguments:  IN
!               zmom -> Momentum level heights (above-ground height)
!               csza -> Cosine of solar zenith angle
!               cldf -> Cloud fraction within a given layer
!               cldw -> Average liquid water content in gridsquare, kg/kg
!
!             IN/OUT
!                jc  -> Correction factor for J values
!============================================================================

  subroutine mach_adom2_jcorr (zmom, cldw, cldf, csza, jc, gni, gnk)
   use chm_consphychm_mod, only: rho_h2o
   implicit none

   integer(kind=4), intent (in) :: gni, gnk
   real(kind=4),    intent (in) :: zmom(gni, gnk)
   real(kind=4),    intent (in) :: cldw(gni, gnk)
   real(kind=4),    intent (in) :: cldf(gni, gnk)
   real(kind=4),    intent (in) :: csza(gni)
   real(kind=4),    intent(out) :: jc  (gni, gnk)
!
! Local variables
!
   real(kind=4)            :: den
   real(kind=4), parameter :: rdrop   = 1.0e-05   !radius of typical cloud drop, m.
   real(kind=4), parameter :: lconmin = 1.0e-04   !minimum liquid water content, g/m3,
   real(kind=4), parameter :: alim    = 1.0e-09   !smallest allowable cloud fraction before layer can be considered to be "cloudless"
   integer(kind=4)         :: i, k
   integer(kind=4)         :: kmin(gni), kmax(gni)
   integer(kind=4)         :: kwat(gni, gnk)
   real(kind=4)            :: chim(gni), tau(gni)
   real(kind=4)            :: fbelow(gni, gnk), fin(gni, gnk), fabove(gni, gnk)
   real(kind=4)            :: aip(gni, gnk), afrac(gni, gnk)

   den = 1.0 / (2.0 * rho_h2o * rdrop)
!  maximum solar zenith angle will be less than or equal to 60 degrees
!  (1.047 radians); cos(chi) > 0.5
   chim = max(0.5, csza)

! (2) avoid division by zero error for cloud fraction:
!  afrac is 1 if the area fraction is greater than the minimum
!  allowed, zero otherwise.
   do k = 1, gnk
      do i = 1, gni
         aip(i, k) = 1.0 / max(cldf(i, k), alim)
         afrac(i, k) = max(sign(1.0, cldf(i, k) - alim), 0.0)
      end do
   end do
!  calculate total cloud optical depth:
   do i = 1, gni
      tau(i) = 0.0
   end do
   do k = 1, gnk
      do i = 1, gni
!  jc temporarily used to store cloud liquid water content
!  in g/m3. multiply by zero if area fraction of cloud is
!  less than lower limit.
         jc(i, k) = max(cldw(i, k) * afrac(i, k), 0.0)
      end do
!  Liquid water content was average for entire grid square.
!  Divide this by the area fraction to get the liquid water
!  content, assuming that all of the liquid water is in
!  clouds.  Note that ami was set earlier to ensure that
!  division by zero error does not occur.  Presumably, jc should
!  be zero if nu is zero anyway; this is just in case there's
!  different round off levels for the two variables.
      do i = 1, gni
         jc(i, k) = jc(i, k) * aip(i, k)
      end do
   end do

!  calculate the total optical depth in the vertical direction
!  due to clouds.  Note that the momentum height surfaces are
!  assumed to be the boundaries of each model layer, except for
!  the ground level.

!  all layers except surface
   do k = 1, gnk-1
      do i = 1, gni
         tau(i) = tau(i) + 3.0 * jc(i, k) * (zmom(i, k) - zmom(i, k+1)) * den
      end do
   end do
!  bottom layer:
   do i = 1, gni
         tau(i) = tau(i) + 3.0 * jc(i, gnk) * zmom(i, gnk) * den
   end do

!  calculate cloud transmission coefficient; replace optical
!  depth with transmission coefficient:
   do i = 1, gni
      tau(i) = (5.0 - exp(-tau(i))) / (4.0 + 0.42 * tau(i))
   end do

!  determine maximum, minimum layer indices for which
!  liquid water is present.  kwat will be the 1.0 value
!  if cloud exists in the given model layer, 0 if no cloud
!  exists in the layer.
   do k = 1, gnk
      do i = 1, gni
         kwat(i, k) = int(max(sign(1.0, jc(i, k) - lconmin), 0.0))
      end do
   end do
   do i = 1, gni
      kmin(i) = 0
      kmax(i) = 0
   end do

!  again, code is changed due to change in direction
!  of vertical index:
   do k = 1, gnk
      do i = 1, gni
      kmin(i) = max(kmin(i), k * kwat(i, k))
      kmax(i) = max(kmax(i), (gnk - k + 1) * kwat(i, k))
      end do
   end do
   do i=1, gni
      kmax(i) = (gnk - kmax(i) + 1)
   end do
!  kmin(i) now contains the vertical index of the lowest altitude
!  layer which contains cloud, kmax(i) contains the vertical index of
!  the highest altitude layer index which contains cloud.
!
!  we now have a "cloud", extending from
!  layer kmin(i) up to layer kmax(i), with transmission coefficient
!  tau(i).  calculate fcld parameter (in jc array, temporarily):
   jc = 1.0

!  calculate below, within, above cloud values:
   do k = 1, gnk
      do i = 1, gni
         fbelow(i, k) = 1.6 * tau(i) * chim(i)
         fin(i, k) = 1.4 * chim(i)
         fabove(i, k) = 1.0 + (1.0 - tau(i)) * chim(i)
      end do
   end do

!  assign these to the appropriate layers:
!  note modification for layer structure:
   do i = 1, gni
      do k = kmin(i) + 1, gnk
         jc(i, k) = fbelow(i, k)
      end do
      do k = kmax(i), kmin(i)
         jc(i, k) = fin(i, k)
      end do
      do k = 1, kmax(i) - 1
         jc(i, k) = fabove(i, k)
      end do
   end do
!  calculate cloud fraction, place in tau array.  the maximum value
!  of any given layer is used as the area fraction for the entire
!  average cloud.  this is a conservative estimate, making the
!  assumption that all layers line up vertically.
   tau = 0.0
   do k = 1, gnk
      do i = 1, gni
         tau(i) = max(cldf(i, k), tau(i))
      end do
   end do

!  calculate the final clear sky correction factor:
   do k = 1, gnk
      do i = 1, gni
         jc(i, k) = 1.0 + tau(i) * (jc(i, k) - 1.0)
      end do
   end do

   return
  end subroutine mach_adom2_jcorr

!============================================================================
! Arguments:  IN
!
!                 zmount   -> layer midpoint elevation - height above ground (m)
!                 skyo     -> sky clearness ratio (0. to 1.)
!                 csza     -> cosine solar zenith angle
!                 gni      -> number of grid point along X
!                 gnk      -> number of vertical levels
!             OUT
!                 rg       -> photolysis reaction rates
!============================================================================
  subroutine mach_gas_lookup_jvals(zmount, skyo, csza, rg, gni, gnk, npt)
   use chm_consphychm_mod, only: pi
   implicit none
!
   integer(kind=4), intent (in) :: gni, gnk, npt
   real(kind=4),    intent (in) :: zmount(gni, gnk)
   real(kind=4),    intent (in) :: skyo  (gni, gnk)
   real(kind=4),    intent (in) :: csza  (gni)
   real(kind=8),    intent(out) :: rg    (npt, ntype2+ntype4)
!
! Local variables
!
   integer(kind=4) :: nn, i, k, jpt, iz1
   real(kind=4)    :: riz1(npt), riz2(npt), sky(npt), znode(npt), zen(npt)
   integer(kind=4) :: iza1(npt), iza2(npt), izen1(npt), izen2(npt)
!  Air concentration (constant) in molecules/cc (=> a(:, 1) in mach_adom2_uprate
   real(kind=4), parameter :: air_conc = 60.0
!------------------------------------------------------------------
!-----type 2 rate constants

! --- retrieve photolytic rates stored as a function of
! --- solar zenith angle and elevation.
!------------------------------------------------------------------------
!
   rg = 0.0

   nn = ((nk_start - 1) * gni) + 1

   jpt = nn
   do k = nk_start, gnk
      do i = 1, gni
         sky(jpt) = skyo(i, k)
         zen(jpt) = acos(csza(i)) * 180.0 / pi ! convert radians to degrees
         znode(jpt) = zmount(i, k)
         jpt = jpt + 1
      end do
   end do

   do iz1 = 1, nzhite - 1
      do jpt = nn, npt
         if(znode(jpt) >= zhite(iz1) .and. znode(jpt) < zhite(iz1 + 1)) then
            iza1(jpt) = iz1
            iza2(jpt) = iz1 + 1
         end if
      end do
   end do

!  use top of lookup table as upper cutoff on heights ;
!  prevent extrapolation (and code errors) above lookup table
!  top.

   do jpt = nn, npt
      if(znode(jpt) >= zhite(nzhite)) then
         iza1(jpt) = nzhite
         iza2(jpt) = nzhite
      end if
   end do

!------------------------------------------------------------------------
! --- determine upper and lower bounds of zenith angle
!-------------------------------------------------------------------------

   do iz1 = 1, nzenth - 1
      do jpt = nn, npt
         if(zen(jpt) >= zenith(iz1) .and. zen(jpt) <= zenith(iz1+1)) then
            izen1(jpt) = iz1
            izen2(jpt) = iz1 + 1
         end if
      end do
   end do

! --- interpolate absolute rates to zenith angle

   do i = 1, ntype2
      do jpt = nn, npt
         riz1(jpt) = ratz2(izen1(jpt), iza1(jpt), i) +   &
                     (ratz2(izen2(jpt), iza1(jpt), i) -  &
                     ratz2(izen1(jpt), iza1(jpt), i)) *  &
                     (zen(jpt) - zenith(izen1(jpt))) /   &
                     (zenith(izen2(jpt)) - zenith(izen1(jpt)))

         riz2(jpt) = ratz2(izen1(jpt), iza2(jpt), i) +   &
                     (ratz2(izen2(jpt), iza2(jpt), i) -  &
                     ratz2(izen1(jpt), iza2(jpt), i)) *  &
                     (zen(jpt) - zenith(izen1(jpt))) /   &
                     (zenith(izen2(jpt)) - zenith(izen1(jpt)))
      end do

! --- interpolate to node elevation & adjust for uv attenuation

      do jpt = nn, npt

         if(znode(jpt) < zhite(nzhite) .and. znode(jpt) >= zhite(1) ) then
!  height of model gridpoint is less than upper limit in lookup table;

            rg(jpt, i) =  dble(air_conc * sky(jpt) *              &
                          (riz1(jpt) + (riz2(jpt) - riz1(jpt)) *  &
                          (znode(jpt) - zhite(iza1(jpt))) /       &
                          (zhite(iza2(jpt)) - zhite(iza1(jpt)))))
         else
!  height of model gridpoint is greater than upper limit in lookup table;
!  use uppermost lookup table value for all model gridpoint heights greater
!  than top of lookup table:
            rg(jpt, i) = dble(air_conc * sky(jpt) * riz1(jpt))
         end if
      end do
   end do
!
!------- type 4 rate constants------------------------------------------
!
!       proportional to type2 based on zenith angle
!       all type 4 are 1st order
!       all type 4 are proportional to rate for reaction no. 1
! --- interpolate proportionality factor to the zenith angle
!-----------------------------------------------------------------------
   do i = 1, ntype4
      k = i + ntype2
      do jpt = nn, npt
         rg(jpt, k) = rg(jpt, 1) * dble((ratz4(izen1(jpt), i) +         &
                        (ratz4(izen2(jpt), i) - ratz4(izen1(jpt), i)) * &
                        (zen(jpt) - zenith(izen1(jpt))) /               &
                        (zenith(izen2(jpt)) - zenith(izen1(jpt)))) )
      end do
   end do

   return
  end subroutine mach_gas_lookup_jvals

!============================================================================
!-------------------------------------------------------------------------
!
! --- calculate pressure/temp. dependent alkane coefficients
! --- prevent extrapolation outside pressure and temperature ranges
!     px = min(atpres, pres(npres)) and px = max(px, pres(1))
!     tx = min(t, temp(ntemp)) and tx = max(tx, temp(1))
!-------------------------------------------------------------------------
! Arguments:  IN
!
!                 p2d      -> atmospheric pressuure (Pa)
!                 tp       -> temperature (deg k)
!                 npt      -> number of grid points (passed through argument list of subr. call)
!                 gni      -> number of grid point along X
!             OUT
!                 bcs      -> product coefficients
!============================================================================
  subroutine mach_adom2_bcs_coeffs(tplus, p2d, bcs, npt, gni, gnk)
   implicit none
!
   integer(kind=4), intent (in) :: npt
   integer(kind=4), intent (in) :: gni, gnk
   real(kind=4),    intent (in) :: tplus(gni, gnk)
   real(kind=4),    intent (in) :: p2d  (gni, gnk)
   real(kind=8),    intent(out) :: bcs  (npt, nprcf)
!
! Local variables
!

   real(kind=4), dimension(npres) :: pres = (/ 0.2,  0.4,  0.6,  0.8,  1.0 /)
   real(kind=4), dimension(ntemp) :: temp = (/ 210.0, 230.0, 250.0, 270.0, 290.0, 310.0 /)

   real(kind=4)    :: bt1(npt, ncoeff), bt2(npt, ncoeff)
   real(kind=4)    :: px(npt), pf(npt), tx(npt), tf(npt), atpres

   integer(kind=4) :: kp1(npt), kp2(npt), jt1(npt), jt2(npt), isw(npt)
   integer(kind=4) :: nn, i, k, jpt, kp, jt

   bcs = 0.0d0

   nn = ((nk_start - 1) * gni) + 1

   jpt = nn
   do k = nk_start, gnk
      do i = 1, gni
         atpres = p2d(i, k) * inv_norm_press

         px(jpt) = max(min(atpres, pres(npres)), pres(1))
         tx(jpt) = max(min(tplus(i, k), temp(ntemp)), temp(1))
         jpt = jpt + 1
      end do
   end do

!-------------------------------------------------------------------------
! --- find pressure indices
!-------------------------------------------------------------------------

   isw(1:npt) = 0

   do kp = npm1, 1, -1
      do jpt = nn, npt
         if (isw(jpt) /= 1) kp1(jpt) = kp
      end do
      do jpt = nn, npt
         if (px(jpt) >= pres(kp)) then
            isw(jpt) = 1
         else
            isw(jpt) = 0
         end if
      end do
   end do
   do jpt = nn, npt
      kp1(jpt) = max(kp1(jpt), 1)
      kp2(jpt) = kp1(jpt) + 1
      pf(jpt) = (px(jpt) - pres(kp1(jpt))) /  &
               (pres(kp2(jpt)) - pres(kp1(jpt)))
   end do

!--------------------------------------------------------------------------
! --- find temperature indices
!-------------------------------------------------------------------------

   isw(1:npt) = 0

   do jt = ntm1, 1, -1
      do jpt = nn, npt
         if (isw(jpt) /= 1) jt1(jpt) = jt
      end do
      do jpt = nn, npt
         if (tx(jpt) >= temp(jt)) then
            isw(jpt) = 1
         else
            isw(jpt) = 0
         end if
      end do
   end do

   do jpt = nn, npt
      jt1(jpt) = max(jt1(jpt), 1)
      jt2(jpt) = jt1(jpt) + 1
      tf(jpt) = (tx(jpt) - temp(jt1(jpt))) / (temp(jt2(jpt)) - temp(jt1(jpt)))
   end do

! --- interpolate product coefficient data to desired pressure
   do i = 1, ncoeff - 1
      do jpt = nn, npt
         bt1(jpt, i) = balka(i, jt1(jpt), kp1(jpt)) +                     &
                       pf(jpt) * (balka(i, jt1(jpt), kp2(jpt)) - balka(i, jt1(jpt), kp1(jpt)))

         bt2(jpt, i) = balka(i, jt2(jpt), kp1(jpt)) +                     &
                       pf(jpt) * (balka(i, jt2(jpt), kp2(jpt)) - balka(i, jt2(jpt), kp1(jpt)))
      end do

!     interpolate product coefficient data to desired temperature
      do jpt = nn, npt
         bcs(jpt, i) = dble(bt1(jpt, i) + tf(jpt) * (bt2(jpt, i) - bt1(jpt, i)))
      end do
   end do

!  remaining coefficients (alkenes and aromatics) are constant and
   do jpt = nn, npt
      bcs(jpt, 08) = dble(y)
      bcs(jpt, 09) = dble(y + 2.00 * yc)
      bcs(jpt, 10) = dble(0.64 * y)
      bcs(jpt, 11) = dble(0.50 * y + yc)
      bcs(jpt, 12) = dble(0.13 * y + 0.27 * yc)
      bcs(jpt, 13) = dble(0.17 * y + 0.21 * yc)
      bcs(jpt, 14) = dble(0.06 * y + 0.12 * yc)
      bcs(jpt, 15) = dble(0.28 * y)
      bcs(jpt, 16) = dble(0.40 * y)
      bcs(jpt, 17) = dble(yc)
      bcs(jpt, 18) = dble(0.40 * y)
      bcs(jpt, 19) = dble(0.20 * y)
      bcs(jpt, 20) = dble(0.20 * y + 0.40 * yc)
      bcs(jpt, 21) = dble(0.60 * y)
      bcs(jpt, 25) = dble(0.20 * y)
      bcs(jpt, 26) = dble(0.10 * y + 0.20 * yc)
! --- aromatics
      bcs(jpt, 22) = dble(0.65 * z + 0.49 * zc)
      bcs(jpt, 23) = dble(0.316 * z + 0.86 * zc)
      bcs(jpt, 24) = dble(0.095 * z)
   end do

   do jpt = nn, npt
      bcs(jpt, ncoeff) = bcs(jpt, 4) + bcs(jpt, 5) + bcs(jpt, 6)
   end do

   return
  end subroutine mach_adom2_bcs_coeffs
!
!============================================================================
! Arguments:  IN
!
!                 p2d      -> atmospheric pressuure (Pa)
!                 tp       -> temperature (deg k)
!                 cno      -> no concentration (ppm)
!                 rg_jvals -> photolysis reaction rates
!                 npt      -> number of grid points (passed through argument list of subr. call)
!                 gni      -> number of grid point along X
!             OUT
!                 bg     -> product coefficients
!          IN/OUT
!                 rg     -> reaction rates
!============================================================================
  subroutine mach_adom2_uprate(tp, p2d, cno, rg_jvals, rg, bg, npt, gni, gnk)
   implicit none
!
   integer(kind=4), intent (in) :: npt
   integer(kind=4), intent (in) :: gni, gnk
   real(kind=4),    intent (in) :: tp      (gni, gnk)
   real(kind=4),    intent (in) :: p2d     (gni, gnk)
   real(kind=4),    intent (in) :: cno     (npt)
   real(kind=8),    intent (in) :: rg_jvals(npt, ntype2+ntype4)
   real(kind=8),    intent(out) :: rg      (npt, nreac)
   real(kind=8),    intent(out) :: bg      (npt, nprcf)
!
! Local variables
!
   integer(kind=4) :: nn, i, j, k, jpt
   real(kind=8)    :: ttmp
   real(kind=4)    :: a(npt, 3)
   real(kind=4)    :: airm(npt), atpres(npt), ee(npt)
   real(kind=4)    :: t(npt), ti(npt), pt
   real(kind=4)    :: rkif(npt), rkm(npt), rkmok(npt), rk0(npt)

!!!#  TODO transformer tableau en vrai 2d et eliminer npt et remplacer par (chm_ni, chm_nk)...

   rg = 0.0d0

   nn = ((nk_start - 1) * gni) + 1

   jpt = nn
   do k = nk_start, gnk
      do i = 1, gni
         t(jpt) = tp(i, k)
         atpres(jpt) = p2d(i, k) * inv_norm_press
         ti(jpt) = 1.0 / t(jpt)
         pt = atpres(jpt) * ti(jpt)

!     compute concentration of air in molecules/cc
         airm(jpt) = 7.3387e21 * pt
         a(jpt, 1) = 60.0
         a(jpt, 2) = 4.40e17  * pt
         a(jpt, 3) = 3.23e33  * (pt * pt)

         jpt = jpt + 1
      end do
   end do


!__type 0 rate constants: k =  arat0 (in cc/(molecule-sec)
!  convert rate constants to ppm-min units

   do i = 1, ntype0
      j = jord0(i)
      k = indx0(i)

      do jpt = nn, npt
         rg(jpt, k) = dble(a(jpt, j) * arat0(i))
      end do
   end do

   do i = 1, ntype1
      j = jord1(i)
      k = indx1(i)
      do jpt = nn, npt
         rg(jpt, k) = dble(a(jpt, j)) * arat1(i) * dble(exp(brat1(i) * ti(jpt)))
      end do
   end do

!---- TYPE 2 and TYPE 4 are photolysis rates skip if MESSY is invoked
   if (.not. chm_messy_jval_l) then

      do i = 1, ntype2
         k = indx2(i)
         rg(:, k) = rg_jvals(:, i)
      end do
      do i = 1, ntype4
         k = indx4(i)
         rg(:, k) = rg_jvals(:, i+ntype2)
      end do

   end if
!-----------------------------------------------------------------------
!-------- type 5 rate constants (temperature & pressure dependent)
!
!               b*(t**c)m                       1
!     type 5 = ----------- a**ee   , ee = -------------------
!              1+b(t**c)m                           b(t**c)m
!                --------              1 + ( log10( -------- ) )**2
!                 d(t**e)                          d(t**e)
!------------------------------------------------------------------------

   do i = 1, ntype5
      k = indx5(i)
      j = jord5(i)

      do jpt = nn, npt
         rk0(jpt) = b5(i) * (t(jpt) ** c5(i))
      end do
      do jpt = nn, npt
         rkif(jpt) = d5(i) * (t(jpt) ** e5(i))
      end do
      do jpt = nn, npt
         rkm(jpt) = rk0(jpt) * airm(jpt)
      end do
      do jpt = nn, npt
         rkmok(jpt) = rkm(jpt) / rkif(jpt)
      end do
      do jpt = nn, npt
         ee(jpt) = 1.0 / (1.0 + (alog10(rkmok(jpt))) ** 2)
      end do
      do jpt = nn, npt
         rg(jpt, k) = dble(a(jpt, j) * (rkm(jpt) / (1.0 + rkmok(jpt))) *  &
                     (a5(i)) ** ee(jpt) )
      end do
   end do

!------------------------------------------------------------------------
!
!------- type 6 (special functions)

   do jpt = nn, npt
      rg(jpt, 2)  = dble(a(jpt, 3) * 3.0e-28 / (t(jpt) ** 2.3))
      rg(jpt, 10) = dble(a(jpt, 1)) * rg(jpt, 9) * dble( 9.09e26 *   &
                    exp(-11200.0 * ti(jpt)) / a(jpt, 2))
      rg(jpt, 25) = dble(a(jpt, 2) * 1.50e-13 * (1.0 + 0.6 * atpres(jpt)))
   end do

   do jpt = nn, npt
      rg(jpt, 29) = dble(a(jpt, 1)) * rg(jpt, 28) * dble( 4.76e26 *  &
                    exp(-10940.0 * ti(jpt)) / a(jpt, 2))

!      rg(jpt, 36) =  a(jpt, 3)*2.7e-54*airm(jpt)*exp(3137*ti(jpt))

      rg(jpt, 36) = dble(a(jpt, 3) * 1.6432e-27 * airm(jpt) * exp(3137.0 * ti(jpt)))
      rg(jpt, 36) = rg(jpt, 36) * 1.6432d-27
      rg(jpt, 42) = rg(jpt, 36)
   end do

   do jpt = nn, npt
!  the original formulation for this rate allows for a difference term
!  severely prone to round-off error.  a revised form will be used
!  here.  p.a. makar, march 2007
!         rg(jpt, 53) = dble(a(jpt, 2) * 1.1e-13 * (1.0 - 20.0 /  &
!                     (20.0 + 4.2e-18 * exp(180.0 / t(jpt)) *  &
!                      cno(jpt) * airm(jpt))) )
      ttmp = 4.2d-18 * exp(180.0d0 * dble(ti(jpt))) * dble(cno(jpt)) * dble(airm(jpt))
      rg(jpt, 53) = dble(a(jpt, 2)) * 1.1d-13 * (ttmp / (20.d0 + ttmp))

      rg(jpt, 69) = dble(a(jpt, 2) * 1.27e-17 * (t(jpt) ** 2) * exp(14.0 * ti(jpt)))

! --- alka + oh    reaction rate
      rg(jpt, 70) = dble(a(jpt, 2) * (1.017e-11 * exp(-354.0 * ti(jpt)) * xx +  &
                     2.312e-11 * exp(-289.0 * ti(jpt)) * xc))

! --- alke + oh    reaction rate
      rg(jpt, 88) = dble(a(jpt, 2) * (5.323e-12 * exp(504.0 * ti(jpt)) * y +    &
                     1.074e-11 * exp(549.0 * ti(jpt)) * yc))

! --- alke + o3    reaction rate
      rg(jpt, 89) = dble(a(jpt, 2) * (1.323e-14 * exp(-2105.0 * ti(jpt)) * y +  &
                     7.333e-15 * exp(-1137.0 * ti(jpt)) * yc))
      rg(jpt, 90) = dble(a(jpt, 2) * (1.18e-11 * exp(-324.0 * ti(jpt)) * y +    &
                     2.26e-11 * exp(10.0 * ti(jpt)) * yc))

! --- alke + no3    reaction rate
      rg(jpt, 91) = dble(a(jpt, 2) * (1.143e-11 * exp(-1935.0 * ti(jpt)) * y +  &
                     3.23e-11 * exp(-975.0 * ti(jpt)) * yc))

! --- arom + oh    reaction rate
      rg(jpt, 101) = dble(a(jpt, 2) *  &
                     (1.407e-11 * exp(116.0 * ti(jpt)) * z + 4.77e-11 * zc))
   end do

! Evaluate the bg product coefficients
   call mach_adom2_bcs_coeffs(tp, p2d, bg, npt, gni, gnk)
!
   return
  end subroutine mach_adom2_uprate


!============================================================================!
!
! Name           : mach_adom2_chemi.ftn90
!
! Description    : Initialize Chemistry Variables. !
!
! Extra info     : Set of values used was obtained
!                  from adom model, version: adomiib(napap)   level: 10/21/89
!
!============================================================================
   subroutine mach_adom2_chemi()
   implicit none
!
!  Local variables
!
   integer (kind=4):: i, j, k
   real (kind=4)   :: cpar, cpar1, cpar2, cole, cole1, cole2, caro, caro1, caro2

   real(kind=4), dimension(ncoeff, ntemp, npres) :: balk4 = reshape((/                    &
               .2253, .603, .487, .060, .940, .434, 1.374, .2179, .574, .518,             &
               .053, .947, .444, 1.391, .2099, .554, .551, .046, .954, .489, 1.443,       &
               .2028, .597, .528, .040, .960, .570, 1.530, .1958, .791, .406, .034, .966, &
               .710, 1.676, .1922, .988, .269, .028, .972, .845, 1.817,                   &
               .2222, .590, .450, .104, .896, .422, 1.318, .2154, .563, .475, .091, .909, &
               .422, 1.331, .2078, .541, .519, .078, .922, .449, 1.372,                   &
               .2006, .552, .524, .065, .935, .514, 1.449, .1939, .683, .453, .052, .948, &
               .630, 1.578, .1887, .898, .314, .042, .958, .781, 1.739,                   &
               .2196, .580, .421, .138, .862, .414, 1.275, .2132, .555, .448, .120, .880, &
               .410, 1.291, .2060, .531, .491, .101, .899, .427, 1.326,                   &
               .1985, .534, .515, .082, .918, .482, 1.400, .1925, .629, .472, .065, .935, &
               .586, 1.521, .1869, .836, .345, .051, .949, .736, 1.685,                   &
               .2172, .571, .398, .166, .834, .406, 1.240, .2112, .549, .427, .143, .857, &
               .403, 1.259, .2045, .526, .471, .119, .881, .415, 1.296,                   &
               .1974, .523, .507, .096, .904, .460, 1.365, .1915, .597, .480, .075, .925, &
               .555, 1.481, .1859, .790, .367, .058, .942, .702, 1.644,                   &
               .2151, .565, .379, .190, .810, .401, 1.212, .2095, .544, .410, .163, .837, &
               .397, 1.234, .2032, .522, .454, .134, .866, .406, 1.272,                   &
               .1965, .516, .498, .107, .893, .446, 1.339, .1906, .575, .484, .082, .918, &
               .533, 1.450, .1852, .753, .384, .064, .936, .674, 1.610/) , (/ncoeff, ntemp, npres/))

   real(kind=4), dimension(ncoeff, ntemp, npres) :: balk7 = reshape(                      &
            (/.0002, .318, .769, .272, .728, .600, 1.329, .0016, .320, .849,              &
               .231, .769, .660, 1.428, .0066, .331, .917, .189, .811, .750, 1.561,       &
               .0166, .365, .958, .151, .849, .851, 1.700, .0319, .439, .953, .120, .880, &
               .959, 1.839, .0523, .522, .934, .095, .905, 1.056, 1.961,                  &
               .0001, .302, .665, .347, .653, .544, 1.197, .0007, .307, .755, .296, .704, &
               .595, 1.299, .0040, .317, .841, .245, .755, .675, 1.430,                   &
               .0109, .344, .901, .198, .802, .774, 1.576, .0245, .402, .920, .158, .842, &
               .883, 1.725, .0396, .484, .911, .125, .875, .988, 1.863,                   &
               .0000, .291, .609, .388, .612, .514, 1.126, .0004, .299, .701, .334, .666, &
               .562, 1.228, .0028, .310, .794, .279, .721, .633, 1.354,                   &
               .0083, .333, .864, .228, .772, .728, 1.501, .0192, .383, .897, .183, .817, &
               .836, 1.653, .0351, .462, .898, .144, .856, .949, 1.804,                   &
               .0000, .284, .571, .416, .584, .494, 1.078, .0003, .292, .664, .361, .639, &
               .541, 1.180, .0021, .304, .760, .304, .696, .606, 1.302,                   &
               .0067, .326, .836, .250, .750, .697, 1.447, .0166, .371, .879, .201, .799, &
               .805, 1.603, .0320, .446, .889, .158, .842, .920, 1.762,                   &
               .0000, .278, .543, .438, .562, .480, 1.043, .0002, .287, .636, .381, .619, &
               .525, 1.143, .0017, .300, .733, .323, .677, .586, 1.262,                   &
               .0056, .320, .814, .267, .733, .673, 1.406, .0147, .362, .865, .215, .784, &
               .780, 1.565, .0297, .434, .882, .169, .831, .898, 1.729/), (/ncoeff, ntemp, npres/))
!
! --- stoichiometric coefficients

!     cpar = c4-c5 fraction of >c3 alkanes on a ppmc basis
!     cole = terminal alkene fraction of  >c2 alkenes on a ppmc basis
!     caro = di-alkylbenzene  fraction of di- & tri-alkylbenzenes (ppmc)

! --- alkanes (paraffins)
   cpar = 0.375
   cpar1 = cpar / 4.3881163
   cpar2 = (1.0 - cpar) / 7.82713536
   xx = cpar1 / (cpar1 + cpar2)
   xc = 1.0 - xx

! --- alkenes (olefins)
   cole = 0.500
   cole1 = cole / 3.5
   cole2 = (1.0 - cole) / 4.93
   y = cole1 / (cole1 + cole2)
   yc = 1.0 - y

! --- aromatics
   caro = 0.685
   caro1 = caro / 8.7166478
   caro2 = (1.0 - caro) / 9.1243209
   z = caro1 / (caro1 + caro2)
   zc = 1.0 - z

! --- alkanes
   do k = 1, npres
      do j = 1, ntemp
         do i = 1, ncoeff
            balka(i, j, k) = xx * balk4(i, j, k) + xc * balk7(i, j, k)
         end do
      end do
   end do

   return
   end subroutine mach_adom2_chemi

end module mach_adom2_rates_mod
