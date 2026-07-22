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
! Fichier/File   : mach_hetv_hetvcall.ftn90
! Creation       : V. Bouchet, P. Makar, S. Menard, S. Gravel, B. Pabla
! Description    : Interface between AURAMS and hetv
!                  Hetv solves the heterogeneous system for NH4-NO3-SO4, vectorizing
!                  over the gridpoint dimension and utilizing the systems of
!                  equations set out in Nenes and Pandis 1998 ISORROPIA code.
! Extra info     : Split HETV cases into 2 broad areas, based on the metastable
!                  state option (chm_hetchem_metstlb_l).
!                  - (S. Gravel and A. Akingunola, August 2016)
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_hetvcall(npts, ghetf, hetf, t_i, rh_i, rho_i, jlat, kount)
   use mach_hetv_mod,         only: maxhet, maxghet
   use mach_cam_utils_mod,    only: isize
!!if_off
   use mach_hetv_headers_mod, only: mach_hetv_main_2cases, mach_hetv_main_12cases
   use chm_nml_mod,           only: chm_hetchem_metstlb_l
   use mach_hetv_mod,         only: lolimit
   use chm_utils_mod,         only: chm_error_l

   implicit none
!!if_on
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: jlat
   integer(kind=4), intent   (in) :: kount
   real(kind=8),    intent   (in) :: t_i   (npts)
   real(kind=8),    intent   (in) :: rh_i  (npts)
   real(kind=8),    intent   (in) :: rho_i (npts)
   real(kind=8),    intent(inout) :: ghetf (npts, maxghet)
   real(kind=8),    intent(inout) :: hetf  (npts, maxhet, isize)
!!if_off
!
! Local variables
!
   integer(kind=4)              :: i, l, n
   integer(kind=4), parameter   :: isosp = 14
   real(kind=8), parameter, dimension(11) :: atomic_weight = (/&
!           SO4         NO3        NH4       H2SO4       HNO3         NH3
      96.0636d0,  62.0049d0,  18.0385d0,  98.0795d0,  63.0128d0, 17.03056d0, &
!    (NH4)2SO4       NH4NO3    NH4HSO4  (NH4)3H(SO4)2    HSO4
     132.1406d0, 80.04344d0, 115.1104d0, 247.2506d0,  97.0715d0 /)
!

   real(kind=8)                          :: zeros = 0.0d0
!  hetv specific declaration
   real(kind=8), dimension(npts, maxghet):: massb, orgmass, gasin, gasout
   real(kind=8), dimension(npts, maxhet) :: solin, solout
   real(kind=8), dimension(isosp, npts)  :: outp
!  Makar:  test:
   real(kind=8), dimension(isosp, npts)  :: inp
   real(kind=4), dimension(npts)         :: case_number
   real(kind=8), dimension(npts)         :: so4_i, no3_i, nh4_i, hso4_i, amsul_i,    &
                                            ambis_i, amnit_i, leto_i, hno3_i, nh3_i, &
                                            h_i, lwn_i

!  convert to hetv input (in moles/kg of air from kg/kg of air)
!  (pressure expected in atm, convert from mb)

   so4_i   = 0.0d0
   no3_i   = 0.0d0
   nh4_i   = 0.0d0
   hso4_i  = 0.0d0
   hno3_i  = 0.0d0
   nh3_i   = 0.0d0
   h_i     = 0.0d0

   do n = 1, isize
      do i = 1, npts
         so4_i(i) = so4_i(i) + (hetf(i, 1, n) / atomic_weight(1))  ! SO4
         no3_i(i) = no3_i(i) + (hetf(i, 2, n) / atomic_weight(2))  ! NO3
         nh4_i(i) = nh4_i(i) + (hetf(i, 3, n) / atomic_weight(3))  ! NH4
      end do
   end do

   do i = 1, npts
!  SO4
      so4_i(i) = so4_i(i) + (ghetf(i, 1) / atomic_weight(4))
      so4_i(i) = so4_i(i) * 1.0d3
!  NO3 & HNO3
      hno3_i(i) = hno3_i(i) + (ghetf(i, 2) / atomic_weight(5))
      hno3_i(i) = hno3_i(i) * 1.0d3
      no3_i(i)  = no3_i(i)  * 1.0d3
!  NH4 & NH3
      nh3_i(i) = nh3_i(i) + (ghetf(i, 3) / atomic_weight(6))
      nh3_i(i) = nh3_i(i) * 1.0d3
      nh4_i(i) = nh4_i(i) * 1.0d3
   end do

!  Na & Cl (not included in hetv yet - May 2002)
!  temp(K), rh (relative humidity on a 0 to 1 scale)

!  Makar, test:  save input values for later checking:
   do i = 1, npts
      inp(1, i)  = zeros  !(lwn_i)                            !h2o(kg water/kg air)
      inp(2, i)  = h_i(i)     * 1.0d-3                        !h+(aq)
      inp(3, i)  = nh4_i(i)   * atomic_weight(3)  * 1.0d-3    !nh4+(aq)
      inp(4, i)  = so4_i(i)   * atomic_weight(1)  * 1.0d-3    !so4(aq)
      inp(5, i)  = hso4_i(i)  * atomic_weight(11) * 1.0d-3    !hso4(aq)
      inp(6, i)  = no3_i(i)   * atomic_weight(2)  * 1.0d-3    !no3(aq)
      inp(7, i)  = zeros  !(amsul_i)                          !(nh4)2so4(s)
      inp(8, i)  = zeros  !(amnit_i)                          !nh4no3(s)
      inp(9, i)  = zeros  !(notincl)                          !h2so4(aq)
      inp(10, i) = zeros  !(ambis_i)                          !nh4hso4(s)
      inp(11, i) = zeros  !(leto_i)                           !(nh4)3h(so4)2(s)
      inp(12, i) = nh3_i(i)   * atomic_weight(6)  * 1.0d-3    !nh3(g)
      inp(13, i) = hno3_i(i)  * atomic_weight(5)  * 1.0d-3    !hno3(g)
      inp(14, i) = zeros  ! (notincl)                         !h2so4(g)
   end do

   if (chm_hetchem_metstlb_l) then
      call mach_hetv_main_2cases(npts, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                                 h_i, nh3_i, lwn_i, t_i, rh_i, rho_i,       &
                                 case_number )
      if (chm_error_l) return
      
      amsul_i = 0.0d0
      ambis_i = 0.0d0
      amnit_i = 0.0d0
      leto_i  = 0.0d0
   else
      call mach_hetv_main_12cases(npts, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                                  h_i, nh3_i, amsul_i, ambis_i, amnit_i,     &
                                  leto_i, lwn_i, t_i, rh_i, rho_i, case_number)
      if (chm_error_l) return
   endif
!
!  convert outputs
   do i = 1, npts
      outp(1, i)  = lwn_i(i)                                ! h2o(kg water/kg air)
      outp(2, i)  = h_i(i)     * 1.0d-3                     ! h+(aq)
      outp(3, i)  = nh4_i(i)   * atomic_weight(3)  * 1.0d-3 ! nh4+(aq)
      outp(4, i)  = so4_i(i)   * atomic_weight(1)  * 1.0d-3 ! so4(aq)
      outp(5, i)  = hso4_i(i)  * atomic_weight(11) * 1.0d-3 ! hso4(aq)
      outp(6, i)  = no3_i(i)   * atomic_weight(2)  * 1.0d-3 ! no3(aq)
      outp(7, i)  = amsul_i(i) * atomic_weight(7)  * 1.0d-3 ! (nh4)2so4(s)
      outp(8, i)  = amnit_i(i) * atomic_weight(8)  * 1.0d-3 ! nh4no3(s)
      outp(9, i)  = zeros !(notincl)                        ! h2so4(aq)
      outp(10, i) = ambis_i(i) * atomic_weight(9)  * 1.0d-3 ! nh4hso4(s)
      outp(11, i) = leto_i(i)  * atomic_weight(10) * 1.0d-3 ! (nh4)3h(so4)2(s)
      outp(12, i) = nh3_i(i)   * atomic_weight(6)  * 1.0d-3 ! nh3(g)
      outp(13, i) = hno3_i(i)  * atomic_weight(5)  * 1.0d-3 ! hno3(g)
      outp(14, i) = zeros !(notincl)                        ! h2so4(g)
   end do

   do l = 1, isosp
      do i = 1, npts
         if (outp(l, i) < 0.0d0) then
            write(0, *) '### Error in mach_hetv_hetvcall ###'
            write(0, *) '### Stopping after hetv with negative output detected'
            write(0, *) '# Relative humidity and temperature: ',rh_i(i),t_i(i)
            write(0, *) '# kount jlat i species_number value:',kount, jlat, i, l, outp(l, i)
            write(0, *) '# Density:              '            , rho_i(i)
            write(0, *) '###    Error in case subroutine number : ',case_number(i)
            write(0, *) '# input and output lwn: '             , inp(1,i),outp(1,i)
            write(0, *) '# input and output H+: '              , inp(2,i),outp(2,i)
            write(0, *) '# input and output NH4+: '            , inp(3,i),outp(3,i)
            write(0, *) '# input and output SO4-: '            , inp(4,i),outp(4,i)
            write(0, *) '# input and output HSO4-: '           , inp(5,i),outp(5,i)
            write(0, *) '# input and output NO3-: '            , inp(6,i),outp(6,i)
            write(0, *) '# input and output (NH4)2SO4(s): '    , inp(7,i),outp(7,i)
            write(0, *) '# input and output NH4NO3(s): '       , inp(8,i),outp(8,i)
            write(0, *) '# input and output H2SO4(aq): '       , inp(9,i),outp(9,i)
            write(0, *) '# input and output NH4HSO4(s): '      , inp(10,i),outp(10,i)
            write(0, *) '# input and output (NH4)3H(SO4)2(s): ', inp(11,i),outp(11,i)
            write(0, *) '# input and output NH3(g): '          , inp(12,i),outp(12,i)
            write(0, *) '# input and output HNO3(g): '         , inp(13,i),outp(13,i)
            write(0, *) '# input and output H2SO4(g): '        , inp(14,i),outp(14,i)
            write(0, *) '###         ABORT         ###'
            chm_error_l = .true.
            return
         end if
      end do
   end do

!  calculate some diagnostics
   solin = 0.0d0
   do l = 1, maxhet
      do i = 1, npts
         do n = 1, isize
            solin(i, l) = solin(i, l) + hetf(i, l, n)
         end do
      end do
   end do

   gasin = ghetf

   do i = 1, npts
      solout(i, 1) = outp(4, i) + atomic_weight(1) * (outp(5, i) / atomic_weight(11) +         &
                     outp(7, i) / atomic_weight(7) + outp(10, i) / atomic_weight(9) +          &
                     2.0d0 * outp(11, i) / atomic_weight(10))
      solout(i, 2) = outp(6, i) + atomic_weight(2) * (outp(8, i) / atomic_weight(8))
      solout(i, 3) = outp(3, i) + atomic_weight(3) * (2.0d0 * outp(7, i) / atomic_weight(7) +  &
                     outp(8, i) / atomic_weight(8) + outp(10, i) / atomic_weight(9) +          &
                     3.0d0 * outp(11, i) / atomic_weight(10))

      gasout(i, 1) = outp(14, i)
      gasout(i, 2) = outp(13, i)
      gasout(i, 3) = outp(12, i)
   end do

!  quick fix: put everything in the first bin for particle for now
!  set a min other than zero for mass balance purpose

   do i = 1, npts
      do l = 1, 3

!  Change:  P.A. Makar, Sept 2003:  use the same lower number limit
!  for the output values as the input values; if its less than
!  1E-24 for all, treat the counters as if "no change has occurred"

         if ((solout(i, l) <= lolimit .and. gasout(i, l) <= lolimit) .and.  &
            (solin(i, l) <= lolimit .and. gasin(i, l) <= lolimit)) then
            solout(i, l) = solin(i, l)
            gasout(i, l) = gasin(i, l)
         end if
         hetf(i, l, 1) = solout(i, l)
         ghetf(i, l) = gasout(i, l)
      end do
   end do

   do i = 1, npts
      do l = 1, 3
         massb(i, l)   = (gasin(i, l) - gasout(i, l)) / atomic_weight(l + 3) + (solin(i, l) - solout(i, l)) / atomic_weight(l)
         orgmass(i, l) = gasin(i, l) / atomic_weight(l + 3) + solin(i, l) / atomic_weight(l)
         orgmass(i, l) = orgmass(i, l) * 1.0d-3
      end do
   end do

   do i = 1, npts
      do l = 1, 3
         if (massb(i, l) > orgmass(i, l) .or. massb(i, l) < -1.0d0 * orgmass(i, l)) then
            write(0, *) '### Error in mach_hetv_hetvcall ###'
            write(0, *) '# mass balance pb ', i, l, kount, jlat
            write(0, *) '#', gasin(i, l), gasout(i, l), atomic_weight(l + 3)
            write(0, *) '#', solin(i, l), solout(i, l), atomic_weight(l), orgmass(i, l)
            write(0, *) '# verooutp ', (outp(n, i), n = 1, isosp)
            write(0, *) '###         ABORT         ###'
            write(0, *) '# Relative humidity and temperature: ',rh_i(i),t_i(i)
            write(0, *) '# Density:                           ', rho_i(i)
            write(0, *) '###    Error in case subroutine number : ',case_number(i)
            write(0, *) '# input and output lwn: '             , inp(1,i),outp(1,i)
            write(0, *) '# input and output H+: '              , inp(2,i),outp(2,i)
            write(0, *) '# input and output NH4+: '            , inp(3,i),outp(3,i)
            write(0, *) '# input and output SO4-: '            , inp(4,i),outp(4,i)
            write(0, *) '# input and output HSO4-: '           , inp(5,i),outp(5,i)
            write(0, *) '# input and output NO3-: '            , inp(6,i),outp(6,i)
            write(0, *) '# input and output (NH4)2SO4(s): '    , inp(7,i),outp(7,i)
            write(0, *) '# input and output NH4NO3(s): '       , inp(8,i),outp(8,i)
            write(0, *) '# input and output H2SO4(aq): '       , inp(9,i),outp(9,i)
            write(0, *) '# input and output NH4HSO4(s): '      , inp(10,i),outp(10,i)
            write(0, *) '# input and output (NH4)3H(SO4)2(s): ', inp(11,i),outp(11,i)
            write(0, *) '# input and output NH3(g): '          , inp(12,i),outp(12,i)
            write(0, *) '# input and output HNO3(g): '         , inp(13,i),outp(13,i)
            write(0, *) '# input and output H2SO4(g): '        , inp(14,i),outp(14,i)
            chm_error_l = .true.
            return
         end if
      end do
   end do

   return
end
