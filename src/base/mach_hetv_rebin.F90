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
! Fichier/File   : mach_hetv_rebin.ftn90
! Creation       : V. Bouchet, S. Menard, P. Makar, S. Gravel, B. pabla
! Description    : Rebinning of the heterogeneous chemistry
! Extra info     :
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
! Note           : variables gamma and sum have been renamed to hgamma and hsum

!!if_on
subroutine mach_hetv_rebin(het, dhet_chem, tempk, zpres, binnum, kount, jlat, &
                           pni, pnk)
   use mach_hetv_mod,      only: maxhet
   use mach_cam_utils_mod, only: isize
!!if_off
   use chm_consphychm_mod, only: mwt_air, kboltz, avno, pi
   use mach_cam_utils_mod, only: binrange
   use mach_hetv_mod,      only: itermax
   use chm_utils_mod,      only: chm_error_l
   use mach_hetv_mod,      only: itermax

   implicit none
!!if_on
    integer(kind=4), intent   (in) :: kount
    integer(kind=4), intent   (in) :: jlat
    integer(kind=4), intent   (in) :: pni, pnk
    real(kind=8),    intent(inout) :: het      (pni, pnk, maxhet, isize)
    real(kind=8),    intent(inout) :: dhet_chem(pni, pnk, maxhet)
    real(kind=4),    intent   (in) :: tempk    (pni, pnk)
    real(kind=4),    intent   (in) :: zpres    (pni, pnk)
    real(kind=4),    intent   (in) :: binnum   (pni, pnk, isize)
!!if_off
!
! Local variables
!
   integer(kind=4), parameter      :: nsp = 3
   integer(kind=4)                 :: i, k, l, n, iter
   real(kind=8), parameter         :: sslppa = 101325.0d0   ! standard sea-level pressure (Pa)
   real(kind=8), dimension(nsp)    :: mwt = (/98.0795d0, 63.0128d0, 17.03056d0/) ! H2SO4     HNO3     NH3
   real(kind=8), dimension(nsp)    :: hgamma = (/1.0d0, 0.1d0, 0.4d0/)           ! H2SO4     HNO3     NH3
   real(kind=8), parameter :: one_third = 1.0d0 / 3.0d0, pi_dp = acos(-1.0d0) !dble(pi)
   real(kind=8) :: dg, dg_b, dg_pow, pavn(nsp)
   real(kind=8) :: radius(isize), ki(pni, pnk, nsp, isize)
   real(kind=8) :: kn
   real(kind=8) :: diff(pni, pnk, nsp), hsum(pni, pnk, nsp)

   real(kind=8) :: totav, totap, delta, ffrac, orgmass
!  new variables added, makar dec 2008 revision:
   real(kind=8) :: mass_to_be_removed, mass_available
   real(kind=8) :: totfrac
   real(kind=8) :: mass_lost_total, mass_lost_bin
   real(kind=8) :: mass_to_be_removed_initial

!  Convert to the right units radius (m to cm), press (mb to atm)
 
!  Calculate the diffusivity and ki
   do l = 1, nsp
      dg     = 0.36855513d0 * mwt(l) + 6.2869566d0
      dg_b   = sqrt((dble(mwt_air) + mwt(l)) / (dble(mwt_air) * mwt(l)))
      dg_pow = (20.1d0 ** one_third + dg ** one_third) ** 2
      pavn(l) = 8.0d0 * dble(kboltz) / (pi_dp * mwt(l) / dble(avno))
      do k = 1, pnk
         do i = 1, pni
            diff(i, k, l) = dg_b * 1.0d-3 * (dble(tempk(i, k)) ** 1.75d0) / &
                            dg_pow / (dble(zpres(i, k)) / sslppa)
         end do
      end do
   end do

   do n = 1, isize
!     Average radius of each bin (in cm)
      radius(n) = 0.5d0 * 100.0d0 * dble(binrange(1, n) + binrange(2, n))

      do l = 1, nsp
         do k = 1, pnk
            do i = 1, pni
               kn = 3.0d0 * diff(i, k, l) / (radius(n) *    &
                    sqrt(pavn(l) * dble(tempk(i, k))))
               ki(i, k, l, n) = 4.0d0 * pi_dp * radius(n) * diff(i, k, l) /  &
                                (1.0d0 + (4.0d0 * kn / (3.0d0 * hgamma(l))) *  &
                                (1.0d0 - 0.47d0 * hgamma(l) / (1.0d0 + kn)))
            end do
         end do
      end do
   end do
!  calculate sum of kiNi
   do l = 1, nsp
      do k = 1, pnk
         do i = 1, pni
            hsum(i, k, l) = 0.0d0     ! Initialize hsum
            do n = 1, isize
               hsum(i, k, l) = hsum(i, k, l) + (ki(i, k, l, n) * dble(binnum(i, k, n)))
            end do
         end do
      end do
   end do

!  Chose how to redistribute
!  Redistribute
   do l = 1, nsp
      do k = 1, pnk
         do i = 1, pni

            if (dhet_chem(i, k, l) >= 0.0d0) then
               totap = 0.0d0
               totav = 0.0d0
               do n = 1, isize
                  totav = totav + het(i, k, l, n)
!  Precision problems with summation at single precision
                  ffrac = dhet_chem(i, k, l) * ki(i, k, l, n) * dble(binnum(i, k, n)) / hsum(i, k, l)
                  totap = totap + het(i, k, l, n) + ffrac
                  het(i, k, l, n) = het(i, k, l, n) + ffrac
               end do

               if (totap > 1.001d0 * (totav + dhet_chem(i, k, l)) .or. totap < 0.999d0 * (totav + dhet_chem(i, k, l))) then
                  orgmass = (totav + dhet_chem(i, k, l)) * 1.0d-3
                  delta = totap - (totav + dhet_chem(i, k, l))
                  write(0, *) '### Error in mach_hetv_rebin ###'
                  write(0, *) '# mass balance pb ', i, l, k, kount, jlat
                  write(0, *) '#', totav, totap, dhet_chem(i, k, l)
                  write(0, *) '#', delta, orgmass
                  write(0, *) '# rebin>0 '
                  write(0, *) '###         ABORT         ###'
                  chm_error_l = .true.
                  return
               end if
            else

!  This branch is taken if the effect of the bulk-calculated inorganic 
!  heterogeneous chemistry was to remove particle mass.
!
!  Determine the bulk mass that is to be removed from the particle phase
!  for the current horizontal gridpoint (i), vertical gridpoint (k), and
!  chemical species (l), across all bins (n):
               mass_available = 0.0d0
               do n = 1, isize
                  mass_available = mass_available + het(i, k, l, n)
               end do

!  The mass available for removal is that in the particle phase
!  in the previous time step prior to the application of
!  the inorganic chemistry partitioning code (i.e. "het" is
!  the value before the inorganic chemistry partitioning code
!  was called, and the sum of het over all bins is the total
!  particle mass that may be removed for the given species).
!  The mass checks in the inorganic chemistry partitioning code
!  are assumed to conserve mass, and the bulk amount of mass lost
!  from the particle phase, calculated from that code, is
!  dhet_chem.  mass_available, as calculated above, therefore
!  should be .ge. abs(dhet_chem).  Round off errors
!  may make the mass being removed greater than the mass available:
!  the following min statement prevents this from happening.
               mass_to_be_removed = min(mass_available, -dhet_chem(i, k, l))
               mass_to_be_removed_initial = mass_to_be_removed

!  Iterative loop:  mass is removed from the larger bins first
               iter = 0
               do while (mass_to_be_removed > 0.0d0 .and. iter < itermax)
                  iter = iter + 1
!  For the current iteration, for the particle bins which have non-zero mass
!  (i.e. are available for mass removal), calculate the
!  ratios of mass that should go to each bin based on gas diffusion:
                  totfrac = 0.0d0
                  do n = isize, 1, -1
                     if (het(i, k, l, n) > 0.0d0) then
                        totfrac = totfrac + (ki(i, k, l, n) * dble(binnum(i, k, n)))
                     end if
                  end do

!  Zero the amount of mass to be lost during the current iteration for this species:
                  mass_lost_total = 0.0d0

!  Determine the fraction of mass to be lost from each bin considering diffusive transfer of gas alone:
                  do n = isize, 1, -1
                     mass_lost_bin = 0.0d0
                     if (het(i, k, l, n) > 0.0d0) then
                        mass_lost_bin = min((ki(i, k, l, n) * dble(binnum(i, k, n)) &
                                         / totfrac * mass_to_be_removed), &
                                         het(i, k, l, n))
                     end if

!  Reduce the mass in the bin, applying a zero minimum to avoid round-off error
                     het(i, k, l, n) = max((het(i, k, l, n) - mass_lost_bin), 0.0d0)

!  Increment the total amount of mass lost from the particle phase for the given species:
                     mass_lost_total = mass_lost_total + mass_lost_bin

!  Repeat the process for the next bin.
                  end do

!  If mass_lost_total is < mass_to_be_removed, some of the bins in the above iteration were totally depleted of
!  mass, and the remaining mass must be taken from another pass through the iteration.
!  Increment the amount of mass to be removed, employing a max statement to avoid zero round-off:
                  mass_to_be_removed = max(mass_to_be_removed - mass_lost_total, 0.0d0)
               end do

!  Determine the fraction of mass to be removed that still remains:
               totfrac = mass_to_be_removed / mass_to_be_removed_initial

!  Either more than 5 iterations were required to re-sort the mass, or there's no mass left to sort
!  (a round off error of 1e-6 is allowed). If the iteration limit was reached, stop the code with an error message:
               if (totfrac .gt. 1.0d-6 .and. iter > itermax) then
                  write(0, *) '### Error in mach_hetv_rebin ###'
                  write(0, *) '# Maximum number of iterations exceeded #'
                  write(0, *) '#Mass available in array het: ', mass_available
                  write(0, *) '#Mass to be removed in dhet: ', dhet_chem(i, k, l)
                  write(0, *) '#Mass still to be removed after ', itermax, ' iterations: ', mass_to_be_removed
                  write(0, *) '#Initial mass to be removed and ratio: ', mass_to_be_removed_initial, totfrac
                  write(0, *) '###         ABORT         ###'
                  chm_error_l = .true.
                  return
               end if

!  The mass has been successfully redistributed over the bins.  Continue with
!  the next horizontal level(i), vertical level(k), and bulk species (l)
            end if
         end do
      end do
   end do
return
end
