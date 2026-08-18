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
! Fichier/File   : mach_cam_coagd.ftn90
! Creation       : S. Gong, E. Girard, S. Gravel and for GEM-MACH, June 2008
!
! Description    : Calculate coagulation processes among size bins
!                  Brownian, turbulent and gravitational coagulation are
!                  computed in this module.
!
! Extra info     : - Previous version based on gelbard et al [1980]
!                    scheme of coagulation from gelbard et al. with
!                    pre-calculated coefficients of coagulation was used.
!                    (E.Girard, Mar 01, 1995)
!                  - Recode the coagulation module except E. Girard adapting
!                    coagulation coefficients funtion beta from previous
!                    code. featuring:
!                    A. A unique do loop to handel all the coagulation
!                       processes among size bins and elliminate the
!                       duplication of computing both k(1, 2) and k(2, 1)
!                       which are equal.
!                    B. Vectorization of the code for best efficiency.
!                    C. The coefficients are calculated based on the real
!                       size and density of particles.
!                    D. Simplify the integration procedure to be run for
!                       a 3-d climate model.
!                    E. Tendency calculation is based on the original
!                       dynamical coagulation equations where particle
!                       number concentrations are prognostic variables.
!                    F. Physical properties such as mean free path,
!                       viscosity and setlling velocity are computed
!                       in a consistant way as in other aerosol routines.
!                    G. Coagc can now run every time step for much
!                       smaller overhead.
!                       (S. Gong Sept 23, 1997)
!                  - Modified for multi-component aerosols. (S. Gong, Oct 01, 1997)
!                  - Adapt the jacobson's scheme to conserve volume [mass]
!                    of each species after coagulation.
!                    [Jacobson et.al. atm. env. 1994] (S. Gong, Nov 06, 1998)
!                  - Introduction of igf and igfij to spped up the code.
!                    the percentage of coagulation in gcmiii wend down from
!                    17% to 1% of total cpu time. (S. Gong, May 12, 2000)
!
!
! Arguments:  IN
!               throw   -> Temp
!              roarow   -> Air density (kg/m3)
!              rhsize   -> Wet radius
!               icob    -> size bin number with which coagulation will apply
!              pdepv    -> gravitational settling velocity
!               pdiff   -> diffusion coefficient
!               mae     -> 0
!              rhop     -> Final wet density
!               amu     -> Air's dynamic viscosity
!
!             OUT
!              RTCOA    -> Coagulation rate
!
!             IN/OUT
!              XROW     -> Tracers concentration update from coagulation
!
! Local variables:
!              BETA     -> Coagulation coefficient
!              AERONUM  -> Number conconcenration (#/kg)
!              TOTMAS   -> Total mass mixing ratio for each bin (all components)
!============================================================================
!
!!if_on
subroutine mach_cam_coagd(throw, roarow, rtcoa, rhsize, xrow, pdepv, pdiff, &
                          mae, rhop, amu, pni, pnk)
   use mach_cam_utils_mod, only: isize, ntr
!!if_off
   use chm_utils_mod,      only: chm_timestep
   use mach_cam_utils_mod, only: icob, iae1, icom, igf, igfij, coagfr, pvol, tmin
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: mae
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: throw (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent  (out) :: rtcoa (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
!!if_off
!
!  local variables
!
   integer(kind=4)         :: n, l, il, i, j, k, nn, ij
   integer(kind=4)         :: ik, ip, nt, no
   real(kind=4)            :: dx, v1, dl3, gx0, diffx, diffy, dsum, cbar, gmean
   real(kind=4)            :: rtloss, turb1, beta0, beta1, vkp
   real(kind=4)            :: l1
   real(kind=4)            :: totmas (pni, pnk, isize)
   real(kind=4)            :: aeronum(pni, pnk, isize)
   real(kind=4)            :: oldnum (pni, pnk, isize)
   real(kind=4)            :: beta   (pni, pnk, isize, isize)
   real(kind=4)            :: binloss(pni, pnk, isize)
   real(kind=4)            :: aerop1 (pni, pnk, isize)
   real(kind=4)            :: cbar12 (pni, pnk, isize)
   real(kind=4)            :: gx     (pni, pnk, isize)
   real(kind=4)            :: sgsum  (pni, pnk)
   real(kind=4)            :: dbsum  (pni, pnk)
   real(kind=4), parameter :: turbds = 0.002, stick = 1.0, xiao = 1.0e6

   rtcoa    = 0.0

!  non start-run begins

!  update the current aerosol number

   totmas   = 0.0
   do n = 1, isize
      do l = 1 + mae, pnk
         do i = 1, pni
            do nt = 1, icom
               no = n + isize * (nt - 1) + (iae1 - 1)
               totmas(i, l, n) = totmas(i, l, n) + xrow(i, l, no)
            end do
            aeronum(i, l, n) = totmas(i, l, n) / (pvol(n) * rhop(i, l, n))
            oldnum(i, l, n) = aeronum(i, l, n) * roarow(i, l)
         end do
      end do
   end do

!  In case of 0 totmas after cloud chem, aeronum is assigned with a value so that
!  the number density will be of 0.1*XIAO (# m-3). This is introduced to avoid
!  overflow in the following calculation, and should not affect the final coagulation
!  calculation since bins with number density smaller than XIAO (1.E06 m-3) are
!  avoided in the final calculation. (Wanmin, Mar 26, 2001)
!
   where (totmas <= 0) oldnum = 0.1 * xiao

   do n = 1, icob
      do l = 1 + mae, pnk
         do il = 1, pni
            dx = 2.0 * rhsize(il, l, n)
            v1 = pvol(n) * rhop(il, l, n)
            cbar12(il, l, n) = 3.51568e-23 * throw(il, l) / v1
            l1 = 2.5465 * pdiff(il, l, n) / sqrt(cbar12(il, l, n))
            dl3 = (dx + l1) * (dx + l1) * (dx + l1)
            gx0 = (dl3 - (dx * dx + l1 * l1) ** 1.5) / (3.0 * dx * l1) - dx
            gx(il, l, n) = gx0 * gx0
         end do
      end do
   end do

   do i = 1, icob
      do j = i, icob
         do l = 1 + mae, pnk
            do il = 1, pni
!  diffusion coefficients
               diffx = pdiff(il, l, i)
               diffy = pdiff(il, l, j)

               dsum = 2.0 * (rhsize(il, l, i) + rhsize(il, l, j))

!  brownian coagulation coefficient        [v1, v2, vr - particle mass, kg]
               cbar  = sqrt(cbar12(il, l, i) + cbar12(il, l, j))
               gmean = sqrt(gx(il, l, i) + gx(il, l, j))
               beta0 = 6.2832 * (diffx + diffy) * dsum / (dsum /         &
                       (dsum + 2.0 * gmean) + 8.0 * (diffx + diffy) /    &
                       (cbar * dsum * stick))

!  add gravitational coagulation
               beta1 = beta0 + 0.7854 * dsum ** 2                        &
                       * abs(pdepv(il, l, i) - pdepv(il, l, j))

!  add turbulent coagulation
               turb1 = .1618 * sqrt(turbds * roarow(il, l) / amu(il, l)) &
                       * dsum * dsum * dsum
               beta(il, l, i, j) = beta1 + stick * turb1
               beta(il, l, j, i) = beta(il, l, i, j)
            end do
         end do
      end do
   end do

   aerop1 = 0.0
   do k = 1, icob
      dbsum   = 0.0
      sgsum   = 0.0
      binloss = 0.0
      do j = 1, icob
         do l = 1 + mae, pnk
            do il = 1, pni

!  number loss of k due to collision with j [1-isize]
               if (oldnum(il, l, j) > xiao) then
                  binloss(il, l, j) = (1.0 - coagfr(k, j, k)) * beta(il, l, k, j) * oldnum(il, l, j)
                  sgsum(il, l) = sgsum(il, l) + binloss(il, l, j)
               end if
            end do
         end do
      end do

      do ij = 1, igf(k)     !gathered points for coagfr
         i = igfij(k, ij, 1)
         j = igfij(k, ij, 2)
         do l = 1 + mae, pnk
            do il = 1, pni
!  volume gain of k due to collision of i and j [=<k]-- [m3 s-1]
               dbsum(il, l) = dbsum(il, l) + coagfr(i, j, k) * beta(il, l, i, j) &
                              * pvol(i) * aerop1(il, l, i) * oldnum(il, l, j)
            end do
         end do
      end do

!  total number of k after coagulation [# m-3]
      do l = 1 + mae, pnk
         do il = 1, pni
            aerop1(il, l, k) = (oldnum(il, l, k) + chm_timestep / pvol(k) &
                               * dbsum(il, l)) / (1.0 + chm_timestep * sgsum(il, l))
         end do
      end do

!  Mass balance for each species of k. Loss of k is the sum of mass gained by all j
!  sgsum: dimension less, aerop1: # m-3
      do l = 1 + mae, pnk
         do il = 1, pni

!  lost tendency of k
!  aerop1(il, l, k) * sgsum(il, l) * pvol(k) is the volume loss rate of bin k due to the
!  coagulation by /(oldnum*pvol(k)) *xrow(il, l, ik) the loss tendency of ik is obtained.
            if (oldnum(il, l, k) > xiao) then
               do nn = 1, icom
                  ik = (nn - 1) * isize + k + (iae1 - 1)
                  vkp = aerop1(il, l, k) / oldnum(il, l, k) * xrow(il, l, ik)
                  rtcoa(il, l, ik) = rtcoa(il, l, ik) - vkp * sgsum(il, l)
!
! gain tendency of i due to loss of [k, j]
                  do j = 1, icob
                     do i = max(j, k), icob
                        if (coagfr(k, j, i) > 0.0) then
                           ip = (nn - 1) * isize + i + (iae1 - 1)
                           rtloss = vkp * binloss(il, l, j) * coagfr(k, j, i)
                           rtcoa(il, l, ip) = rtcoa(il, l, ip) + rtloss
                        end if
                     end do
                  end do
               end do      !end of nn loop
            end if      !end of oldnum conditional
         end do
      end do
!
   end do         !end of k  loop

!  update tracer due to coagulations
   do k = 1, isize - 1
      do nn = 1, icom
         ik = (nn - 1) * isize + k + (iae1 - 1)
         do l = 1 + mae, pnk
            do il = 1, pni
               xrow(il, l, ik) = max(tmin, xrow(il, l, ik) + rtcoa(il, l, ik) * chm_timestep)
            end do
         end do
      end do
   end do

   return
end subroutine mach_cam_coagd
