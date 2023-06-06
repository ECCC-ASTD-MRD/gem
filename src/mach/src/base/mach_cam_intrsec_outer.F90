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
!  Modification:  the mach_cam_intrsec.ftn90 algorithm has been modified to
!  allow a subdivision of the input bin structure to improve accuracy for the
!  mass transfer.  6 options are possible, making use of the user input
!  variable "chm_intrsec_ver" on the call list.  The code
!  splits the input particle mass, according to the version number,
!  calculates the corresponding aerosol number in each bin,
!  and then splits the condensation rate by version number.
!
!  The options (chm_intrsec_ver = ) are:
!   1 = mass evenly distributed over sub-bins by radius '
!       condensation rate redistributed using area ratio.'
!   2 = mass evenly distributed over sub-bins by ln radius'
!       condensation rate redistributed using area ratio.'
!   3 = mass distributed according to locally scaled trimodal lognormal distribution'
!       condensation rate redistributed using area ratio.'
!   4 = mass evenly distributed over sub-bins by radius '
!       fuchs-sutugin equation used to redistribute condensation rate'
!   5 = mass evenly distributed over sub-bins by ln radius '
!       fuchs-sutugin equation used to redistribute condensation rate'
!   6 = mass distributed according to locally scaled trimodal lognormal distribution'
!       fuchs-sutugin equation used to redistribute condensation rate.'
!
!  P.A. Makar, W. Gong, S. Gong, February 2009
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_cam_intrsec.ftn90
! Creation       : S. Gong, S. Gravel and B. Pabla for GEM-MACH, June 2008
!
! Description    : This module computes intersectional transport of aerosols
!                  due to condensation or cloud processes
!
! Extra info     : First version created by S. Gong Aug 11 1997 for CAM
!
! Arguments:  IN
!               rtcond   -> Mass transfer rate onto each particle size bin
!               totmas   -> Total mass of aerosol in each bin
!               nn       -> Loop index on aerosol types
!              rtcond    -> Mass transfer rate onto each particle size bin
!              rtcond_n  -> Mass transfer rate onto each particle sub bin
!               aeronum  -> Number concentration of aerosols
!               aeronum_n-> Number concentration of aerosols in sub bins
!                  mae   -> 0
!               xrow_n   -> Tracer concentration with sub-bins for the particles
!
!             IN/OUT
!              XROW      -> Tracer concentration in each bin before/after intersection tranport
!
!             MODULE: CAM_UTILS
!               ntr      -> Total number of trace substances (gases and aerosols)
!               iae1     -> Index of first aerosol in trace substance list (ntr)
!               rhop0    -> dry aerosol density [ug/m^3]
!               isize    -> Number of size bins
!============================================================================
!
!!if_on
 subroutine mach_cam_intrsec_outer(xrow, nn, rtcond, aeronum, mae, pres, &
                                   tp, roarow, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
!!if_off
   use chm_consphychm_mod,   only: mwt_air, pi, rgasi
   use chm_nml_mod,          only: chm_intrsec_ver, chm_intrsec_ndiv
   use mach_cam_headers_mod, only: mach_cam_intrsec_inner
   use mach_cam_utils_mod,   only: nbnd, iae1, icom, nsb, rhop0, &
                                   dvn, sub_pvol, sub_rn, sub_a_n
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: mae
   integer(kind=4), intent   (in) :: nn
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rtcond (pni, pnk, isize)
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent   (in) :: pres   (pni, pnk)
   real(kind=4),    intent   (in) :: tp     (pni, pnk)
   real(kind=4),    intent   (in) :: roarow (pni, pnk)
!!if_off
!
!  local variables
!
   integer(kind=4)                           :: k, np, pos, l, il, n, i
!  New variables:
   real(kind=4) , dimension(pni, pnk, nbnd)  :: aeronum_n, rtcond_n
   real(kind=4) , dimension(pni, pnk, nsb)   :: xrow_n
   real(kind=4)    :: kn, an
   real(kind=4)    :: dfso4, frp, fn
   real(kind=8)    :: suma
   real(kind=4)    :: adv_so4, adv_air, net_vol, mwfactor
   integer(kind=4) :: j, jsub, jk
!
   real(kind=4), parameter :: ace = 0.02, xiao = 1.0e3, cub = 1.0 / 3.0
   real(kind=4), parameter :: so4mw = 98.1e-03
   real(kind=8), parameter :: smf1 = 1.D-25   ! small number limit for total particle area / kg air
!
! -------------------------------------------------------------
!  skip subdivision when the 12-bin structure is uesd
!
   if (isize == 12) then
      do k = 1, isize
         do il = 1, pnk
            do i = 1, pni
               rtcond_n(i, il, k) = rtcond(i, il, k)
               aeronum_n(i, il, k) = aeronum(i, il, k)
!
               do np = 1, icom
                  pos = (np - 1) * isize + k + iae1 - 1
                  jk  = (np - 1) * isize + k
                  xrow_n(i, il, jk) = xrow(i, il, pos)
               end do
            end do
         end do
      end do
!
   else
!
!  Commence subdivision of original bins:
!
!  Assign new particle mass, and new particle number based on individual 
!  particle volume, density,:
      do k = 1, isize
         do j = 1, chm_intrsec_ndiv
            jsub = (k - 1) * chm_intrsec_ndiv + j
            do il = 1, pnk
               do i = 1, pni
                  net_vol = 0.0
                  do np = 1, icom
                     pos = (np - 1) * isize + k + iae1 - 1
                     jk = (np - 1) * nbnd + jsub
                     xrow_n(i, il, jk) = xrow(i, il, pos) * dvn(jsub)
                     net_vol = net_vol + xrow_n(i, il, jk) / rhop0(np)
                  end do
                  aeronum_n(i, il, jsub)  = net_vol / sub_pvol(jsub)
               end do
            end do
         end do
      end do
!
!  Determine revised condensation rate:
      if (chm_intrsec_ver == 1 .or. chm_intrsec_ver == 2 .or. chm_intrsec_ver == 3 ) then
         do k = 1, isize
            do il = 1, pnk
               do i = 1, pni
                  suma = 0.D0
                  do j = 1, chm_intrsec_ndiv
                     jk = (k - 1) * chm_intrsec_ndiv + j
                     suma = suma + dble(sub_a_n(jk)) * dble(aeronum_n(i, il, jk))
                  end do
                  if (suma > 0.0D0) then
                     do j = 1, chm_intrsec_ndiv
                        jk = (k - 1) * chm_intrsec_ndiv + j
                        rtcond_n(i, il, jk) = rtcond(i, il, k) * real(      &
                             dble(sub_a_n(jk)) * dble(aeronum_n(i, il, jk)) &
                             / max(suma, smf1) )
                     end do
                  else
!  Summed aerosol number * area = zero.  Condensation rate therefore set to zero:
                     do j = 1, chm_intrsec_ndiv
                        jk = (k - 1) * chm_intrsec_ndiv + j
                        rtcond_n(i, il, jk) = 0.0
                     end do
                  end if
               end do
            end do
         end do
      end if

!  Alternative condensation rate (based on fuchs-sutugin equation)
!  lifted from Sunling Gong's sulfate condensation code.
!  Note that dry radius is used here instead of wet.
      if (chm_intrsec_ver == 4 .or. chm_intrsec_ver == 5 .or. chm_intrsec_ver == 6) then
!
         rtcond_n = 0.0
!
         !atomic diffusion volume (adv_*) [cm3] Makar et al, 1998
         adv_air = 0.369 * mwt_air + 6.29
         adv_so4 = 0.369 * (so4mw * 1.0e3) + 6.29 
         mwfactor = (adv_air ** cub + adv_so4 ** cub) ** 2
!
         do n = 1, nbnd
            do l = 1 + mae, pnk
               do i = 1, pni
                  if (aeronum_n(i, l, n) * roarow(i, l) > xiao) then
!  diffusion coefficient of h2so4 [gas] [Perry and Green, 1984]
!  [m2 s-1]
!  where 0.21145 = [(ma+mg)/mamg]^1/2, ma=28.97, mg=98.1 [h2so4]
                     dfso4 = 1.0e-7 * tp(i, l) ** 1.75 * 0.21145 / &
                             (pres(i, l) / 1.01325e5 * mwfactor)
                             

!  the mean free path of vapour sulphuric acid
                     frp = 3.0 * sqrt(pi * 98.1 / (8.0e3 * rgasi * tp(i, l))) * dfso4

!  ace --> accommodation coefficient ~0.02
                     kn = frp / sub_rn(n)
                     fn = (1.0 + kn) / (1.0 + 1.71 * kn + 1.33 * kn * kn)
                     an = 1.0 / (1.0 + 1.33 * kn * fn * (1.0 / ace - 1.0))

!  the condensation coefficients for combined computation next.
!  the aeronum*roarow accounts for all the particles in bin n.
                     rtcond_n(i, l, n) = 4.0 * pi * sub_rn(n) * dfso4 * fn *  &
                                         an * aeronum_n(i, l, n) * roarow(i, l)
                  end if
               end do
            end do
         end do
         do k = 1, isize
            do il = 1, pnk
               do i = 1, pni
                  suma = 0.D0
                  do j = 1, chm_intrsec_ndiv
                     jk = (k - 1) * chm_intrsec_ndiv + j
                     suma = suma + dble(rtcond_n(i, il, jk))
                  end do
!  If the sum of condensation rates is > zero, use it to ratio
!  the original condensation rate.  If the sum is zero, there is no condensation,
!  and the rtcond_n values are already zero; no need to reset them to zero.
                  if (suma > 0.0D0) then
                     do j = 1, chm_intrsec_ndiv
                        jk = (k - 1) * chm_intrsec_ndiv + j
                        rtcond_n(i, il, jk) = rtcond(i, il, k) * &
                                             real(dble(rtcond_n(i, il, jk)) / max(suma, smf1))
                     end do
                  end if
               end do
            end do
         end do
      end if
!
   end if

!  Finished subdivision of inputs for condensation.
!  Call the original intrsec code to determine the mass transfer:
!
   call mach_cam_intrsec_inner(xrow_n, nn, rtcond_n, aeronum_n, mae, pni, pnk)

!  Reassign the mass back into the original xrow locations:
   if (isize == 12) then
      do k = 1, isize
         do np = 1, icom
            pos = (np - 1) * isize + k + iae1 - 1
            jk  = (np - 1) * isize + k
            do il = 1, pnk
               do i = 1, pni
                  xrow(i, il, pos) = xrow_n(i, il, jk)
               end do
            end do
         end do
      end do
   else
      do k = 1, isize
         do np = 1, icom
            pos = (np - 1) * isize + k + iae1 - 1
            do il = 1, pnk
               do i = 1, pni
                  xrow(i, il, pos) = 0.0
                  do j = 1, chm_intrsec_ndiv
                     jsub = (k - 1) * chm_intrsec_ndiv + j
                     jk = (np - 1) * nbnd + jsub
                     xrow(i, il, pos) = xrow(i, il, pos) + xrow_n(i, il, jk)
                  end do
               end do
            end do
         end do
      end do
   endif

   return
end subroutine mach_cam_intrsec_outer
