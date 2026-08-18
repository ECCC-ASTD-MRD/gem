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
! Fichier/File   : mach_cam_sulfate.ftn90
! Creation       : S. Gong, B. Pabla and S. Gravel for GEM-MACH, June 2008
! Description    : Sulfate aerosol formation routine
!                  (1) Calculate oh concentration
!                  (2) Oxidation of SO2 by OH
!                  (3) Nucleation, i.e. generation of new sulphate aerosols
! Extra info     : - First version created by S. Gong Aug 14 1996 for CAM
!                  - Vectorization and add chemistry chemical fields of OH, H2O2, O3 and
!                    NO3 from ncar images model are used [daily averaged] (S. Gong, May 26, 1997)
!                  - Combined nucleation and condensation rates into one differential
!                    equation and solved analytically. ((S. Gong, jun 4, 1998)
!                  - New nucleation scheme of Kulmala, [1998] K, Von Salzen, was
!                    introduced into the combined analytical equation. (S. Gong, Dec 15, 1998)
!                  - Combined h2so4 production, condensation and nucleation into
!                    an approximate solution with a new statement function. (S. Gong, May 09, 2000)
!                  - Added double precision for some variables and implemented veronique changes... (B. Pabla)
!                  - Optimized code by reducing a number of 2D arrays, removing unnecessary 
!                    variables and reordering do loops (A. Akingunola, 2018)
!
! Arguments:  IN
!                aeronum -> Number-Concentration of Aerosol
!                ntr     -> total number of trace substances (gases and aerosols)
!                roarow  -> Air density (kg/m3)
!                pres    -> Mid layer pressure [Pa]
!                tp      -> temperature
!                rh      -> Relative humidity
!                rhsize  -> Wet radius
!                rtso2   -> so2 oxidation rate
!                mae     -> 0
!
!             OUT
!
!                RTNUCL  -> Nucleation
!                RTCOND  -> Mass transfer rate onto each particle size bi
!                PCOND   -> fractional condensation to each bin
!
!             IN/OUT
!
!                XROW    -> tracers concentration
!
!============================================================================
!
!!if_on
subroutine mach_cam_sulfate(aeronum, xrow, roarow, pres, tp, rh, rhsize, &
                            rtnucl, rtcond, rtso2, pcond, mae, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
!!if_off
   use chm_utils_mod,        only: chm_timestep
   use chm_consphychm_mod,   only: avno, mwt_air, pi, rgasi
   use mach_cam_headers_mod, only: mach_cam_intrsec_outer
   use mach_cam_utils_mod,   only: binrange, igs_SO2, igs_SO4, iae1, iae_SU, &
                                   rhop0, ursv, condnu
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: mae
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: roarow (pni, pnk)
   real(kind=4),    intent   (in) :: pres   (pni, pnk)
   real(kind=4),    intent   (in) :: tp     (pni, pnk)
   real(kind=4),    intent   (in) :: rh     (pni, pnk)
   real(kind=4),    intent   (in) :: rhsize (pni, pnk, isize)
   real(kind=4),    intent  (out) :: rtnucl (pni, pnk)
   real(kind=4),    intent  (out) :: rtcond (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rtso2  (pni, pnk)
   real(kind=4),    intent  (out) :: pcond  (pni, pnk, isize)
!!if_off
!
! Local variables
!
   integer(kind=4) :: i, l, n, inew, ii
   real(kind=4)    :: dfso4, frp, fn, kn, blg, mwfactor
   real(kind=8)    :: totcdnd, afak, so4td, aexpd, rtso4d, efact, exppp, ucf, condnu_dp(15)
   real(kind=4)    :: rh01
   real(kind=4)    :: adv_air, adv_so4, anw, adelta, aconcr
   real(kind=4)    :: apara, aparb, aparc, apard, apare, aexp, an, atempe
   real(kind=4)    :: so4t0, so4t(15)
   real(kind=4)    :: oldso4, newso4, sumso4, avail_so4, rtso4
   real(kind=4)    :: ratio, rtcondsum, totcd, totcdn, nuclsum, condsum
   real(kind=4), parameter :: ace   = 0.02,    xiao = 1.0e3,    cub = 1.0 / 3.0
   real(kind=4), parameter :: so4mw = 98.1e-03, akb = 1.380662e-23, ara = 1.0

   pcond  = 0.0
   rtcond = 0.0
   rtnucl = 0.0

   condnu_dp = dble(condnu)
   
! unit-change factor:
   ucf = ursv * binrange(1,1) * binrange(1,1) * binrange(1,1) * rhop0(iae_SU) * 1.0d+06

! h2so4 diffusion coefficient factor
   !atomic diffusion volume (adv_*) [cm3] Makar et al, 1998
   adv_air = 0.369 * mwt_air + 6.29
   adv_so4 = 0.369 * (so4mw * 1.0e3) + 6.29 
   mwfactor = (adv_air ** cub + adv_so4 ** cub) ** 2

   do l = 1 + mae, pnk
      do i = 1, pni
!
!  condensation of h2so4 [g] to exsiting particles.
!
!  diffusion coefficient of h2so4 [gas] [Perry and Green, 1984]  [m2 s-1]
!  where 0.21145 = [(ma+mg)/mamg]^1/2, ma=28.97, mg=98.1 [h2so4]
         dfso4 = 1.0e-7 * tp(i, l) ** 1.75 * 0.21145 / (pres(i, l) /  &
                 1.01325e5 * mwfactor)

!  the mean free path of vapour sulphuric acid
         frp = 3.0 * sqrt(pi * 98.1 / (8.0e3 * rgasi * tp(i, l))) * dfso4
         totcdn = 0.0
         do n = 1, isize
            if (aeronum(i, l, n) * roarow(i, l) > xiao) then
!  ace --> accommodation coefficient ~0.02
               kn = frp / rhsize(i, l, n)
               fn = (1.0 + kn) / (1.0 + 1.71 * kn + 1.33 * kn * kn)
               an = 1.0 / (1.0 + 1.33 * kn * fn * (1.0 / ace - 1.0))

!  the condensation coefficients for combined computation next.
!  the aeronum*roarow accounts for all the particles in bin n.
               rtcond(i, l, n) = 4.0 * pi * rhsize(i, l, n) * dfso4 * fn * an *  &
                                 aeronum(i, l, n) * roarow(i, l)
               totcdn = totcdn + rtcond(i, l, n)
            end if
         end do
!
!  fractional condensation to each bin for use in soa
         do n = 1, isize
            if (aeronum(i, l, n) * roarow(i, l) >  xiao) then
               pcond(i, l, n) = min(rtcond(i, l, n) / totcdn, 1.0)
            end if
         end do
!
         rtso4 = max(0.0, -98.1 / 64.6 * rtso2(i, l))
         so4t0 = max(0.0, xrow(i, l, igs_SO4))
         if (so4t0 > 1.0e-25 .and. totcdn > 0.0 .and. xrow(i, l, igs_SO2) > 1.0e-15) then
!  nucleation
!  scheme by kulmala et al [1998]
!  this portion of code was contributed by knut von salzen
!  for all models - gcm, rcm or aurams
            rh01   = max(min(rh(i, l), 1.0), 0.1)
            atempe = max(min(tp(i, l), 298.0), 233.0)

!  water vapour concentration in 1/cm^3
            anw = rh01 * 1.0e-06 * exp(77.34491296 - 7235.42451 / atempe &
                  - 8.2 * log(atempe) + 5.7113e-03 * atempe) / (akb * atempe)
            adelta = atempe / 273.15
            apara = 25.1289 - 4890.8 / atempe - 2.2479 * adelta * rh01
            aparb = 7643.4 / atempe - 1.9712 * adelta / rh01
            aparc = -1743.3 / atempe
            apard = 1.2233 - 0.0154 * ara / (ara + rh01) - 0.0415 * log(anw) &
                    + 0.0016 * atempe
            apare = 0.0102

!  calculation of threshold H2SO4 concentration
            aconcr = -14.5125 + 0.1335 * atempe - 10.5462 * rh01 &
                     + 1958.4 * rh01 / atempe
            aexp = apara + aparb * apare

!  H2SO4 in kg/kg in combined equation and J in per cm^3/s

            exppp = dble(aexp * log(avno * roarow(i, l) / so4mw * 1.0d-06) - apara * aconcr)
            efact = exp(exppp)
            afak = exp(aparb * apard + aparc) * efact
!  to kg/kg/s
            afak = afak * ucf / roarow(i, l)

            oldso4 = so4t0
            totcdn = min(totcdn, 0.025)
            totcdnd = dble(totcdn)
!  nucleation coefficients
            aexpd = dble(aexp)
            rtso4d = dble(rtso4)

            so4t(1) = so4(totcdnd, afak, dble(so4t0), aexpd, rtso4d, condnu_dp(1))
            do ii = 2, 15
               so4td = dble(so4t(ii - 1))
               so4t(ii) = so4(totcdnd, afak, so4td, aexpd, rtso4d, condnu_dp(ii))
            end do

!  Update which occurred here takes place later, after a total mass check (P.A. Makar, March 2009)
!  xrow(i, l, igs_SO4) = so4t(15)
            newso4 = so4t(15)
!
!  for nucleation
            if (afak > 0.0d0) then
               blg = real(log(afak))
               nuclsum = 0.0
               do ii = 1, 15
                  nuclsum = nuclsum + &
                            (condnu(ii)* exp(blg + aexp * log(so4t(ii))))
               end do
               rtnucl(i, l) = nuclsum / chm_timestep
            end if

!  for condensation
            condsum = 0.0
            do ii = 1, 15
               condsum = condsum + (condnu(ii) * so4t(ii))
            end do
            totcd = totcdn * condsum / chm_timestep
!
!  percentage of so4 condensed to each bin
            do n = 1, isize
               if (aeronum(i, l, n) * roarow(i, l) >  xiao) then
                   rtcond(i, l, n) = pcond(i, l, n) * totcd
               else
                   rtcond(i, l, n) = 0.0
               end if
            end do
!
!  At this stage, the calculated values of the change in H2SO4(g) { = newso4 - oldso4},
!  the change of sulphate particle mass due to nucleation { = rtnucl * chm_timestep },
!  and the change in sulphate particle mass due to condensation { = rtcond * chm_timestep },
!  are all scaled to ensure that the total increase in sulphate mass is =
!  the sulphate mass available via gas-phase chemistry { = -98.1 / 64.6 * rtso2(i, l) * chm_timestep).

            rtcondsum = 0.0
            do n = 1, isize
               rtcondsum = rtcondsum + rtcond(i, l, n)
            end do

            sumso4 = max(newso4 - oldso4, 0.0) + rtnucl(i, l) * chm_timestep + rtcondsum * chm_timestep
            avail_so4 = rtso4 * chm_timestep
            if (sumso4 > 1.E-27) then
               ratio = avail_so4 / sumso4
            else
               ratio = 1.0
            end if

! Scale nucleation rate:
            rtnucl(i, l) = rtnucl(i, l) * ratio
!
! Scale condensation rates:
            rtcond(i, l, 1:isize) = rtcond(i, l, 1:isize) * ratio
!
!  Scale mass change of H2SO4 gas:
            newso4 = oldso4 + max(newso4 - oldso4, 0.0) * ratio
!
!  Carry out operations postponed until mass conservation was ensured:
            xrow(i, l, igs_SO4) = newso4
!  production of new particles by nucleation
            inew = isize * (iae_SU - 1) + 1 + (iae1 - 1)
            xrow(i, l, inew) = xrow(i, l, inew) + rtnucl(i, l) * chm_timestep
!
!  don't loose h2so4, if concentration is too low.
         else
            xrow(i, l, igs_SO4) = xrow(i, l, igs_SO4) + rtso4 * chm_timestep
!  The following is necessary to conserve sulphate mass in the above step:
            rtcond(i, l, 1:isize) = 0.0

         end if
      end do
   end do

!  call to compute the intersectional transport due to condensation process.
   call mach_cam_intrsec_outer(xrow, iae_SU, rtcond, aeronum, mae, pres, &
                               tp, roarow, pni, pnk)
   return
   
   contains
   
!  Function to compute so4 concentration due to combined
!  nucleation and condensation processes as a function of time.
   real(kind=4) function so4(a, b, c, d, pp, t0)
      implicit none
      real(kind=8), intent(in) :: a, b, c, d, pp, t0
      real(kind=8)             :: fomega, expp, so4tt
      real(kind=4), parameter  :: akb = 1.380662e-23
      
      fomega = a + b * c ** (d - 1.0d0)
      expp   = exp(-fomega * t0)
      so4tt  = ((1.0d0 - expp) * pp + expp * fomega * c) / fomega
      so4    = max(akb, real(so4tt))
      
      return
   end function so4
   
end subroutine mach_cam_sulfate
