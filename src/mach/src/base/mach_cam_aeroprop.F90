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

! !============================================================================!
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
! Fichier/File   : mach_cam_aeroprop.ftn90
! Creation       : S. Gong, G. Lesins, V. Bouchet, C. Stroud, S. Gravel, and
!                  B. Pabla for GEM-MACH, June 2008
! Description    : Aerosol property calculation
!
! Extra info     : - First version created by S. Gong Dec 05 1994 for CAM
!                  - Vectorized the whole program and add working spaces.
!                    (S. Gong, Jan 19, 1996)
!                  - Combine previous aeroden and growth routines to compute
!                    the aerosol ambient properties. (s. gong, jun 25, 1997)
!                  - Modified parameterization valid for all sizes and rh
!                    uses empirical osmotic coefficients (g. lesins, mar 13, 1998)
!                  - Added nitrate, bc, dust added deliquescence & recryst data
!                    (G. Lesins, Jun 19, 1998)
!                  - Added organic aerosol type - treated as sulphate
!                    (G. Lesins, May 05, 1999)
!                  - Add the aerosol gravitational settling velocity calculation
!                    which was done several times in other places.
!                    (S. Gong, Sept 04, 1999)
!                  - Split aerosol oc into primary and secondary components
!                    (C. Stroud, Jul 2004)
!                  - Changed OC Water Uptake Parameterization (C. Stroud, Jun 2005)
!                  - Per suggestion of Paul Makar, Verica Savic-Jovcic modified
!                    calculation of particle settling velocity to avoid supersonic
!                    velocities for large particles (copied Paul's code from v1 in Oct 2015)
!                  - Re-order aerosol components as used in other parts of CAM, 
!                    as specified in mach_cam_utils_mod (aeroname). Adjust the
!                    parameters accordingly (A. Akingunola 2017)
!                  - TOTMAS (total mass mixing ratio for each bin, all components)
!                    and TRWTROW (aerosol liquid water content for each bin) are
!                    now local/optional variables. (A. Akingunola 2017)
!
! Arguments:  IN
!               ntr     -> Total number of trace substances (gases and aerosols)
!             rhrow     -> Relative humidity
!             throw     -> Temp
!             roarow    -> Air density (kg/m3)
!             rgrid     -> Mass mixing ratio for each trace substance
!               amu     -> Air's dynamic viscosity
!              amfp     -> Mean molecular free path
!
!            OUT
!              RHSIZE   -> Wet radius
!              RHOP     -> Final wet density
!              AERONUM  -> Number conconcenration (#/kg)
!              PDIFF    -> Diffusion coefficients
!              PDEPV    -> gravitational settling velocity
!              TRWTROW  -> Aerosol-bound water concentration
!
!            The remaining arguments are working arrays added by lesins
!                  bmix -> solute coefficient in kohler equation
!                  fr1  -> estimated radius ratio at relative humidity=1
!                  fc   -> critical relative humidity (supersaturated)
!
!============================================================================
!
!!if_on
subroutine mach_cam_aeroprop(rhsize, rhop, rhrow, throw, rgrid, aeronum,  &
                             pni, pnk, amu, amfp, roarow, pdiff, pdepv, trwtrow)
   use mach_cam_utils_mod, only: isize, ntr
!!if_off
   use chm_consphychm_mod, only: avno, grav, pi, rgasi, rgasv
   use mach_cam_utils_mod, only: icom, iae1, rhop0, pvol, binrange, &
                                 numsol, nu_sol, mw_sol, deliq, recry, phit
   implicit none
!!if_on
   integer(kind=4), intent (in)           :: pni, pnk
   real(kind=4),    intent(out)           :: rhsize (pni, pnk, isize)
   real(kind=4),    intent(out)           :: rhop   (pni, pnk, isize)
   real(kind=4),    intent (in)           :: rhrow  (pni, pnk)
   real(kind=4),    intent (in)           :: throw  (pni, pnk)
   real(kind=4),    intent (in)           :: rgrid  (pni, pnk, ntr)
   real(kind=4),    intent(out)           :: aeronum(pni, pnk, isize)
   real(kind=4),    intent (in), optional :: amu    (pni, pnk)
   real(kind=4),    intent (in), optional :: amfp   (pni, pnk)
   real(kind=4),    intent (in), optional :: roarow (pni, pnk)
   real(kind=4),    intent(out), optional :: pdiff  (pni, pnk, isize)
   real(kind=4),    intent(out), optional :: pdepv  (pni, pnk, isize)
   real(kind=4),    intent(out), optional :: trwtrow(pni, pnk, isize)
!!if_off
!
!  local variables
   integer(kind=4)            :: nt, n, no, l, i, iter
!  OC factor is a factor used to multiply the organic carbon osmotic
!  coefficients to force a difference from sulphate
   real(kind=4), parameter    :: aa1 = 1.257, aa2 = 0.4, aa3 = 1.1,  &
                                 cub = 1.0 / 3.0
   real(kind=4)               :: tramass, theta, dsr, frctest, awx1, &
                                 frx1, frx3
   real(kind=4)               :: bpr, vv, dd, ff, prii, priiv, cfac, &
                                 taurel, amob, boltzk

! following section contains data arrays for all aerosol types
   real(kind=4)               :: phik(numsol)

! mww is the molecular weight of water
! denw is the density of water                        ( kg / m^3 )
! sfcten is the surface tension between water and air ( j / m^2 )
   real(kind=4), parameter  :: mww = 18.015, denw = 1000.0, sfcten = 0.076
   real(kind=8), parameter  :: sixth = 1.0d0 / 6.0d0
   real(kind=4)     :: avesize(isize)
   real(kind=8)     :: r, q, d
   real(kind=8)     :: xx, yy, np
   real(kind=4)     :: sig, delrho, delp
! numiter is the number of iterations performed for kohler equation
! numiter=2 should give better than 1% accuracy
! numiter must not exceed the declared size of the last index in awx
   integer(kind=4) :: numiter = 2
!  deliqcry is a logical switch to handle size if rh is between
!  crystallization and deliquescence.  if =.true. then final size is
!  weighted average between dry and deliquesced sizes.  if =.false.
!  then final size is the fully deliquesced size.
   logical(kind=4) :: deliqcry =.true.
   real(kind=4)    :: phix  (pni, pnk, isize)
   real(kind=4)    :: fmso  (pni, pnk, isize)
   real(kind=4)    :: phiaw1(pni, pnk, isize)
   real(kind=4)    :: amw   (pni, pnk, isize)
   real(kind=4)    :: bmix  (pni, pnk, isize)
   real(kind=4)    :: fr1   (pni, pnk, isize)
   real(kind=4)    :: fc    (pni, pnk, isize)
   real(kind=4)    :: anu   (pni, pnk, isize)
   real(kind=4)    :: a     (pni, pnk, isize)
   real(kind=4)    :: frc   (pni, pnk, isize)
   real(kind=4)    :: deliqs(pni, pnk, isize)
   real(kind=4)    :: recrys(pni, pnk, isize)
   real(kind=4)    :: fmo   (pni, pnk, numsol, isize)
   real(kind=4)    :: awx   (pni, pnk, isize, 2)
   real(kind=4)    :: totvol(pni, pnk, isize)
   real(kind=4)    :: totmas(pni, pnk, isize)
   real(kind=4)    :: rh    (pni, pnk)
   logical(kind=4) :: crt_frc_fc = .false.

   amw      = 0.0
   anu      = 0.0
   deliqs   = 0.0
   fmso     = 0.0
   phiaw1   = 0.0
   recrys   = 0.0
   totvol   = 0.0
   totmas   = 0.0

!
! Set bounds for relative humidity (note: this is already done in mach_pm_chem)
   rh = min(1.0, max(0.0, rhrow))
!
!  total dry mass mixing ratio & dry aerosol composite density of aerosol in each bin
!
!  rho = (m1 + m2 + m3) / (m1 / rho1 + m2 / rho2 + m3 / rho3)
!  in this do loop, rhop(*, *, *) is only holding the denominator of this equation

   do n = 1, isize
      avesize(n) = (binrange(1, n) + binrange(2, n)) * 0.5
      do nt = 1, icom
         no = n + isize * (nt - 1) + (iae1 - 1)
         do l = 1, pnk
            do i = 1, pni
               tramass = max(1.0e-33, rgrid(i, l, no))
               totmas(i, l, n) = totmas(i, l, n) + tramass
               totvol(i, l, n) = totvol(i, l, n) + tramass / rhop0(nt)
            end do
         end do
      end do
!            
! compute the mass fraction of each dry aerosol component, fmo
      do nt = 1, numsol
         no = n + isize * (nt - 1) + (iae1 - 1)
         do l = 1, pnk
            do i = 1, pni
               tramass = max(1.0e-33, rgrid(i, l, no))
               fmo(i, l, nt, n) = tramass / totmas(i, l, n)
            end do
         end do
      end do
   end do

!
!  compute phi at aw=1 for each solute
!  compute kohler b factor for each solute
   do nt = 1, numsol
      phik(nt) = phit(3, nt) + phit(4, nt) + phit(5, nt) + phit(6, nt)
   end do

! compute the mass fraction of each dry aerosol component, fmo
! compute the solute mass fraction, fmso
! compute average nu and average molecular weight of soluble part
! compute phi at aw=1 for the mixed aerosol
! compute average deliquescence and recrystallization points
   do n = 1, isize
      do nt = 1, numsol
         do l = 1, pnk
            do i = 1, pni
               fmso(i, l, n) = fmso(i, l, n) + fmo(i, l, nt, n)
               amw(i, l, n) = amw(i, l, n) + fmo(i, l, nt, n) / mw_sol(nt)
               anu(i, l, n) = anu(i, l, n) + &
                              nu_sol(nt) * fmo(i, l, nt, n) / mw_sol(nt)
               phiaw1(i, l, n) = phiaw1(i, l, n) + phik(nt) * &
                                 nu_sol(nt) * fmo(i, l, nt, n) / mw_sol(nt)
               deliqs(i, l, n) = deliqs(i, l, n) + fmo(i, l, nt, n) * deliq(nt)
               recrys(i, l, n) = recrys(i, l, n) + fmo(i, l, nt, n) * recry(nt)
            end do
         end do
      end do
   end do

   where (fmso /= 0.0)
      amw = fmso / amw
      anu = anu * amw / fmso
      phiaw1 = phiaw1 * amw / (fmso * anu)
      deliqs = deliqs / fmso
      recrys = recrys / fmso
   elsewhere
      deliqs = 1.1
      recrys = 1.1
   end where
!
!  compute aerosol number concentration (#/kg_air)
!  compute kohler a' factor and b factor for mixture
!  compute radius ratio at rh=1 (fr1)
   do n = 1, isize
      do l = 1, pnk
         do i = 1, pni
            rhop(i, l, n) = totmas(i, l, n) / totvol(i, l, n)
            aeronum(i, l, n) = totvol(i, l, n) / pvol(n)

            a(i, l, n) = 2.0 * sfcten / (denw * rgasv * throw(i, l) * avesize(n))
            if (amw(i, l, n) /= 0.0) then
               bmix(i, l, n) = anu(i, l, n) * phiaw1(i, l, n) * mww * rhop(i, l, n) *  &
                               fmso(i, l, n) / (amw(i, l, n) * denw)
            else
               bmix(i, l, n) = 0.0
            end if
            q = -bmix(i, l, n) / (3.0 * a(i, l, n))
            d = q * q * q + 0.25d0
            if (d < 0.0 .and. q < 0.0) then
               theta = real(acos(0.5 / sqrt(-q * q * q)))
               fr1(i, l, n) = real(2.0 * sqrt(-q) * cos(theta / 3.0))
            else
               dsr = real(sqrt(d))
               fr1(i, l, n) = real((0.5 + dsr) ** cub + (0.5 - dsr) ** cub)
            end if
         end do
      end do
   end do

!  this section can be commented out if you are not interested
!  in computing the critical radius and critical relative humidity
   if (crt_frc_fc) then
!  compute critical radus ratio (frc)
!  compute critical relative humidity (fc)
      do n = 1, isize
         do l = 1, pnk
            do i = 1, pni
               q = -bmix(i, l, n) / a(i, l, n)
               d = q * q * q + 1.0d0
               if (d < 0.0 .and. q < 0.0) then
                  theta = real(acos(1.0 / sqrt(-q * q * q)))
                  frc(i, l, n) = real(2.0 * sqrt(-q) * cos(theta / 3.0))
               else
                  dsr = real(sqrt(d))
                  frc(i, l, n) = real((1.0 + dsr) ** cub + (1.0 - dsr) ** cub)
               end if
               frctest = frc(i, l, n) ** 3.0
               if (frctest == 1.0) then
                  fc(i, l, n) = 1.0
               else
                  fc(i, l, n) = exp(a(i, l, n) / frc(i, l, n) - bmix(i, l, n) /  &
                                (frctest - 1.0))
               end if
            end do
         end do
      end do
   end if

!
! solve for the wet radius after water uptake
   do iter = 1, numiter
      phix = 0.0
      do n = 1, isize
         do nt = 1, numsol
            do l = 1, pnk
               do i = 1, pni
                  if (iter == 1) then
                     awx(i, l, n, iter) = rh(i, l) * exp(-a(i, l, n) / (0.8 * fr1(i, l, n)))
                     if (awx(i, l, n, iter) > 1.0) awx(i, l, n, iter) = 1.0
                     awx1 = awx(i, l, n, iter)
                  else
                     awx1 = awx(i, l, n, iter - 1)
                  end if
                  if (awx1 > phit(2, nt)) then
                     phix(i, l, n) = phix(i, l, n) + (phit(3, nt) * awx1 * awx1 * awx1 +  &
                                     phit(4, nt) * awx1 * awx1 + phit(5, nt) * awx1 +     &
                                     phit(6, nt)) * nu_sol(nt) * fmo(i, l, nt, n) / mw_sol(nt)
                  else
                     phix(i, l, n) = phix(i, l, n) +                            &
                                     (phit(7, nt) * awx1 * awx1 * awx1 * awx1 + &
                                      phit(8, nt) * awx1 * awx1 * awx1 +        &
                                      phit(9, nt) * awx1 * awx1 +               &
                                      phit(10, nt) * awx1 + phit(11, nt)) *     &
                                      nu_sol(nt) * fmo(i, l, nt, n) / mw_sol(nt)
                  end if
               end do
            end do
         end do

         do l = 1, pnk
            do i = 1, pni
               if (iter == 1) then
                  awx1 = awx(i, l, n, iter)
               else
                  awx1 = awx(i, l, n, iter - 1)
               end if
               if (fmso(i, l, n) /= 0.0 .and. awx1 > 0.0) then
                  phix(i, l, n) = phix(i, l, n) * amw(i, l, n) / (fmso(i, l, n) * anu(i, l, n))
                  frx1 = (1.0 - anu(i, l, n) * phix(i, l, n) * mww * rhop(i, l, n) *  &
                          fmso(i, l, n) / (amw(i, l, n) * denw * log(awx1))) ** cub
               else
                  phix(i, l, n) = 0.0
                  frx1 = 1.0
               end if

!  this section greatly improves accuracy for rh ~ 1
!  it can be commented out if speed is more important than accuracy

               if (rh(i, l) > 0.98 .and. frx1 > 1.0) then
                  frx3 = frx1 * frx1 * frx1
                  bpr = anu(i, l, n) * phix(i, l, n) * mww * rhop(i, l, n) *  &
                        fmso(i, l, n) * frx3 / (amw(i, l, n) * denw * (frx3 - 1.0))
                  q = dble(-a(i, l, n) / (3.0 * bpr))
                  r = dble(-log(rh(i, l)) / (2.0 * bpr))
                  d = q * q * q + r * r
                  if (d < 0.0 .and. q < 0.0) then
                     theta = real(acos(r / sqrt(-q * q * q)))
                     frx1 = real(1.0 / (2.0 * sqrt(-q) * cos(theta / 3.0)))
                  else
                     vv = real(abs(r + sqrt(d)))
                     dd = real(abs(r - sqrt(d)))
                     frx1 = real(1.0 / ((vv) ** cub + (dd) ** cub))
                  end if
               end if
               awx(i, l, n, iter) = rh(i, l) * exp(-a(i, l, n) / frx1)
               rhsize(i, l, n) = avesize(n) * frx1
            end do
         end do
      end do
   end do

! adjust size if rh < deliquescence point
! compute the mean density of the wet aerosol

   do n = 1, isize
      do l = 1, pnk
         do i = 1, pni
            if (rh(i, l) < deliqs(i, l, n)) then
               if (rh(i, l) < recrys(i, l, n)) then
                  rhsize(i, l, n) = avesize(n)
               else
                  if (deliqcry) then
                     rhsize(i, l, n) = avesize(n) + (rhsize(i, l, n) - avesize(n)) *  &
                     (rh(i, l) - recrys(i, l, n)) / (deliqs(i, l, n) - recrys(i, l, n))
                  end if
               end if
            end if
            ff = avesize(n) / rhsize(i, l, n)
            rhop(i, l, n) = denw + max(0.0, ff * ff * ff * (rhop(i, l, n) - denw))
         end do
      end do
   end do

! Evaluate water-bound aerosol for estimating aerosol optical properties and direct feedback effects
!
   if (present(trwtrow)) then
! compute the aerosol liquid water content of each size bin [kg/kg air]
      trwtrow = max(0.0, ((4.189 * rhsize**3 * aeronum * rhop) - totmas))
      return
   end if
   
   boltzk    = rgasi / avno
   do n = 1, isize
      do l = 1, pnk
         do i = 1, pni
! aerosol gravitational settling velocity and diffusion coefficient
            prii = 2.0 / 9.0 * grav / amu(i, l)
            priiv = prii * (rhop(i, l, n) - roarow(i, l))

! cunningham slip correction factor and relaxation time = vg/grav.
            cfac = 1.0 + amfp(i, l) / rhsize(i, l, n) * (aa1 + aa2 *  &
                   exp(-aa3 * rhsize(i, l, n) / amfp(i, l)))

! stokes friction and diffusion coefficients.
            amob = 6.0 * pi * amu(i, l) * rhsize(i, l, n) / cfac
            pdiff(i, l, n) = boltzk * throw(i, l) / amob

! gravitational settling velocity.
!
!   The formula for the settling or terminal velocity in the following section of code is
!   based on Beard, K.V., "Terminal Velocity and Shape of Cloud and Precipitation Drops Aloft",
!   Journal of the Atmospheric Sciences, vol 33, pp 851-864, 1976, formulae in Table 1, with the
!   exception of the wet diameter size < 19.E-06, wherein the full Cunningham slip correction
!   was used (Beard's equation (3)) rather than his approximation in Table 1.
!
            delp = 2.0 * rhsize(i, l, n)
            delrho = rhop(i, l, n) - roarow(i, l)
            if (delp <= 19.e-06) then
!  Small particles, diameters less than 19 um:
               taurel = max(priiv * rhsize(i, l, n) ** 2 * cfac / grav, 0.0)
               pdepv(i, l, n) = taurel * grav
            else
               if (delp <= 1.07e-03) then
!  Larger particles (e.g. with significant water uptake), diameters between 19 um and 1.07 mm:
                  xx = log(dble( 4. * roarow(i, l) * delrho * grav &
                         * delp ** 3 / (3. * amu(i, l) * amu(i, l)) ))
                  yy = -3.18657d+00 + ((((((-0.327815d-05 * xx) + 0.855176d-04) * xx &
                                    - 0.578878d-03 ) * xx - 0.987059d-03 ) * xx &
                                    - 0.153193d-02 ) * xx + 0.992696d+00 ) * xx
                  pdepv(i, l, n) = amu(i, l) * cfac / (roarow(i, l) * delp ) * real(exp(yy))
               else
!
!  Particles with diameters > 1.07 mm:  droplets; surface tension effects must be considered:
!  SIG = surface tension of water, N/m
!  Surface tension temperature dependance based on data from:
!  http://www.engineeringtoolbox.com/water-surface-tension-d_597.html, accessed June 20, 2011.
!
                  sig = -1.66912980e-04 * throw(i, l) + 1.21570339e-01
!
                  np = dble(sig ** 3 * roarow(i, l) ** 2 / (amu(i, l) ** 4 * delrho * grav))
                  np = np ** sixth
                  xx = log(dble(4. * delrho * grav * (delp ** 2) / (3. * sig)) * np)
                  yy = -5.00015d+00 + ((((( 0.238449d-02 * xx) - 0.542819d-01) * xx &
                                + 0.475294d+00) * xx - 0.204914d+01) * xx           &
                                + 0.523778d+01) * xx
                  pdepv(i, l, n) = amu(i, l) * real(np * exp(yy)) / (roarow(i, l) * delp)
               end if
            end if

         end do
      end do
   end do

   return
end subroutine mach_cam_aeroprop
