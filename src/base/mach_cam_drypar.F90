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
! Fichier/File   : mach_cam_drypar.ftn90
! Creation       : S. Gong, L. Zhang, S. Gravel and B. Pabla for GEM-MACH, June 2008
! Description    : Calculation of dry deposition of aerosol particles.
!
! Extra info     : - First version created by S. Gong Aug 04 1994 for CAM
!                    Method:
!                       parameterization of dry deposition scheme by Filippo Giorgi
!                       [JGR, 91(D9), 1986] was used.
!                       Two layer model.
!                  - Vectorized the whole program and added working spaces.(S. Gong, Dec 19, 1996)
!                  - Modify dry deposition parameterization (L. Zhang, Mar 19, 1998)
!                  - Pass all control variables for 15 land-use from the control file (S. Gong, Aug 01, 1999)
!                  - The gravitational velocity is now done in aeroprop (S. Gong, Sep 04, 1999)
!                  - ra, usi, seasn are now calculated in drygas (L. Zhang, Sep 10, 2000)
!
! Arguments:  IN
!               thlev   -> Layer thickness [m]
!               roarow  -> Air density (kg/m3)
!               pdiff   -> diffusion coefficient
!               rhsize  -> Wet radius
!               fland   -> Landuse
!               xrow    -> tracer concentration
!               iseasn  -> Assigned season descriptors
!               Ra      -> Aerodynamic resistance for species "species" (s/m)
!               usi     -> friction velocity
!               amu     -> Air's dynamic viscosity
!
!             IN/OUT
!              PDEPV    -> gravitational settling velocity
!
!============================================================================
!
!!if_on
subroutine mach_cam_drypar(thlev, roarow, pdiff, rhsize, amu, fland, &
                           xrow, pdepv, iseasn, ra, usi, pni, pnk)
   use mach_cam_utils_mod, only: isize, ntr
   use mach_drydep_mod,    only: lucprm
!!if_off
   use chm_utils_mod,      only: chm_timestep, global_debug
   use chm_consphychm_mod, only: grav
   use chm_nml_mod,        only: chm_intrsec_ndiv
   use mach_cam_utils_mod, only: nsb, icom, iae1, nbnd, binrange, dvn, sub_rn
   use mach_drydep_mod,    only: aest, pgamma, pllp
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: thlev (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
   real(kind=4),    intent   (in) :: fland (pni, lucprm)
   real(kind=4),    intent   (in) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent(inout) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: ra    (pni, lucprm)
   real(kind=4),    intent   (in) :: usi   (pni)
   integer(kind=4), intent   (in) :: iseasn(pni)
!!if_off
!
! local variables
!
!SLG Sub V
   logical (kind=4) :: local_dbg
   integer (kind=4) :: l, n, i, ic, isea, np
   integer (kind=4) :: j, jsub, k, il, jk, pos
   real(kind=4)     :: st, eb, eim, ein, r1, rs, vdf
   real(kind=4)     :: vdp      (pni)
   real(kind=4)     :: anu      (pni)
   real(kind=4)     :: schm     (pni)
   real(kind=8)     :: tm       (pni, pnk, isize)
   real(kind=8)     :: tmr      (pni, pnk, nbnd)
   real(kind=4)     :: sub_Pdepv(pni, pnk, nbnd)
   real(kind=8)     :: vx       (pni, pnk, isize)

   real(kind=4)     :: aa, bb, v1, v2, cc
   real(kind=8)     :: rr1(isize)
   real(kind=8)     :: rdb, fcn1, fcn2
   real(kind=8)     :: smf1 = 1.D-20     !  minimum particle mass mixing ratio
   real(kind=8)     :: smf2 = 1.D-20     !  minimum value of vx, where vx = exp(- adt2/th * dep_vel)

   real(kind=4)     :: a1       (pni)
   real(kind=4)     :: b1       (pni)
   real(kind=4)     :: c1       (pni, pnk)
   real(kind=4)     :: sub_xrow (pni, pnk, nsb)
!SLG end

   local_dbg = (.false. .or. global_debug)

   sub_Pdepv = 0.0
!
!  loop for the number of aerosols (size bins)
!  for the surface layer, dry deposition is calculated.
!  note: the gravitational velocity was computed in mach_cam_aeroprop
!
! SLG - Sub dry
!
   if (nbnd > isize) then  ! case with redistribution
!
!  On input, pdepv contains the original bin structure settling velocities.
!  Here, these will be broken up over sub-bins using a simple scaling function.
!  settling velocity only (same as non-urban for the moment, sunling to provide
!  new function
!  The function cc is a polynomial fit to typical log10(settling velocity) as a
!  function of log10(radius):
      c1 = 0.0
      do n = 1, isize
!  Median size of the aerosols
         rr1(n) = dlog10(0.5d0 * dble(binrange(1, n) + binrange(2, n)))
!
         cc = real(-0.02614934d0 * rr1(n)**3 -0.37534441d0 * rr1(n)**2 + &
                 0.22485540d0 * rr1(n) + 5.77001978d0)
         do l = 1, pnk
            do i = 1, pni
               c1(i, l) = c1(i, l) + log10(pdepv(i, l, n)) - cc        !settling velocity
            end do
         end do
      end do

      c1 = c1 / real(isize)
!
!  c1 now contains the average offset between the idealized
!  simple function for {log10(settling velocity) = (cc)} and
!  the precalculated more detailed settling velocity, for
!  the original bin structure.  This offset will now be used
!  to improve the fit for the simple scaling function, for
!  all vertical levels
!  Note that this is applied only for the levels above the
!  ground level (pnk), since the latter is corrected using a
!  separate function.
      do n = 1, nbnd
         rdb = dble(log10(sub_rn(n)))
         fcn1 = -0.02614934d0 * rdb**3 -0.37534441d0 *rdb**2 + &
                 0.22485540d0 * rdb + 5.77001978d0
         do l = 1, pnk
            do i = 1, pni
! settling velocity, calculated using simple function, then offset
! according to more detailed calculation at the original bin distribution:
               sub_pdepv(i, l, n) = 10**(real( fcn1 + dble(c1(i, l))))
            end do
         end do
      end do

   end if ! nbnd > isize
!
! sub_Pdepv now contains the scaled settling velocities at the sub-bins
! Calculate the deposition velocities at the original bins:
   l = pnk
!
!  air's kinematic viscosity
   do i = 1, pni
      anu(i) = amu(i, l) / roarow(i, l)
   end do
!
   do n = 1, isize
      do i = 1, pni
         vdp(i) = 0.0
      end do
      do ic = 1, lucprm
         do i = 1, pni
            isea = iseasn(i)
            if (fland(i, ic) >= 0.005) then
               schm(i) = anu(i) / pdiff(i, l, n)

!  calculate middle-variles needed for eb, ein, em, etc., then ra, rs, vdf
               if (pllp(ic, isea) <= 0.0) then
                  st = pdepv(i, l, n) / grav * usi(i) * usi(i) / anu(i)
                  ein = 0.0
               else
                  st = pdepv(i, l, n) / grav * usi(i) / pllp(ic, isea) * 1000.0
                  ein = (1000.0 * 2.0 * rhsize(i, l, n) / pllp(ic, isea)) ** 2 * 0.5
                  ein = min(ein, 0.6)
               end if
               eb = schm(i) ** pgamma(ic)
               eim = (st / (st + aest(ic))) ** 2
               eim = min(eim, 0.6)
               r1 = max(0.5, exp(-st ** 0.5))
               rs = 1.0 / 3.0 / usi(i) / (eb + eim + ein) / r1
               vdf = pdepv(i, l, n) + 1.0 / (ra(i, ic) + rs)
               vdp(i) = vdp(i) + vdf * fland(i, ic)
            end if
         end do
      end do

      do i = 1, pni
         pdepv(i, l, n) = vdp(i)
      end do
   end do
!
!SLG Sub-dry calculations
!  a1 and b1 are the average displacement between the
!  3rd order polynomial describing
!  log10 (deposition velocity + settling velocity) as a
!  function of radius for the sub-bins, and the deposition
!  velocity for the full bin, at the two bin locations,
!  for all non-urban land uses and urban land uses, respectively.
!
   if (nbnd > isize) then  ! case with redistribution
      a1 = 0.0
      b1 = 0.0
      do n = 1, isize
!  non-urban
         aa = real(0.11325186d0 * rr1(n)**3 + 2.85005039d0 * rr1(n)**2 + &
                   22.72885201d0 * rr1(n) + 55.053510890d0)
!  urban
         bb = real(0.21069888d0 * rr1(n)**3 + 4.69614907d0 * rr1(n)**2 + &
                   33.94123621d0 * rr1(n) + 77.28509785d0)
         do i = 1, pni
            a1(i) = a1(i) + log10(pdepv(i, l, n)) - aa   !for all other landuse
            b1(i) = b1(i) + log10(pdepv(i, l, n)) - bb   !for urban category
         end do
      end do
!
      a1 = a1 / real(isize)
      b1 = b1 / real(isize)
!
! a1, b1 are now the average offset between the simple functions for
! log10(deposition velocity+settling velocity) as a function of size, and the
! detailed calculation at the original bin sizes.  The same functions will now
! be used to calculate the deposition velocity at the sub-bins, and the
! offset will be used to correct the simple parameterized curve for
! deposition velocity + settling velocity as a function of radius to
! the original bin values rough location at the original isize bins.
!
!
! calculate the deposition velocity + settling velocity for each sub-bin,
! adding in the offset to correct the values to the
! curves calculated in more detail for the full bins:
!
      l = pnk
      do n = 1, nbnd
         rdb = dble(log10(sub_rn(n)))
         fcn1 = 0.11325186d0 * rdb**3 + 2.85005039d0 * rdb**2 + &
                22.72885201d0 * rdb + 55.053510890d0
         fcn2 = 0.21069888d0 * rdb**3 + 4.69614907d0 * rdb**2 + &
                33.94123621d0 * rdb + 77.28509785d0
         do i = 1, pni
! non-urban surfaces
            v1 = 10**real(fcn1 + dble(a1(i)))
! urban surfaces
            v2 = 10**real(fcn2 + dble(b1(i)))
!  Weighting:  deposition velocity is assumed to dominate
!  in lowest model layer and is a function of the relative
!  urban versus non-urban land use categories.  For all
!  upper layers, settling velocity function is used
! weight by land use to get sub-bin deposition velocity:
            sub_Pdepv(i, l, n) = v1 * (1.0 - fland(i, 15)) + v2 * fland(i, 15)
         end do
      end do
!
!  Create deposition and fall velocities for original bin structure which
!  preserve the net total removal during the time step of the
!  sub-bin structure.  This is done to reduce subsequent computational load.
!
!  First, assign new particle mass into the sub-bins:
      do k = 1, isize
         do j = 1, chm_intrsec_ndiv
            jsub = (k - 1) * chm_intrsec_ndiv + j
            do np = 1, icom
               pos = (np - 1) * isize + k + iae1 - 1
               jk = (np - 1) * nbnd + jsub
               do il = 1, pnk
                  do i = 1, pni
                     sub_xrow(i, il, jk) = xrow(i, il, pos) * dvn(jsub)
                  end do
               end do
            end do
         end do
      end do
!
! Second, determine the total mass ratio at each subbin-size at each gridpoint:
!
      tm  = 0.d0  !  total mass at gridpoint across all species for each original bin size
      tmr = 0.d0  !  total mass ratio for each subbin size
      do k = 1, nbnd
         do np = 1, icom
            jk = (np - 1) * nbnd + k
            do il = 1, pnk
               do i = 1, pni
                  tmr(i, il, k) = tmr(i, il, k) + dble(sub_xrow(i, il, jk))
               end do
            end do
         end do
      end do
      do k = 1, isize
         do j = 1, chm_intrsec_ndiv
            jsub = (k - 1) * chm_intrsec_ndiv + j
            do il = 1, pnk
               do i = 1, pni
                  tm(i, il, k) = tm(i, il, k) + tmr(i, il, jsub)
               end do
            end do
         end do
      end do
!  Calculate total mass ratio at each gridpoint as a function of
!  size.  A lower cutoff is employed to prevent division by zero errors.
      do k = 1, isize
         do j = 1, chm_intrsec_ndiv
            jsub = (k - 1) * chm_intrsec_ndiv + j
            do il = 1, pnk
               do i = 1, pni
                  tmr(i, il, jsub) = tmr(i, il, jsub) / max(smf1, tm(i, il, k))
               end do
            end do
         end do
      end do
!  Determine mass weighted deposition velocities
!  Note that the weighting preserves the total mass
!  loss due to fall and/or deposition relative to the
!  sub-bin distribution:
      vx = 0.d0
      do k = 1, isize
         do j = 1, chm_intrsec_ndiv
            jsub = (k - 1) * chm_intrsec_ndiv + j
            do il = 1, pnk
               do i = 1, pni
                  vx(i, il, k) = vx(i, il, k) + tmr(i, il, jsub) * dexp(dble(-sub_pdepv(i, il, jsub) * chm_timestep / thlev(i, il)) )
               end do
            end do
         end do
      end do

!  Apply limiters: vx must be greater than zero and less than unity
!  in order for the resulting deposition velocity to be bounded:
      do k = 1, isize
         do il = 1, pnk
            do i = 1, pni
               vx(i, il, k) = min(vx(i, il, k), 1.0d0)
               vx(i, il, k) = max(vx(i, il, k), smf2)
            end do
         end do
      end do
! Replace original deposition velocities with recalculated ones:
      do k = 1, isize
         do il = 1, pnk
            do i = 1, pni
               if (vx(i, il, k) > smf2 .and. tm(i, il, k) > smf1) then
                  pdepv(i, il, k) = -thlev(i, il) / chm_timestep * real(dlog(vx(i, il, k)))
               end if
            end do
         end do
      end do
!   
   end if ! nbnd > isize
!
  return
end subroutine mach_cam_drypar
