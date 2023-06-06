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
! Fichier/File   : mach_cam_sfss.ftn90
! Creation       : S. Gong, P. Huang, S. Gravel and B. Pabla for GEM-MACH, Sept. 2008
!
! Description    : A subroutine to evaluate the generation rate of sea salt aerosols
!                  as a function of meteorologic conditions.
!
! Extra info     : - First version created by S. Gong Nov 21 1994 for CAM
!                    Aerosol surface flux computation for sea-salt
!                    Method:
!                      1.  Monahan et al (Oceanic Whitecaps, E.C. Monahan
!                          and G.Mac Niocail (eds), 1986)
!                      2.  Smith et al.
!                      3.  Gong, Global Biogeochemical Cycles 2003
!                  - As in previous SFSALT but optimized by vectorization
!                    (S. Gong, Dec 19, 1996)
!                  - New landuse catagory fland array is used to replace gc.
!                    the fland is based upon eos 1x1 km resolution data sets.
!                    (S. Gong, Dec 12, 1997)
!                  - Flux integration over each size bin is done now when model
!                    starts or restarts. the integrated efficient is saved in
!                    fintrow. (S. Gong Jun 11, 1998)
!                  - Add the Smith's sea-salt surface flux function.
!                    (S. Gong, Sep 02, 1999)
!                  - The surf zone production of sea spray function of de Leeuw
!                     et al (2001) was added with ocean cover less than 20% grid.
!                    (S. Gong, Aug 01, 2001)
!                  - The particle mass of each sub-devided bin is now computed
!                    in the do pos=1,n*35 loop: pmass. it was done on an averaged
!                    mass in a bin. Version 5C. (S. Gong, Oct 31, 2001)
!                  - A new scheme - GONG-MONAHAN [Gong, Global Biogeochemical
!                    Cycles,2003] was used re replace the old Monahan's scheme so
!                    that the appliable range of sea-salt production is expanded
!                    to less than 0.1 miron. (S. Gong, Nov 04, 2002)
!
! Arguments:  IN
!                tdiag      -> Diagnostic level temperature
!                surfwd     -> surface wind module
!                fland      -> Landuse
!
!             OUT
!                RSFROW     -> sea salt surface flux
!
!============================================================================
!
!!if_on
subroutine mach_cam_sfss (tdiag, surfwd, rsfrow, fland, pni)
   use mach_cam_utils_mod, only: isize
   use mach_drydep_mod,    only: lucprm
!!if_off
   use chm_nml_mod,        only: chm_nosurfzoneflux_l, chm_seaflux_s
   use mach_cam_utils_mod, only: rhop0, iae_SS, binrange, ursv
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: pni
   real(kind=4),    intent   (in) :: tdiag (pni)
   real(kind=4),    intent   (in) :: surfwd(pni)
   real(kind=4),    intent  (out) :: rsfrow(pni, isize)
   real(kind=4),    intent   (in) :: fland (pni, lucprm)
!!if_off
!
! Local variables
!
   integer(kind=4) :: n, i, pos
   real(kind=4)    :: rwi, pmass, rr1, rr2, cc3, del, delr0
   real(kind=4)    :: a1, a2, b1, b2
   real(kind=4), parameter :: log08 = log10(0.8), cub = 1.0 / 3.0
   real(kind=4), dimension(pni) :: r1, r2, df0, df1, df2
   real(kind=4), dimension(pni, isize, 3) :: fintrow
   real(kind=4), dimension(4, 4), parameter :: cc = reshape((/ &
                         0.7674, 3.079, 2.572e-11, -1.424, &
                         0.3926, 3.101, 4.190e-11, -1.404, &
                         0.2789, 3.115, 5.415e-11, -1.399, &
                         0.4809, 3.082, 3.110e-11, -1.428 /), (/4, 4/))

   rsfrow = 0.0

!  otherwise compute the integration of sea-salt for this model run - fintrow
   fintrow = 0.0

   do n = 1, isize
      rr1 = binrange(1, n) * 1.0e6     !from m to um
      rr2 = binrange(2, n) * 1.0e6
      a1 = cc(1, 1) * rr1 ** cc(2, 1)
      a2 = cc(1, 1) * rr2 ** cc(2, 1)
      b1 = rr1 ** cc(4, 1)
      b2 = rr2 ** cc(4, 1)
      do i = 1, pni
!
!  convert the dry size of sea-salt aerosols into wet size by taking rh into account.
!  The wet size (80% rh) is used in the surface flux calculation.
!
!  temperature effect

         cc3   = cc(3, 1) * (1.0 + 0.004 * (298.0 - tdiag(i)))
         r1(i) = (rr1 ** 3 + a1 / (cc3 * b1 - log08)) ** cub
         r2(i) = (rr2 ** 3 + a2 / (cc3 * b2 - log08)) ** cub

!  df0 is for the indirect mechanism (via bubbles)

         df0(i) = 0.0
         df1(i) = 0.0
         df2(i) = 0.0
      end do
!
      do pos = 1, n * 35
         do i = 1, pni
            del = (rr2 - rr1) / real(n * 35)           !del in um
            rwi = (rr1 + del * real(pos - 1)) * 1.e-6  !from um to m

!  dry volume and mass of an aerosol particle of size rwi

            pmass = ursv * (rwi * rwi * rwi) * rhop0(iae_SS)
            delr0 = (r2(i) - r1(i)) / real(n * 35)

!  gc(pni) = 1 for sea-ice excluded from sea-salt production
!  NOTE: Sea-ice presence accounted for already in fland in mach_landuse
            if (fland(i, 14) > 0.0) then
               if (chm_seaflux_s(1:12) == 'GONG_MONAHAN') then  !gong's scheme
                  df0(i) = df0(i) + delr0 * pmass * flx0(r1(i) + delr0 * real(pos - 1))

               else                     !smith's scheme
                  df0(i) = df0(i) + delr0 * pmass * flx1(r1(i) + delr0 * real(pos - 1))
                  df1(i) = df1(i) + delr0 * pmass * flx2(r1(i) + delr0 * real(pos - 1))
               end if

!  surf zone production

               if ((.not.chm_nosurfzoneflux_l) .and. fland(i, 14) < 1.0) then
                  df2(i) = df2(i) + delr0 * 2.0 * pmass * flx3(r1(i) * 2.0 &
                         + delr0 * 2.0 * real(pos - 1))
               end if
            end if
         end do
      end do
      do i = 1, pni
         fintrow(i, n, 1) = df0(i)
         fintrow(i, n, 2) = df1(i)
         fintrow(i, n, 3) = df2(i)
      end do
   end do
!
!  insertion du flux de surface
!
   do i = 1, pni
      if (surfwd(i) > 0.0 .and. fland(i, 14) > 0.0) then
         cc3 = exp(0.23 * surfwd(i))
         if (chm_seaflux_s(1:12) == 'GONG_MONAHAN') then  !gong-monahan's scheme
            a1 = surfwd(i) ** 3.41
            a2 = 0.0
         else !smith's scheme
            a1 = 10 ** (0.0676 * surfwd(i) + 2.430)
            a2 = 10 ** (0.9590 * sqrt(surfwd(i)) - 1.476)
         end if
         do n = 1, isize
            rsfrow(i, n) = (a1 * fintrow(i, n, 1) + a2 * fintrow(i, n, 2)) * &
                              fland(i, 14)
!  calculating surfing production flux
            rsfrow(i, n) = rsfrow(i, n) + fintrow(i, n, 3) * cc3 * fland(i, 14)
         end do
      end if
   end do
   return
   
   contains
   
!  (1) Gong-Monahan's source
   real(kind=4) function flx0(x)
      implicit none
      real(kind=4), intent(in) :: x
      real(kind=4)             :: b
      
      b  = (0.433333-log10(x)) / 0.433333
      flx0 = 1.373 * x ** (-4.7 * (1.+35.0 * x) ** (-0.017 * x ** (-1.44))) * &
             (1.0 + 0.057 * x ** 3.45) * 10.0 ** (1.607 * exp( - b * b))
      return
   end function flx0
      
!  (2) Smith's source
   real(kind=4) function flx1(x)
      implicit none
      real(kind=4), intent(in) :: x
       
      flx1 = exp(-3.1 * (log(x / 2.1)) ** 2)
      return
   end function flx1
      
   real(kind=4) function flx2(x)
      implicit none
      real(kind=4), intent(in) :: x

      flx2 = exp(-3.3 * (log(x / 9.2)) ** 2)
      return
   end function flx2
      
!  (3)  Surf zone production
   real(kind=4) function flx3(x)
      implicit none
      real(kind=4), intent(in) :: x

      flx3 = 1.1 * x ** (-1.65)
      return
   end function flx3
      
end subroutine mach_cam_sfss
