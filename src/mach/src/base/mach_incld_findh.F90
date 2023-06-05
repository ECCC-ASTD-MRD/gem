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
! Fichier/File   : mach_incld_findh.ftn90
! Creation       : S. Menard,  W. Gong, S. Gravel, V. Bouchet, B. Pabla, GEM-MACH,  June 2008.
!
! Description    : Compute Array of Aq. and Gas species, aqueous rate constant array,
!                  aqueous variable coefficients, charge balance function.
!
! Extra info     :
!
! Arguments:
!           IN

!              nptsnz -->    Total number of grids to integrate
!
!           INOUT
!              T      -->    Initial component totals
!                            1-Total S(IV)
!                            2-Total NO3
!                            3-Total Ammonia
!                            4-Total Carbon Dioxide
!                            5-Net charge due to SO4= and CAT1
!              C      -->    Aq. proton concentration
!                            C= 0.1*[H+] or 10.0 * [H+]
!                            (units :  M/L)
!              R      -->    Aqueous rate constant array
!              B      -->    Aqueous variable coefficients
!
!=============================================================================
!
!!if_on
subroutine mach_incld_findh (T, C, R, B, nptsnz)
!!if_off
   use mach_incld_headers_mod, only: mach_incld_fff
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nptsnz
   real(kind=4),    intent(inout) :: t(nptsnz,5)
   real(kind=4),    intent(inout) :: c(nptsnz,12)
   real(kind=4),    intent(inout) :: r(nptsnz,25)
   real(kind=4),    intent(inout) :: b(nptsnz,5)
!!if_off
!
! local variables
!
   integer(kind=4), parameter             :: maxit = 100, maxitr = 10
   integer(kind=4)                        :: ivec, i, nleft1, isum, itd, itnd, ipos, ii, jk
   real(kind=4),    parameter             :: eps = 0.001
   real(kind=4)                           :: f2f0, dx21, dx20, dx
   integer(kind=4), dimension(nptsnz)     :: igrid, jgrid, iges, icheck, itmp, ichecktmp
   real(kind=4),    dimension(nptsnz)     :: x0tmp, x1tmp, x2tmp, f0tmp, f1tmp, check1,   &
                                             x0, x1, x2, x00, x10, f0, f1, f2, f00
   real(kind=4),    dimension(nptsnz, 5)  :: ttmp, btmp
   real(kind=4),    dimension(nptsnz, 25) :: rtmp

   do ivec = 1, nptsnz
      icheck(ivec) = 0
      iges(ivec) = 1
      check1(ivec) = 0.0
      x0(ivec) = 0.1 * c(ivec, 5)
      x1(ivec) = 10.0 * c(ivec, 5)
!  modifications to solver : 12.19.85 - pk
      x00(ivec) = x0(ivec)
      x10(ivec) = x1(ivec)
!  modifications end
   end do

   do i = 1, 5
      call mach_incld_fff(F0, T, X0, R, B, nptsnz, 1)
      call mach_incld_fff(F1, T, X1, R, B, nptsnz, 1)
!  modifications
      if (i == 1) then
         f00 = f0
      end if
!  modifications end
      do ivec = 1, nptsnz
         if (f0(ivec) * f1(ivec) >= 0.0) then
            iges(ivec) = iges(ivec) + 1
            x0(ivec) = 0.1 * x0(ivec)
            x1(ivec) = 10.0 * x1(ivec)
         end if
      end do
   end do

   where ((f0 * f1) >= 0.0) check1 = 1.0

!  modifications
   do ivec = 1, nptsnz
      if (iges(ivec) > 1) then
         if (f0(ivec) * f00(ivec) < 0.0) then
            x1(ivec) = x10(ivec)
         else
         x0(ivec) = x00(ivec)
         end if
      end if
   end do
!  modifications end

   nleft1 = 1
   do ivec = 1, nptsnz
      igrid(ivec) = ivec
   end do
   i = 1

   do while ((nleft1 - 1) < nptsnz .and. i <= maxit)
 
!  modifications
      if (i < maxitr) then
!  regula-falsi method
         do ivec = nleft1, nptsnz
            x2(ivec) = (x0(ivec) * f1(ivec) - x1(ivec) * f0(ivec)) / (f1(ivec) - f0(ivec))
!  setting a lower limit on trial proton concentration
            x2(ivec) = max(x2(ivec), 1.0e-12)
         end do
!  regula-falsi method
      else
!  half-interval method
         do ivec = nleft1, nptsnz
            x2(ivec) = 0.5 * (x0(ivec) + x1(ivec))
         end do
!  half-interval method
      end if
!  modifications end

      call mach_incld_fff(F2, T, X2, R, B, nptsnz, nleft1)
      do ivec = nleft1, nptsnz
         f2f0 = f2(ivec) * f0(ivec)
         if ((x2(ivec) > 0.0) .and. (f2(ivec) == 0.0)) icheck(ivec) = 1
         dx21 = abs((x2(ivec) - x1(ivec)) / x2(ivec))
         dx20 = abs((x2(ivec) - x0(ivec)) / x2(ivec))
         dx = min(dx21, dx20)
         if (dx <= eps) icheck(ivec) = 1
         if (f2f0 <= 0.0) then
            x1(ivec) = x2(ivec)
            f1(ivec) = f2(ivec)
         else
            x0(ivec) = x2(ivec)
            f0(ivec) = f2(ivec)
         end if
      end do
      isum = 0
      do ivec = nleft1, nptsnz
         isum = isum + icheck(ivec)
      end do
      if (isum > 0) then
!  rearrange array index
         do ivec = 1, nptsnz
            jgrid(ivec) = ivec
            itmp(ivec) = igrid(ivec)
         end do

         itd = nleft1 - 1
         itnd = nptsnz + 1
         do ivec = nleft1, nptsnz
            itd = itd + icheck(ivec)
            itnd = itnd + (icheck(ivec) - 1)
            ipos = itd * icheck(ivec) + itnd * (1 - icheck(ivec))
            jgrid(ipos) = ivec
         end do
!
!  shifting arrays
!
         do ii = 1, 5
            do ivec = nleft1, nptsnz
               jk = jgrid(ivec)
               btmp(ivec, ii) = b(jk, ii)
               ttmp(ivec, ii) = t(jk, ii)
            end do
            do ivec = nleft1, nptsnz
               b(ivec, ii) = btmp(ivec, ii)
               t(ivec, ii) = ttmp(ivec, ii)
            end do
         end do
         do ii = 1, 25
            do ivec = nleft1, nptsnz
               jk = jgrid(ivec)
               rtmp(ivec, ii) = r(jk, ii)
            end do
            do ivec = nleft1, nptsnz
               r(ivec, ii) = rtmp(ivec, ii)
            end do
         end do
         do ivec = nleft1, nptsnz
            jk = jgrid(ivec)
            ichecktmp(ivec) = icheck(jk)
            x1tmp(ivec) = x1(jk)
            x0tmp(ivec) = x0(jk)
            f0tmp(ivec) = f0(jk)
            f1tmp(ivec) = f1(jk)
            x2tmp(ivec) = x2(jk)
         end do
         do ivec = nleft1, nptsnz
            icheck(ivec) = ichecktmp(ivec)
            x1(ivec) = x1tmp(ivec)
            x0(ivec) = x0tmp(ivec)
            f0(ivec) = f0tmp(ivec)
            f1(ivec) = f1tmp(ivec)
            x2(ivec) = x2tmp(ivec)
         end do
         nleft1 = nleft1 + isum
         do ivec = 1, nptsnz
            igrid(ivec) = itmp(jgrid(ivec))
         end do
      end if
      i = i + 1
   end do
!  all grids have converged, or maxit is reached
!  rearrange array (x2) back to the original grids

   do ivec = 1, nptsnz
      jk = igrid(ivec)
      x2tmp(jk) = x2(ivec)
   end do

   do ivec = 1, nptsnz
      x2(ivec) = x2tmp(ivec)
   end do

!  more  array need to be brought back to their original format
!  they are:  b, t, r .    s.m./dec99.

   do ii = 1, 5
      do ivec = 1, nptsnz
         jk = igrid(ivec)
         btmp(jk, ii) = b(ivec, ii)
         ttmp(jk, ii) = t(ivec, ii)
      end do
      do ivec = 1, nptsnz
         b(ivec, ii) = btmp(ivec, ii)
         t(ivec, ii) = ttmp(ivec, ii)
      end do
   end do

   do ii = 1, 25
      do ivec = 1, nptsnz
         jk = igrid(ivec)
         rtmp(jk, ii) = r(ivec, ii)
      end do
      do ivec = 1, nptsnz
         r(ivec, ii) = rtmp(ivec, ii)
      end do
   end do

   do ivec = 1, nptsnz
      c(ivec, 5) = c(ivec, 5) * check1(ivec) + x2(ivec) * (1.0 - check1(ivec))
   end do

   return
end
