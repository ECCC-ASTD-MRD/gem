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
!
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
! Fichier/File   : mach_aurams_cldcv2d.ftn90
! Creation       : A. Dastoor, W. Gong, S. Gravel, GEM-MACH, July 2008
! Description    : Compute atmospheric pressure
!
! Extra info     :
!
! Arguments  IN
!                 fnuage -->  3d cloud fraction
!            OUT
!                 cc2d   -->  2d cloud cover above that level
!
!=============================================================================
!
!!if_on
subroutine mach_aurams_cldcv2d (cc2d, fnuage, pni, pnk)
!!if_off
   implicit none
!!if_on
   integer(kind=4), intent (in) :: pni, pnk
   real(kind=4),    intent (in) :: fnuage(pni, pnk)
   real(kind=4),    intent(out) :: cc2d  (pni, pnk)
!!if_off
!
! local variables
!
   integer(kind=4)  :: i, l, jind, j
   real(kind=4)     :: xnu
   real(kind=4)     :: f(pni, pnk), tmem(pni), trmin(pni)

! note: vertical levels are from model top to surface.
   do i = 1, pnk
      do l = 1, pni
         f(l, i) = 1.0
         cc2d(l, i) = fnuage(l, i)
         tmem(l) = 1.0
         trmin(l) = 1.0
      end do
   end do

   do j = 2, pnk
      jind = j - 2
      jind = max(jind, 1)
      do l = 1, pni

! transmission through cloud layer xnu
         xnu = 1.0 - fnuage(l, j - 1)
         if (fnuage(l, jind) < 0.01) then
! entering a new cloud isolated from upper one
! keep in memory transmission down to top of new cloud tmem
! trmin is minimum transmission in cloud layer processed
! basic idea is random overlap (hence tmem x xnu = cloud transmittance)
! for cloud layers separated by clear ones; but maximum overlap
! (hence cloud transmittance is min (tmem, xnu))for adjacent cloud layers.
               tmem(l) = f(l, j - 1)
               trmin(l) = xnu
            else
! inside a cloud use maximum overlap
! compute minimum transmission between adjacent layers
               trmin(l) = min(trmin(l), xnu)
            end if

            f(l, j) = tmem(l) * trmin(l)
            cc2d(l, j) = 1.0 - f(l, j)

      end do
   end do
   return
end subroutine mach_aurams_cldcv2d
