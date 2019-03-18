!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

module coriolis
   implicit none
   public
   save

   !  CORIOLIS FACTORS FOR U and V grids                                  |
   !______________________________________________________________________|
   !                    |                                                 |
   ! NAME               | DESCRIPTION                                     |
   !--------------------|-------------------------------------------------|
   ! Cori_fcoru_8       | coriolis factor computed on the U grid          |
   ! Cori_fcorv_8       | coriolis factor computed on the V grid          |
   !----------------------------------------------------------------------

   real*8, dimension(:,:), allocatable :: cori_fcoru_8, cori_fcorv_8

   public :: set_coriolis

contains

   subroutine set_coriolis(F_x_8, F_y_8, F_xu_8, F_yv_8, F_rot_8, Minx,Maxx,Miny,Maxy)

      ! coriolis - compute coriolis factor
      use dcst
      use ctrl
      use gem_options
      use HORgrid_options
      use glb_ld
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx,Maxx,Miny,Maxy
      real*8, intent(in) :: F_x_8 (Minx:Maxx), & ! longitudes in radians PHI grid
                            F_y_8 (Miny:Maxy), & ! latitudes in radians PHI grid
                            F_xu_8(Minx:Maxx), & ! longitudes in radians U grid
                            F_yv_8(Miny:Maxy), & ! latitudes in radians V grid
                            F_rot_8(3,3)         ! rotation matrix of the grid

      integer :: i, j
      real*8  :: s0, ang, c0, sa, ca
      !
      !-------------------------------------------------------------------
      !
      allocate( cori_fcoru_8(l_minx:l_maxx,l_miny:l_maxy),&
                cori_fcorv_8(l_minx:l_maxx,l_miny:l_maxy) )

      if ( Ctrl_theoc_L ) then
         cori_fcoru_8 = 0.d0
         cori_fcorv_8 = 0.d0
         return
      end if

      if ( Ctrl_canonical_williamson_L ) then
         call canonical_coriolis ( cori_fcoru_8, cori_fcorv_8, F_x_8, F_y_8,  &
                                   F_xu_8, F_yv_8, F_rot_8, Minx, Maxx, Miny, Maxy)
         return
      end if

      s0 = F_rot_8(3,3)

      if ( abs( (abs(s0)-1.d0) ) > 1.0e-10 ) then
         ang = atan2( F_rot_8(2,3), F_rot_8(1,3) )
      else
         s0 = sign( 1.d0, s0 )
         ang = 0.d0
      end if

      c0 = sqrt( max( 0.d0, 1.d0 - s0**2 ) )

      !	processing coriolis FACTOR on V grid

      do j=1-G_haloy,l_nj+G_haloy
         sa = ( 2.d0 * Dcst_omega_8 ) * s0 * sin(F_yv_8(j))
         ca = ( 2.d0 * Dcst_omega_8 ) * c0 * cos(F_yv_8(j))
         do i=1-G_halox,l_ni+G_halox
            cori_fcorv_8(i,j) = ca * cos(F_x_8(i)-ang) + sa
         end do
      end do

      !	processing coriolis FACTOR on U grid

      do j=1-G_haloy,l_nj+G_haloy
         sa = ( 2.d0 * Dcst_omega_8 ) * s0 * sin(F_y_8(j))
         ca = ( 2.d0 * Dcst_omega_8 ) * c0 * cos(F_y_8(j))
         do i=1-G_halox,l_ni+G_halox
            cori_fcoru_8(i,j) = ca * cos(F_xu_8(i) - ang) + sa
         end do
      end do

      return
   end subroutine set_coriolis

end module coriolis
