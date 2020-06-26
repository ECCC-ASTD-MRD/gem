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

!**s/r pw_update_GW - Update physical quantities WZ and GZ

      subroutine pw_update_GW()
      use dynkernel_options
      use gem_timing
      use glb_ld
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use metric
      use tdpack
      implicit none
#include <arch_specific.hf>

      integer :: i,j,k
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk+1) :: fi
!     ________________________________________________________________

      if (Schm_autobar_L) return

      call gemtime_start ( 5, 'PW_UPDATE', 0)

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H') then
         pw_wz_plus = 0. ; pw_gz_plus = 0.
         pw_me_plus = 0.
         return ! Not yet implemented
      end if

      do k=1,l_nk
         do j=1,l_nj
            do i=1,l_ni
               pw_wz_plus(i,j,k) = wt1(i,j,k)
            end do
         end do
      end do
      do j=1,l_nj
         do i=1,l_ni
            pw_me_plus(i,j)= fis0(i,j)
         end do
      end do

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

         do k=1,l_nk
            do j=1,l_nj
            do i=1,l_ni
               pw_gz_plus(i,j,k)= grav_8*zmom_8(i,j,k)
            end do
            end do
         end do

      else

         call diag_fi (fi, st1, tt1, qt1, l_minx,l_maxx,l_miny,l_maxy,&
                       G_nk, 1, l_ni, 1, l_nj)

         do k=1,l_nk
            do j=1,l_nj
               do i=1,l_ni
                  pw_gz_plus(i,j,k)= fi(i,j,k)
               end do
            end do
         end do

      end if

      call gemtime_stop (5)
!     ________________________________________________________________
!
      return
      end subroutine pw_update_GW
