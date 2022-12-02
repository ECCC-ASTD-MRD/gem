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

      subroutine pw_update_GW_hlt ()
      use dynkernel_options
      use gem_options
      use omp_timing
      use omp_timing
      use glb_ld
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use mem_tstp
      use metric
      use tdpack
      implicit none

      integer i,j,k
!
!-------------------------------------------------------------------
!
      if (Schm_autobar_L) return

      call gtmg_start (5, 'PW_UPDATE', 0)

!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy, l_nj+G_haloy
            do i=1-G_halox, l_ni+G_halox
               pw_wz_plus(i,j,k)= wt1(i,j,k)
               pw_gz_plus(i,j,k)= GVM%zmom_8(i,j,k)*grav_8
            end do
         end do
      end do
!$omp end do
!$omp do
      do j=1-G_haloy, l_nj+G_haloy
         do i=1-G_halox, l_ni+G_halox
            pw_me_plus(i,j)= fis0(i,j)
         end do
      end do
!$omp end do

      call gtmg_stop (5)
!     ________________________________________________________________
!
      return
      end subroutine pw_update_GW_hlt
