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

!**s/r pw_update_T - Update physical temperature from virtual temperature tt1

      subroutine pw_update_T_hlt ()
      use glb_ld
      use gem_options
      use gmm_pw
      use gmm_vt1
      use tr3d
      use mem_tracers
      use mem_tstp
      use tdpack
      use omp_timing
      implicit none

      integer i,j,k
      real, dimension(:,:,:), pointer :: wk
!
!     ________________________________________________________________
!
      call gtmg_start (5, 'PW_UPDATE', 0)

      wk(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(1:)
      call sumhydro_hlt (wk, l_minx,l_maxx,l_miny,l_maxy, l_nk,Tr3d_ntr, trt1)

!$omp do collapse(2)
      do k= 1,l_nk
         do j=1-g_haloy, l_nj+g_haloy
            do i=1-g_halox, l_ni+g_halox
               pw_tt_plus(i,j,k) = fottvh(tt1(i,j,k),tracers_p(tr3d_hu)%pntr(i,j,k),wk(i,j,k))
            end do
         end do
      end do
!$omp end do

      call gtmg_stop (5)
!
!     ________________________________________________________________
!
      return
      end subroutine pw_update_T_hlt
