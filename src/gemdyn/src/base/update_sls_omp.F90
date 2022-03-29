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

!**s/r get_topo - update sls with F_src

      subroutine update_sls_omp (F_src,F_dest,minx,maxx,miny,maxy)
      use dynkernel_options
      use gem_options
      use glb_ld
      use gmm_geof
      use tdpack
      implicit none
      
      integer, intent(in) :: Minx, Maxx, Miny, Maxy
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in ) :: F_src
      real, dimension(Minx:Maxx,Miny:Maxy), intent(out) :: F_dest

      integer :: i,j
      real(kind=REAL64) :: oneoRT
!
!---------------------------------------------------------------------
!
      oneoRT=1.d0 / (rgasd_8 * Tcdk_8)
      if ( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P' ) then
!$omp do
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               F_dest(i,j) = -F_src(i,j) * oneoRT
            end do
         end do
!$omp enddo
      else
!$omp do
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               F_dest(i,j) = F_src(i,j)
            end do
         end do
!$omp enddo
      end if
!     
!---------------------------------------------------------------------
!
      return
      end subroutine update_sls_omp
