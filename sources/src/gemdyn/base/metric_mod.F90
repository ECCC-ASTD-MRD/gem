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

module metric
   implicit none
   public
   save

   real*8, dimension(:,:,:), allocatable :: mc_g1u, mc_g1c, mc_g2v, mc_g2c, mc_g3w, mc_g3c

   public :: set_metric
contains

!**   s/r set_metric - calculate metric coefficients
   subroutine set_metric
      use HORgrid_options
      use gem_options
      use gmm_geof
      use geomh
      use tdpack
      use gmm_itf_mod
      use glb_ld
      use lun
      use cstv
      use ver
      implicit none
#include <arch_specific.hf>

      integer k,istat,i,j
      real*8, parameter :: zero=0.d0, one=1.d0, half=.5d0
     !real temp(l_minx:l_maxx,l_miny:l_maxy,G_nk)
!
!     ---------------------------------------------------------------

      allocate ( mc_g1u(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 mc_g1c(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 mc_g2v(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 mc_g2c(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 mc_g3w(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 mc_g3c(l_minx:l_maxx,l_miny:l_maxy,G_nk) )

      istat = gmm_get (gmmk_fis0_s,fis0)

      mc_g1u=0.d0; mc_g2v=0.d0; mc_g1c=0.d0; mc_g2c=0.d0; mc_g3w=1.d0; mc_g3c=1.d0
      do k=1,G_nk
         do j=1,l_nj
         do i=1,l_ni
            mc_g1u(i,j,k)=Ver_b_8%t(k)*(fis0(i+1,j)-fis0(i-1,j))*half*geomh_invDX_8(j) &
                         *Ver_idz_8%t(k)/(grav_8+Ver_dbdz_8%t(k)*fis0(i,j))
            mc_g1c(i,j,k)=Ver_b_8%t(k)*(fis0(i+1,j)-fis0(i,j))*geomh_invDXMu_8(j) &
                         *Ver_idz_8%t(k)/(grav_8+Ver_dbdz_8%t(k)*half*(fis0(i+1,j)+fis0(i,j)))
            mc_g2v(i,j,k)=Ver_b_8%t(k)*(fis0(i,j+1)-fis0(i,j-1))*half*geomh_invDY_8 &
                         *Ver_idz_8%t(k)/(grav_8+Ver_dbdz_8%t(k)*fis0(i,j))
            mc_g2c(i,j,k)=Ver_b_8%t(k)*(fis0(i,j+1)-fis0(i,j))*geomh_invDY_8 &
                         *Ver_idz_8%t(k)/(grav_8+Ver_dbdz_8%t(k)*half*(fis0(i,j+1)+fis0(i,j)))
            mc_g3w(i,j,k)=one/(one+Ver_dbdz_8%t(k)*fis0(i,j)/grav_8)
            mc_g3c(i,j,k)=one/(one+Ver_dbdz_8%m(k)*fis0(i,j)/grav_8)
         end do
         end do
      end do

     !temp=mc_g1c
     !call glbstat2 ( temp,'G1C',"metric",l_minx,l_maxx,l_miny,l_maxy, &
     !                1,G_nk, 1,G_ni,1,G_nj,1,G_nk )

      if (Grd_yinyang_L) then
         print*,'things to do in set_metric for yinyang'
         stop
      end if
!     ---------------------------------------------------------------
!
      return
   end subroutine set_metric

end module metric
