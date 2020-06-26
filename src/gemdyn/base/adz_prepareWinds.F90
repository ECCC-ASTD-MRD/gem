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

      subroutine adz_prepareWinds (uut1, vvt1, zzt1, &
                                   uut0, vvt0, zzt0, &
                               Minx,Maxx,Miny,Maxy,NK)
      use adz_mem
      use dcst
      use gem_options
      use glb_ld
      use ver
      use gem_timing
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer,intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,NK),intent(in)  :: &
                            uut1, vvt1, zzt1, uut0, vvt0, zzt0
      integer :: i,j,k,km2,nrow
      real(kind=REAL64), parameter :: alpha1=-1.d0/16.d0, alpha2=9.d0/16.d0
      real(kind=REAL64) :: xx, x1, x2, x3, x4, w1, w2, w3, w4

!
!     ---------------------------------------------------------------
!
      nrow = 0
      call rpn_comm_xch_halox( uut1, l_minx,l_maxx,l_miny,l_maxy      ,&
                               l_ni, l_nj, l_nk, Adz_halox, Adz_haloy ,&
                               .false.,.false., Adz_uu_ext, Adz_lminx ,&
                               Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, nrow)
      call rpn_comm_xch_halox( vvt1, l_minx,l_maxx,l_miny,l_maxy      ,&
                               l_ni, l_nj, l_nk, Adz_halox, Adz_haloy ,&
                               .false.,.false., Adz_vv_ext, Adz_lminx ,&
                               Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, nrow)
      call rpn_comm_xch_halox( zzt1, l_minx,l_maxx,l_miny,l_maxy      ,&
                               l_ni, l_nj, l_nk, Adz_halox, Adz_haloy ,&
                               .false.,.false., Adz_ww_ext, Adz_lminx ,&
                               Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, nrow)
      call rpn_comm_xch_halo( uut0,l_minx,l_maxx,l_miny,l_maxy,&
                               l_niu,l_nj,l_nk,G_halox,G_haloy,&
                               G_periodx,G_periody,G_niu,0)
      call rpn_comm_xch_halo( vvt0,l_minx,l_maxx,l_miny,l_maxy,&
                               l_ni,l_njv,l_nk,G_halox,G_haloy,&
                               G_periodx,G_periody,G_ni,0)

      call time_trace_barr(gem_time_trace,10000+Adz_icn,Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

      do k =1, l_nk
         do j =1-Adz_haloy, l_nj+Adz_haloy
            do i =1-Adz_halox, l_ni+Adz_halox
               Adz_uvw_d(1,i,j,k)= Dcst_inv_rayt_8 * Adz_uu_ext(i,j,k) * Adz_cy_8(j)
               Adz_uvw_d(2,i,j,k)= Dcst_inv_rayt_8 * Adz_vv_ext(i,j,k)
            end do
         end do
      end do

      do k =1, l_nk
         do j =1, l_nj
            do i =1, l_ni
               Adz_uu_arr(i,j,k) = ((uut0(i-2,j,k) + uut0(i+1,j,k))*alpha1 &
                                   +(uut0(i  ,j,k) + uut0(i-1,j,k))*alpha2)&
                                   * Dcst_inv_rayt_8 * Adz_cy_8(j)
               Adz_vv_arr(i,j,k) = ((vvt0(i,j-2,k) + vvt0(i,j+1,k))*alpha1 &
                                   +(vvt0(i  ,j,k) + vvt0(i,j-1,k))*alpha2)&
                                   * Dcst_inv_rayt_8
            end do
         end do
      end do

      do k=2, l_nk-1
         xx = Ver_z_8%m(k)
         x1 = Ver_z_8%x(k-2)
         x2 = Ver_z_8%x(k-1)
         x3 = Ver_z_8%x(k)
         x4 = Ver_z_8%x(k+1)
         w1 = lag3(xx, x1, x2, x3, x4)
         w2 = lag3(xx, x2, x1, x3, x4)
         w3 = lag3(xx, x3, x1, x2, x4)
         w4 = lag3(xx, x4, x1, x2, x3)
         km2=max(1,k-2)
         if (k == 2) w1=0.d0
         do j = 1,l_nj
            do i = 1,l_ni
               Adz_ww_arr(i,j,k) = w1 * zzt0(i,j,km2) + w2 * zzt0(i,j,k-1) &
                                 + w3 * zzt0(i,j,k  ) + w4 * zzt0(i,j,k+1)
            end do
         end do
         do j =1-Adz_haloy, l_nj+Adz_haloy
            do i =1-Adz_halox, l_ni+Adz_halox
               Adz_uvw_d(3,i,j,k)= w1 * Adz_ww_ext(i,j,km2) + &
                                   w2 * Adz_ww_ext(i,j,k-1) &
                                 + w3 * Adz_ww_ext(i,j,k  ) + &
                                   w4 * Adz_ww_ext(i,j,k+1)
            end do
         end do
      end do

      w1 = (Ver_z_8%x(0)-Ver_z_8%m(1)) / (Ver_z_8%x(0)-Ver_z_8%x(1))
      w2 = (Ver_z_8%m(l_nk  )-Ver_z_8%x(l_nk)) / &
          (Ver_z_8%x(l_nk-1)-Ver_z_8%x(l_nk))
      do j =1-Adz_haloy, l_nj+Adz_haloy
         do i =1-Adz_halox, l_ni+Adz_halox
            Adz_uvw_d(3,i,j,   1) = w1 * Adz_ww_ext(i,j,     1)
            Adz_uvw_d(3,i,j,l_nk) = w2 * Adz_ww_ext(i,j,l_nk-1)
         end do
      end do
      do j= 1, l_nj
         do i= 1, l_ni
            Adz_ww_arr(i,j,   1)= w1 * zzt0(i,j,     1)
            Adz_ww_arr(i,j,l_nk)= w2 * zzt0(i,j,l_nk-1)
         end do
      end do

      call adz_perturbWinds(alpha1, alpha2)
      
      call time_trace_barr(gem_time_trace,10100+Adz_icn,Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)
!
!     ---------------------------------------------------------------
!
      return

      contains

         real(kind=REAL64) function lag3(zz, z1, z2, z3, z4)
            implicit none
            real(kind=REAL64), intent(in)  :: zz, z1, z2, z3, z4

            lag3 = ( ((zz - z2) * (zz - z3) * (zz - z4) ) / &
                     ((z1 - z2) * (z1 - z3) * (z1 - z4) ) )
         end function lag3

      end subroutine adz_prepareWinds
