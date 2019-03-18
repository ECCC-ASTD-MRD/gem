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

      subroutine adz_traject (uut1, vvt1, zzt1, &
                              uut0, vvt0, zzt0, &
                          Minx,Maxx,Miny,Maxy,NK)
      use glb_ld
      use cstv
      use geomh
      use ver
      use adz_options
      use adz_mem
      implicit none
#include <arch_specific.hf>

      integer,intent(in) :: Minx,Maxx,Miny,Maxy,NK
      real, dimension(Minx:Maxx,Miny:Maxy,NK),intent(in)  :: &
                            uut1, vvt1, zzt1, uut0, vvt0, zzt0

      integer :: iter,i,j,k,kk1,nnk
      integer,dimension(l_ni) :: kk
      real*8, dimension(l_ni) :: xm,ym,zm
      real*8 pos
!
!     ---------------------------------------------------------------
!
      call adz_prepareWinds (uut1, vvt1, zzt1, &
                             uut0, vvt0, zzt0, &
                         Minx,Maxx,Miny,Maxy,NK)

      nnk = Adz_2dnh *(l_nk-Adz_k0m+1)

      do  iter = 1, Adz_niter

         call adz_tricub_wind ( Adz_uu_dep(1,1,Adz_k0m),&
                                Adz_vv_dep(1,1,Adz_k0m),&
                                Adz_ww_dep(1,1,Adz_k0m),&
                     Adz_uvw_d,Adz_pxyzm(1,1,1,Adz_k0m),nnk)

!- Compute departure positions

          do k= Adz_k0m, l_nk
             do j= 1, l_nj
                do i= 1, l_ni
                  xm(i) = dble(i+l_i0-1) - &
                    ( Cstv_dtD_8 * Adz_uu_dep(i,j,k)  &
                    + Cstv_dtA_8 * Adz_uu_arr(i,j,k) )&
                    * geomh_inv_hx_8
                  ym(i) = dble(j+l_j0-1) - &
                    ( Cstv_dtD_8 * Adz_vv_dep(i,j,k)  &
                    + Cstv_dtA_8 * Adz_vv_arr(i,j,k) )&
                    * geomh_inv_hy_8
                  pos = Ver_z_8%m(k) - Cstv_dtzD_8* Adz_ww_dep(i,j,k) &
                                     - Cstv_dtzA_8* Adz_ww_arr(i,j,k)
                  zm(i) = min(max(pos,Ver_zmin_8),Ver_zmax_8)

                  kk1 = (zm(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                  kk1 = adz_search_m(kk1)
                  if ( zm(i) > ver_z_8%m(min(kk1+1,l_nk+1))) kk1= kk1 + 1
                  if ( zm(i) < ver_z_8%m(kk1             ) ) kk1= kk1 - 1
                  kk(i) = kk1
               end do
!DIR$ IVDEP
               do i= 1, l_ni
                  Adz_wpxyz(i,j,k,1) = xm(i)
                  Adz_pxyzm(1,i,j,k) = min(max(xm(i),Adz_iminposx),&
                                                     Adz_imaxposx)
                  Adz_wpxyz(i,j,k,2) = ym(i)
                  Adz_pxyzm(2,i,j,k) = min(max(ym(i),Adz_iminposy),&
                                                     Adz_imaxposy)
                  kk1 = min(l_nk+1,max(0,kk(i)))
                  Adz_wpxyz(i,j,k,3) = (zm(i)-ver_z_8%m(kk1))&
                                     *Adz_odelz_m(kk1) + dble(kk1)
                  Adz_pxyzm(3,i,j,k) = Adz_wpxyz(i,j,k,3)
                end do
            end do
         end do
      end do
! Interpolate Momentun positions on Thermo Level and U V  grid

      call adz_interp_traj (Adz_k0,Adz_k0t)

      Adz_niter = Adz_itraj
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_traject
