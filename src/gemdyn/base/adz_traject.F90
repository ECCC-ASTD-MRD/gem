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
      use gem_timing
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer,intent(in) :: Minx,Maxx,Miny,Maxy,NK
      real, dimension(Minx:Maxx,Miny:Maxy,NK),intent(in)  :: &
                            uut1, vvt1, zzt1, uut0, vvt0, zzt0

      include "tricublin_f90.inc"
      integer :: iter,i,j,k,kk1,nnk,nb
      integer,dimension(l_ni) :: kk
      real(kind=REAL64), dimension(l_ni) :: xm,ym,zm
      real(kind=REAL64) pos
!
!     ---------------------------------------------------------------
!
      call time_trace_barr(gem_time_trace, 1, Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)
      call adz_prepareWinds (uut1, vvt1, zzt1, &
                             uut0, vvt0, zzt0, &
                         Minx,Maxx,Miny,Maxy,NK)

      nnk = Adz_2dnh *(l_nk-Adz_k0m+1)

      do  iter = 1, Adz_niter

         call adz_ondemand (Adz_expq,Adz_i0,Adz_in,Adz_j0,Adz_jn    ,&
               Adz_k0,l_nk,Adz_wpxyz(-1,-1,1,1),Adz_wpxyz(-1,-1,1,2),&
               Adz_wpxyz(-1,-1,1,3),l_ni,l_nj,l_nk)

         call tricublin_zyx3_n ( Adz_uvw_dep,Adz_uvw_d(1,1,1,1), &
                                 Adz_pxyzm,Adz_cpntr_q,Adz_3dnh )

         Adz_cnt_traj = Adz_cnt_traj+1
         call time_trace_barr(gem_time_trace, 10200+Adz_cnt_traj,&
                     Gem_trace_barr, Ptopo_intracomm, MPI_BARRIER)
!- Compute departure positions

          do k= Adz_k0m, l_nk
             do j= 1, l_nj
                do i= 1, l_ni
                  xm(i) = dble(i+l_i0-1) - &
                    ( Cstv_dtD_8 * Adz_uvw_dep(1,i,j,k)  &
                    + Cstv_dtA_8 * Adz_uu_arr(i,j,k) )&
                    * geomh_inv_hx_8
                  ym(i) = dble(j+l_j0-1) - &
                    ( Cstv_dtD_8 * Adz_uvw_dep(2,i,j,k)  &
                    + Cstv_dtA_8 * Adz_vv_arr(i,j,k) )&
                    * geomh_inv_hy_8
                  pos = Ver_z_8%m(k)- Cstv_dtzD_8* Adz_uvw_dep(3,i,j,k) &
                                    - Cstv_dtzA_8* Adz_ww_arr(i,j,k)
                  zm(i) = min(max(pos,Ver_zmin_8),Ver_zmax_8)

                  kk1 = (zm(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                  kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                  kk1 = adz_search_m(kk1)
                  if ( sig * zm(i) > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                  if ( sig * zm(i) < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
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
                  nb  = max(min(kk1,G_nk-1),1)
                  Adz_wpxyz(i,j,k,3) = (zm(i)-ver_z_8%m(nb))&
                                     *Adz_odelz_m(nb) + dble(nb)
                  Adz_pxyzm(3,i,j,k) = Adz_wpxyz(i,j,k,3)
                end do
            end do
         end do
         call time_trace_barr(gem_time_trace, 10300+Adz_cnt_traj,&
                     Gem_trace_barr, Ptopo_intracomm, MPI_BARRIER)

      end do

      Adz_pm   (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, Adz_k0:l_nk)=&
      Adz_pxyzm(:,Adz_i0:Adz_in, Adz_j0:Adz_jn, Adz_k0:l_nk)

      call adz_ondemand (Adz_expq,Adz_i0,Adz_in,Adz_j0,Adz_jn       ,&
               Adz_k0,l_nk,Adz_wpxyz(-1,-1,1,1),Adz_wpxyz(-1,-1,1,2),&
               Adz_wpxyz(-1,-1,1,3),l_ni,l_nj,l_nk)

! Interpolate Momentun positions on Thermo Level and U V  grid

      call adz_interp_traj ()
      call time_trace_barr(gem_time_trace,10400+Adz_icn,Gem_trace_barr,&
                           Ptopo_intracomm, MPI_BARRIER)

      Adz_niter = Adz_itraj
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_traject
