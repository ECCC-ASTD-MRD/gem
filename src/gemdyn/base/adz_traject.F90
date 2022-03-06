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

      subroutine adz_traject (F_dtA_8, F_dtzA_8, F_dtD_8, F_dtzD_8)
      use glb_ld
      use geomh
      use ver
      use adz_options
      use adz_mem
      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dtA_8,F_dtzA_8,F_dtD_8,F_dtzD_8

      include "tricublin_f90.inc"
      integer :: iter,i,j,k,kk1,nb,k00
      integer,dimension(l_ni) :: kk
      real(kind=REAL64), dimension(l_ni) :: xm,ym,zm
      real(kind=REAL64) pos
!
!     ---------------------------------------------------------------
!
      call adz_prepareWinds ()

      k00=Adz_k0m
      if (Adz_k0>1) k00=1

      do  iter = 1, Adz_niter
!$omp do
         do k= k00, l_nk
             call tricublin_zyx3_n ( Adz_uvw_dep(1,1,1,k),Adz_uvw_d(1,1,1,1), &
                                     Adz_pxyzm(1:3,:,:,k),Adz_cpntr_q,Adz_2dnh)
             do j= 1, l_nj
                do i= 1, l_ni
                  xm(i) = dble(i+l_i0-1) - &
                    ( F_dtD_8 * Adz_uvw_dep(1,i,j,k)  &
                    + F_dtA_8 * Adz_uu_arr(i,j,k) )&
                    * geomh_inv_hx_8
                  ym(i) = dble(j+l_j0-1) - &
                    ( F_dtD_8 * Adz_uvw_dep(2,i,j,k)  &
                    + F_dtA_8 * Adz_vv_arr(i,j,k) )&
                    * geomh_inv_hy_8
                  pos = Ver_z_8%m(k)- F_dtzD_8* Adz_uvw_dep(3,i,j,k) &
                                    - F_dtzA_8* Adz_ww_arr(i,j,k)
                  zm(i) = min(max(pos,Ver_zmin_8),Ver_zmax_8)

                  kk1 = (zm(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                  kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                  kk1 = adz_search_m(kk1)
                  if ( sig * zm(i) > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                  if ( sig * zm(i) < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
                  kk(i) = kk1
               end do
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
         enddo
!$omp enddo
      end do
      
!!$OMP BARRIER

!$omp do
      do k= Adz_k0, l_nk
         Adz_pm   (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)=&
         Adz_pxyzm(:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)
      end do
!$omp enddo nowait
      
!$omp single
      Adz_niter = Adz_itraj
!$omp end single

!     ---------------------------------------------------------------
!
      return
      end subroutine adz_traject
