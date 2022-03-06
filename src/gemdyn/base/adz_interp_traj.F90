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

      subroutine adz_interp_traj (F_dtzD_8, F_dtzA_8, F_dt_8)
      use adz_mem
      use adz_options
      use glb_ld
      use dcst
      use gmm_vt0
      use ver
      use geomh
      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64), intent(IN) :: F_dtzD_8, F_dtzA_8, F_dt_8

      integer :: i,j,k,k00
      real(kind=REAL64) :: xt,yt,zt, wdt,ww,wp,wm
      real(kind=REAL64) :: zmin_bound, zmax_bound
      real(kind=REAL64), parameter :: &
      half= 0.5d0, aa= -0.0625d0, bb= +0.5625d0
!
!---------------------------------------------------------------------
!
      zmin_bound = dble(0)
      zmax_bound = dble(l_nk+1)
      k00=max(Adz_k0t-1,1)
      if ((adz_BC_LAM_flux/=0).and.(Adz_k0>1)) k00=1
      ww=Ver_wmstar_8(l_nk)
      wp=(Ver_z_8%t(l_nk  )-Ver_z_8%m(l_nk-1))*Ver_idz_8%t(l_nk-1)
      wm=1.d0-wp

!$omp do collapse(2)
      do k= Adz_k0, l_nk
         do j= Adz_j0, Adz_jn
            do i= Adz_i0u, Adz_inu
               xt =  aa * (Adz_wpxyz(i-1,j,k,1) + Adz_wpxyz(i+2,j,k,1)) &
                   + bb * (Adz_wpxyz(i  ,j,k,1) + Adz_wpxyz(i+1,j,k,1)) - 0.5d0
               yt =  aa * (Adz_wpxyz(i-1,j,k,2) + Adz_wpxyz(i+2,j,k,2)) &
                   + bb * (Adz_wpxyz(i  ,j,k,2) + Adz_wpxyz(i+1,j,k,2))
               zt =  aa * (Adz_wpxyz(i-1,j,k,3) + Adz_wpxyz(i+2,j,k,3)) &
                   + bb * (Adz_wpxyz(i  ,j,k,3) + Adz_wpxyz(i+1,j,k,3))
               Adz_pmu(1,i,j,k)= min(max(xt,Adz_iminposx),Adz_imaxposx)
               Adz_pmu(2,i,j,k)= min(max(yt,Adz_iminposy),Adz_imaxposy)
               Adz_pmu(3,i,j,k)= min(max(zt,zmin_bound  ),zmax_bound  )
            end do
         end do
      end do
!$omp enddo nowait
!$omp do collapse(2)
      do k= Adz_k0, l_nk
         do j= Adz_j0v, Adz_jnv
            do i= Adz_i0, Adz_in
               xt =  aa * (Adz_wpxyz(i,j-1,k,1) + Adz_wpxyz(i,j+2,k,1)) &
                   + bb * (Adz_wpxyz(i,j  ,k,1) + Adz_wpxyz(i,j+1,k,1))
               yt =  aa * (Adz_wpxyz(i,j-1,k,2) + Adz_wpxyz(i,j+2,k,2)) &
                   + bb * (Adz_wpxyz(i,j  ,k,2) + Adz_wpxyz(i,j+1,k,2)) - 0.5d0
               zt =  aa * (Adz_wpxyz(i,j-1,k,3) + Adz_wpxyz(i,j+2,k,3)) &
                   + bb * (Adz_wpxyz(i,j  ,k,3) + Adz_wpxyz(i,j+1,k,3))
               Adz_pmv(1,i,j,k)= min(max(xt,Adz_iminposx),Adz_imaxposx)
               Adz_pmv(2,i,j,k)= min(max(yt,Adz_iminposy),Adz_imaxposy)
               Adz_pmv(3,i,j,k)= min(max(zt,zmin_bound  ),zmax_bound  )
           end do
         end do
      end do
!$omp enddo nowait

!$omp do collapse(2)
      do k= max(k00,2), l_nk-2
         do j= Adz_j0b, Adz_jnb
            do i= Adz_i0b, Adz_inb
               xt = adz_vw1t(k)*Adz_wpxyz(i,j,k-1,1)+ &
                    adz_vw2t(k)*Adz_wpxyz(i,j,k  ,1)+ &
                    adz_vw3t(k)*Adz_wpxyz(i,j,k+1,1)+ &
                    adz_vw4t(k)*Adz_wpxyz(i,j,k+2,1)
               yt = adz_vw1t(k)*Adz_wpxyz(i,j,k-1,2)+ &
                    adz_vw2t(k)*Adz_wpxyz(i,j,k  ,2)+ &
                    adz_vw3t(k)*Adz_wpxyz(i,j,k+1,2)+ &
                    adz_vw4t(k)*Adz_wpxyz(i,j,k+2,2)
               wdt= adz_vw1t(k)*Adz_uvw_dep(3,i,j,k-1)+ &
                    adz_vw2t(k)*Adz_uvw_dep(3,i,j,k  )+ &
                    adz_vw3t(k)*Adz_uvw_dep(3,i,j,k+1)+ &
                    adz_vw4t(k)*Adz_uvw_dep(3,i,j,k+2)
               zt = Ver_z_8%t(k) - F_dtzD_8*  wdt &
                                 - F_dtzA_8* zdt0(i,j,k)
               Adz_pb(1,i,j,k)= min(max(xt,Adz_iminposx),Adz_imaxposx)
               Adz_pb(2,i,j,k)= min(max(yt,Adz_iminposy),Adz_imaxposy)
               Adz_pb(3,i,j,k)= zindx(zt)
            end do
         end do
      end do
!$omp enddo nowait

      k= l_nk-1
!$omp do
      do j= Adz_j0b, Adz_jnb
         do i= Adz_i0b, Adz_inb
            xt = (Adz_wpxyz(i,j,k,1)+Adz_wpxyz(i,j,k+1,1))*half
            yt = (Adz_wpxyz(i,j,k,2)+Adz_wpxyz(i,j,k+1,2))*half
            wdt= (Adz_uvw_dep(3,i,j,k)+Adz_uvw_dep(3,i,j,k+1))*half
            zt =Ver_z_8%t(k) - F_dtzD_8*  wdt &
                             - F_dtzA_8*zdt0(i,j,k)
            Adz_pb(1,i,j,k)= min(max(xt,Adz_iminposx),Adz_imaxposx)
            Adz_pb(2,i,j,k)= min(max(yt,Adz_iminposy),Adz_imaxposy)
            Adz_pb(3,i,j,k)= zindx(zt)
         end do
      end do
!$omp enddo nowait

      k= l_nk
      if (Adz_slt_winds) then
!$omp do
         do j= Adz_j0b, Adz_jnb
            do i= Adz_i0b, Adz_inb
               xt= dble(i+l_i0-1) - Dcst_inv_rayt_8 * F_dt_8 * &
                      Adz_uslt(i,j)*Adz_cy_8(j) * geomh_inv_hx_8
               yt= dble(j+l_j0-1) - Dcst_inv_rayt_8 * F_dt_8 * &
                      Adz_vslt(i,j)             * geomh_inv_hy_8
               Adz_pb(1,i,j,k)= min(max(xt,Adz_iminposx),Adz_imaxposx)
               Adz_pb(2,i,j,k)= min(max(yt,Adz_iminposy),Adz_imaxposy)
            end do
         end do
!$omp enddo nowait
      else
!$omp do collapse(2)
         do j= Adz_j0b, Adz_jnb
            do i= Adz_i0b, Adz_inb
               xt = wp*Adz_wpxyz(i,j,k,1)+wm*Adz_wpxyz(i,j,k-1,1)
               yt = wp*Adz_wpxyz(i,j,k,2)+wm*Adz_wpxyz(i,j,k-1,2)
               Adz_pb(1,i,j,k)= min(max(xt,Adz_iminposx),Adz_imaxposx)
               Adz_pb(2,i,j,k)= min(max(yt,Adz_iminposy),Adz_imaxposy)
            end do
         end do
!$omp enddo nowait
      end if
!$omp do
      do j= Adz_j0b, Adz_jnb
         do i= Adz_i0b, Adz_inb
            wdt= (Adz_uvw_dep(3,i,j,k-1)+Adz_uvw_dep(3,i,j,k))*half
            zt = Ver_z_8%t(k) - ww*(F_dtzD_8*  wdt &
                                  + F_dtzA_8*zdt0(i,j,k-1))
            Adz_pb(3,i,j,k)= zindx(zt)
         end do
      end do
!$omp enddo
      
      if (k00==1) then
!$omp do
         do j= Adz_j0b, Adz_jnb
            do i= Adz_i0b, Adz_inb
               xt = (Adz_wpxyz(i,j,1,1 )+Adz_wpxyz (i,j,2,1))*half
               yt = (Adz_wpxyz(i,j,1,2 )+Adz_wpxyz (i,j,2,2))*half
               wdt= (Adz_uvw_dep(3,i,j,1)+Adz_uvw_dep(3,i,j,2))*half
               zt = Ver_z_8%t(1) - F_dtzD_8*wdt &
                                 - F_dtzA_8*zdt0(i,j,1)
               Adz_pb(1,i,j,1)= min(max(xt,Adz_iminposx),Adz_imaxposx)
               Adz_pb(2,i,j,1)= min(max(yt,Adz_iminposy),Adz_imaxposy)
               Adz_pb(3,i,j,1)= zindx(zt)
            end do
         end do
!$omp enddo
      endif
      
!$omp do
      do k= Adz_k0t, l_nk
         Adz_pt  (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)=&
         Adz_pb  (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)
      end do
!$omp enddo
!
!---------------------------------------------------------------------
!
      return

contains

      real(kind=REAL64) function zindx (posz)
      use, intrinsic :: iso_fortran_env
      implicit none
      real(kind=REAL64), intent(inout) :: posz
      integer :: kk1,nb
      posz = min(max(posz,Ver_zmin_8),Ver_zmax_8)
      kk1 = (posz - ver_z_8%t(0)  ) * adz_ovdzt_8 + 1.d0
      kk1 = min(max(1,kk1),ubound(adz_search_t,1))
      kk1 = adz_search_t(kk1)
      if ( sig*posz >  sig*ver_z_8%t(min(kk1+1,l_nk+1))) kk1= kk1 + 1
      if ( sig*posz <  sig*ver_z_8%t(kk1)              ) kk1= kk1 - 1
      kk1 = min(l_nk+1,max(0,kk1))
      nb  = max(min(kk1,G_nk-1),1)
      zindx= (posz-ver_z_8%t(nb))*Adz_odelz_t(nb) + dble(nb)
      return
      end function zindx

      end subroutine adz_interp_traj

