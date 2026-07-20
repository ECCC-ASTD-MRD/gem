! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numeriquefdg1
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

!**s/r - 3D_diffusion operator   computation for GEM_H
!
      subroutine  hzd_scal_alh ( F_Sol1,HzdlnR, Minx, Maxx, Miny, Maxy,Nk)
      use gem_options
      use gmm_vt1
      use geomh
      use glb_ld
      use cstv
      use ver
      use metric
      use hzd_mod
      use dcst
      use hvdif_options
      use step_options
      use gmm_geof
      use tdpack
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, NK
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent (inout) :: F_Sol1
      real  HzdlnR
      logical, save :: done=.false.

!author
!       Abdessamad Qaddouri -  2018
!
!revision
! v5.0 - Qaddouri A.       - initial version

      integer j,i,k
      integer kd0 , k00, k01
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)

      real(kind=REAL64)   C1, C2 
      real(kind=REAL64) Jz,qkm,qkp
      real(kind=REAL64) C1_8,C2_8,C,ski,skpi,skip,skpip
      real(kind=REAL64) dcoef


! kd0 given by user      
      kd0=hzd_hyb_top
      k00=kd0
      k01=kd0+1
      if (kd0.ne.1) k01=kd0

      if (Hzd_pwr_z==2) then 
         dcoef = 0.25*HzdlnR*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      else
         dcoef = 0.25*sqrt(HzdlnR)*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      endif

! Apply Horizontal diffusion along z

       Afdg1 = 0.0d0
       Bfdg1 = 0.0d0
       add_v8 =0.0d0
       bdd_v8=0.0d0
       cdd_v81=0.d0
       cdd_v82=0.d0
       fdg2_4 =0.0
       cflux=0.0
      do k = 1, nk+1
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w-1, l_ni-pil_e+1
               Afdg1  (i,j,k) = zero
               Bfdg1  (i,j,k) = zero
               Afdg2  (i,j,k) = zero
               Bfdg2  (i,j,k) = zero
               add_v8 (i,j,k) = zero
               bdd_v8 (i,j,k) = zero
               cdd_v8 (i,j,k) = zero
               cdd_v81(i,j,k) = zero
               cdd_v82(i,j,k) = zero
               cflux  (i,j,k) = 0.0
               fdg2_4 (i,j,k) = 0.0 
            enddo
         enddo
      enddo

!Field  before diffusion on T-level K  on phii,j
!$omp do collapse(2)
         do k = 1, nk
            do j=1+pil_s-1, l_nj-pil_n+1
               do i=1+pil_w-1, l_ni-pil_e+1
                  fdg2_4(i,j,k )=F_Sol1(i,j,k)
               enddo
            enddo
         enddo
!$omp enddo

!$omp single
         call rpn_comm_xch_halo(fdg2_4,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,Nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
!$omp end single

         k=k00
!$omp do 
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
               Afdg1(i,j,k) = dcoef*((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j))
               Afdg2(i,j,k) = ((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j))
            enddo
         enddo
!$omp enddo

         k= NK
!$omp do 
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
               C1_8=half*(fdg2_4(i,j,k)*Jxt(i,j,k,1)-fdg2_4(i,j,k-1)*Jxt(i,j,k,2) + &
                         fdg2_4(i+1,j,k )*Jxt(i,j,k,3)-fdg2_4(i+1,j,k-1)*Jxt(i,j,k,4))
               Afdg1(i,j,k) = dcoef* ((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j) - half*C1_8)
               Afdg2(i,j,k) = ((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j)) 
            enddo
         enddo
!$omp enddo
!$omp do collapse(2)
         do k = k01,Nk-1
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w-1, l_ni-pil_e
               qkm =  half* (fdg2_4(i,j,k    )*Jxt(i,j,k,1) - fdg2_4(i,j,k-1  )*Jxt(i,j,k,2) +&
                             fdg2_4(i+1,j,k  )*Jxt(i,j,k,3) - fdg2_4(i+1,j,k-1)*Jxt(i,j,k,4))
               qkp = half*(fdg2_4(i,j,k+1  )*Jxt(i,j,k,5) - fdg2_4(i,j,k )*Jxt(i,j,k,6) +&
                           fdg2_4(i+1,j,k+1)*Jxt(i,j,k,7) - fdg2_4(i+1,j,k )*Jxt(i,j,k,8))
               C=   half*(qkp+qkm)

               Afdg1(i,j,k) = dcoef*((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j) - C)
               Afdg2(i,j,k) =((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j))
               enddo
            enddo
         enddo
!$omp enddo
         k=k00
!$omp do
         do j=1+pil_s-1, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               Bfdg1(i,j,k) = dcoef*(Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k)) * geomh_invDYMv_8(j) & 
                         * geomh_cyv_8(j)
               Bfdg2(i,j,k) =(Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
            enddo
         enddo
!$omp enddo

         k= NK
!$omp do
         do j=1+pil_s-1, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               C1_8=half*(fdg2_4(i,j,k  )*Jyt(i,j,k,1)-fdg2_4(i,j,k-1)*Jyt(i,j,k,2) +&
                          fdg2_4(i,j+1,  k)*Jyt(i,j,k,3)-fdg2_4(i,j+1,k-1)*Jyt(i,j,k,4))
               Bfdg1(i,j,k) = dcoef*((Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j) - half*C1_8)* &
                             geomh_cyv_8(j)
               Bfdg2(i,j,k) = (Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
            enddo
         enddo
!$omp enddo
!$omp do collapse(2)
         do k = k01,Nk-1
            do j=1+pil_s-1, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  qkm=half*(fdg2_4(i,j,k   )*Jyt(i,j,k,1)-fdg2_4(i,j,k-1)*Jyt(i,j,k,2  )+&
                            fdg2_4(i,j+1,k )*Jyt(i,j,k,3)-fdg2_4(i,j+1,k-1)*Jyt(i,j,k,4))
                  qkp=half*(fdg2_4(i,j,k+1)*Jyt(i,j,k,5)-fdg2_4(i,j,k )*Jyt(i,j,k,6)+&
                            fdg2_4(i,j+1,k+1)*Jyt(i,j,k,7)-fdg2_4(i,j+1,k )*Jyt(i,j,k,8))
                  C = half*(qkp+qkm)

                  Bfdg1(i,j,k) =dcoef* ((Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j) - C) &
                        * geomh_cyv_8(j)
                  Bfdg2(i,j,k) =(Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
               enddo
            enddo
         enddo
!$omp enddo

! Apply divergence
!$omp do collapse(2)
         do k = k00, nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  add_v8(i,j,k) = Jm(i,j,k,1)*(Afdg1 (i,j,k)-Afdg1 (i-1,j,k))
                  bdd_v8(i,j,k) = Jm(i,j,k,2)*(Bfdg1 (i,j,k)-Bfdg1 (i,j-1,k))
               enddo
            enddo
         enddo
!$omp enddo
!  flux
!$omp do collapse(2)
      do k=k00+1,Nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
                  !Jzz(i,j,k,1)  = (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))
                  cdd_v81(i,j,k) = (Ver_wp_8%m(k)*(half*(Afdg2(i-1,j,k)+Afdg2(i,j,k)))+ &
                          Ver_wm_8%m(k)*(half*(Afdg2(i-1,j,k-1)+Afdg2(i,j,k-1))))
! put zero if using stencil
                  cdd_v81(i,j,k) = half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i,j,k))*Jzz(i,j,k,1)* dcoef*cdd_v81(i,j,k)
                  cdd_v82(i,j,k) = Ver_wp_8%m(k)*(half*(Bfdg2(i,j-1,k)+Bfdg2(i,j,k)))+ &
                                   Ver_wm_8%m(k)*(half*(Bfdg2(i,j-1,k-1)+Bfdg2(i,j,k-1)))
! put zero if using stencil
                  cdd_v82(i,j,k)= half*(GVM%mc_Jy_8(i,j-1,k)+GVM%mc_Jy_8(i,j,k))* Jzz(i,j,k,1)*dcoef*geomh_cy_8(j)*&
                                  cdd_v82(i,j,k)
               enddo
            enddo
         enddo
!$omp enddo

!$omp do collapse(2)
         do k=k00,Nk-1
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  cflux(i,j,k)=Jzz(i,j,k,2)*Ver_idz_8%t(k)* &
                  ((cdd_v81(i,j,k+1)-cdd_v81(i,j,k)) + &
                    (cdd_v82(i,j,k+1)-cdd_v82(i,j,k))*geomh_invcy_8(j) )
                  F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-cflux(i,j,k))
               enddo
            enddo
         enddo
!$omp enddo
         k=Nk
!$omp do
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               cflux(i,j,k)= (one-(ver_z_8%t(Nk)-ver_z_8%t(Nk-1))/(ver_z_8%t(Nk+1)-ver_z_8%t(Nk-1)))*&
                                    cflux(i,j,Nk-1)
               F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-cflux(i,j,k))
            enddo
         enddo
!$omp enddo

!implicit
!$omp do
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               do k = k00+1 , Nk
                  F_sol1(i,j,k) = F_sol1(i,j,k) - W_zdt(i,j,k) * F_sol1(i,j,k-1)
               enddo
               F_sol1(i,j,Nk) = F_sol1(i,j,Nk) / b_zdt(i,j,Nk)
               do k = Nk-1, k00, -1
                  F_sol1(i,j,k) = (F_sol1(i,j,k) - c_zdt(i,j,k) * F_sol1(i,j,k+1)) / b_zdt(i,j,k)
               enddo
            enddo
         enddo
!$omp enddo

      ! Hybrid diffusion if hzd_hyb_bot >0
      if(hzd_hyb_bot >0) then
!$omp do collapse(2) 
         do k = nk-hzd_hyb_bot+1, nk
            do j=1+pil_s-1, l_nj-pil_n+1
               do i=1+pil_w-1, l_ni-pil_e+1
                  F_sol1(i,j,k) = fdg2_4(i,j,k )
               enddo
            enddo
         enddo
!$omp enddo
      endif

! hybrid difusion for first kd0+1 level
      if (kd0.ne.1) then
!$omp do collapse(2) 
          do k = 1, kd0 
            do j=1+pil_s-1, l_nj-pil_n+1
               do i=1+pil_w-1, l_ni-pil_e+1
                  F_sol1(i,j,k) = fdg2_4(i,j,k )
               enddo
            enddo
         enddo
!$omp enddo
      endif

      return
      end

