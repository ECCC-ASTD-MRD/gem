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
      subroutine  hzd_u_alh ( F_Sol1,HzdlnR, Minx, Maxx, Miny, Maxy,Nk) 

      use step_options
      use gem_options
      use geomh
      use glb_ld
      use cstv
      use ver
      use metric
      use hzd_mod
      use hvdif_options
      use dcst
      use tdpack
      use gmm_geof
      use ptopo
      use step_options
      use omp_timing
      use stat_mpi, only: statf_dm
!
      use, intrinsic :: iso_fortran_env
      implicit none
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, NK
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent (inout) :: F_Sol1
      real  HzdlnR


!author
!       Abdessamad Qaddouri -  2018
!
!revision
! v5.0 - Qaddouri A.       - initial version


      integer j,i,k,kp,km
      integer kd0 , k00, k01
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)
      real   fdg2_4(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64)   Afdg1(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64)   Bfdg1(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64)   Afdg2(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64)   Bfdg2(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) C1_8,C2_8,C,C3_8
      real(kind=REAL64) bdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) add_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) cdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) cdd_v82(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) cdd_v81(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) crit_coef,base_coefT, dcoef

      real(kind=REAL64)  ztht_8(l_minx:l_maxx, l_miny:l_maxy,0:Nk+1)
!     ------Ui-1,j+1-------------phii,j+1------------Ui,j+1----------------------------
!                                 Vi,j                Bi,j 
!     ------Ui-1,j---------------Ai,jandphii,j--------Ui,j-------------------------
!
!     ------Ui-1,j-1---------------phii,j-1-----------Ui,j-1-----------------------------

!  kd0 given by user      
      kd0=hzd_hyb_top

      k00=kd0
      k01=kd0+1 
      if (kd0.ne.1) k01=kd0

      call gtmg_start (69, 'HZD_u_alh', 65)

      if (Hzd_pwr_z==2) then 
         dcoef = 0.25*HzdlnR*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      else
         dcoef = 0.25*sqrt(HzdlnR)*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      endif

      Afdg1 = .0d0
      Bfdg1 = .0d0
      Afdg2= zero
      Bfdg2= zero
      add_v8 =0.0d0
      bdd_v8=0.0d0
      cdd_v8=0.0d0
      cdd_v81=0.d0
      cdd_v82=0.0d0
      fdg2_4 =0.0 
!
      do k = 1, nk
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w-1, l_ni-pil_e+1
               fdg2_4(i,j,k )=F_Sol1(i,j,k)
            enddo
         enddo
      enddo

      call rpn_comm_xch_halo(fdg2_4,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,Nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
! aplly gradient
! gradient component along X
      k=k00
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e+1
             Afdg1(i,j,k) = dcoef*((skpu(i,j,k,1)*fdg2_4(i,j,k) &
                         -sku(i,j,k,1)*fdg2_4(i-1,j,k))*geomh_invDXM_8(j)) ! on M-level k and phi-ij position 
             Afdg2(i,j,k) = Afdg1(i,j,k) 
          enddo
        enddo
!
      k= NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e+1
               C2_8= half*Ver_wm_8%m(k)*((fdg2_4(i-1,j,k)*GVM%mc_Jx_8(i-1,j,k)-&
                      fdg2_4(i-1,j,k-1)*GVM%mc_Jx_8(i-1,j,k-1))+ &
                      (fdg2_4(i,j,k)*GVM%mc_Jx_8(i,j,k)-&
                      fdg2_4(i,j,k-1)*GVM%mc_Jx_8(i,j,k-1)))*Ver_idz_8%t(k-1)
              Afdg1(i,j,k) =dcoef*((skpu(i,j,k,1)*fdg2_4(i,j,k)&
                            -sku(i,j,k,1)*fdg2_4(i-1,j,k)) *geomh_invDXM_8(j)-C2_8)
              Afdg2(i,j,k) = dcoef*((skpu(i,j,k,1)*fdg2_4(i,j,k) &
                          -sku(i,j,k,1)*fdg2_4(i-1,j,k))*geomh_invDXM_8(j)) ! on M-level k and phi-ij position 
          enddo
        enddo

      do k = k01,Nk-1
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e+1
               C=half* ((fdg2_4(i,j,k+1)*GVM%mc_Jx_8(i,j,k+1)-&
                     fdg2_4(i,j,k  )*GVM%mc_Jx_8(i,j,k  ))+&
                     (fdg2_4(i-1,j,k+1)*GVM%mc_Jx_8(i-1,j,k+1)-&
                     fdg2_4(i-1,j,k  )*GVM%mc_Jx_8(i-1,j,k  )))*Ver_idz_8%t(k)
               C3_8=half* ((fdg2_4(i,j,k)*GVM%mc_Jx_8(i,j,k)-&
                     fdg2_4(i,j,k-1)*GVM%mc_Jx_8(i,j,k-1)) +&
                     (fdg2_4(i-1,j,k)*GVM%mc_Jx_8(i-1,j,k)-&
                     fdg2_4(i-1,j,k-1)*GVM%mc_Jx_8(i-1,j,k-1)))*Ver_idz_8%t(k-1)
! (D(jx*U)/Dzeta)  on M-level k and Ai,j position
               C2_8  = Ver_wp_8%m(k)*C+Ver_wm_8%m(k)* C3_8
! D(J_zeta*U)/Dx-(D(jx*U)/Dzeta) on  M-level k and phi,j position
              Afdg1(i,j,k) =dcoef*((skpu(i,j,k,1)*fdg2_4(i,j,k) & 
                         -sku(i,j,k,1)*fdg2_4(i-1,j,k))*geomh_invDXM_8(j)-C2_8)
              Afdg2(i,j,k) =dcoef*((skpu(i,j,k,1)*fdg2_4(i,j,k) &
                         -sku(i,j,k,1)*fdg2_4(i-1,j,k))*geomh_invDXM_8(j)) ! on M-level k and phi-ij position 
            enddo
         enddo
         enddo
! Gradient component Along Y
         k=k00
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
            Bfdg1(i,j,k) =dcoef*((sku(i,j,k,2)*fdg2_4(i,j+1,k) -skpu(i,j,k,2)*fdg2_4(i,j,k))* geomh_invDY_8)  ! M-level B_ij position 
            Bfdg2(i,j,k) =Bfdg1(i,j,k)

          enddo
        enddo
!
      k= NK
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
! D(jyU)/Dzeta on k T-level on ui,j position
             C2_8= half*Ver_wm_8%m(k)*((Jyp(i,j,k,1)*fdg2_4(i,j,k)-Jy(i,j,k,1)*fdg2_4(i,j,k-1)) + &
                   (Jyp(i,j,k,2)*fdg2_4(i,j+1,k)-Jy(i,j,k,2)*fdg2_4(i,j+1,k-1)))*Ver_idz_8%t(k-1)
               Bfdg1(i,j,k) = dcoef*((sku(i,j,k,2)*fdg2_4(i,j+1,k) -skpu(i,j,k,2)*fdg2_4(i,j,k))* geomh_invDY_8 - C2_8)
               Bfdg2(i,j,k) =dcoef*((sku(i,j,k,2)*fdg2_4(i,j+1,k) -skpu(i,j,k,2)*fdg2_4(i,j,k))* geomh_invDY_8)  ! M-level B_ij position 

          enddo
        enddo
!
      do k = k01,Nk-1
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
! D(jyU)/Dzeta on k-1 T-level on ui,j position
             C2_8=half*( (Ver_wm_8%m(k)*((Jyp(i,j,k,1)*fdg2_4(i,j,k)-Jy(i,j,k,1)*fdg2_4(i,j,k-1)) + &
                        (Jyp(i,j,k,2)*fdg2_4(i,j+1,k)-Jy(i,j,k,2)*fdg2_4(i,j+1,k-1)))*Ver_idz_8%t(k-1)) +&
                       (Ver_wp_8%m(k)*( (Jyp(i,j,k,3)*fdg2_4(i,j,k+1)-Jy(i,j,k,3)*fdg2_4(i,j,k)) + &
                         (Jyp(i,j,k,4)*fdg2_4(i,j+1,k+1)-Jy(i,j,k,4)*fdg2_4(i,j+1,k)))*Ver_idz_8%t(k)))
! D(J_zeta*U)/Dy-(D(jy*U)/Dzeta) on  M-level k and B_ij position
               Bfdg1(i,j,k) =  dcoef*((sku(i,j,k,2)*fdg2_4(i,j+1,k) -skpu(i,j,k,2)*fdg2_4(i,j,k))* geomh_invDY_8 - C2_8)
               Bfdg2(i,j,k) =  dcoef*((sku(i,j,k,2)*fdg2_4(i,j+1,k) -skpu(i,j,k,2)*fdg2_4(i,j,k))* geomh_invDY_8 )  ! M-level B_ij position 

            enddo
         enddo
         enddo
! Apply divergence
      do k = k00, nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
!inverse(D(z)/D(zeta)) on k M-level and Ui,j position
!X-divergence (Dz/Dzeta)^-1*DAfdg/Dx on  k M-level and Ui,j position
         add_v8(i,j,k) = xfactu(i,j,k)* (Afdg1 (i+1,j,k)-Afdg1 (i,j,k))*geomh_invDXM_8(j)
         bdd_v8(i,j,k) = xfactu(i,j,k)* (Bfdg1 (i,j,k)*geomh_cyv_8(j)-Bfdg1 (i,j-1,k)*geomh_cyv_8(j-1))*geomh_invDYM_8(j)
            enddo
         enddo
      enddo
!  flux
      do k = k00+1,Nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
                !cdd_v81(i,j,k) =half*Ver_wm_8%m(k)* ((GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i-1,j,k)) *Afdg2(i,j,k)/Jzx(i,j,k,1) -&
                cdd_v81(i,j,k) = half*Ver_wm_8%m(k)*(Afdg2(i,j,k  )*Jizx(i,j,k,1 ) - &
                                                     Afdg2(i,j,k-1)*Jizxm(i,j,k,1))*Ver_idz_8%t(k-1) + &
                                 half*Ver_wp_8%m(k)*(Afdg2(i,j,k+1)*Jizpx(i,j,k,1) - & 
                                                     Afdg2(i,j,k  )*Jizx(i,j,k,1 ))*Ver_idz_8%t(k)

                cdd_v82(i,j,k) = half*Ver_wm_8%m(k)*(Bfdg2(i,j,k  )*Jizpx(i,j,k,2) -&
                                                     Bfdg2(i,j,k-1)*Jizxm(i,j,k,2))*Ver_idz_8%t(k-1) + &
                                 half*Ver_wp_8%m(k)*(Bfdg2(i,j,k+1)*Jizx(i,j,k,3 ) -&
                                                     Bfdg2(i,j,k  )*Jizpx(i,j,k,2))*Ver_idz_8%t(k)
           enddo
         enddo
      enddo

      k=Nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
!Dz(Jx/Jz A) on T-level K-1 and position Ai,j
                cdd_v81(i,j,k)=  half*Ver_wm_8%m(k)*(Afdg2(i,j,k  )*Jizx (i,j,k,1) -&
                                                     Afdg2(i,j,k-1)*Jizxm(i,j,k,1))*Ver_idz_8%t(k-1)
                cdd_v82(i,j,k)=  half*Ver_wm_8%m(k)*(Bfdg2(i,j,k  )*jizpx(i,j,k,2) -&
                                                     Bfdg2(i,j,k-1)*Jizxm(i,j,k,2))*Ver_idz_8%t(k-1)
           enddo
         enddo

      do k=k00,NK-1
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
! xfact* cdd on k M-level and Ui,j position
               cdd_v8(i,j,k)= xfactu(i,j,k) *(half*(cdd_v81(i,j,k)+ cdd_v81(i+1,j,k))+ half*(cdd_v82(i,j,k)+ cdd_v82(i,j-1,k))/geomh_cy_8(j))
! boundary condition
               F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-cdd_v8(i,j,k))

            enddo
         enddo
      enddo

      k=NK
      do j=1+pil_s, l_nj-pil_n
         do i=1+pil_w, l_ni-pil_e
            F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-ru*cdd_v8(i,j,k-1))
         enddo
      enddo

      do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
              do k = k00+1 , Nk
                F_sol1(i,j,k) = F_sol1(i,j,k) - W_u(i,j,k) * F_sol1(i,j,k- 1)
              enddo
                 F_sol1(i,j,Nk) = F_sol1(i,j,Nk) / b_u(i,j,Nk)
                  do  k = Nk-1, k00, -1
                  F_sol1(i,j,k) = (F_sol1(i,j,k) - c_u(i,j,k) * F_sol1(i,j,k + 1)) / b_u(i,j,k)
                  enddo
             enddo
          enddo

      ! Hybrid diffusion if hzd_hyb_bot >0
      if(hzd_hyb_bot > 0) then
         do k = nk-hzd_hyb_bot+1, nk
            do j=1+pil_s-1, l_nj-pil_n+1
               do i=1+pil_w-1, l_ni-pil_e+1
                  F_sol1(i,j,k) = fdg2_4(i,j,k )
               enddo
            enddo
         enddo
      endif
! hybrid difusion for first kd0+1 level
      if (kd0.ne.1) then
          do k = 1, kd0 
            do j=1+pil_s-1, l_nj-pil_n+1
               do i=1+pil_w-1, l_ni-pil_e+1
                  F_sol1(i,j,k) = fdg2_4(i,j,k )
               enddo
            enddo
         enddo
      endif


        call gtmg_stop (69)

      return
      end

