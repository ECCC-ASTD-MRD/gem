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
      subroutine  hzd_v_alh ( F_Sol1,HzdlnR , Minx, Maxx, Miny, Maxy,Nk)
      use gem_options
      use geomh
      use glb_ld
      use cstv
      use ver
      use metric
      use hzd_mod
      use hvdif_options
      use dcst
!
      use tdpack
      use gmm_geof
      use ptopo
      use step_options
      use stat_mpi, only: statf_dm
      use omp_timing
!
      use, intrinsic :: iso_fortran_env
      implicit none
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, NK
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent (inout) :: F_Sol1

!author
!       Abdessamad Qaddouri -  2018
!
!revision
! v5.0 - Qaddouri A.       - initial version


      integer j,i,k
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)
      integer  km, kp
      integer kd0 , k00, k01
      real(kind=REAL64)   Afdg1(l_minx:l_maxx, l_miny:l_maxy,Nk),Afdg2(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64)   Bfdg1(l_minx:l_maxx, l_miny:l_maxy,Nk),Bfdg2(l_minx:l_maxx, l_miny:l_maxy,Nk)

      real(kind=REAL64) Jzpi,Jz,Jzm,Jzmpi
      real(kind=REAL64) C1_8,C2_8,C,ski,skpi,C3_8
      real(kind=REAL64) bdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) add_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) cdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) cdd_v81(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) cdd_v82(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) dcoef
      real  HzdlnR
      real  fdg2_4(l_minx:l_maxx, l_miny:l_maxy,Nk+1)

!
!     ---------------------Vi,j--------Afdgi,j ----------------------------------
!
!     ---------------------Bfdgi,j---------Ui,j--------------------------------
!
!     ---------------------------------------------------------------

      call gtmg_start (68, 'HZD_v_alh', 65)
!  kd0 given by user      
      kd0=hzd_hyb_top

      k00=kd0
      k01=kd0+1
      if (kd0.ne.1) k01=kd0

      if (Hzd_pwr_z==2) then 
         dcoef = 0.25*HzdlnR*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      else
         dcoef = 0.25*sqrt(HzdlnR)*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      endif

      Afdg1 = .0d0
      Bfdg1 = .0d0
      Afdg2 = .0d0
      Bfdg2 = .0d0
      add_v8 =0.0d0
      bdd_v8=0.0d0
      cdd_v8=0.0d0
      cdd_v81=0.0d0
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
            do i=1+pil_w-1, l_ni-pil_e
! D(J_zeta*V)/Dx-(D(jx*V)/Dzeta) on  M-level k and phi,j position
           Afdg1(i,j,k) =dcoef*((skpv(i,j,k,1)*fdg2_4(i+1,j,k) -skv(i,j,k,1)*fdg2_4(i,j,k))*geomh_invDXv_8(j))
           Afdg2(i,j,k) =dcoef*((skpv(i,j,k,1)*fdg2_4(i+1,j,k) -skv(i,j,k,1)*fdg2_4(i,j,k))*geomh_invDXv_8(j))
          enddo
        enddo
!
      k= NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
! (D(jx*V)/Dzeta)  on T-level k and Vij position
               C= (Jxp(i,j,k,1)*fdg2_4(i,j,k)-Jx(i,j,k,1)*fdg2_4(i,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! Jx on M K level and Vi+1,j 
! (D(jx*V)/Dzeta)  on T-level k and Vi+1j position
               C2_8= (Jxp(i,j,k,2)*fdg2_4(i+1,j,k)-Jx(i,j,k,2)*fdg2_4(i+1,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! (D(jx*V)/Dzeta)  on T-level k-1 and Aij position
               C2_8= half*( C+ C2_8)
! (D(jx*V)/Dzeta)  on M-level k and Aij position! zero lower condition
               C2_8= Ver_wp_8%m(k)*zero + Ver_wm_8%m(k) * C2_8          
! D(J_zeta*V)/Dx-(D(jx*V)/Dzeta) on  M-level k and phi,j position
              Afdg1(i,j,k) =dcoef*((skpv(i,j,k,1)*fdg2_4(i+1,j,k) -skv(i,j,k,1)*fdg2_4(i,j,k)) *geomh_invDXv_8(j)-C2_8)
              Afdg2(i,j,k) =dcoef*((skpv(i,j,k,1)*fdg2_4(i+1,j,k) -skv(i,j,k,1)*fdg2_4(i,j,k)) *geomh_invDXv_8(j))
          enddo
        enddo
!       
      do k = k01,Nk-1
         km=min(k-1,1)
         kp = max(k+1,NK)
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
! (D(jx*V)/Dzeta)  on T-level k and Vij position
               C= (Jxp(i,j,k,1)*fdg2_4(i,j,k)-Jx(i,j,k,1)*fdg2_4(i,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! (D(jx*V)/Dzeta)  on T-level k and Vi+1j position
               C2_8= (Jxp(i,j,k,2)*fdg2_4(i+1,j,k)-Jx(i,j,k,2)*fdg2_4(i+1,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! (D(jx*V)/Dzeta)  on T-level k-1 and Aij position
               C1_8= half*( C+ C2_8)
! (D(jx*V)/Dzeta)  on T-level k and Vij position
               C= (Jxp(i,j,k,3)*fdg2_4(i,j,k+1)-Jx(i,j,k,3)*fdg2_4(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! (D(jx*V)/Dzeta)  on T-level k and Vi+1j position
               C2_8= (Jxp(i,j,k,4)*fdg2_4(i+1,j,k+1)-Jx(i,j,k,4)*fdg2_4(i+1,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! (D(jx*V)/Dzeta)  on T-level k and Aij position
               C2_8= half*( C+ C2_8)

! (D(jx*V)/Dzeta)  on M-level k and Vij position
               C2_8  = Ver_wp_8%m(k)*C2_8+Ver_wm_8%m(k)* C1_8
! D(J_zeta*V)/Dx-(D(jx*V)/Dzeta) on  M-level k and phi,j position
              Afdg1(i,j,k) =dcoef*((skpv(i,j,k,1)*fdg2_4(i+1,j,k) -skv(i,j,k,1)*fdg2_4(i,j,k))*geomh_invDXv_8(j)-C2_8)
              Afdg2(i,j,k) =dcoef*((skpv(i,j,k,1)*fdg2_4(i+1,j,k) -skv(i,j,k,1)*fdg2_4(i,j,k))*geomh_invDXv_8(j))
            enddo
         enddo
         enddo
! Gradient component Along Y
         k=k00
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
            Bfdg1(i,j,k) =dcoef*((skpv(i,j,k,2)*fdg2_4(i,j,k) -skv(i,j,k,2)*fdg2_4(i,j-1,k))* geomh_invDY_8)
            Bfdg2(i,j,k) =dcoef*((skpv(i,j,k,2)*fdg2_4(i,j,k) -skv(i,j,k,2)*fdg2_4(i,j-1,k))* geomh_invDY_8)
          enddo
        enddo
      k= NK
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
! (D(jy*V)/Dzeta)  on T-level k-1 and Bij position
               C2_8= (GVM%mc_Jy_8(i,j,k)*fdg2_4(i,j,k)-GVM%mc_Jy_8(i,j,k-1)*fdg2_4(i,j,k-1))*Ver_idz_8%t(k-1)
               C3_8= (GVM%mc_Jy_8(i,j-1,k)*fdg2_4(i,j-1,k)-GVM%mc_Jy_8(i,j-1,k-1)*fdg2_4(i,j-1,k-1))*Ver_idz_8%t(k-1)
               C= half*(C2_8+C3_8)
! (D(jy*V)/Dzeta)  on M-level k and Bij position. lower BD condition
               C2_8= Ver_wp_8%m(k)*zero +Ver_wm_8%m(k)*C
! D(J_zeta*V)/Dy-(D(jy*V)/Dzeta) on  M-level k and Bi,j position
            Bfdg1(i,j,k) =dcoef*((skpv(i,j,k,2)*fdg2_4(i,j,k) -skv(i,j,k,2)*fdg2_4(i,j-1,k))* geomh_invDY_8 - C2_8)
            Bfdg2(i,j,k) =dcoef*((skpv(i,j,k,2)*fdg2_4(i,j,k) -skv(i,j,k,2)*fdg2_4(i,j-1,k))* geomh_invDY_8 )
          enddo
        enddo


      do k = k01,Nk-1
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
! (D(jy*V)/Dzeta)  on T-level k-1 and Bij position
               C= half*((GVM%mc_Jy_8(i,j,k  )*fdg2_4(i,j,k  )-GVM%mc_Jy_8(i,j,k-1)*fdg2_4(i,j,k-1   )) +&
                        (GVM%mc_Jy_8(i,j-1,k)*fdg2_4(i,j-1,k)-GVM%mc_Jy_8(i,j-1,k-1)*fdg2_4(i,j-1,k-1)))&
                            *Ver_idz_8%t(k-1)
! (D(jy*V)/Dzeta)  on T-level k and Bij position
               C3_8= half*( (GVM%mc_Jy_8(i,j,k+1)*fdg2_4(i,j,k+1)-GVM%mc_Jy_8(i,j,k)*fdg2_4(i,j,k)) +&
                            (GVM%mc_Jy_8(i,j-1,k+1)*fdg2_4(i,j-1,k+1)-GVM%mc_Jy_8(i,j-1,k)*fdg2_4(i,j-1,k)) )&
                            *Ver_idz_8%t(k)
! (D(jy*V)/Dzeta)  on M-level k and Bij position
               C2_8= Ver_wp_8%m(k)*C3_8+Ver_wm_8%m(k)*C
! D(J_zeta*V)/Dy-(D(jy*V)/Dzeta) on  M-level k and Bi,j position
            Bfdg1(i,j,k) =dcoef*((skpv(i,j,k,2)*fdg2_4(i,j,k) -skv(i,j,k,2)*fdg2_4(i,j-1,k))* geomh_invDY_8 - C2_8)
            Bfdg2(i,j,k) =dcoef*((skpv(i,j,k,2)*fdg2_4(i,j,k) -skv(i,j,k,2)*fdg2_4(i,j-1,k))* geomh_invDY_8)
            enddo
         enddo
         enddo
! Apply divergence
      do k = k00, nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
!X-divergence (Dz/Dzeta)^-1*DAfdg/Dx on  k M-level and Vi,j position
               add_v8(i,j,k) = xfactv(i,j,k)* (Afdg1 (i,j,k)-Afdg1 (i-1,j,k))*geomh_invDXv_8(j)
            enddo
         enddo
      enddo
      do k = k00, nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
          bdd_v8(i,j,k) = xfactv(i,j,k)*(Bfdg1 (i,j+1,k)*geomh_cy_8(j+1)-Bfdg1 (i,j,k)*geomh_cy_8(j))*geomh_invDYM_8 (j)
            enddo
         enddo
      enddo
!  flux
! start at pil_w and finish at l_nj-pil_n+1 ! important pour le calcul de sol a la fin
      do k = k00+1,Nk
         do j=1+pil_s, l_nj-pil_n+1
            do i=pil_w, l_ni-pil_e
               cdd_v81(i,j,k-1)=  half*((GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i,j+1,k)) *Afdg2(i,j,k)*Jzv(i,j,k,1) -&
                      (GVM%mc_Jx_8(i,j,k-1) +GVM%mc_Jx_8(i,j+1,k-1)) *Afdg2(i,j,k-1)*Jzvm(i,j,k,1))*Ver_idz_8%t(k-1)
!   Dz(Jy/Jz B) on T-level K-1 and position Bi,j
               cdd_v82(i,j,k-1)=  half* ((GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i,j-1,k)) *Bfdg2(i,j,k)*Jzv(i,j,k,2) -&
                      (GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i,j-1,k-1)) *Bfdg2(i,j,k-1)*Jzvm(i,j,k,2))*Ver_idz_8%t(k-1)
           enddo
         enddo
        enddo

      k=k00  
      do j=1+pil_s, l_nj-pil_n
         do i=1+pil_w, l_ni-pil_e
               F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k))
         enddo
      enddo
      do k=k00+1,NK-1
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
! Dz/Dzeta^-1* Dz(Jx/Jz A+ Jy/Jz B) on M-level K and position Bi,j
               cdd_v8(i,j,k)=Jzv(i,j,k,3)*(  &
                           half *( ( (Ver_wp_8%m(k)*cdd_v81(i,j,k)+Ver_wm_8%m(k)* cdd_v81(i,j,k-1))+&  
                          (Ver_wp_8%m(k)*cdd_v81(i+1,j,k)+Ver_wm_8%m(k)* cdd_v81(i+1,j,k-1)) ) + &
                          ((Ver_wp_8%m(k)*cdd_v82(i,j,k)+Ver_wm_8%m(k)* cdd_v82(i,j,k-1)) +&
                          (Ver_wp_8%m(k)*cdd_v82(i,j+1,k)+Ver_wm_8%m(k)* cdd_v82(i,j+1,k-1)))/geomh_cyv_8(j)) )
               F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-cdd_v8(i,j,k))
            enddo
         enddo
      enddo

      k=Nk  
      do j=1+pil_s, l_nj-pil_n
         do i=1+pil_w, l_ni-pil_e
            !cdd_v8(i,j,k)= (one-(ver_z_8%m(Nk)-ver_z_8%m(Nk-1))/(ver_z_8%m(Nk+1)-ver_z_8%m(Nk-1)))*&
            !cdd_v8(i,j,k)= rv*&
            !            cdd_v8(i,j,Nk-1)
            !F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-cdd_v8(i,j,k))
            F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-rv*cdd_v8(i,j,k-1))
         enddo
      enddo

! Implicit v
      do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
              do k = k00+1 , Nk
                F_sol1(i,j,k) = F_sol1(i,j,k) - W_v(i,j,k) * F_sol1(i,j,k- 1)
              enddo
                 F_sol1(i,j,Nk) = F_sol1(i,j,Nk) / b_v(i,j,Nk)
                  do  k = Nk-1, k00, -1
                  F_sol1(i,j,k) = (F_sol1(i,j,k) - c_v(i,j,k) * F_sol1(i,j,k + 1)) / b_v(i,j,k)
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


       call gtmg_stop (68)

      return
      end

