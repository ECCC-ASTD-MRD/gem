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
! Ajout
      subroutine  hzd_u1_alh ( F_Sol1,HzdlnR,Minx, Maxx, Miny, Maxy,Nk,Niter) 
! fin Ajout
!
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
      use stat_mpi, only: statf_dm
!
      use, intrinsic :: iso_fortran_env
      implicit none
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, NK
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent (inout) :: F_Sol1
      real(kind=REAL64) xfact
      real  HzdlnR
  


!author
!       Abdessamad Qaddouri -  2018
!
!revision
! v5.0 - Qaddouri A.       - initial version


      integer j,i,k,halox,haloy
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)
      real(kind=REAL64)   Afdg1(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64)   Bfdg1(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
!
      real(kind=REAL64)   Afdg2(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64)   Bfdg2(l_minx:l_maxx, l_miny:l_maxy,Nk+1)

!
      real   fdg2_4(l_minx:l_maxx, l_miny:l_maxy,Nk+1)

      real(kind=REAL64) Jzpi,Jz,Jzm,Jzmpi
      real(kind=REAL64) C1_8,C2_8,C,ski,skpi,C3_8
      real(kind=REAL64) bdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) add_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) cdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      integer iter ,Niter
      real(kind=REAL64) F_coef_8(1:NK),Jy, Jyp
      real(kind=REAL64) cdd_v82(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) cdd_v81(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) crit_coef,base_coefT
      integer kp,km 

!     ------Ui-1,j+1-------------phii,j+1------------Ui,j+1----------------------------
!                                 Vi,j                Bi,j 
!     ------Ui-1,j---------------Ai,jandphii,j--------Ui,j-------------------------
!
!     ------Ui-1,j-1---------------phii,j-1-----------Ui,j-1-----------------------------

      real(kind=REAL64)  ztht_8(l_minx:l_maxx, l_miny:l_maxy,0:Nk+1)
      real(kind=REAL64), dimension (:,:,:,:), allocatable :: stencil_V
      real(kind=REAL64) a(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) b(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) d(l_minx:l_maxx, l_miny:l_maxy,Nk),W,Jxx

!
      do j=1-G_haloy,l_nj+G_haloy
        do i=1-G_halox,l_ni+G_halox
          do k=1 ,NK
           ztht_8(i,j,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*fis0(i,j)+Ver_c_8%t(k)*sls(i,j))/grav_8
          enddo
            ztht_8(i,j,0)   =   GVM%zmom_8(i,j,0)
            ztht_8(i,j,Nk+1)=   GVM%zmom_8(i,j,Nk+1)
        enddo
      enddo

      crit_coef =(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      if (Hzd_pwr_z==2)  then
              base_coefT= 0.25d0*HzdlnR*crit_coef
      else       
               base_coefT=0.25d0*sqrt(HzdlnR)*crit_coef 
      endif            

      F_coef_8(1:NK)= base_coefT 
      halox=1
      haloy=halox
     do iter =1,Niter 
      Afdg1 = .0d0
      Bfdg1 = .0d0
!
      Afdg2= zero
      Bfdg2= zero
!
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
! aplly gradient
! gradient component along X
      k=1
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e+1
!Dz/Dzeta M-level k on Ui,j position
               Jzpi =(ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))!M-level k on phi+1,j 
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))!M-level k on phi,j
               skpi= (Jzpi+Jz)*half !Dz/Dzeta on M-level k and  Ui,j position
!Dz/Dzeta M-level k on Ui-1,j position!
               Jzpi =(ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))!M-level k on phi,j
               Jz   =(ztht_8(i-1  ,j,k)-ztht_8(i-1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))!M-level k on phi-1,j
               ski = (Jzpi+Jz)*half   !Dz/Dzeta on M-level k  and  Ui-1,j position 
! (D(jx*U)/Dzeta)  on T-level k and Uij position
               C2_8= (fdg2_4(i,j,k+1)*GVM%mc_Jx_8(i,j,k+1)-&
                     fdg2_4(i,j,k  )*GVM%mc_Jx_8(i,j,k  ))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! (D(jx*U)/Dzeta)  on T-level k and Ui-1,j position
               C= (fdg2_4(i-1,j,k+1)*GVM%mc_Jx_8(i-1,j,k+1)-&
                     fdg2_4(i-1,j,k  )*GVM%mc_Jx_8(i-1,j,k  ))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! (D(jx*U)/Dzeta)  on T-level k and phi,j=Ai,j position
               C= half*(C2_8+C)
!               
! (D(jx*U)/Dzeta)  on T-level k-1 and Ai,j position = upper Bd condition            
               C3_8=zero
! (D(jx*U)/Dzeta)  on M-level k and phi,j position
               C2_8  = Ver_wp_8%m(k)*C+Ver_wm_8%m(k)* C3_8
! D(J_zeta*U)/Dx-(D(jx*U)/Dzeta) on  M-level k and phi,j position
           Afdg1(i,j,k) = F_coef_8(k)*((skpi*fdg2_4(i,j,k) &
-ski*fdg2_4(i-1,j,k))*geomh_invDXM_8(j) - C2_8) ! on M-level k and phi-ij position 
!
           Afdg2(i,j,k) = F_coef_8(k)*((skpi*fdg2_4(i,j,k) &
-ski*fdg2_4(i-1,j,k))*geomh_invDXM_8(j) - zero) ! on M-level k and phi-ij position 
          enddo
        enddo
!
      k= NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e+1
!Dz/Dzeta M-level k on Ui,j position
               Jzpi =(ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               skpi= (Jzpi+Jz)*half  !Dz/Dzeta on M-level k  and  Ui,j position
!Dz/Dzeta M-level k on Ui-1,j position
               Jzpi =(ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i-1  ,j,k)-ztht_8(i-1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               ski = (Jzpi+Jz)*half !Dz/Dzeta on M-level k  and  Ui-1,j position
! (D(jx*U)/Dzeta)  on T-level k-1 and Uij position
               C2_8= (fdg2_4(i,j,k)*GVM%mc_Jx_8(i,j,k)-&
                     fdg2_4(i,j,k-1  )*GVM%mc_Jx_8(i,j,k-1  ))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! (D(jx*U)/Dzeta)  on T-level k-1 and Ui-1,j position
               C3_8= (fdg2_4(i-1,j,k)*GVM%mc_Jx_8(i-1,j,k)-&
                     fdg2_4(i-1,j,k-1  )*GVM%mc_Jx_8(i-1,j,k-1  ))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! (D(jx*U)/Dzeta)  on T-level k-1 and Ai,j position
               C3_8= half*(C2_8+C3_8)
! (D(jx*U)/Dzeta)  on T-level k and phi,j position: lower BD condition
               C=zero
! (D(jx*U)/Dzeta)  on M-level k and phi,j position
               C2_8  = Ver_wp_8%m(k)*C+Ver_wm_8%m(k)* C3_8
! D(J_zeta*U)/Dx-(D(jx*U)/Dzeta) on  M-level k and phi,j position
              Afdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k)&
-ski*fdg2_4(i-1,j,k)) *geomh_invDXM_8(j)-C2_8)
!
           Afdg2(i,j,k) = F_coef_8(k)*((skpi*fdg2_4(i,j,k) &
-ski*fdg2_4(i-1,j,k))*geomh_invDXM_8(j) - zero) ! on M-level k and phi-ij position 
          enddo
        enddo
!       
      do k = 2,Nk-1
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e+1
               Jzpi =(ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!Dz/Dzeta on M-level k  and  Ui,j position
               skpi= (Jzpi+Jz)*half
!
               Jzpi =(ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i-1  ,j,k)-ztht_8(i-1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!Dz/Dzeta on M-level k  and  Ui-1,j position
               ski = (Jzpi+Jz)*half
!
               C2_8= (fdg2_4(i,j,k+1)*GVM%mc_Jx_8(i,j,k+1)-&
                     fdg2_4(i,j,k  )*GVM%mc_Jx_8(i,j,k  ))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C= (fdg2_4(i-1,j,k+1)*GVM%mc_Jx_8(i-1,j,k+1)-&
                     fdg2_4(i-1,j,k  )*GVM%mc_Jx_8(i-1,j,k  ))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! (D(jx*U)/Dzeta)  on T-level k and Ai,j position
               C= half*(C2_8+C)
!
               C2_8= (fdg2_4(i,j,k)*GVM%mc_Jx_8(i,j,k)-&
                     fdg2_4(i,j,k-1  )*GVM%mc_Jx_8(i,j,k-1  ))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
               C3_8= (fdg2_4(i-1,j,k)*GVM%mc_Jx_8(i-1,j,k)-&
                     fdg2_4(i-1,j,k-1  )*GVM%mc_Jx_8(i-1,j,k-1  ))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! (D(jx*U)/Dzeta)  on T-level k-1 and Ai,j position
               C3_8= half*(C2_8+C3_8)
! (D(jx*U)/Dzeta)  on M-level k and Ai,j position
               C2_8  = Ver_wp_8%m(k)*C+Ver_wm_8%m(k)* C3_8
! D(J_zeta*U)/Dx-(D(jx*U)/Dzeta) on  M-level k and phi,j position
              Afdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k) -ski*fdg2_4(i-1,j,k))*geomh_invDXM_8(j)-C2_8)
!
           Afdg2(i,j,k) = F_coef_8(k)*((skpi*fdg2_4(i,j,k) &
-ski*fdg2_4(i-1,j,k))*geomh_invDXM_8(j) - zero) ! on M-level k and phi-ij position 
            enddo
         enddo
         enddo
! Gradient component Along Y
         k=1
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
               Jzpi =(ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!Dz/Dzeta on M-level and  Ui,j position
               skpi= (Jzpi+Jz)*half 
!
               Jzpi =(ztht_8(i+1,j+1,k)-ztht_8(i+1,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j+1,k)-ztht_8(i  ,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!  Dz/Dzeta on M-level and  Ui,j+1 position
               ski= (Jzpi+Jz)*half  
! jy k M-level  U_i,j position
               Jy  = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) & 
                   +half * (GVM%mc_Jy_8(i,j-1,k) +GVM%mc_Jy_8(i+1,j-1,k)))

! jy k+1 M-level  U_i,j position
               Jyp = half*(half * (GVM%mc_Jy_8(i,j,k+1) +GVM%mc_Jy_8(i+1,j,k+1))  &
                       +half * (GVM%mc_Jy_8(i,j-1,k+1) +GVM%mc_Jy_8(i+1,j-1,k+1)))
! D(jyU)/Dzeta on k T-level on ui,j position
             C2_8= (Jyp*fdg2_4(i,j,k+1)-Jy*fdg2_4(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! jy k M-level  U_i,j+1 position
               Jy  = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) &
                   +half * (GVM%mc_Jy_8(i,j+1,k) +GVM%mc_Jy_8(i+1,j+1,k)))

! jy k+1 M-level  U_i,j position
               Jyp = half*(half * (GVM%mc_Jy_8(i,j,k+1) +GVM%mc_Jy_8(i+1,j,k+1)) &
                    +   half * (GVM%mc_Jy_8(i,j+1,k+1) +GVM%mc_Jy_8(i+1,j+1,k+1)))
! D(jyU)/Dzeta on k T-level on ui,j+1 position
             C= (Jyp*fdg2_4(i,j+1,k+1)-Jy*fdg2_4(i,j+1,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! D(jyU)/Dzeta on k T-level on Bij
              C2_8= half*(C2_8+C)
! D(jyU)/Dzeta on k M-level  on Bi,j position ! condition zero en haut
               C2_8  = Ver_wp_8%m(k)*C2_8 +Ver_wm_8%m(k)* zero 
! D(J_zeta*U)/Dy-(D(jy*U)/Dzeta) on  M-level k and B_ij position
            Bfdg1(i,j,k) =F_coef_8(k)*((ski*fdg2_4(i,j+1,k) -skpi*fdg2_4(i,j,k))* geomh_invDY_8 - c2_8)  ! M-level B_ij position 
!
          Bfdg2(i,j,k) =F_coef_8(k)*((ski*fdg2_4(i,j+1,k) -skpi*fdg2_4(i,j,k))* geomh_invDY_8 - zero)  ! M-level B_ij position 

          enddo
        enddo
!
      k= NK
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
               Jzpi =(ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!Dz/Dzeta on k M-level and  Ui,j position
               skpi= (Jzpi+Jz)*half ! ui,j,k
!
               Jzpi =(ztht_8(i+1,j+1,k)-ztht_8(i+1,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j+1,k)-ztht_8(i  ,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!Dz/Dzeta on k M-level and  Ui,j+1 position
               ski= (Jzpi+Jz)*half ! ui,j+1,k
!
! jy k-1 M-level  U_i,j position
               Jy  = half*(half * (GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i+1,j,k-1))&
                  + half * (GVM%mc_Jy_8(i,j-1,k-1) +GVM%mc_Jy_8(i+1,j-1,k-1)))

! jy k M-level  U_i,j position
               Jyp = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k))&
                      + half * (GVM%mc_Jy_8(i,j-1,k) +GVM%mc_Jy_8(i+1,j-1,k)))
! D(jyU)/Dzeta on k T-level on ui,j position
             C2_8= (Jyp*fdg2_4(i,j,k)-Jy*fdg2_4(i,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))

! jy k-1 M-level  U_i,j+1 position
               Jy  = half*(half * (GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i+1,j,k-1)) &
                  + half * (GVM%mc_Jy_8(i,j+1,k-1) +GVM%mc_Jy_8(i+1,j+1,k-1)))

! jy k M-level  U_i,j position
               Jyp = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k))&
                      + half * (GVM%mc_Jy_8(i,j+1,k) +GVM%mc_Jy_8(i+1,j+1,k)))
! D(jyU)/Dzeta on k-1 T-level on ui,j+1 position
             C= (Jyp*fdg2_4(i,j+1,k)-Jy*fdg2_4(i,j+1,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! D(jyU)/Dzeta on k-1 T-level on Bij
              C2_8= half*(C2_8+C)
! D(jyU)/Dzeta on k M-level  on Bi,j position ! condition zero enbas 
               C2_8  = Ver_wp_8%m(k)*zero +Ver_wm_8%m(k)* C2_8

! D(J_zeta*U)/Dy-(D(jy*U)/Dzeta) on  M-level k and B_ij position
               Bfdg1(i,j,k) = F_coef_8(k)*((ski*fdg2_4(i,j+1,k) -skpi*fdg2_4(i,j,k))* geomh_invDY_8 - C2_8)
!
         Bfdg2(i,j,k) =F_coef_8(k)*((ski*fdg2_4(i,j+1,k) -skpi*fdg2_4(i,j,k))* geomh_invDY_8 - zero)  ! M-level B_ij position 


          enddo
        enddo
!
      do k = 2,Nk-1
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e

               Jzpi =(ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!Dz/Dzeta on k  M-level and  Ui,j position
               skpi= (Jzpi+Jz)*half ! Dz/Dzeta on M-level k on Uij position 
!
               Jzpi =(ztht_8(i+1,j+1,k)-ztht_8(i+1,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j+1,k)-ztht_8(i  ,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!Dz/Dzeta on k M-level and  Ui,j+1 position
               ski= (Jzpi+Jz)*half ! Dz/Dzeta on M-level k on Uij+1 position
!
! jy k-1 M-level  U_i,j position
               Jy  = half*(half * (GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i+1,j,k-1))&
                  + half * (GVM%mc_Jy_8(i,j-1,k-1) +GVM%mc_Jy_8(i+1,j-1,k-1)))

! jy k M-level  U_i,j position
               Jyp = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k))&
                      + half * (GVM%mc_Jy_8(i,j-1,k) +GVM%mc_Jy_8(i+1,j-1,k)))
! D(jyU)/Dzeta on k-1 T-level on ui,j position
             C2_8= (Jyp*fdg2_4(i,j,k)-Jy*fdg2_4(i,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! jy k-1 M-level  B_i,j+1 position
               Jy  = half*(half * (GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i+1,j,k-1)) &
                  + half * (GVM%mc_Jy_8(i,j+1,k-1) +GVM%mc_Jy_8(i+1,j+1,k-1)))

! jy k M-level  B_i,j+1 position
               Jyp = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k))&
                      + half * (GVM%mc_Jy_8(i,j+1,k) +GVM%mc_Jy_8(i+1,j+1,k)))
! D(jyU)/Dzeta on k-1 T-level on ui,j+1 position
             C= (Jyp*fdg2_4(i,j+1,k)-Jy*fdg2_4(i,j+1,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! D(jyU)/Dzeta on k-1 T-level on Bij
              C1_8= half*(C2_8+C)
!
! jy k M-level  U_i,j position
               Jy  = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) &
                   +half * (GVM%mc_Jy_8(i,j-1,k) +GVM%mc_Jy_8(i+1,j-1,k)))

! jy k+1 M-level  U_i,j position
               Jyp = half*(half * (GVM%mc_Jy_8(i,j,k+1) +GVM%mc_Jy_8(i+1,j,k+1))  &
                       +half * (GVM%mc_Jy_8(i,j-1,k+1) +GVM%mc_Jy_8(i+1,j-1,k+1)))
! D(jyU)/Dzeta on k T-level on ui,j position
             C2_8= (Jyp*fdg2_4(i,j,k+1)-Jy*fdg2_4(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! jy k M-level  U_i,j+1 position
               Jy  = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) &
                   +half * (GVM%mc_Jy_8(i,j+1,k) +GVM%mc_Jy_8(i+1,j+1,k)))

! jy k+1 M-level  U_i,j position
               Jyp = half*(half * (GVM%mc_Jy_8(i,j,k+1) +GVM%mc_Jy_8(i+1,j,k+1)) &
                    +   half * (GVM%mc_Jy_8(i,j+1,k+1) +GVM%mc_Jy_8(i+1,j+1,k+1)))
! D(jyU)/Dzeta on k T-level on ui,j+1 position
             C= (Jyp*fdg2_4(i,j+1,k+1)-Jy*fdg2_4(i,j+1,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! D(jyU)/Dzeta on k T-level on Bij
              C2_8= half*(C2_8+C)
! D(jyU)/Dzeta on k M-level  on Bi,j position ! condition zero enbas 
               C2_8  = Ver_wp_8%m(k)*C2_8 +Ver_wm_8%m(k)* C1_8
!
! D(J_zeta*U)/Dy-(D(jy*U)/Dzeta) on  M-level k and B_ij position
               Bfdg1(i,j,k) =  F_coef_8(k)*((ski*fdg2_4(i,j+1,k) -skpi*fdg2_4(i,j,k))* geomh_invDY_8 - C2_8)
!
         Bfdg2(i,j,k) =F_coef_8(k)*((ski*fdg2_4(i,j+1,k) -skpi*fdg2_4(i,j,k))* geomh_invDY_8 - zero)  ! M-level B_ij position 

            enddo
         enddo
         enddo
! Apply divergence
      do k = 1, nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
!inverse(D(z)/D(zeta)) on k M-level and Ui,j position
         xfact= half*(((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i+1,j,k  )-ztht_8(i+1,j,k-1)))+&
                      ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1)))) 
!X-divergence (Dz/Dzeta)^-1*DAfdg/Dx on  k M-level and Ui,j position
         add_v8(i,j,k) = xfact* (Afdg1 (i+1,j,k)-Afdg1 (i,j,k))*geomh_invDXM_8(j)
            enddo
         enddo
      enddo
      do k = 1, nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
!inverse(D(z)/D(zeta)) on k M-level and Ui,j position
         xfact= half*(((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i+1,j,k  )-ztht_8(i+1,j,k-1)))+&
                      ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1))))
!Y-divergence (Dz/Dzeta)^-1*(costheta*Bfdg)/DY on  k M-level and Ui,j position
         bdd_v8(i,j,k) = xfact* (Bfdg1 (i,j,k)*geomh_cyv_8(j)-Bfdg1 (i,j-1,k)*geomh_cyv_8(j-1))*geomh_invDYM_8(j)
            enddo
         enddo
      enddo
!  flux
      do k = 2,Nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
! Dz(jx/Jz A)
! Dz/Dzeta^-1 on M-level k on Afdg_i,j position
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jzm  =(ztht_8(i  ,j,k-1)-ztht_8(i  ,j,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
!Dz(Jx/Jz A) on T-level K-1 and position Ai,j
               C=  half*(GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i-1,j,k)) *Afdg2(i,j,k)/Jz -&
                      half*(GVM%mc_Jx_8(i,j,k-1) +GVM%mc_Jx_8(i-1,j,k-1)) *Afdg2(i,j,k-1)/Jzm 
!               C=  half*(GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i-1,j,k)) *Afdg1(i,j,k)/Jz -&
!                      half*(GVM%mc_Jx_8(i,j,k-1) +GVM%mc_Jx_8(i-1,j,k-1)) *Afdg1(i,j,k-1)/Jzm

               C=  C/(Ver_z_8%m(k)-Ver_z_8%m(k-1)) !Dz(Jx/Jz A) on T-level K-1 and position Ai,j
!
! Dz/Dzeta^-1 on M-level k+1 on Afdg_i,j position
               Jzpi  =(ztht_8(i  ,j,k+1)-ztht_8(i  ,j,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
    
!Dz(Jx/Jz A) on T-level K and position Ai,j
               if (k.eq.NK) then
               C1_8=zero ! lower BD condition
               else
               C1_8=  half*(GVM%mc_Jx_8(i,j,k+1) +GVM%mc_Jx_8(i-1,j,k+1)) *Afdg2(i,j,k+1)/Jzpi -&
                     half* (GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i-1,j,k)) *Afdg2(i,j,k)/Jz
!               C1_8=  half*(GVM%mc_Jx_8(i,j,k+1) +GVM%mc_Jx_8(i-1,j,k+1)) *Afdg1(i,j,k+1)/Jzpi -&
!                     half* (GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i-1,j,k)) *Afdg1(i,j,k)/Jz

               C1_8=  C1_8/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) !Dz(Jx/Jz A) on T-level K and position Ai,j
               endif
!Dz(Jx/Jz A) on M-level K and position Ai,j
               cdd_v81(i,j,k)= Ver_wp_8%m(k)*C1_8+Ver_wm_8%m(k)* C  !Dz(Jx/Jz A) on M-level K and position Ai,j
!               
! D(z)/D(zeta) on k M_level on Vi,j position 
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jzpi  =(ztht_8(i  ,j+1,k)-ztht_8(i  ,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz= half*(Jz+Jzpi)
! D(z)/D(zeta) on k M_level on Vi+1,j position
               Jzm   =(ztht_8(i+1  ,j,k)-ztht_8(i+1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jzpi  =(ztht_8(i+1  ,j+1,k)-ztht_8(i+1  ,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jzpi= half*(Jzm+Jzpi)
! D(z)/D(zeta) on k M_level on Bi,j position
               Jzpi= half*(Jz+Jzpi)  ! Jz on M_level k  and position Bij
! D(z)/D(zeta) on k-1 M_level on Vi,j position
               Jz   =(ztht_8(i  ,j,k-1)-ztht_8(i  ,j,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jzm  =(ztht_8(i  ,j+1,k-1)-ztht_8(i  ,j+1,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jz= half*(Jz+Jzm)
! D(z)/D(zeta) on k-1 M_level on Vi+1,j position
               Jzmpi   =(ztht_8(i+1  ,j,k-1)-ztht_8(i+1  ,j,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jzm  =(ztht_8(i+1  ,j+1,k-1)-ztht_8(i+1  ,j+1,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jzm= half*(Jzmpi+Jzm)
! D(z)/D(zeta) on k-1 M_level on Bi,j position
               Jzm= half*(Jz+Jzm)  ! Jz on M_level k-1  and position Bij

               C=  half*(GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) *Bfdg2(i,j,k)/jzpi -&
                      half*(GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i+1,j,k-1)) *Bfdg2(i,j,k-1)/Jzm
!               C=  half*(GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) *Bfdg1(i,j,k)/jzpi -&
!                      half*(GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i+1,j,k-1)) *Bfdg1(i,j,k-1)/Jzm

!               endif
               C=  C/(Ver_z_8%m(k)-Ver_z_8%m(k-1)) !Dz(Jy/Jz B) on T-level K-1 and position Bi,j
! D(z)/D(zeta) on k+1 M_level on Vi,j position
               Jz   =(ztht_8(i  ,j,k+1)-ztht_8(i  ,j,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jzm  =(ztht_8(i  ,j+1,k+1)-ztht_8(i  ,j+1,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jz= half*(Jz+Jzm)
! D(z)/D(zeta) on k+1 M_level on Vi+1,j position
               Jzmpi   =(ztht_8(i+1  ,j,k+1)-ztht_8(i+1  ,j,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jzm  =(ztht_8(i+1  ,j+1,k+1)-ztht_8(i+1  ,j+1,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jzm= half*(Jzmpi+Jzm)
!D(z)/D(zeta) on k+1 M_level on  Bi,j position
               Jz= half*(Jz+Jzm) !Jz on M_level k  and position Bij
!Dz(Jy/Jz B) on k T_level on Bi,j position
               if (k.eq.NK) then 
               C2_8=zero ! lower BD condition        
               else
               C2_8=  half*(GVM%mc_Jy_8(i,j,k+1) +GVM%mc_Jy_8(i+1,j,k+1)) *Bfdg2(i,j,k+1)/jz -&
                      half*(GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) *Bfdg2(i,j,k)/Jzpi
!               C2_8=  half*(GVM%mc_Jy_8(i,j,k+1) +GVM%mc_Jy_8(i+1,j,k+1)) *Bfdg1(i,j,k+1)/jz -&
!                      half*(GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) *Bfdg1(i,j,k)/Jzpi

!
               C2_8=  C2_8/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) !Dz(Jy/Jz B) on T-level K and position Bi,j
               endif
!Dz(Jy/Jz B) on k M_level on Bi,j position
               C2_8= Ver_wp_8%m(k)*C2_8+Ver_wm_8%m(k)* C  !Dz(Jy/Jz B) on M-level K and position Bi,j
               cdd_v82(i,j,k)=C2_8

           enddo
         enddo
      enddo
          do k=1,NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
!inverse(D(z)/D(zeta)) on k M-level and Ui,j position
            xfact= half*(((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i+1,j,k  )-ztht_8(i+1,j,k-1)))+&
                      ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1))))
! xfact* cdd on k M-level and Ui,j position
            cdd_v8(i,j,k)= xfact *(half*(cdd_v81(i,j,k)+ cdd_v81(i+1,j,k))+ half*(cdd_v82(i,j,k)+ cdd_v82(i,j-1,k))/geomh_cy_8(j))
! boundary condition
      if (k==NK) then
! extrapolation lineaire
       cdd_v8(i,j,k)= (one-(ver_z_8%m(Nk)-ver_z_8%m(Nk-1))/(ver_z_8%m(Nk+1)-ver_z_8%m(Nk-1)))*&
                        cdd_v8(i,j,Nk-1)
!extrapolation cubique
!        cdd_v8(i,j,k)= cdd_v8(i,j,k-2)+(GVM%zmom_8(i  ,j,k) -GVM%zmom_8(i  ,j,k-2))/(GVM%zmom_8(i  ,j,k-1) -GVM%zmom_8(i  ,j,k-2))&
!        *(cdd_v8(i,j,k-1)-cdd_v8(i,j,k-2))
      endif
! solution  on k M-level and Ui,j position
! Ajout
            F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-cdd_v8(i,j,k))
! fin Ajout

             enddo
             enddo
          enddo

! termes implicit restant        
! stencilV
      allocate(stencil_V(l_minx:l_maxx, l_miny:l_maxy,3,Nk))
      stencil_V=zero
        do k=1,NK
         km=max(k-1,1)
         kp=min(NK,k+1)
         do j=1+pil_s, l_nj-pil_n
          do i=1+pil_w, l_ni-pil_e
! Dz/Dzeta^-1 at T-levels K and K-1  and Ui,j position
            Jz= half*(((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k))) +&
              ((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i+1,j,k+1  )-GVM%zmom_8(i+1,j,k))) )
            Jzm= half*(((Ver_z_8%m(k)-Ver_z_8%m(k-1))/(GVM%zmom_8(i,j,k  )-GVM%zmom_8(i,j,k-1))) +&
              ((Ver_z_8%m(k)-Ver_z_8%m(k-1))/(GVM%zmom_8(i+1,j,k  )-GVM%zmom_8(i+1,j,k-1))) )
! JYT ate T_levels K and K-1 and Ui,j position
           C1_8= half*( half*(GVM%mc_Jyt_8(i,j,k)+GVM%mc_Jyt_8(i,j-1,k))+&
                  half*(GVM%mc_Jyt_8(i+1,j,k)+GVM%mc_Jyt_8(i+1,j-1,k)))   
           C2_8=zero  
           if (k.ne.1) C2_8= half*( half*(GVM%mc_Jyt_8(i,j,km)+GVM%mc_Jyt_8(i,j-1,km))+&
                  half*(GVM%mc_Jyt_8(i+1,j,km)+GVM%mc_Jyt_8(i+1,j-1,km)))
!JY ate  M_levels K+1 and K and k-1
            C3_8=zero
            ski=zero
            if(k.ne.NK) C3_8= half*( half*(GVM%mc_Jy_8(i,j,kp)+GVM%mc_Jy_8(i,j-1,kp))+&
                  half*(GVM%mc_Jy_8(i+1,j,kp)+GVM%mc_Jy_8(i+1,j-1,kp)))
            C= half*( half*(GVM%mc_Jy_8(i,j,k)+GVM%mc_Jy_8(i,j-1,k))+&
                  half*(GVM%mc_Jy_8(i+1,j,k)+GVM%mc_Jy_8(i+1,j-1,k)))
            if (k.ne.1)  ski= half*( half*(GVM%mc_Jy_8(i,j,km)+GVM%mc_Jy_8(i,j-1,km))+&
                  half*(GVM%mc_Jy_8(i+1,j,km)+GVM%mc_Jy_8(i+1,j-1,km)))
!
          if (k.ne.NK) &
          stencil_V(i,j,3,k)= one* F_coef_8(k)*Jz*C1_8*&
          C3_8 /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))+&
          one* F_coef_8(k)*Jz*GVM%mc_JxT_8(i,j,k)*&
          GVM%mc_Jx_8(i,j,kp) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))
!
          if (k .ne.1) &
          stencil_V(i,j,2,k)= one* F_coef_8(km)*Jzm*C2_8*&
          ski /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))+&
        one* F_coef_8(km)*Jzm*GVM%mc_Jxt_8(i,j,km)*&
          GVM%mc_Jx_8(i,j,km) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))
!
          Jxx=zero
          if (k.ne.1) Jxx= GVM%mc_Jxt_8(i,j,km)
          stencil_V(i,j,1,k)= -one *F_coef_8(k)*Jz*C1_8*&
              C/((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))-&
          one* F_coef_8(k)*Jz*GVM%mc_Jxt_8(i,j,k)*&
          GVM%mc_Jx_8(i,j,k) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))) &
          -one *F_coef_8(km)*Jzm*C2_8*&
              C/((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))-&
          one* F_coef_8(km)*Jzm*Jxx*&
          GVM%mc_Jx_8(i,j,k) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))
!
            Jz= (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1))

            stencil_V(i,j,3,k )=Jz*  stencil_V(i,j,3,k )
            stencil_V(i,j,2,k )=Jz*  stencil_V(i,j,2,k )
            stencil_V(i,j,1,k )=Jz*  stencil_V(i,j,1,k )
!
!            if ((k.ne.1).and.(k.ne.NK))   F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*(stencil_V(i,j,2,k) *fdg2_4(i,j,k-1)+&
!                    stencil_V(i,j,1,k) *fdg2_4(i,j,k)+stencil_V(i,j,3,k) *fdg2_4(i,j,k+1))
!           if (k.eq.1)   F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*(stencil_V(i,j,3,k) *fdg2_4(i,j,k+1)+&
!                    stencil_V(i,j,1,k) *fdg2_4(i,j,k))
!             if (k==NK) C= (one-(ver_z_8%m(Nk)-ver_z_8%m(Nk-1))/(ver_z_8%m(Nk+1)-ver_z_8%m(Nk-1)))*&
!                 (-stencil_V(i,j,1,Nk-1)* fdg2_4(i,j,Nk-1)-  stencil_V(i,j,2,Nk-1)* fdg2_4(i,j,Nk-2)&
!                        -stencil_V(i,j,3,Nk-1)* fdg2_4(i,j,Nk))
!           if (k.eq.NK)   F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*(stencil_V(i,j,2,k) *fdg2_4(i,j,k-1)+&
!                    stencil_V(i,j,1,k) *fdg2_4(i,j,k))
!             if (k.eq.NK) F_sol1(i,j,k)= F_sol1(i,j,k) - Cstv_dt_8* C

     if (k==NK) then
          stencil_V(i,j,2,k )=  (one-(ver_z_8%m(Nk)-ver_z_8%m(Nk-1))/(ver_z_8%m(Nk+1)-ver_z_8%m(Nk-1)))* stencil_V(i,j,2,k-1 )
          stencil_V(i,j,1,k )= (one-(ver_z_8%m(Nk)-ver_z_8%m(Nk-1))/(ver_z_8%m(Nk+1)-ver_z_8%m(Nk-1)))*(stencil_V(i,j,1,k-1 )+&
                 stencil_V(i,j,3,k-1 ))
      endif
             enddo
         enddo
     enddo

      k=1
        do j=1+pil_s, l_nj-pil_n
           do i=1+pil_w, l_ni-pil_e
            d(i,j,k)= -Cstv_dt_8*stencil_V(i,j,3,k )
            b(i,j,k)=one-Cstv_dt_8*stencil_V(i,j,1,k)

            enddo
         enddo
        do k=2,Nk-1
           do j=1+pil_s, l_nj-pil_n
             do i=1+pil_w, l_ni-pil_e
             a(i,j,k)= -Cstv_dt_8*stencil_V(i,j,2,k)
             b(i,j,k)=one-Cstv_dt_8*stencil_V(i,j,1,k)
             d(i,j,k)=-Cstv_dt_8*stencil_V(i,j,3,k)

             enddo
           enddo
        enddo
        k=Nk
           do j=1+pil_s, l_nj-pil_n
             do i=1+pil_w, l_ni-pil_e
             a(i,j,k)= -Cstv_dt_8*stencil_V(i,j,2,k)
             b(i,j,k)= one-Cstv_dt_8*stencil_V(i,j,1,k)

             enddo
           enddo
!
         deallocate (stencil_V)

      do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
              do k = 2 , Nk
                W = a(i,j,k) / b(i,j,k - 1)
                b(i,j,k) = b(i,j,k) - W * d(i,j,k - 1)
                F_sol1(i,j,k) = F_sol1(i,j,k) - W * F_sol1(i,j,k- 1)
              enddo
                 F_sol1(i,j,Nk) = F_sol1(i,j,Nk) / b(i,j,Nk)
                  do  k = Nk-1, 1, -1
                  F_sol1(i,j,k) = (F_sol1(i,j,k) - d(i,j,k) * F_sol1(i,j,k + 1)) / b(i,j,k)
                  enddo
             enddo
          enddo
   enddo
      ! Hybrid diffusion if hzd_hyb_nk >0
      if(hzd_hyb_nk > 0) then
         do k = nk-hzd_hyb_nk+1, nk
            do j=1+pil_s-1, l_nj-pil_n+1
               do i=1+pil_w-1, l_ni-pil_e+1
                  F_sol1(i,j,k) = fdg2_4(i,j,k )
               enddo
            enddo
         enddo
      endif

      return
      end

