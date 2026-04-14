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
! Ajout
      subroutine  hzd_v1_alh ( F_Sol1,HzdlnR,Minx, Maxx, Miny, Maxy,Nk,Niter)
! fin Ajout

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
!
      use, intrinsic :: iso_fortran_env
      implicit none
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, NK, Niter
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent (inout) :: F_Sol1
      real(kind=REAL64) deno


!author
!       Abdessamad Qaddouri -  2018
!
!revision
! v5.0 - Qaddouri A.       - initial version


      integer j,i,k,halox,haloy
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)
      integer  km, kp
      real(kind=REAL64)   Afdg1(l_minx:l_maxx, l_miny:l_maxy,Nk),Afdg2(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64)   Bfdg1(l_minx:l_maxx, l_miny:l_maxy,Nk),Bfdg2(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real   fdg2_4(l_minx:l_maxx, l_miny:l_maxy,Nk+1)

       real(kind=REAL64) Jzpi,Jz,Jzm,Jzmpi
       real(kind=REAL64) C1_8,C2_8,C,ski,skpi,C3_8
       real(kind=REAL64) bdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
       real(kind=REAL64) add_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
       real(kind=REAL64) cdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
       integer iter
       real(kind=REAL64) F_coef_8(1:NK) ,Jx, Jxp
       real(kind=REAL64) cdd_v81(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
       real(kind=REAL64) cdd_v82(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
       real  HzdlnR

!
!     ---------------------Vi,j--------Afdgi,j ----------------------------------
!
!     ---------------------Bfdgi,j---------Ui,j--------------------------------
!
!     ---------------------------------------------------------------
      real(kind=REAL64)  ztht_8(l_minx:l_maxx, l_miny:l_maxy,0:Nk+1)
      real(kind=REAL64) crit_coef,base_coefT
      real(kind=REAL64), dimension (:,:,:,:), allocatable :: stencil_V
      real(kind=REAL64) cwest,ceast
      real(kind=REAL64) a(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) b(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) d(l_minx:l_maxx, l_miny:l_maxy,Nk),W,Jxx,Jyy

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
      k=1
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
!Dz/Dzeta M-level k on Vi,j position
               Jzpi =(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1)) !Dz/Dzeta M-level k
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1)) !Dz/Dzeta M-level k
               ski= (Jzpi+Jz)*half                                                             !Dz/Dzeta M-level k
!Dz/Dzeta M-level k on Vi+1,j position
               Jzpi =(ztht_8(i+1,j+1,k)-ztht_8(i+1,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1)) !Dz/Dzeta M-level k
               Jz   =(ztht_8(i+1  ,j,k)-ztht_8(i+1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1)) !Dz/Dzeta M-level k
               skpi= (Jzpi+Jz)*half                                                                !Dz/Dzeta M-level k  
! Jx on M K+1 level and Vi,j 
               Jxp= half*(half*(GVM%mc_Jx_8(i,j,k+1)+GVM%mc_Jx_8(i,j+1,k+1)) &
                  + half*(GVM%mc_Jx_8(i-1,j,k+1)+GVM%mc_Jx_8(i-1,j+1,k+1)))
! Jx on M K level and Vi,j 
               Jx= half*(half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)) &
                  + half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i-1,j+1,k)))

! (D(jx*V)/Dzeta)  on T-level k and Vij position
               C= (Jxp*fdg2_4(i,j,k+1)-Jx*fdg2_4(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))         
! Jx on M K+1 level and Vi+1,j 
               Jxp= half*(half*(GVM%mc_Jx_8(i+1,j,k+1)+GVM%mc_Jx_8(i+1,j+1,k+1)) &
                  + half*(GVM%mc_Jx_8(i,j,k+1)+GVM%mc_Jx_8(i,j+1,k+1)))
! Jx on M K level and Vi,j 
               Jx= half*(half*(GVM%mc_Jx_8(i+1,j,k)+GVM%mc_Jx_8(i+1,j+1,k)) &
                  + half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)))
! (D(jx*V)/Dzeta)  on T-level k and Vi+1j position
               C2_8= (Jxp*fdg2_4(i+1,j,k+1)-Jx*fdg2_4(i+1,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))     
! (D(jx*V)/Dzeta)  on T-level k and Aij position
               C2_8= half*( C+ C2_8)
! (D(jx*V)/Dzeta)  on M-level k and Aij position! zero upper condition
               C2_8= Ver_wp_8%m(k)*C2_8 + Ver_wm_8%m(k) *zero                    
! D(J_zeta*V)/Dx-(D(jx*V)/Dzeta) on  M-level k and phi,j position
           Afdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i+1,j,k) -ski*fdg2_4(i,j,k))*geomh_invDXv_8(j)- C2_8 )
           Afdg2(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i+1,j,k) -ski*fdg2_4(i,j,k))*geomh_invDXv_8(j)- zero )

          enddo
        enddo
!
      k= NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
!Dz/Dzeta M-level k on Vi,j position
               Jzpi =(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               ski= (Jzpi+Jz)*half
!!Dz/Dzeta M-level k on Vi+1,j position
               Jzpi =(ztht_8(i+1,j+1,k)-ztht_8(i+1,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i+1  ,j,k)-ztht_8(i+1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               skpi= (Jzpi+Jz)*half
! Jx on M K level and Vi,j 
               Jxp= half*(half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)) &
                   +half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i-1,j+1,k)))
! Jx on M K-1 level and Vi,j 
               Jx= half*(half*(GVM%mc_Jx_8(i,j,k-1)+GVM%mc_Jx_8(i,j+1,k-1))&
                  + half*(GVM%mc_Jx_8(i-1,j,k-1)+GVM%mc_Jx_8(i-1,j+1,k-1)))

! (D(jx*V)/Dzeta)  on T-level k and Vij position
               C= (Jxp*fdg2_4(i,j,k)-Jx*fdg2_4(i,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! Jx on M K level and Vi+1,j 
               Jxp= half*(half*(GVM%mc_Jx_8(i+1,j,k)+GVM%mc_Jx_8(i+1,j+1,k)) &
                  +  half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)))
! Jx on M K-1 level and Vi,j 
               Jx= half*(half*(GVM%mc_Jx_8(i+1,j,k-1)+GVM%mc_Jx_8(i+1,j+1,k-1)) &
                  + half*(GVM%mc_Jx_8(i,j,k-1)+GVM%mc_Jx_8(i,j+1,k-1)))
! (D(jx*V)/Dzeta)  on T-level k and Vi+1j position
               C2_8= (Jxp*fdg2_4(i+1,j,k)-Jx*fdg2_4(i+1,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! (D(jx*V)/Dzeta)  on T-level k-1 and Aij position
               C2_8= half*( C+ C2_8)
! (D(jx*V)/Dzeta)  on M-level k and Aij position! zero lower condition
               C2_8= Ver_wp_8%m(k)*zero + Ver_wm_8%m(k) * C2_8          
! D(J_zeta*V)/Dx-(D(jx*V)/Dzeta) on  M-level k and phi,j position
              Afdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i+1,j,k) -ski*fdg2_4(i,j,k)) *geomh_invDXv_8(j)-C2_8)
              Afdg2(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i+1,j,k) -ski*fdg2_4(i,j,k)) *geomh_invDXv_8(j)-zero)

          enddo
        enddo
!       
      do k = 2,Nk-1
         km=min(k-1,1)
         kp = max(k+1,NK)
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
!Dz/Dzeta M-level k on Vi,j position
               Jzpi =(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               ski= (Jzpi+Jz)*half
!Dz/Dzeta M-level k on Vi+1,j position
               Jzpi =(ztht_8(i+1,j+1,k)-ztht_8(i+1,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i+1  ,j,k)-ztht_8(i+1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               skpi= (Jzpi+Jz)*half
! Jx on M K level and Vi,j 
               Jxp= half*(half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)) &
                   +half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i-1,j+1,k)))
! Jx on M K-1 level and Vi,j 
               Jx= half*(half*(GVM%mc_Jx_8(i,j,k-1)+GVM%mc_Jx_8(i,j+1,k-1))&
                  + half*(GVM%mc_Jx_8(i-1,j,k-1)+GVM%mc_Jx_8(i-1,j+1,k-1)))

! (D(jx*V)/Dzeta)  on T-level k and Vij position
               C= (Jxp*fdg2_4(i,j,k)-Jx*fdg2_4(i,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! Jx on M K level and Vi+1,j 
               Jxp= half*(half*(GVM%mc_Jx_8(i+1,j,k)+GVM%mc_Jx_8(i+1,j+1,k)) &
                  +  half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)))
! Jx on M K-1 level and Vi,j 
               Jx= half*(half*(GVM%mc_Jx_8(i+1,j,k-1)+GVM%mc_Jx_8(i+1,j+1,k-1)) &
                  + half*(GVM%mc_Jx_8(i,j,k-1)+GVM%mc_Jx_8(i,j+1,k-1)))
! (D(jx*V)/Dzeta)  on T-level k and Vi+1j position
               C2_8= (Jxp*fdg2_4(i+1,j,k)-Jx*fdg2_4(i+1,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
! (D(jx*V)/Dzeta)  on T-level k-1 and Aij position
               C1_8= half*( C+ C2_8)
! Jx on M K+1 level and Vi,j 
               Jxp= half*(half*(GVM%mc_Jx_8(i,j,k+1)+GVM%mc_Jx_8(i,j+1,k+1)) &
                  + half*(GVM%mc_Jx_8(i-1,j,k+1)+GVM%mc_Jx_8(i-1,j+1,k+1)))
! Jx on M K level and Vi,j 
               Jx= half*(half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)) &
                  + half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i-1,j+1,k)))
! (D(jx*V)/Dzeta)  on T-level k and Vij position
               C= (Jxp*fdg2_4(i,j,k+1)-Jx*fdg2_4(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! Jx on M K+1 level and Vi+1,j 
               Jxp= half*(half*(GVM%mc_Jx_8(i+1,j,k+1)+GVM%mc_Jx_8(i+1,j+1,k+1)) &
                  + half*(GVM%mc_Jx_8(i,j,k+1)+GVM%mc_Jx_8(i,j+1,k+1)))
! Jx on M K level and Vi,j 
               Jx= half*(half*(GVM%mc_Jx_8(i+1,j,k)+GVM%mc_Jx_8(i+1,j+1,k)) &
                  + half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)))
! (D(jx*V)/Dzeta)  on T-level k and Vi+1j position
               C2_8= (Jxp*fdg2_4(i+1,j,k+1)-Jx*fdg2_4(i+1,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
! (D(jx*V)/Dzeta)  on T-level k and Aij position
               C2_8= half*( C+ C2_8)

! (D(jx*V)/Dzeta)  on M-level k and Vij position
               C2_8  = Ver_wp_8%m(k)*C2_8+Ver_wm_8%m(k)* C1_8
! D(J_zeta*V)/Dx-(D(jx*V)/Dzeta) on  M-level k and phi,j position
              Afdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i+1,j,k) -ski*fdg2_4(i,j,k))*geomh_invDXv_8(j)-C2_8)
              Afdg2(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i+1,j,k) -ski*fdg2_4(i,j,k))*geomh_invDXv_8(j)-zero)
            enddo
         enddo
         enddo
! Gradient component Along Y
         k=1
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
!Dz/Dzeta M-level k on Vi,j position
               Jzpi =(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               skpi= (Jzpi+Jz)*half
!Dz/Dzeta M-level k on Vi,j-1 position
               Jzpi =(ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j-1,k)-ztht_8(i  ,j-1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               ski = (Jzpi+Jz)*half
! (D(jy*V)/Dzeta)  on T-level k and Vij position and Vi,j-1 position >>>> Bij position
               C2_8= (GVM%mc_Jy_8(i,j,k+1)*fdg2_4(i,j,k+1)-GVM%mc_Jy_8(i,j,k)*fdg2_4(i,j,k))&
                            /(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C3_8= (GVM%mc_Jy_8(i,j-1,k+1)*fdg2_4(i,j-1,k+1)-GVM%mc_Jy_8(i,j-1,k)*fdg2_4(i,j-1,k))&
                            /(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C= half*(C2_8+C3_8)
! (D(jy*V)/Dzeta)  on M-level k and Bij position. upper BD condition
               C2_8= Ver_wp_8%m(k)*C + Ver_wm_8%m(k)* zero
! D(J_zeta*V)/Dy-(D(jy*V)/Dzeta) on  M-level k and Bi,j=phi,j position
 !          Bfdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k) -ski*fdg2_4(i,j-1,k))* geomh_invDYM_8(j) - C2_8)
            Bfdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k) -ski*fdg2_4(i,j-1,k))* geomh_invDY_8 - C2_8)
            Bfdg2(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k) -ski*fdg2_4(i,j-1,k))* geomh_invDY_8 - zero)
          enddo
        enddo
!
      k= NK
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
!Dz/Dzeta M-level k on Vi,j position
               Jzpi =(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               skpi= (Jzpi+Jz)*half
!Dz/Dzeta M-level k on Vi,j-1 position
               Jzpi =(ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j-1,k)-ztht_8(i  ,j-1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               ski = (Jzpi+Jz)*half
! (D(jy*V)/Dzeta)  on T-level k-1 and Bij position
               C2_8= (GVM%mc_Jy_8(i,j,k)*fdg2_4(i,j,k)-GVM%mc_Jy_8(i,j,k-1)*fdg2_4(i,j,k-1))&
                            /(Ver_z_8%m(k)-Ver_z_8%m(k-1))
               C3_8= (GVM%mc_Jy_8(i,j-1,k)*fdg2_4(i,j-1,k)-GVM%mc_Jy_8(i,j-1,k-1)*fdg2_4(i,j-1,k-1))&
                            /(Ver_z_8%m(k)-Ver_z_8%m(k-1))
               C= half*(C2_8+C3_8)
! (D(jy*V)/Dzeta)  on M-level k and Bij position. lower BD condition
               C2_8= Ver_wp_8%m(k)*zero +Ver_wm_8%m(k)*C
! D(J_zeta*V)/Dy-(D(jy*V)/Dzeta) on  M-level k and Bi,j position
           Bfdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k) -ski*fdg2_4(i,j-1,k))* geomh_invDY_8 - C2_8)
           Bfdg2(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k) -ski*fdg2_4(i,j-1,k))* geomh_invDY_8 - zero)

!
          enddo
        enddo


      do k = 2,Nk-1
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
!Dz/Dzeta M-level k on Vi,j position
               Jzpi =(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               skpi= (Jzpi+Jz)*half
!
!Dz/Dzeta M-level k on Vi,j-1 position
               Jzpi =(ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j-1,k)-ztht_8(i  ,j-1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               ski = (Jzpi+Jz)*half
! (D(jy*V)/Dzeta)  on T-level k-1 and Bij position
               C2_8= (GVM%mc_Jy_8(i,j,k)*fdg2_4(i,j,k)-GVM%mc_Jy_8(i,j,k-1)*fdg2_4(i,j,k-1))&
                            /(Ver_z_8%m(k)-Ver_z_8%m(k-1))
               C3_8= (GVM%mc_Jy_8(i,j-1,k)*fdg2_4(i,j-1,k)-GVM%mc_Jy_8(i,j-1,k-1)*fdg2_4(i,j-1,k-1))&
                            /(Ver_z_8%m(k)-Ver_z_8%m(k-1))
               C= half*(C2_8+C3_8)

! (D(jy*V)/Dzeta)  on T-level k and Bij position
               C2_8= (GVM%mc_Jy_8(i,j,k+1)*fdg2_4(i,j,k+1)-GVM%mc_Jy_8(i,j,k)*fdg2_4(i,j,k))&
                            /(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C3_8= (GVM%mc_Jy_8(i,j-1,k+1)*fdg2_4(i,j-1,k+1)-GVM%mc_Jy_8(i,j-1,k)*fdg2_4(i,j-1,k))&
                            /(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C3_8= half*(C2_8+C3_8)
! (D(jy*V)/Dzeta)  on M-level k and Bij position
               C2_8= Ver_wp_8%m(k)*C3_8+Ver_wm_8%m(k)*C
! D(J_zeta*V)/Dy-(D(jy*V)/Dzeta) on  M-level k and Bi,j position
           Bfdg1(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k) -ski*fdg2_4(i,j-1,k))* geomh_invDY_8- C2_8)
           Bfdg2(i,j,k) =F_coef_8(k)*((skpi*fdg2_4(i,j,k) -ski*fdg2_4(i,j-1,k))* geomh_invDY_8- zero)
!
            enddo
         enddo
         enddo
!

! Apply divergence
      do k = 1, nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
! Dzeta/Dz^-1  on  k M-level and Vi,j position
         C    =  ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j+1,k  )-ztht_8(i,j+1,k-1)))
         deno =  ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1)))
         deno =  half*(deno+C)               ! Dzeta/Dz 
!X-divergence (Dz/Dzeta)^-1*DAfdg/Dx on  k M-level and Vi,j position
         add_v8(i,j,k) =    deno*  (Afdg1 (i,j,k)-Afdg1 (i-1,j,k))*geomh_invDXv_8(j)
            enddo
         enddo
      enddo
      do k = 1, nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
! Dzeta/Dz^-1  on  k M-level and Vi,j position
         C    =  ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j+1,k  )-ztht_8(i,j+1,k-1)))
         deno =  ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1)))
         deno =  half*(deno+C)               ! Dzeta/Dz
!Y-divergence (Dz/Dzeta)^-1*DBfdg/Dy on  k M-level and Vi,j position
!         bdd_v8(i,j,k) =    deno* (Bfdg1 (i,j+1,k)*geomh_cy_8(j+1)-Bfdg1 (i,j,k)*geomh_cy_8(j))*geomh_invDY_8*geomh_invcy_8(j)
          bdd_v8(i,j,k) =    deno* (Bfdg1 (i,j+1,k)*geomh_cy_8(j+1)-Bfdg1 (i,j,k)*geomh_cy_8(j))*geomh_invDYM_8 (j)
            enddo
         enddo
      enddo
!  flux
! start at pil_w and finish at l_nj-pil_n+1 ! important pour le calcul de sol a la fin
      do k = 2,Nk
         do j=1+pil_s, l_nj-pil_n+1
            do i=pil_w, l_ni-pil_e
! Dz(jx/Jz A)
!Dz/Dzeta M-level k on Vi,j position
               Jz   = (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))
               Jzpi  =(Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j+1,k)-ztht_8(i  ,j+1,k-1))
               Jz= half*(Jz+Jzpi)
!Dz/Dzeta M-level k on Vi+1,j position
               Jzm   =(Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i+1  ,j,k)-ztht_8(i+1  ,j,k-1))
               Jzpi  =(Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i+1  ,j+1,k)-ztht_8(i+1  ,j+1,k-1))
               Jzpi= half*(Jzm+Jzpi)
!Dz/Dzeta M-level k on Ai,j position
               Jz = half*(Jz + Jzpi) ! Jz on M-level k on Aij point
! Dz/Dzeta M-level k-1 on Vi,j position
               Jzm   =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i  ,j,k-1)-ztht_8(i  ,j,k-2))
               Jzpi  =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i  ,j+1,k-1)-ztht_8(i  ,j+1,k-2))
               Jzm= half*(Jzm+Jzpi)
!Dz/Dzeta M-level k-1 on Vi+1,j position
               Jzmpi =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i+1  ,j,k-1)-ztht_8(i+1  ,j,k-2))
               Jzpi  =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i+1  ,j+1,k-1)-ztht_8(i+1  ,j+1,k-2))
               Jzpi= half*(Jzmpi+Jzpi)
!Dz/Dzeta M-level k-1 on Ai,j position
               Jzm = half*(Jzm + Jzpi) ! Jz on M-level k-1 on Aij point

               C=  half*(GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i,j+1,k)) *Afdg2(i,j,k)*Jz -&
                      half*(GVM%mc_Jx_8(i,j,k-1) +GVM%mc_Jx_8(i,j+1,k-1)) *Afdg2(i,j,k-1)*Jzm
!               C=  half*(GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i,j+1,k)) *Afdg1(i,j,k)*Jz -&
!                      half*(GVM%mc_Jx_8(i,j,k-1) +GVM%mc_Jx_8(i,j+1,k-1)) *Afdg1(i,j,k-1)*Jzm

!   Dz(Jx/Jz A) on T-level K-1 and position Ai,j
               C=  C/(Ver_z_8%m(k)-Ver_z_8%m(k-1)) !Dz(Jx/Jz A) on T-level K-1 and position Ai,j
                cdd_v81(i,j,k-1)=C
!Dz(Jy/Jz B)
!Jz on M-level k on Bij point
               Jz   =(Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1)) 
!Jz on M-level k-1 on Bij point
               Jzm   =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i  ,j,k-1)-ztht_8(i  ,j,k-2))
!   Dz(Jy/Jz B) on T-level K-1 and position Bi,j
               C=  half*(GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i,j-1,k)) *Bfdg2(i,j,k)*Jz -&
                      half*(GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i,j-1,k-1)) *Bfdg2(i,j,k-1)*Jzm

!               C=  half*(GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i,j-1,k)) *Bfdg1(i,j,k)*Jz -&
!                      half*(GVM%mc_Jy_8(i,j,k-1) +GVM%mc_Jy_8(i,j-1,k-1)) *Bfdg1(i,j,k-1)*Jzm
               C=  C/(Ver_z_8%m(k)-Ver_z_8%m(k-1)) 
               cdd_v82(i,j,k-1)=C
!
           enddo
         enddo
        enddo
          do k=1,NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
         if (k==1) then 
             cdd_v8(i,j,k) = zero 
         else
! Dz/Dzeta M-level k on Vi,j position
            Jz= (Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i,j,k)-ztht_8(i,j,k-1))
            Jzpi= (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))
            Jz = half*(Jz+Jzpi)
! Dz/Dzeta^-1* Dz(Jx/Jz A+ Jy/Jz B) on M-level K and position Bi,j
            C1_8= half *( (Ver_wp_8%m(k)*cdd_v81(i,j,k)+Ver_wm_8%m(k)* cdd_v81(i,j,k-1))+&  
                   (Ver_wp_8%m(k)*cdd_v81(i+1,j,k)+Ver_wm_8%m(k)* cdd_v81(i+1,j,k-1)) )
            C2_8= half*((Ver_wp_8%m(k)*cdd_v82(i,j,k)+Ver_wm_8%m(k)* cdd_v82(i,j,k-1)) +&
                   (Ver_wp_8%m(k)*cdd_v82(i,j+1,k)+Ver_wm_8%m(k)* cdd_v82(i,j+1,k-1)))
            C= (C1_8+C2_8/geomh_cyv_8(j)) *Jz
            cdd_v8(i,j,k)= C
! boundary condition
      if (k==NK) then
       cdd_v8(i,j,k)= (one-(ver_z_8%m(Nk)-ver_z_8%m(Nk-1))/(ver_z_8%m(Nk+1)-ver_z_8%m(Nk-1)))*&
                        cdd_v8(i,j,Nk-1)
      endif
      endif
! Ajout
            F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-cdd_v8(i,j,k))
! fin Ajout
             enddo
             enddo
          enddo
!
! termes implicit restant        
! stencilV
      allocate(stencil_V(l_minx:l_maxx, l_miny:l_maxy,3,Nk))
      stencil_V=zero
        do k=1,NK
           km = max(k-1,1)
           kp = min(k+1, Nk)
         do j=1+pil_s, l_nj-pil_n
          do i=1+pil_w, l_ni-pil_e
! Dz/Dzeta^-1 at T-levels K and K-1  and Ui,j position
            Jz= half*(((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k))) +&
              ((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j+1,k+1  )-GVM%zmom_8(i,j+1,k))) )
            Jzm= half*(((Ver_z_8%m(k)-Ver_z_8%m(k-1))/(GVM%zmom_8(i,j,k  )-GVM%zmom_8(i,j,k-1))) +&
              ((Ver_z_8%m(k)-Ver_z_8%m(k-1))/(GVM%zmom_8(i,j+1,k  )-GVM%zmom_8(i,j+1,k-1))) )
! JXT ate T_levels K and K-1 and Ui,j position
            C1_8= half*( half*(GVM%mc_Jxt_8(i,j,k)+GVM%mc_Jxt_8(i-1,j,k))+&
                  half*(GVM%mc_Jxt_8(i,j+1,k)+GVM%mc_Jxt_8(i-1,j+1,k)))
            C2_8=zero
            if(k.ne.1)  C2_8= half*( half*(GVM%mc_Jxt_8(i,j,km)+GVM%mc_Jxt_8(i-1,j,km))+&
                  half*(GVM%mc_Jxt_8(i,j+1,km)+GVM%mc_Jxt_8(i-1,j+1,km)))
!JX ate  M_levels K+1 and K and k-1
            C3_8=zero
            if(k.ne.NK)  C3_8= half*( half*(GVM%mc_Jx_8(i,j,kp)+GVM%mc_Jx_8(i-1,j,kp))+&
                  half*(GVM%mc_Jx_8(i,j+1,kp)+GVM%mc_Jx_8(i-1,j+1,kp)))
            C= half*( half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i-1,j,k))+&
                  half*(GVM%mc_Jx_8(i,j+1,k)+GVM%mc_Jx_8(i-1,j+1,k)))
            ski=zero
            if(k.ne.1)  ski= half*( half*(GVM%mc_Jx_8(i,j,km)+GVM%mc_Jx_8(i-1,j,km))+&
                  half*(GVM%mc_Jx_8(i,j+1,km)+GVM%mc_Jx_8(i-1,j+1,km)))
!
          if (k.ne.NK) &
          stencil_V(i,j,3,k)= one* F_coef_8(k)*Jz*C1_8*&
          C3_8 /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))+&
        one* F_coef_8(k)*Jz*GVM%mc_JyT_8(i,j,k)*&
          GVM%mc_Jy_8(i,j,kp) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))
!
          if (k .ne.1) &
          stencil_V(i,j,2,k)= one* F_coef_8(km)*Jzm*C2_8*&
          ski /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))+&
        one* F_coef_8(km)*Jzm*GVM%mc_Jyt_8(i,j,km)*&
          GVM%mc_Jy_8(i,j,km) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))
!
          Jyy=zero
          if(k.ne.1) Jyy= GVM%mc_Jyt_8(i,j,km)    
          stencil_V(i,j,1,k)= -one *F_coef_8(k)*Jz*C1_8*&
              C/((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))-&
          one* F_coef_8(k)*Jz*GVM%mc_Jyt_8(i,j,k)*&
          GVM%mc_Jy_8(i,j,k) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))) &
          -one *F_coef_8(km)*Jzm*C2_8*&
              C/((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))-&
          one* F_coef_8(km)*Jzm*Jyy*&
          GVM%mc_Jy_8(i,j,k) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))
!
        Jz= (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1))

        stencil_V(i,j,3,k)=Jz *stencil_V(i,j,3,k)
        stencil_V(i,j,2,k)=Jz *stencil_V(i,j,2,k)
        stencil_V(i,j,1,k)=Jz *stencil_V(i,j,1,k)

!
    if (k==NK) then
          stencil_V(i,j,2,k )=  (one-(ver_z_8%m(Nk)-ver_z_8%m(Nk-1))/(ver_z_8%m(Nk+1)-ver_z_8%m(Nk-1)))* stencil_V(i,j,2,k-1 )
          stencil_V(i,j,1,k )= (one-(ver_z_8%m(Nk)-ver_z_8%m(Nk-1))/(ver_z_8%m(Nk+1)-ver_z_8%m(Nk-1)))*(stencil_V(i,j,1,k-1 )+&
                 stencil_V(i,j,3,k-1 ))
      endif

!         ceast=zero
!          cwest= zero
!           if (k.ne.1) cwest = stencil_V(i,j,2,k) *fdg2_4(i,j,km)
!            if (k.ne.Nk) ceast= stencil_V(i,j,3,k)* fdg2_4(i,j,kp)
!            C1_8= stencil_V(i,j,1,k)*fdg2_4(i,j,k)+ cwest+ ceast
!            F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( zero+ C1_8)
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

!
! enddo iter

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

