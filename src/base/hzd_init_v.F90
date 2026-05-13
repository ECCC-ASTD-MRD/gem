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
      subroutine  hzd_init_v ()
      use gem_options
      use geomh
      use glb_ld
      use cstv
      use dcst
      use ver
      use metric
      use hzd_mod
      use hvdif_options
      use gmm_geof
      use tdpack
      use ptopo
!
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer j,i,k,km,kp
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)

      real(kind=REAL64) Jzpi,Jz,Jzm,Jzmi,qkm,qkp,Jzmpi
      real(kind=REAL64) C1_8,C2_8,C,ski,skpi,C3_8
      real(kind=REAL64)  ztht_8(l_minx:l_maxx, l_miny:l_maxy,0:l_nk+1),Jxx,Jyy
      real(kind=REAL64) dcoef,beta_imp,beta_exp
      real(kind=REAL64) deno
      real(kind=REAL64), dimension (:,:,:,:), allocatable :: vsten
!
      allocate(skpv(l_minx:l_maxx, l_miny:l_maxy,l_nk,2) , &
               skv (l_minx:l_maxx, l_miny:l_maxy,l_nk,2))

      allocate(jx(l_minx:l_maxx, l_miny:l_maxy,l_nk,4) , &
              jxp(l_minx:l_maxx, l_miny:l_maxy,l_nk,4)) 
      allocate(jzv(l_minx:l_maxx, l_miny:l_maxy,l_nk,3), &
              jzvm(l_minx:l_maxx, l_miny:l_maxy,l_nk,2)) 
      !        jzpix(l_minx:l_maxx, l_miny:l_maxy,l_nk,3))

      allocate(xfactv(l_minx:l_maxx, l_miny:l_maxy,l_nk)) 

      allocate(vsten(l_minx:l_maxx, l_miny:l_maxy,3,l_nk))
      allocate( a_v(l_minx:l_maxx, l_miny:l_maxy,l_nk),&
               b_v(l_minx:l_maxx, l_miny:l_maxy,l_nk), &
               c_v(l_minx:l_maxx, l_miny:l_maxy,l_nk), &
               W_v(l_minx:l_maxx, l_miny:l_maxy,l_nk) )

      beta_imp = one

      skpv=zero; skv=zero;
      jx=zero; jxp=zero
      jzv=zero; jzvm=zero
      xfactv=zero

      vsten=zero
      a_v=zero
      b_v=zero
      c_v=zero
      W_v=zero
      
      if (Hzd_pwr_z==2) then 
         dcoef = 0.25*Hzd_lnr_z*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      else
         dcoef = 0.25*sqrt(Hzd_lnr_z)*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      endif

      rv= (one-(ver_z_8%m(l_nk)-ver_z_8%m(l_nk-1))/(ver_z_8%m(l_nk+1)-ver_z_8%m(l_nk-1)))

      do j=1-G_haloy,l_nj+G_haloy
        do i=1-G_halox,l_ni+G_halox
          do k=1 ,l_nk
           ztht_8(i,j,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*fis0(i,j)+Ver_c_8%t(k)*sls(i,j))/grav_8
          enddo
            ztht_8(i,j,0)   =   GVM%zmom_8(i,j,0)
            ztht_8(i,j,l_nk+1)=   GVM%zmom_8(i,j,l_nk+1)
        enddo
      enddo

      do k=1,l_nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
!Dz/Dzeta M-level k on Vi,j position
               Jzpi =(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1)) !Dz/Dzeta M-level k
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1)) !Dz/Dzeta M-level k
               skv(i,j,k,1)= (Jzpi+Jz)*half                                                             !Dz/Dzeta M-level k
!Dz/Dzeta M-level k on Vi+1,j position
               Jzpi =(ztht_8(i+1,j+1,k)-ztht_8(i+1,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1)) !Dz/Dzeta M-level k
               Jz   =(ztht_8(i+1  ,j,k)-ztht_8(i+1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1)) !Dz/Dzeta M-level k
               skpv(i,j,k,1)= (Jzpi+Jz)*half                                                                !Dz/Dzeta M-level k  
! Jx on M K+1 level and Vi,j 
          enddo
        enddo
      enddo

! Gradient component Along Y
      do k = 1,l_nk
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
!Dz/Dzeta M-level k on Vi,j position
               Jzpi =(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               skpv(i,j,k,2)= (Jzpi+Jz)*half
!
!Dz/Dzeta M-level k on Vi,j-1 position
               Jzpi =(ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j-1,k)-ztht_8(i  ,j-1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               skv(i,j,k,2) = (Jzpi+Jz)*half
            enddo
         enddo
      enddo    

      do k = 1,l_nk
         km=max(k-1,1)
         kp=min(k+1,l_nk)
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
! Jx on M K level and Vi,j 
               Jxp(i,j,k,1)= half*(half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)) &
                   +half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i-1,j+1,k)))
! Jx on M K-1 level and Vi,j 
               Jx(i,j,k,1)= half*(half*(GVM%mc_Jx_8(i,j,km)+GVM%mc_Jx_8(i,j+1,km))&
                  + half*(GVM%mc_Jx_8(i-1,j,km)+GVM%mc_Jx_8(i-1,j+1,km)))
! Jx on M K level and Vi+1,j 
               Jxp(i,j,k,2)= half*(half*(GVM%mc_Jx_8(i+1,j,k)+GVM%mc_Jx_8(i+1,j+1,k)) &
                  +  half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)))
! Jx on M K-1 level and Vi,j 
               Jx(i,j,k,2)= half*(half*(GVM%mc_Jx_8(i+1,j,km)+GVM%mc_Jx_8(i+1,j+1,km)) &
                  + half*(GVM%mc_Jx_8(i,j,km)+GVM%mc_Jx_8(i,j+1,km)))
! Jx on M K+1 level and Vi,j 
               Jxp(i,j,k,3)= half*(half*(GVM%mc_Jx_8(i,j,kp)+GVM%mc_Jx_8(i,j+1,kp)) &
                  + half*(GVM%mc_Jx_8(i-1,j,kp)+GVM%mc_Jx_8(i-1,j+1,kp)))
! Jx on M K level and Vi,j 
               Jx(i,j,k,3)= half*(half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)) &
                  + half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i-1,j+1,k)))
! Jx on M K+1 level and Vi+1,j 
               Jxp(i,j,k,4)= half*(half*(GVM%mc_Jx_8(i+1,j,kp)+GVM%mc_Jx_8(i+1,j+1,kp)) &
                  + half*(GVM%mc_Jx_8(i,j,kp)+GVM%mc_Jx_8(i,j+1,kp)))
! Jx on M K level and Vi,j 
               Jx(i,j,k,4)= half*(half*(GVM%mc_Jx_8(i+1,j,k)+GVM%mc_Jx_8(i+1,j+1,k)) &
                  + half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i,j+1,k)))
          enddo
        enddo
      enddo
!       

!  flux
! start at pil_w and finish at l_nj-pil_n+1 ! important pour le calcul de sol a la fin
      do k = 2,l_nk
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
               Jzv(i,j,k,1) = half*(Jz + Jzpi) ! Jz on M-level k on Aij point
! Dz/Dzeta M-level k-1 on Vi,j position
               Jzm   =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i  ,j,k-1)-ztht_8(i  ,j,k-2))
               Jzpi  =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i  ,j+1,k-1)-ztht_8(i  ,j+1,k-2))
               Jzm= half*(Jzm+Jzpi)
!Dz/Dzeta M-level k-1 on Vi+1,j position
               Jzmpi =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i+1  ,j,k-1)-ztht_8(i+1  ,j,k-2))
               Jzpi  =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i+1  ,j+1,k-1)-ztht_8(i+1  ,j+1,k-2))
               Jzpi= half*(Jzmpi+Jzpi)
!Dz/Dzeta M-level k-1 on Ai,j position
               Jzvm(i,j,k,1) = half*(Jzm + Jzpi) ! Jz on M-level k-1 on Aij point

!Jz on M-level k on Bij point
               Jzv(i,j,k,2)   =(Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1)) 
!Jz on M-level k-1 on Bij point
               Jzvm(i,j,k,2)   =(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i  ,j,k-1)-ztht_8(i  ,j,k-2))

               Jz= (Ver_z_8%t(k-1)-Ver_z_8%t(k-2))/(ztht_8(i,j,k)-ztht_8(i,j,k-1))
               Jzpi= (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))
               Jzv(i,j,k,3) = half*(Jz+Jzpi)
!
           enddo
         enddo
        enddo
! Apply divergence
      do k = 1, l_nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
! Dzeta/Dz^-1  on  k M-level and Vi,j position
         C    =  ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j+1,k  )-ztht_8(i,j+1,k-1)))
         deno =  ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1)))
         xfactv(i,j,k) =  half*(deno+C)               ! Dzeta/Dz 
!X-divergence (Dz/Dzeta)^-1*DAfdg/Dx on  k M-level and Vi,j position
         !add_v8(i,j,k) =    deno*  (Afdg1 (i,j,k)-Afdg1 (i-1,j,k))*geomh_invDXv_8(j)
            enddo
         enddo
      enddo

        do k=1,l_nk
           km = max(k-1,1)
           kp = min(k+1, l_nk)
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
            if(k.ne.l_nk)  C3_8= half*( half*(GVM%mc_Jx_8(i,j,kp)+GVM%mc_Jx_8(i-1,j,kp))+&
                  half*(GVM%mc_Jx_8(i,j+1,kp)+GVM%mc_Jx_8(i-1,j+1,kp)))
            C= half*( half*(GVM%mc_Jx_8(i,j,k)+GVM%mc_Jx_8(i-1,j,k))+&
                  half*(GVM%mc_Jx_8(i,j+1,k)+GVM%mc_Jx_8(i-1,j+1,k)))
            ski=zero
            if(k.ne.1)  ski= half*( half*(GVM%mc_Jx_8(i,j,km)+GVM%mc_Jx_8(i-1,j,km))+&
                  half*(GVM%mc_Jx_8(i,j+1,km)+GVM%mc_Jx_8(i-1,j+1,km)))
!
          if (k.ne.l_nk) &
          vsten(i,j,3,k)= one*dcoef*Jz*C1_8*&
          C3_8 /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))+&
        one*dcoef*Jz*GVM%mc_JyT_8(i,j,k)*&
          GVM%mc_Jy_8(i,j,kp) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))
!
          if (k .ne.1) &
          vsten(i,j,2,k)= one* dcoef*Jzm*C2_8*&
          ski /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))+&
        one* dcoef*Jzm*GVM%mc_Jyt_8(i,j,km)*&
          GVM%mc_Jy_8(i,j,km) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))
!
          Jyy=zero
          if(k.ne.1) Jyy= GVM%mc_Jyt_8(i,j,km)    
          vsten(i,j,1,k)= -one *dcoef*Jz*C1_8*&
              C/((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))-&
          one*dcoef*Jz*GVM%mc_Jyt_8(i,j,k)*&
          GVM%mc_Jy_8(i,j,k) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))) &
          -one *dcoef*Jzm*C2_8*&
              C/((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))-&
          one* dcoef*Jzm*Jyy*&
          GVM%mc_Jy_8(i,j,k) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))
!
        Jz= (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1))

        vsten(i,j,3,k)=Jz *vsten(i,j,3,k)
        vsten(i,j,2,k)=Jz *vsten(i,j,2,k)
        vsten(i,j,1,k)=Jz *vsten(i,j,1,k)

!
    	if (k==l_nk) then
        	  vsten(i,j,2,k )=  (one-(ver_z_8%m(l_nk)-ver_z_8%m(l_nk-1))/(ver_z_8%m(l_nk+1)-ver_z_8%m(l_nk-1)))* vsten(i,j,2,k-1 )
          	vsten(i,j,1,k )= (one-(ver_z_8%m(l_nk)-ver_z_8%m(l_nk-1))/(ver_z_8%m(l_nk+1)-ver_z_8%m(l_nk-1)))*(vsten(i,j,1,k-1 )+&
                 vsten(i,j,3,k-1 ))

     	endif
             enddo
         enddo
      enddo
      

!goto 100
         beta_imp=beta_imp*Cstv_dt_8
         k=1
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
                c_v(i,j,k)= -beta_imp*vsten(i,j,3,k )
                b_v(i,j,k)=one-beta_imp*vsten(i,j,1,k)
            enddo
         enddo
         do k=2,l_nk-1
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  a_v(i,j,k)= -beta_imp*vsten(i,j,2,k)
                  b_v(i,j,k)= one-beta_imp*vsten(i,j,1,k)
                  c_v(i,j,k)= -beta_imp*vsten(i,j,3,k)
               enddo
            enddo
         enddo
         k=l_nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  a_v(i,j,k)= -beta_imp*vsten(i,j,2,k)
                  b_v(i,j,k)= one-beta_imp*vsten(i,j,1,k)

               enddo
            enddo
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               do k = 2 , l_nk
                  W_v(i,j,k) = a_v(i,j,k) / b_v(i,j,k - 1)
                  b_v(i,j,k) = b_v(i,j,k) - W_v(i,j,k) * c_v(i,j,k - 1)
               enddo
            enddo
         enddo

      deallocate(vsten)

      return
      end
