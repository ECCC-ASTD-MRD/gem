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
      subroutine  hzd_init_u ()
      use gem_options
      use geomh
      use glb_ld
      use cstv
      use dcst
      use ver
      use metric
      use hzd_mod
      use hvdif_options
!
      use gmm_geof
      use tdpack
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      integer j,i,k,km,kp
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)

      real(kind=REAL64) Jzx,Jzxm,Jzpx
      real(kind=REAL64) Jzpi,Jz,Jzm,Jzmi,qkm,qkp,Jzmpi
      real(kind=REAL64) C1_8,C2_8,C,ski,skpi,C3_8
      real(kind=REAL64)  ztht_8(l_minx:l_maxx, l_miny:l_maxy,0:l_nk+1),Jxx,Jyy
      real(kind=REAL64) dcoef,beta_imp
      real(kind=REAL64), dimension (:,:,:,:), allocatable :: vsten
!

      allocate(skpu(l_minx:l_maxx, l_miny:l_maxy,l_nk,2) , &
               sku (l_minx:l_maxx, l_miny:l_maxy,l_nk,2))

      allocate(jy(l_minx:l_maxx, l_miny:l_maxy,l_nk,4) , &
              jyp(l_minx:l_maxx, l_miny:l_maxy,l_nk,4)) 
      allocate(jizx(l_minx:l_maxx, l_miny:l_maxy,l_nk,3), &
              jizxm(l_minx:l_maxx, l_miny:l_maxy,l_nk,3) , &
              jizpx(l_minx:l_maxx, l_miny:l_maxy,l_nk,3))

      allocate(xfactu(l_minx:l_maxx, l_miny:l_maxy,l_nk)) 

      allocate(vsten(l_minx:l_maxx, l_miny:l_maxy,3,l_nk))

      allocate( a_u(l_minx:l_maxx, l_miny:l_maxy,l_nk),&
               b_u(l_minx:l_maxx, l_miny:l_maxy,l_nk), &
               c_u(l_minx:l_maxx, l_miny:l_maxy,l_nk), &
               W_u(l_minx:l_maxx, l_miny:l_maxy,l_nk) )

      beta_imp = one

      skpu=0.; sku=0. 

      jy=0.; jyp=0.;  jizx=0.; jizxm=0.; jizpx=0.

      vsten=zero
      a_u=zero
      b_u=zero
      c_u=zero
      W_u=zero
      

      if (Hzd_pwr_z==2) then 
         dcoef = 0.25*Hzd_lnr_z*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      else
         dcoef = 0.25*sqrt(Hzd_lnr_z)*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      endif

      ru= (one-(ver_z_8%m(l_nk)-ver_z_8%m(l_nk-1))/(ver_z_8%m(l_nk+1)-ver_z_8%m(l_nk-1)))

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
            do i=1+pil_w-1, l_ni-pil_e+1
!Dz/Dzeta M-level k on Ui,j position
               Jzpi =(ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))!M-level k on phi+1,j 
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))!M-level k on phi,j
               skpu(i,j,k,1)= (Jzpi+Jz)*half !Dz/Dzeta on M-level k and  Ui,j position
!Dz/Dzeta M-level k on Ui-1,j position!
               Jzpi =(ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))!M-level k on phi,j
               Jz   =(ztht_8(i-1  ,j,k)-ztht_8(i-1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))!M-level k on phi-1,j
               sku(i,j,k,1) = (Jzpi+Jz)*half   !Dz/Dzeta on M-level k  and  Ui-1,j position 

! Gradient component Along Y
            enddo   
         enddo
      enddo

      do k=1,l_nk
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
               Jzpi =(ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!Dz/Dzeta on M-level and  Ui,j position
               skpu(i,j,k,2)= (Jzpi+Jz)*half 
!
               Jzpi =(ztht_8(i+1,j+1,k)-ztht_8(i+1,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j+1,k)-ztht_8(i  ,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!  Dz/Dzeta on M-level and  Ui,j+1 position
              sku(i,j,k,2)= (Jzpi+Jz)*half  
            enddo   
         enddo
      enddo

      do k = 1,l_nk
         km=max(1,k-1)
         kp=min(l_nk,k+1)
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w, l_ni-pil_e
! jy k-1 M-level  U_i,j position
               Jy(i,j,k,1) = half*(half * (GVM%mc_Jy_8(i,j,km) +GVM%mc_Jy_8(i+1,j,km))&
                  + half * (GVM%mc_Jy_8(i,j-1,km) +GVM%mc_Jy_8(i+1,j-1,km)))
! jy k M-level  U_i,j position
               Jyp(i,j,k,1)= half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k))&
                      + half * (GVM%mc_Jy_8(i,j-1,k) +GVM%mc_Jy_8(i+1,j-1,k)))
! jy km M-level  B_i,j+1 position
               Jy(i,j,k,2) = half*(half * (GVM%mc_Jy_8(i,j,km) +GVM%mc_Jy_8(i+1,j,km)) &
                  + half * (GVM%mc_Jy_8(i,j+1,km) +GVM%mc_Jy_8(i+1,j+1,km)))

! jy k M-level  B_i,j+1 position
               Jyp(i,j,k,2)= half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k))&
                      + half * (GVM%mc_Jy_8(i,j+1,k) +GVM%mc_Jy_8(i+1,j+1,k)))
!
! jy k M-level  U_i,j position
               Jy(i,j,k,3) = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) &
                   +half * (GVM%mc_Jy_8(i,j-1,k) +GVM%mc_Jy_8(i+1,j-1,k)))

! jy k+1 M-level  U_i,j position
               Jyp(i,j,k,3)= half*(half * (GVM%mc_Jy_8(i,j,kp) +GVM%mc_Jy_8(i+1,j,kp))  &
                       +half * (GVM%mc_Jy_8(i,j-1,kp) +GVM%mc_Jy_8(i+1,j-1,kp)))
! jy k M-level  U_i,j+1 position
               Jy(i,j,k,4) = half*(half * (GVM%mc_Jy_8(i,j,k) +GVM%mc_Jy_8(i+1,j,k)) &
                   +half * (GVM%mc_Jy_8(i,j+1,k) +GVM%mc_Jy_8(i+1,j+1,k)))
! jy kp M-level  U_i,j position
               Jyp(i,j,k,4)= half*(half * (GVM%mc_Jy_8(i,j,kp) +GVM%mc_Jy_8(i+1,j,kp)) &
                    +   half * (GVM%mc_Jy_8(i,j+1,kp) +GVM%mc_Jy_8(i+1,j+1,kp)))

            enddo
         enddo
         enddo

!  flux
      do k = 2,l_nk
         kp=min(l_nk,k+1)
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
! Dz(jx/Jz A)
! Dz/Dzeta^-1 on M-level k on Afdg_i,j position
               Jzx   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jizx(i,j,k,1)  =(GVM%mc_Jx_8(i,j,k) +GVM%mc_Jx_8(i-1,j,k))/Jzx

               Jzxm  =(ztht_8(i  ,j,k-1)-ztht_8(i  ,j,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jizxm(i,j,k,1) =(GVM%mc_Jx_8(i,j,k-1) +GVM%mc_Jx_8(i-1,j,k-1))/Jzxm
! Dz/Dzeta^-1 on M-level k+1 on Afdg_i,j position
               Jzpx  =(ztht_8(i  ,j,k+1)-ztht_8(i  ,j,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jizpx(i,j,k,1)  =(GVM%mc_Jx_8(i,j,kp) +GVM%mc_Jx_8(i-1,j,kp)) /Jzpx
    
! D(z)/D(zeta) on k M_level on Vi,j position 
               !Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz   =(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jzpi  =(ztht_8(i  ,j+1,k)-ztht_8(i  ,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jz= half*(Jz+Jzpi)
! D(z)/D(zeta) on k M_level on Vi+1,j position
               Jzm   =(ztht_8(i+1  ,j,k)-ztht_8(i+1  ,j,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jzpi  =(ztht_8(i+1  ,j+1,k)-ztht_8(i+1  ,j+1,k-1))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
               Jzpi= half*(Jzm+Jzpi)
! D(z)/D(zeta) on k M_level on Bi,j position
               Jzpx= half*(Jz+Jzpi)  ! Jz on M_level k  and position Bij
               Jizpx(i,j,k,2)=((GVM%mc_Jy_8(i,j,k  ) +GVM%mc_Jy_8(i+1,j,k  )) )/jzpx
! D(z)/D(zeta) on k-1 M_level on Vi,j position
               Jz   =(ztht_8(i  ,j,k-1)-ztht_8(i  ,j,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jzm  =(ztht_8(i  ,j+1,k-1)-ztht_8(i  ,j+1,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jz= half*(Jz+Jzm)
! D(z)/D(zeta) on k-1 M_level on Vi+1,j position
               Jzmpi   =(ztht_8(i+1  ,j,k-1)-ztht_8(i+1  ,j,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jzm  =(ztht_8(i+1  ,j+1,k-1)-ztht_8(i+1  ,j+1,k-2))/(Ver_z_8%t(k-1)-Ver_z_8%t(k-2))
               Jzm= half*(Jzmpi+Jzm)
! D(z)/D(zeta) on k-1 M_level on Bi,j position
               Jzxm= half*(Jz+Jzm)  ! Jz on M_level k-1  and position Bij
               Jizxm(i,j,k,2)= ((GVM%mc_Jy_8(i,j,k-1)  +GVM%mc_Jy_8(i+1,j,k-1)) )/Jzxm

! D(z)/D(zeta) on k+1 M_level on Vi,j position
               Jz   =(ztht_8(i  ,j,k+1)-ztht_8(i  ,j,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jzm  =(ztht_8(i  ,j+1,k+1)-ztht_8(i  ,j+1,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jz= half*(Jz+Jzm)
! D(z)/D(zeta) on k+1 M_level on Vi+1,j position
               Jzmpi   =(ztht_8(i+1  ,j,k+1)-ztht_8(i+1  ,j,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jzm  =(ztht_8(i+1  ,j+1,k+1)-ztht_8(i+1  ,j+1,k))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               Jzm= half*(Jzmpi+Jzm)
!D(z)/D(zeta) on k+1 M_level on  Bi,j position
               Jzx= half*(Jz+Jzm) !Jz on M_level k  and position Bij
               Jizx(i,j,k,3)=(GVM%mc_Jy_8(i,j,kp) +GVM%mc_Jy_8(i+1,j,kp)) /Jzx 

           enddo
         enddo
      enddo
! Apply divergence
      do k = 1, l_nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
!inverse(D(z)/D(zeta)) on k M-level and Ui,j position
         xfactu(i,j,k)= half*(((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i+1,j,k  )-ztht_8(i+1,j,k-1)))+&
                      ((Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1)))) 
            enddo
         enddo
      enddo

! stencilV
        do k=1,l_nk
         km=max(k-1,1)
         kp=min(l_nk,k+1)
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
            if(k.ne.l_nk) C3_8= half*( half*(GVM%mc_Jy_8(i,j,kp)+GVM%mc_Jy_8(i,j-1,kp))+&
                  half*(GVM%mc_Jy_8(i+1,j,kp)+GVM%mc_Jy_8(i+1,j-1,kp)))
            C= half*( half*(GVM%mc_Jy_8(i,j,k)+GVM%mc_Jy_8(i,j-1,k))+&
                  half*(GVM%mc_Jy_8(i+1,j,k)+GVM%mc_Jy_8(i+1,j-1,k)))
            if (k.ne.1)  ski= half*( half*(GVM%mc_Jy_8(i,j,km)+GVM%mc_Jy_8(i,j-1,km))+&
                  half*(GVM%mc_Jy_8(i+1,j,km)+GVM%mc_Jy_8(i+1,j-1,km)))
!
          if (k.ne.l_nk) &
          vsten(i,j,3,k)= one* dcoef*Jz*C1_8*&
          C3_8 /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))+&
          one* dcoef*Jz*GVM%mc_JxT_8(i,j,k)*&
          GVM%mc_Jx_8(i,j,kp) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))
!
          if (k .ne.1) &
          vsten(i,j,2,k)= one* dcoef*Jzm*C2_8*&
          ski /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))+&
        one* dcoef*Jzm*GVM%mc_Jxt_8(i,j,km)*&
          GVM%mc_Jx_8(i,j,km) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))
!
          Jxx=zero
          if (k.ne.1) Jxx= GVM%mc_Jxt_8(i,j,km)
          vsten(i,j,1,k)= -one *dcoef*Jz*C1_8*&
              C/((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))-&
          one* dcoef*Jz*GVM%mc_Jxt_8(i,j,k)*&
          GVM%mc_Jx_8(i,j,k) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))) &
          -one *dcoef*Jzm*C2_8*&
              C/((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))-&
          one* dcoef*Jzm*Jxx*&
          GVM%mc_Jx_8(i,j,k) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)))
!
            Jz= (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i,j,k  )-ztht_8(i,j,k-1))

            vsten(i,j,3,k )=Jz*  vsten(i,j,3,k )
            vsten(i,j,2,k )=Jz*  vsten(i,j,2,k )
            vsten(i,j,1,k )=Jz*  vsten(i,j,1,k )
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
                c_u(i,j,k)= -beta_imp*vsten(i,j,3,k )
                b_u(i,j,k)=one-beta_imp*vsten(i,j,1,k)
            enddo
         enddo
         do k=2,l_nk-1
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  a_u(i,j,k)= -beta_imp*vsten(i,j,2,k)
                  b_u(i,j,k)= one-beta_imp*vsten(i,j,1,k)
                  c_u(i,j,k)= -beta_imp*vsten(i,j,3,k)
               enddo
            enddo
         enddo
         k=l_nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  a_u(i,j,k)= -beta_imp*vsten(i,j,2,k)
                  b_u(i,j,k)= one-beta_imp*vsten(i,j,1,k)

               enddo
            enddo
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               do k = 2 , l_nk
                  W_u(i,j,k) = a_u(i,j,k) / b_u(i,j,k - 1)
                  b_u(i,j,k) = b_u(i,j,k) - W_u(i,j,k) * c_u(i,j,k - 1)
               enddo
            enddo
         enddo

       deallocate(vsten)

      return
      end
