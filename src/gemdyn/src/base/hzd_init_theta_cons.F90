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
      subroutine  hzd_init_theta_cons ()
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

      use, intrinsic :: iso_fortran_env
      implicit none

      integer j,i,k,km,kp
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)
      real(kind=REAL64) Jzp1,Jz1,Jzm1
      real(kind=REAL64)  ztht_8(l_minx:l_maxx, l_miny:l_maxy,0:l_nk+1)
      real dcoef
!
      allocate(jzpt(l_minx:l_maxx, l_miny:l_maxy,l_nk,2), &
               jzt(l_minx:l_maxx, l_miny:l_maxy,l_nk,2))
      allocate(jxt(l_minx:l_maxx, l_miny:l_maxy,l_nk,8))
      allocate(jyt(l_minx:l_maxx, l_miny:l_maxy,l_nk,8))
      allocate(jm(l_minx:l_maxx, l_miny:l_maxy,l_nk,2))
      allocate(jzz(l_minx:l_maxx, l_miny:l_maxy,l_nk,2))

      jzpt=zero ; jzt=zero 
      jxt=zero  ; jyt=zero
      jm=zero; jzz=zero

      do j=1-G_haloy,l_nj+G_haloy
        do i=1-G_halox,l_ni+G_halox
          do k=1 ,l_nk
           ztht_8(i,j,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*fis0(i,j)+Ver_c_8%t(k)*sls(i,j))/grav_8
          enddo
            ztht_8(i,j,0)   =   GVM%zmom_8(i,j,0)
            ztht_8(i,j,l_nk+1)=   GVM%zmom_8(i,j,l_nk+1)
        enddo
      enddo

      if (Hzd_pwr_z==2) then 
         dcoef = 0.25*Hzd_lnr_theta_z*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      else
         dcoef = 0.25*sqrt(Hzd_lnr_theta_z)*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      endif

         do k = 1,l_nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w-1, l_ni-pil_e
                  Jzpt(i,j,k,1)= (GVM%zmom_8(i+1,j,k+1)-GVM%zmom_8(i+1,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
                  Jzt(i,j,k,1)  = (GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k  ))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               enddo
            enddo
         enddo
         do k = 1,l_nk
            do j=1+pil_s-1, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  Jzpt(i,j,k,2)= (GVM%zmom_8(i,j+1,k+1)-GVM%zmom_8(i,j+1,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
                  Jzt(i,j,k,2)  = (GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)) /(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               enddo
            enddo
         enddo

         do k = 1,l_nk
            km=max(k-1,1)
            kp=min(k+1,l_nk)
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w-1, l_ni-pil_e
               Jxt(i,j,k,1)=half*(GVM%mc_Jxt_8(i-1,j,k  )+GVM%mc_Jxt_8(i,j  ,k   ))*Ver_idz_8%m(k)
               Jxt(i,j,k,2)=half*(GVM%mc_Jxt_8(i-1,j,km) +GVM%mc_Jxt_8(i,j  ,km  ))*Ver_idz_8%m(k)
               Jxt(i,j,k,3)=half*(GVM%mc_Jxt_8(i,j,k    )+GVM%mc_Jxt_8(i+1,j,k   ))*Ver_idz_8%m(k)
               Jxt(i,j,k,4)=half*(GVM%mc_Jxt_8(i,j,km  ) +GVM%mc_Jxt_8(i+1,j,km  ))*Ver_idz_8%m(k)

               Jxt(i,j,k,5)=half*(GVM%mc_Jxt_8(i-1,j,kp)+GVM%mc_Jxt_8(i,j,kp  ))*Ver_idz_8%m(kp)
               Jxt(i,j,k,6)=half*(GVM%mc_Jxt_8(i-1,j,k  )+GVM%mc_Jxt_8(i,j,k    ))*Ver_idz_8%m(kp)
               Jxt(i,j,k,7)=half*(GVM%mc_Jxt_8(i,j,kp  )+GVM%mc_Jxt_8(i+1,j,kp))*Ver_idz_8%m(kp)
               Jxt(i,j,k,8)=half*(GVM%mc_Jxt_8(i,j,k    )+GVM%mc_Jxt_8(i+1,j,k  ))*Ver_idz_8%m(kp)
               enddo
            enddo
         enddo
         do k = 1,l_nk
            km=max(k-1,1)
            kp=min(k+1,l_nk)
            do j=1+pil_s-1, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  Jyt(i,j,k,1)= half*(GVM%mc_Jyt_8(i,j-1,k  )+GVM%mc_Jyt_8(i,j,k    ))*Ver_idz_8%m(k)
                  Jyt(i,j,k,2)= half*(GVM%mc_Jyt_8(i,j-1,km)+GVM%mc_Jyt_8(i,j,km  ))*Ver_idz_8%m(k)
                  Jyt(i,j,k,3)= half*(GVM%mc_Jyt_8(i,j,k    )+GVM%mc_Jyt_8(i,j+1,k  ))*Ver_idz_8%m(k)
                  Jyt(i,j,k,4)= half*(GVM%mc_Jyt_8(i,j,km  )+GVM%mc_Jyt_8(i,j+1,km))*Ver_idz_8%m(k)
                  Jyt(i,j,k,5)= half*(GVM%mc_Jyt_8(i,j-1,kp)+GVM%mc_Jyt_8(i,j,kp  ))*Ver_idz_8%m(kp)
                  Jyt(i,j,k,6)= half*(GVM%mc_Jyt_8(i,j-1,k  )+GVM%mc_Jyt_8(i,j,k    ))*Ver_idz_8%m(kp)
                  Jyt(i,j,k,7)= half*(GVM%mc_Jyt_8(i,j,kp  )+GVM%mc_Jyt_8(i,j+1,kp))*Ver_idz_8%m(kp)
                  Jyt(i,j,k,8)= half*(GVM%mc_Jyt_8(i,j,k    )+GVM%mc_Jyt_8(i,j+1,k  ))*Ver_idz_8%m(kp)
               enddo
            enddo
         enddo
         do k = 1, l_nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  Jm(i,j,k,1) = ((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k))) &
                                   *geomh_invDXM_8(j)
                  Jm(i,j,k,2) = ((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k))) &
                                   *geomh_invDYMv_8(j) /geomh_cy_8(j)
               enddo
            enddo
         enddo

      do k=1,l_nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
                  Jzz(i,j,k,1)  = (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))
                  Jzz(i,j,k,2)=((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)))
            enddo
         enddo
      enddo

         allocate(stencil_V(l_minx:l_maxx, l_miny:l_maxy,3,l_nk))
         allocate(stencil_V1(l_minx:l_maxx, l_miny:l_maxy,l_nk))
         allocate(stencil_V2(l_minx:l_maxx, l_miny:l_maxy,l_nk))
! stencil_V
         stencil_V=zero
         stencil_V1=zero
         stencil_V2=zero

         do k=2,l_nk-1
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  Jzp1=((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)))
                  Jz1  =(Ver_z_8%t(k+1)-Ver_z_8%t(k))/(ztht_8(i  ,j,k+1  )-ztht_8(i  ,j,k))
                  Jzm1 =(Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k  )-ztht_8(i  ,j,k-1))

                     stencil_V(i,j,3,k)= Jzp1*(&
        half* dcoef*Jz1*half*(GVM%mc_Jy_8(i,j-1,k+1)+GVM%mc_Jy_8(i,j,k+1))*(&
           half*(GVM%mc_Jyt_8(i,j,k+1)+GVM%mc_Jyt_8(i,j-1,k+1)) /((Ver_z_8%t(k+1)-Ver_z_8%t(k))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))))+&
        half* dcoef*Jz1*half*(GVM%mc_Jx_8(i-1,j,k+1)+GVM%mc_Jx_8(i,j,k+1))*(&
         half*(GVM%mc_Jxt_8(i,j,k+1)+GVM%mc_Jxt_8(i-1,j,k+1)) /((Ver_z_8%t(k+1)-Ver_z_8%t(k))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))) )
!
                     stencil_V(i,j,2,k)= Jzp1* (&
      half* dcoef*Jzm1*half*(GVM%mc_Jy_8(i,j-1,k)+GVM%mc_Jy_8(i,j,k))*(&
             half*(GVM%mc_Jyt_8(i,j,k-1)+GVM%mc_Jyt_8(i,j-1,k-1)) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))))+&
      half* dcoef*Jzm1*half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i,j,k))*(&
      half*(GVM%mc_Jxt_8(i,j,k-1)+GVM%mc_Jxt_8(i-1,j,k-1)) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))) )
!
                   stencil_V1(i,j,k)= Jzp1*(&
     - half* dcoef*Jzm1*half*(GVM%mc_Jy_8(i,j-1,k)+GVM%mc_Jy_8(i,j,k))*(&
               half*(GVM%mc_Jyt_8(i,j,k)+GVM%mc_Jyt_8(i,j-1,k)) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))))&
     -half* dcoef*Jzm1*half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i,j,k))*(&
              half*(GVM%mc_Jxt_8(i,j,k)+GVM%mc_Jxt_8(i-1,j,k)) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))) )!&

                   stencil_V2(i,j,k)= Jzp1*(&
       -half*one* dcoef*Jz1*half*(GVM%mc_Jy_8(i,j-1,k+1)+GVM%mc_Jy_8(i,j,k+1))*(&
             half*(GVM%mc_Jyt_8(i,j,k)+GVM%mc_Jyt_8(i,j-1,k)) /((Ver_z_8%t(k+1)-Ver_z_8%t(k))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))))&
  -half* dcoef*Jz1*half*(GVM%mc_Jx_8(i-1,j,k+1)+GVM%mc_Jx_8(i,j,k+1))*(&
                half*(GVM%mc_Jxt_8(i,j,k)+GVM%mc_Jxt_8(i-1,j,k)) /((Ver_z_8%t(k+1)-Ver_z_8%t(k))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))) )

               enddo
            enddo
         enddo

         k=l_nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
                  Jzp1=((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)))
                  Jz1  =(Ver_z_8%t(k+1)-Ver_z_8%t(k))/(ztht_8(i  ,j,k+1  )-ztht_8(i  ,j,k))
                  Jzm1 =(Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k  )-ztht_8(i  ,j,k-1))
             stencil_V(i,j,2,k)= Jzp1 *(&
      half* dcoef*Jzm1*half*(GVM%mc_Jy_8(i,j-1,k)+GVM%mc_Jy_8(i,j,k))*(&
      half*(GVM%mc_Jyt_8(i,j,k-1)+GVM%mc_Jyt_8(i,j-1,k-1)) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))))+&
      half* dcoef*Jzm1*half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i,j,k))*(&
      half*(GVM%mc_Jxt_8(i,j,k-1)+GVM%mc_Jxt_8(i-1,j,k-1)) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))) )
!
                stencil_V1(i,j,k)= Jzp1* ( &
     - half* dcoef*Jzm1*half*(GVM%mc_Jy_8(i,j-1,k)+GVM%mc_Jy_8(i,j,k))*(&
           half*(GVM%mc_Jyt_8(i,j,k)+GVM%mc_Jyt_8(i,j-1,k)) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k))))&
     - half* dcoef*Jzm1*half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i,j,k))*(&
            half*(GVM%mc_Jxt_8(i,j,k)+GVM%mc_Jxt_8(i-1,j,k)) /((Ver_z_8%t(k)-Ver_z_8%t(k-1))*(Ver_z_8%m(k+1)-Ver_z_8%m(k)))) ) !&

              stencil_V(i,j,2,k)=0.d0

           zfact=(one-(ver_z_8%t(l_nk)-ver_z_8%t(l_nk-1))/(ver_z_8%t(l_nk+1)-ver_z_8%t(l_nk-1)))  
         enddo
      enddo   

      return
      end
