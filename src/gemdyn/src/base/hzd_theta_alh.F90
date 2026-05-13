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
      subroutine  hzd_theta_alh ( F_Sol1,HzdlnR, Minx, Maxx, Miny, Maxy,Nk)
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
!
      use ptopo
      use stat_mpi, only: statf_dm
      use gmm_geof
      use tdpack


!
      use, intrinsic :: iso_fortran_env
      implicit none
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


      integer j,i,k,halox,haloy
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)
      real(kind=REAL64)   Afdg1(l_minx:l_maxx, l_miny:l_maxy,Nk),Afdg2(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64)   Bfdg1(l_minx:l_maxx, l_miny:l_maxy,Nk),Bfdg2(l_minx:l_maxx, l_miny:l_maxy,Nk)

      real(kind=REAL64) Jzpi,Jz,Jzm,Jzmi,qkm,qkp,Jzmpi
      real(kind=REAL64) C1_8,C2_8,C,ski,skpi,skip,skpip
      real(kind=REAL64) bdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) add_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) cdd1_v8(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) cdd2_v8(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) Cflux_8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) C1flux_8(l_minx:l_maxx, l_miny:l_maxy,Nk),dsten(l_minx:l_maxx, l_miny:l_maxy)
      real(kind=REAL64)  ztht_8(l_minx:l_maxx, l_miny:l_maxy,0:Nk+1),Jxx,Jyy
      real(kind=REAL64) dcoef
      real cflux(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real fdg2_4(l_minx:l_maxx, l_miny:l_maxy,Nk+1)


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

      if (Hzd_pwr_z==2) then 
         dcoef = 0.25*HzdlnR*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      else
         dcoef = 0.25*sqrt(HzdlnR)*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      endif

      Afdg1 = .0d0
      Bfdg1 = .0d0
      add_v8 =0.0d0
      bdd_v8=0.0d0
      cdd1_v8=0.0d0
      cdd2_v8=0.0d0
      fdg2_4 =0.0
      cflux=0.0
      Cflux_8=zero

      do k = 1, nk
         do j=1+pil_s-1, l_nj-pil_n+1
            do i=1+pil_w-1, l_ni-pil_e+1
               fdg2_4(i,j,k )=F_Sol1(i,j,k)
            enddo
         enddo
      enddo

      call rpn_comm_xch_halo(fdg2_4,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,Nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         k=1
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
               Jzpi= (GVM%zmom_8(i+1,j,k+1)-GVM%zmom_8(i+1,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               Jz  = (GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)) /(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C1_8= (fdg2_4(i,j,k+1)*(half*(GVM%mc_Jxt_8(i-1,j,k+1)+GVM%mc_Jxt_8(i,j,k+1)))-&
                     fdg2_4(i,j,k  )*(half*(GVM%mc_Jxt_8(i-1,j,k  )+GVM%mc_Jxt_8(i,j,k  )))) /(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               C2_8= (fdg2_4(i+1,j,k+1)*(half*(GVM%mc_Jxt_8(i,j,k+1)+GVM%mc_Jxt_8(i+1,j,k+1)))-&
                     fdg2_4(i+1,j,k  )*(half*(GVM%mc_Jxt_8(i,j,k  )+GVM%mc_Jxt_8(i+1,j,k  )))) /(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               C1_8= half*(C1_8+C2_8)
! Boundary condition at k-1
               C=zero
               C= half*(C+C1_8)
               Afdg1(i,j,k) = dcoef*((Jzpi*fdg2_4(i+1,j,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDX_8(j) -C)
               Afdg2(i,j,k) = ((Jzpi*fdg2_4(i+1,j,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDX_8(j))

            enddo
         enddo
!
         k= NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
               Jzpi= (GVM%zmom_8(i+1,j,k+1)-GVM%zmom_8(i+1,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               Jz  = (GVM%zmom_8(i  ,j,k+1)-GVM%zmom_8(i  ,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C1_8=(fdg2_4(i,j,k  )*(half*(GVM%mc_Jxt_8(i-1,j,k  )+GVM%mc_Jxt_8(i,j,k  )))-&
                    fdg2_4(i,j,k-1)*(half*(GVM%mc_Jxt_8(i-1,j,k-1)+GVM%mc_Jxt_8(i,j,k-1)))) / (Ver_z_8%t(k)-Ver_z_8%t(k-1))

               C2_8=(fdg2_4(i+1,j,k  )*(half*(GVM%mc_Jxt_8(i,j,k  )+GVM%mc_Jxt_8(i+1,j,k  )))-&
                    fdg2_4(i+1,j,k-1)*(half*(GVM%mc_Jxt_8(i,j,k-1)+GVM%mc_Jxt_8(i+1,j,k-1)))) / (Ver_z_8%t(k)-Ver_z_8%t(k-1))

            C1_8= half*(C1_8+C2_8)
! Boundary condtion at k+1
            C=zero
            C= half*(C+C1_8)
            Afdg1(i,j,k) = dcoef* ((Jzpi*fdg2_4(i+1,j,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDX_8(j) - C)

            Jzmi = (GVM%zmom_8(i+1,j,k)-GVM%zmom_8(i+1,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
            Jzm  = (GVM%zmom_8(i  ,j,k)-GVM%zmom_8(i  ,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
!
            Afdg2(i,j,k) = ((Jzpi*fdg2_4(i+1,j,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDX_8(j))

            enddo
         enddo
!
         do k = 2,Nk-1
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w-1, l_ni-pil_e
                  ski=(fdg2_4(i,j,k )*(half*(GVM%mc_Jxt_8(i-1,j,k )+GVM%mc_Jxt_8(i,j,k )))-&
                       fdg2_4(i,j,k-1)*(half*(GVM%mc_Jxt_8(i-1,j,k-1)+GVM%mc_Jxt_8(i,j,k-1))))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
                  skip=(fdg2_4(i+1,j,k) *(half*(GVM%mc_Jxt_8(i,j,k )+GVM%mc_Jxt_8(i+1,j,k )))-&
                       fdg2_4(i+1,j,k-1)*(half*(GVM%mc_Jxt_8(i,j,k-1)+GVM%mc_Jxt_8(i+1,j,k-1))))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
                  skpi=(fdg2_4(i,j,k+1)*(half*(GVM%mc_Jxt_8(i-1,j,k+1)+GVM%mc_Jxt_8(i,j,k+1)))-&
                       fdg2_4(i,j,k )*(half*(GVM%mc_Jxt_8(i-1,j,k )+GVM%mc_Jxt_8(i,j,k ))))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
                  skpip=(fdg2_4(i+1,j,k+1)*(half*(GVM%mc_Jxt_8(i,j,k+1)+GVM%mc_Jxt_8(i+1,j,k+1)))-&
                       fdg2_4(i+1,j,k )*(half*(GVM%mc_Jxt_8(i,j,k )+GVM%mc_Jxt_8(i+1,j,k ))))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
                  qkm=  half*(ski+skip)
                  qkp=  half*(skpi+skpip)
                  C=   half*(qkp+qkm)
                  Jzpi= (GVM%zmom_8(i+1,j,k+1)-GVM%zmom_8(i+1,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
                  Jz  = (GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k  ))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))

                  Afdg1(i,j,k) = dcoef*((Jzpi*fdg2_4(i+1,j,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDX_8(j) - C)
                  Jzmi= (GVM%zmom_8(i+1,j,k)-GVM%zmom_8(i+1,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
                  Jzm  = (GVM%zmom_8(i  ,j,k)-GVM%zmom_8(i  ,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
                  Afdg2(i,j,k) =((Jzpi*fdg2_4(i+1,j,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDX_8(j))
               enddo
            enddo
         enddo
         k=1
         do j=1+pil_s-1, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               Jzpi= (GVM%zmom_8(i,j+1,k+1)-GVM%zmom_8(i,j+1,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               Jz  = (GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k  ))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C1_8=(fdg2_4(i,j,k+1)*(half*(GVM%mc_Jyt_8(i,j-1,k+1)+GVM%mc_Jyt_8(i,j,k+1)))-&
                   fdg2_4(i,j,k  )*(half*(GVM%mc_Jyt_8(i,j-1,k  )+GVM%mc_Jyt_8(i,j,k  )))) /(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               C2_8=(fdg2_4(i,j+1,k+1)*(half*(GVM%mc_Jyt_8(i,j,k+1)+GVM%mc_Jyt_8(i,j+1,k+1)))-&
                   fdg2_4(i,j+1,k  )*(half*(GVM%mc_Jyt_8(i,j,k  )+GVM%mc_Jyt_8(i,j+1,k  )))) /(Ver_z_8%t(k+1)-Ver_z_8%t(k))
               C1_8= half*(C1_8+C2_8)
! Boundary condition at k-1
               C=zero
               C= half*(C+C1_8)
               Bfdg1(i,j,k) = (Jzpi*fdg2_4(i,j+1,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j) - C
               Bfdg1(i,j,k) =  dcoef*Bfdg1(i,j,k) * geomh_cyv_8(j)
               Bfdg2(i,j,k) =(Jzpi*fdg2_4(i,j+1,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
            enddo
         enddo
!
         k= NK
         do j=1+pil_s-1, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               Jzpi= (GVM%zmom_8(i,j+1,k+1)-GVM%zmom_8(i,j+1,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               Jz  = (GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)) /(Ver_z_8%m(k+1)-Ver_z_8%m(k))
               C1_8=(fdg2_4(i,j,k  )*(half*(GVM%mc_Jyt_8(i,j-1,k  )+GVM%mc_Jyt_8(i,j,k)))-&
                   fdg2_4(i,j,k-1)*(half*(GVM%mc_Jyt_8(i,j-1,k-1)+GVM%mc_Jyt_8(i,j,k-1))))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
!D(jyt f)/Dzeta on M-level k on phiij+1
               C2_8=(fdg2_4(i,j+1,  k)*(half*(GVM%mc_Jyt_8(i,j,k  )+GVM%mc_Jyt_8(i,j+1,k  )))-&
                   fdg2_4(i,j+1,k-1)*(half*(GVM%mc_Jyt_8(i,j,k-1)+GVM%mc_Jyt_8(i,j+1,k-1))))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))

               C1_8= half*(C1_8+C2_8)
! Boundary condition at k+1
               C=zero
               C= half*(C+C1_8)
               Bfdg1(i,j,k) =(Jzpi*fdg2_4(i,j+1,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j) - C
               Bfdg1(i,j,k) =  dcoef*Bfdg1(i,j,k) * geomh_cyv_8(j)
               Jzmi= (GVM%zmom_8(i,j+1,k)-GVM%zmom_8(i,j+1,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
               Jzm  = (GVM%zmom_8(i  ,j,k)-GVM%zmom_8(i  ,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
               Bfdg2(i,j,k) =(Jzpi*fdg2_4(i,j+1,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
            enddo
         enddo
         do k = 2,Nk-1
            do j=1+pil_s-1, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  ski=(fdg2_4(i,j,k )*(half*(GVM%mc_Jyt_8(i,j-1,k )+GVM%mc_Jyt_8(i,j,k )))-&
                       fdg2_4(i,j,k-1)*(half*(GVM%mc_Jyt_8(i,j-1,k-1)+GVM%mc_Jyt_8(i,j,k-1))))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
                  skip=(fdg2_4(i,j+1,k )*(half*(GVM%mc_Jyt_8(i,j,k )+GVM%mc_Jyt_8(i,j+1,k )))-&
                        fdg2_4(i,j+1,k-1)*(half*(GVM%mc_Jyt_8(i,j,k-1)+GVM%mc_Jyt_8(i,j+1,k-1))))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
                  skpi=(fdg2_4(i,j,k+1)*(half*(GVM%mc_Jyt_8(i,j-1,k+1)+GVM%mc_Jyt_8(i,j,k+1)))-&
                        fdg2_4(i,j,k )*(half*(GVM%mc_Jyt_8(i,j-1,k )+GVM%mc_Jyt_8(i,j,k ))))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
                  skpip=(fdg2_4(i,j+1,k+1)*(half*(GVM%mc_Jyt_8(i,j,k+1)+GVM%mc_Jyt_8(i,j+1,k+1)))-&
                         fdg2_4(i,j+1,k )*(half*(GVM%mc_Jyt_8(i,j,k )+GVM%mc_Jyt_8(i,j+1,k ))))/(Ver_z_8%t(k+1)-Ver_z_8%t(k))
                  qkm = half*(ski+skip)
                  qkp = half*(skpi+skpip)
                  C = half*(qkp+qkm)
                  Jzpi= (GVM%zmom_8(i,j+1,k+1)-GVM%zmom_8(i,j+1,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
                  Jz  = (GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)) /(Ver_z_8%m(k+1)-Ver_z_8%m(k))
                  Bfdg1(i,j,k) = (Jzpi*fdg2_4(i,j+1,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j) - C
         	  Bfdg1(i,j,k) = dcoef*Bfdg1(i,j,k) * geomh_cyv_8(j)
                  Jzmi = (GVM%zmom_8(i,j+1,k)-GVM%zmom_8(i,j+1,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))
                  Jzm  = (GVM%zmom_8(i  ,j,k)-GVM%zmom_8(i  ,j,k-1))/(Ver_z_8%m(k)-Ver_z_8%m(k-1))

         	  Bfdg2(i,j,k) =(Jzpi*fdg2_4(i,j+1,k) - Jz*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
               enddo
            enddo
         enddo

! Apply divergence
         do k = 1, nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  add_v8(i,j,k) = ((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)))*&
                                   (Afdg1 (i,j,k)-Afdg1 (i-1,j,k))*geomh_invDXM_8(j)
               enddo
            enddo
         enddo
         do k = 1, nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  bdd_v8(i,j,k) = ((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)))*&
                     (Bfdg1 (i,j,k)-Bfdg1 (i,j-1,k))*geomh_invDYMv_8(j) /geomh_cy_8(j)
               enddo
            enddo
         enddo

!  flux
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               do k=2,Nk
                  Jz   = (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))
                  C2_8 = half*(GVM%mc_Jxt_8(i-1,j,k)+GVM%mc_Jxt_8(i,j,k))*&
                  fdg2_4(i,j,k)-half*(GVM%mc_Jxt_8(i-1,j,k-1)+GVM%mc_Jxt_8(i,j,k-1))* fdg2_4(i,j,k-1)
                  cdd1_v8(i,j,k) = (Ver_wp_8%m(k)*(half*(Afdg2(i-1,j,k)+Afdg2(i,j,k)))+ &
                          Ver_wm_8%m(k)*(half*(Afdg2(i-1,j,k-1)+Afdg2(i,j,k-1))))
! put zero if using stencil
                  c2_8=zero
                  cdd1_v8(i,j,k) = half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i,j,k))*Jz* dcoef*(cdd1_v8(i,j,k)&
                                             - C2_8 /(Ver_z_8%t(k)-Ver_z_8%t(k-1)))
                  C2_8 = half*(GVM%mc_Jyt_8(i,j-1,k)+GVM%mc_Jyt_8(i,j,k))*&
                  fdg2_4(i,j,k)-half*(GVM%mc_Jyt_8(i,j-1,k-1)+GVM%mc_Jyt_8(i,j,k-1))* fdg2_4(i,j,k-1)
                  cdd2_v8(i,j,k) = Ver_wp_8%m(k)*(half*(Bfdg2(i,j-1,k)+Bfdg2(i,j,k)))+ &
                                   Ver_wm_8%m(k)*(half*(Bfdg2(i,j-1,k-1)+Bfdg2(i,j,k-1)))
! put zero if using stencil
                  c2_8=zero
                  cdd2_v8(i,j,k)= half*(GVM%mc_Jy_8(i,j-1,k)+GVM%mc_Jy_8(i,j,k))* Jz*dcoef*geomh_cy_8(j)*&
                           (cdd2_v8(i,j,k)-C2_8 /(Ver_z_8%t(k)-Ver_z_8%t(k-1)))
               enddo
            enddo
         enddo

         do k=1,Nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  Jz=((Ver_z_8%m(k+1)-Ver_z_8%m(k))/(GVM%zmom_8(i,j,k+1  )-GVM%zmom_8(i,j,k)))
                  C1_8 = (cdd1_v8(i,j,k+1)-cdd1_v8(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
                  C2_8 = (cdd2_v8(i,j,k+1)-cdd2_v8(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k))
                  C2_8 = c2_8/geomh_cy_8(j)
! imprimer pour k=NK ce C
                  C=Jz*(C1_8+C2_8)
                  cflux_8(i,j,k)=C
                  if (k==NK) then
                     cflux_8(i,j,k)= (one-(ver_z_8%t(Nk)-ver_z_8%t(Nk-1))/(ver_z_8%t(Nk+1)-ver_z_8%t(Nk-1)))*&
                                    cflux_8(i,j,Nk-1)
                     C= cflux_8(i,j,k)
                  endif
                  cflux(i,j,k)= cflux_8(i,j,k)

!Field after diffusion on T-level K  on phii,j
                  F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8*( add_v8(i,j,k)+bdd_v8(i,j,k)-Cflux_8(i,j,k))
               enddo
            enddo
         enddo

         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               do k = 2 , Nk
                  F_sol1(i,j,k) = F_sol1(i,j,k) - W_th(i,j,k) * F_sol1(i,j,k-1)
               enddo
               F_sol1(i,j,Nk) = F_sol1(i,j,Nk) / b_th(i,j,Nk)
               do k = Nk-1, 1, -1
                  F_sol1(i,j,k) = (F_sol1(i,j,k) - c_th(i,j,k) * F_sol1(i,j,k+1)) / b_th(i,j,k)
               enddo
            enddo
         enddo

      ! Hybrid diffusion if hzd_hyb_bot >0
      if(hzd_hyb_bot >0) then
         do k = nk-hzd_hyb_bot+1, nk
            do j=1+pil_s-1, l_nj-pil_n+1
               do i=1+pil_w-1, l_ni-pil_e+1
                  F_sol1(i,j,k) = fdg2_4(i,j,k )
               enddo
            enddo
         enddo
      endif

      return
      end

