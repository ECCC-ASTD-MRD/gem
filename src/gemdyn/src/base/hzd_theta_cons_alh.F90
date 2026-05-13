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
      subroutine  hzd_theta_cons_alh ( F_Sol1,HzdlnR, Minx, Maxx, Miny, Maxy,Nk)
      use gem_options
      use gmm_vt1
      use gmm_hzd
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
      integer kd0 , k00, k01
      real(kind=REAL64)    one,half,zero
      parameter( one=1.0d0,half=0.5d0,zero=0.d0)
      real   fdg2_4(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real  crit_coef, base_coefT,cdelta2,cdelta,Creal
      real F_s (l_minx:l_maxx, l_miny:l_maxy,Nk),cs
      real cflux(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real dcoef 

      real(kind=REAL64)   Afdg1(l_minx:l_maxx, l_miny:l_maxy,Nk),Afdg2(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64)   Bfdg1(l_minx:l_maxx, l_miny:l_maxy,Nk),Bfdg2(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64)   C1,C2 
      real(kind=REAL64) Jz, qkm,qkp

      real(kind=REAL64) C1_8,C2_8,C,ski,skpi,skip,skpip
      real(kind=REAL64) bdd_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) add_v8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) cdd1_v8(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64) cdd2_v8(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
      real(kind=REAL64), dimension (:,:,:,:), allocatable :: stencilV
      real(kind=REAL64) a(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) b(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) d(l_minx:l_maxx, l_miny:l_maxy,Nk),W
      real(kind=REAL64) Cflux_8(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real(kind=REAL64) C1flux_8(l_minx:l_maxx, l_miny:l_maxy,Nk),dsten(l_minx:l_maxx, l_miny:l_maxy)
      real(kind=REAL64)  ztht_8(l_minx:l_maxx, l_miny:l_maxy,0:Nk+1),Jxx,Jyy


! kd0 given by user      
      kd0=hzd_hyb_top
      k00=kd0
      k01=kd0+1
      if (kd0.ne.1) k01=kd0

      if (Hzd_pwr_z==2)  then
          dcoef= 0.25*HzdlnR*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      else
          dcoef= 0.25*sqrt(HzdlnR)*(Dcst_rayt_8*geomh_hy_8)**2/Cstv_dt_8
      endif

! Apply Horizontal diffusion along z

         Afdg1 = .0d0
         Bfdg1 = .0d0
         add_v8 =0.0d0
         bdd_v8=0.0d0
         cdd1_v8=0.0d0
         cdd2_v8=0.0d0
         fdg2_4 =0.0
         a=zero
         b=zero
         d=zero
         cflux=0.0
         Cflux_8=zero

!Field  before diffusion on T-level K  on phii,j
!$omp do
         do k = 1, nk
            do j=1+pil_s-1, l_nj-pil_n+1
               do i=1+pil_w-1, l_ni-pil_e+1
                  fdg2_4(i,j,k )=F_Sol1(i,j,k)
               enddo
            enddo
         enddo
!$omp enddo

         call rpn_comm_xch_halo(fdg2_4,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,Nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         k=k00
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
               C1 = dcoef*((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j))
               C2 = ((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j))
               Afdg1(i,j,k)= half*(air_dens(i,j,k)+air_dens(i+1,j,k))*C1
               Afdg2(i,j,k)= half*(air_dens(i,j,k)+air_dens(i+1,j,k))*C2
            enddo
         enddo
!
         k= NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w-1, l_ni-pil_e
               C1_8=half*(fdg2_4(i,j,k)*Jxt(i,j,k,1)-fdg2_4(i,j,k-1)*Jxt(i,j,k,2) + &
                         fdg2_4(i+1,j,k )*Jxt(i,j,k,3)-fdg2_4(i+1,j,k-1)*Jxt(i,j,k,4))
               C1 = dcoef* ((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j) - half*C1_8)
               C2 = ((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j)) 
               Afdg1(i,j,k)= half*(air_dens(i,j,k)+air_dens(i+1,j,k))*C1
               Afdg2(i,j,k)= half*(air_dens(i,j,k)+air_dens(i+1,j,k))*C2
            enddo
         enddo
!
         do k = k01,Nk-1
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w-1, l_ni-pil_e
               qkm =  half* (fdg2_4(i,j,k    )*Jxt(i,j,k,1) - fdg2_4(i,j,k-1  )*Jxt(i,j,k,2) +&
                             fdg2_4(i+1,j,k  )*Jxt(i,j,k,3) - fdg2_4(i+1,j,k-1)*Jxt(i,j,k,4))
               qkp = half*(fdg2_4(i,j,k+1  )*Jxt(i,j,k,5) - fdg2_4(i,j,k )*Jxt(i,j,k,6) +&
                           fdg2_4(i+1,j,k+1)*Jxt(i,j,k,7) - fdg2_4(i+1,j,k )*Jxt(i,j,k,8))
               C=   half*(qkp+qkm)

               C1 = dcoef*((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j) - C)
               C2 =((Jzpt(i,j,k,1)*fdg2_4(i+1,j,k) - Jzt(i,j,k,1)*fdg2_4(i,j,k) ) * geomh_invDX_8(j))
! conservative
               Afdg1(i,j,k)= half*(air_dens(i,j,k)+air_dens(i+1,j,k))*C1
               Afdg2(i,j,k)= half*(air_dens(i,j,k)+air_dens(i+1,j,k))*C2
               enddo
            enddo
         enddo
         k=k00
         do j=1+pil_s-1, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               C1 = dcoef*(Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k)) * geomh_invDYMv_8(j) & 
                         * geomh_cyv_8(j)
               C2 =(Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
! conservative
               Bfdg1(i,j,k)= half*(air_dens(i,j,k)+air_dens(i,j+1,k))* C1
               Bfdg2(i,j,k)= half*(air_dens(i,j,k)+air_dens(i,j+1,k))* C2

            enddo
         enddo

         k= NK
         do j=1+pil_s-1, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               C1_8=half*(fdg2_4(i,j,k  )*Jyt(i,j,k,1)-fdg2_4(i,j,k-1)*Jyt(i,j,k,2) +&
                          fdg2_4(i,j+1,  k)*Jyt(i,j,k,3)-fdg2_4(i,j+1,k-1)*Jyt(i,j,k,4))
               C1 =dcoef*((Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j) - half*C1_8)* &
                             geomh_cyv_8(j)
               C2 =(Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
! conservative
               Bfdg1(i,j,k)= half*(air_dens(i,j,k)+air_dens(i,j+1,k))*C1
               Bfdg2(i,j,k)= half*(air_dens(i,j,k)+air_dens(i,j+1,k))*C2

            enddo
         enddo
         do k = k01,Nk-1
            do j=1+pil_s-1, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  qkm=half*(fdg2_4(i,j,k   )*Jyt(i,j,k,1)-fdg2_4(i,j,k-1)*Jyt(i,j,k,2  )+&
                            fdg2_4(i,j+1,k )*Jyt(i,j,k,3)-fdg2_4(i,j+1,k-1)*Jyt(i,j,k,4))
                  qkp=half*(fdg2_4(i,j,k+1)*Jyt(i,j,k,5)-fdg2_4(i,j,k )*Jyt(i,j,k,6)+&
                            fdg2_4(i,j+1,k+1)*Jyt(i,j,k,7)-fdg2_4(i,j+1,k )*Jyt(i,j,k,8))
                  C = half*(qkp+qkm)

                  C1 =dcoef* ((Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j) - C) &
                        * geomh_cyv_8(j)
                  C2 =(Jzpt(i,j,k,2)*fdg2_4(i,j+1,k) - Jzt(i,j,k,2)*fdg2_4(i,j,k) ) * geomh_invDYMv_8(j)
! conservative
                  Bfdg1(i,j,k)= half*(air_dens(i,j,k)+air_dens(i,j+1,k))*C1
                  Bfdg2(i,j,k)= half*(air_dens(i,j,k)+air_dens(i,j+1,k))*C2
               enddo
            enddo
         enddo

! Apply divergence
         do k = k00, nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  add_v8(i,j,k) = Jm(i,j,k,1)*(Afdg1 (i,j,k)-Afdg1 (i-1,j,k))
                  bdd_v8(i,j,k) = Jm(i,j,k,2)*(Bfdg1 (i,j,k)-Bfdg1 (i,j-1,k))
               enddo
            enddo
         enddo
!  flux
      do k=k00+1,Nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
                  !Jzz(i,j,k,1)  = (Ver_z_8%t(k)-Ver_z_8%t(k-1))/(ztht_8(i  ,j,k)-ztht_8(i  ,j,k-1))
                  cdd1_v8(i,j,k) = (Ver_wp_8%m(k)*(half*(Afdg2(i-1,j,k)+Afdg2(i,j,k)))+ &
                          Ver_wm_8%m(k)*(half*(Afdg2(i-1,j,k-1)+Afdg2(i,j,k-1))))
! put zero if using stencil
                  cdd1_v8(i,j,k) = half*(GVM%mc_Jx_8(i-1,j,k)+GVM%mc_Jx_8(i,j,k))*Jzz(i,j,k,1)* dcoef*cdd1_v8(i,j,k)
                  cdd2_v8(i,j,k) = Ver_wp_8%m(k)*(half*(Bfdg2(i,j-1,k)+Bfdg2(i,j,k)))+ &
                                   Ver_wm_8%m(k)*(half*(Bfdg2(i,j-1,k-1)+Bfdg2(i,j,k-1)))
! put zero if using stencil
                  cdd2_v8(i,j,k)= half*(GVM%mc_Jy_8(i,j-1,k)+GVM%mc_Jy_8(i,j,k))* Jzz(i,j,k,1)*dcoef*geomh_cy_8(j)*&
                                  cdd2_v8(i,j,k)
               enddo
            enddo
         enddo

         do k=k00,Nk-1
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  cflux_8(i,j,k)=Jzz(i,j,k,2)*Ver_idz_8%t(k)* &
                  ((cdd1_v8(i,j,k+1)-cdd1_v8(i,j,k)) + &
                    (cdd2_v8(i,j,k+1)-cdd2_v8(i,j,k))*geomh_invcy_8(j) )
                  F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8/air_dens(i,j,k)*( add_v8(i,j,k)+bdd_v8(i,j,k)-Cflux_8(i,j,k))
               enddo
            enddo
         enddo
         k=Nk
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               cflux_8(i,j,k)= (one-(ver_z_8%t(Nk)-ver_z_8%t(Nk-1))/(ver_z_8%t(Nk+1)-ver_z_8%t(Nk-1)))*&
                                    cflux_8(i,j,Nk-1)
               F_sol1(i,j,k)= F_sol1(i,j,k) + Cstv_dt_8/air_dens(i,j,k)*( add_v8(i,j,k)+bdd_v8(i,j,k)-Cflux_8(i,j,k))
            enddo
         enddo

         allocate(stencilV(l_minx:l_maxx, l_miny:l_maxy,3,Nk))
         stencilV=zero

! conservation 
      do k=k00+1,NK
         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               stencilV(i,j,1,k)= (air_dens_m(i,j,k-1)*stencil_V1(i,j,k) +air_dens_m(i,j,k)* stencil_V2(i,j,k))
               stencilV(i,j,2,k)= air_dens_m(i,j,k-1)*stencil_V(i,j,2,k)
               stencilV(i,j,3,k)= air_dens_m(i,j,k  )*stencil_V(i,j,3,k)
            enddo
         enddo
      enddo

       k=NK
        do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
               stencilV(i,j,2,k )=zfact  * stencilV(i,j,2,k-1 )
               stencilV(i,j,1,k )=zfact  *(stencilV(i,j,1,k-1 )+stencilV(i,j,3,k-1 ))
            enddo
        enddo

         do k=k00+1,Nk
            do j=1+pil_s, l_nj-pil_n
               do i=1+pil_w, l_ni-pil_e
                  a(i,j,k)=-Cstv_dt_8*stencilV(i,j,2,k) /air_dens(i,j,k) 
                  b(i,j,k)=one-Cstv_dt_8*stencilV(i,j,1,k) /air_dens(i,j,k) 
                  d(i,j,k)=-Cstv_dt_8*stencilV(i,j,3,k) /air_dens(i,j,k) 
               enddo
            enddo
         enddo

         deallocate (stencilV)

         do j=1+pil_s, l_nj-pil_n
            do i=1+pil_w, l_ni-pil_e
                  b(i,j,k00)=one 
               do k = k00+1 , Nk
!                  b(i,j,k00)=one
                  W = a(i,j,k) / b(i,j,k - 1)
                  b(i,j,k) = b(i,j,k) - W * d(i,j,k - 1)
                  F_sol1(i,j,k) = F_sol1(i,j,k) - W * F_sol1(i,j,k- 1)
               enddo
               F_sol1(i,j,Nk) = F_sol1(i,j,Nk) / b(i,j,Nk)
               F_sol1(i,j,1) = (F_sol1(i,j,1) - d(i,j,1) * F_sol1(i,j,2)) 
               do k = Nk-1, k00+1, -1
                  F_sol1(i,j,k) = (F_sol1(i,j,k) - d(i,j,k) * F_sol1(i,j,k + 1)) / b(i,j,k)
               enddo
            enddo
         enddo

      ! Hybrid diffusion if hzd_hyb_th_nk >0
      if(hzd_hyb_bot >0) then
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


      return
      end

