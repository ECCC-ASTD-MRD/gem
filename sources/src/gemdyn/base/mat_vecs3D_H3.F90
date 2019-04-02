!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
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

!**s/r mat_vecs3D - 3D_elliptic matrix_vector's computation for GEM_H
!
      subroutine mat_vecs3D_H3 ( F_Sol, F_Rhs, Minx, Maxx, Miny, Maxy, nil, njl, Nk)
      use geomh
      use gem_options
      use tdpack
      use glb_ld
      use cstv
      use ver
      use sol
      use opr
      use metric
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy,nil, njl, NK
      real*8, intent(out) :: F_Rhs(Minx:Maxx,Miny:Maxy,Nk)
      real*8, intent(in) :: F_Sol(Minx:Maxx,Miny:Maxy,Nk)
!author
!       Abdessamad Qaddouri -  2017
!
      integer j,i,k,halox,haloy
      real*8, parameter :: one=1.d0, zero=0.d0, half=0.5d0
      real*8  c_k(Minx:Maxx,Miny:Maxy,Nk),A_k(Minx:Maxx,Miny:Maxy,Nk), B_k(Minx:Maxx,Miny:Maxy,Nk)
      real*8  cdd_v8(Minx:Maxx,Miny:Maxy,Nk)
      integer  km, kp
      real*8  add_v8(Minx:Maxx,Miny:Maxy,Nk)
      real*8  bdd_v8(Minx:Maxx,Miny:Maxy,Nk)
      real   Afdg1(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real   Bfdg1(l_minx:l_maxx, l_miny:l_maxy,Nk)
      real   fdg2(l_minx:l_maxx, l_miny:l_maxy,Nk+1)
!     ---------------------------------------------------------------
!
!     ---------------------------------------------------------------
!
       F_Rhs= 0.0d0

      halox=1
      haloy=halox

      Afdg1 = .0
      Bfdg1 = .0
      A_k=0.0
      B_k=0.0
      C_k=0.0
      fdg2 =0.0
!
!
      do k = 1,Nk-1
         do j=1+sol_pil_s, njl-sol_pil_n
            do i=1+sol_pil_w, nil-sol_pil_e
               c_k(i,j,k) =  gama_8*((F_Sol(i,j,k+1)-F_Sol(i,j,k))*mc_iJz(i  ,j,k )&
                   -mu_8*half*(F_Sol(i,j,k+1)+F_Sol(i,j,k)))

            end do
         end do
      end do
!
!
      k=1
      do j=1+sol_pil_s, njl-sol_pil_n
         do i=1+sol_pil_w, nil-sol_pil_e
            cdd_v8(i,j,k)=(c_k(i,j,k))*Ver_idz_8%m(k)&
                  +(mc_Iz(i,j,k)-epsi_8)*(Ver_wp_8%m(k)*c_k(i,j,k)) -gg_8* F_Sol(i,j,k)
         end do
      end do

      do k = 2,Nk-1
         do j=1+sol_pil_s, njl-sol_pil_n
            do i=1+sol_pil_w, nil-sol_pil_e
               cdd_v8(i,j,k)=(c_k(i,j,k)-c_k(i,j,k-1))*Ver_idz_8%m(k)&
                     +(mc_Iz(i,j,k)-epsi_8)*(Ver_wp_8%m(k)*c_k(i,j,k)+Ver_wm_8%m(k)*c_k(i,j,k-1)) &
                     -gg_8* F_Sol(i,j,k)

            end do
         end do
      end do

      k=Nk
      do j=1+sol_pil_s, njl-sol_pil_n
         do i=1+sol_pil_w, nil-sol_pil_e
            cdd_v8(i,j,k)=(-c_k(i,j,k-1))*Ver_idz_8%m(k)&
               +(mc_Iz(i,j,k)-epsi_8)*(ver_wmA_8(k)*c_k(i,j,k-1)) -gg_8* F_Sol(i,j,k)

         end do
      end do
!

!$omp do

      do k = 1, nk
         fdg2(:,:,k) = .0d0
         do j=1+sol_pil_s, njl-sol_pil_n
            do i=1+sol_pil_w, nil-sol_pil_e
               fdg2(i,j,k)=F_Sol(i,j,k)

            end do
         end do
      end do
!$omp end do

!
      do j=1+sol_pil_s, njl-sol_pil_n
         do i=1+sol_pil_w, nil-sol_pil_e
               fdg2(i,j,Nk+1) = mc_alfas_H_8(i,j) * F_sol(i,j,NK)   &
                              - mc_betas_H_8(i,j) * F_sol(i,j,NK-1)

         end do
      end do
!

      call rpn_comm_xch_halo(fdg2,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )


      do k = 1,Nk
         km=max(k-1,1)
         kp=min(k+1,Nk+1)

         do j=1+sol_pil_s, njl-sol_pil_n
            do i=1+sol_pil_w, l_niu-sol_pil_e
!
               A_k(i,j,k) =(fdg2(i+1,j,k) - fdg2(i,j,k) ) * geomh_invDX_8(j)   - mc_Jx(i,j,k) * (  &
                          Ver_wp_8%m(k)*half*( (fdg2(i+1,j,kp)-fdg2(i+1,j,k ))*mc_iJz(i+1,j,k )   &
                                           +(fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*mc_iJz(i  ,j,k ) ) &
                         +Ver_wm_8%m(k)*half*( (fdg2(i+1,j,k  )-fdg2(i+1,j,km))*mc_iJz(i+1,j,km)   &
                                           +(fdg2(i  ,j,k  )-fdg2(i  ,j,km))*mc_iJz(i  ,j,km) ) )
               Afdg1(i,j,k)= real(A_k(i,j,k))

            end do
         end do
!
         do j=1+sol_pil_s, l_njv-sol_pil_n
            do i=1+sol_pil_w, nil-sol_pil_e

               B_k(i,j,k) =(fdg2(i,j+1,k) - fdg2(i,j,k) ) * geomh_invDYMv_8(j)  - mc_Jy(i,j,k) * ( &
                          Ver_wp_8%m(k)*half*( (fdg2(i,j+1,kp)-fdg2(i,j+1,k ))*mc_iJz(i,j+1,k )   &
                                           +(fdg2(i,j  ,kp)-fdg2(i,j  ,k ))*mc_iJz(i,j  ,k ) ) &
                         +Ver_wm_8%m(k)*half*( (fdg2(i,j+1,k  )-fdg2(i,j+1,km))*mc_iJz(i,j+1,km)   &
                                           +(fdg2(i,j  ,k  )-fdg2(i,j  ,km))*mc_iJz(i,j  ,km) ) )
               Bfdg1(i,j,k)=real(B_k(i,j,k))

            end do
         end do
      end do


      call rpn_comm_xch_halo( Afdg1 , l_minx,l_maxx,l_miny,l_maxy,l_niu,l_nj ,nk  , &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( Bfdg1 , l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_njv,nk  , &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )


!
      do k = 1, nk
         do j=1+sol_pil_s, njl-sol_pil_n
            do i=1+sol_pil_w, nil-sol_pil_e

               add_v8(i,j,k) =   (Afdg1 (i,j,k)-Afdg1 (i-1,j,k))*geomh_invDXM_8(j) &
                         + half * ( mc_Ix(i,j,k)*(Afdg1(i,j,k)+Afdg1(i-1,j,k)))
!
               bdd_v8(i,j,k) = (Bfdg1 (i,j,k)*geomh_cyM_8(j)-Bfdg1 (i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                         + half * (mc_Iy(i,j,k)*(Bfdg1(i,j,k)+Bfdg1(i,j-1,k)) )


            end do
         end do
      end do

!$omp do
      do k=1, NK
         do j=1+sol_pil_s, njl-sol_pil_n
            do i=1+sol_pil_w, nil-sol_pil_e
               F_Rhs(i,j,k)=Cstv_hco0_8*(add_v8(i,j,k)+bdd_v8(i,j,k)+cdd_v8(i,j,k))
            end do
         end do
      end do
!$omp end do


!     ---------------------------------------------------------------
!
      return
      end

