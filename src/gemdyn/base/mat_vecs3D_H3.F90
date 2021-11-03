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
      subroutine mat_vecs3D_H3 (F_Sol, F_Rhs, Minx, Maxx, Miny, Maxy, nil, njl, Nk)
      use geomh
      use gem_options
      use HORgrid_options
      use tdpack
      use glb_ld
      use cstv
      use ver
      use sol
      use opr
      use ptopo
      use metric
      use lam_options
      use dyn_fisl_options
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy,nil, njl, NK
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: F_Rhs
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: F_Sol
!author
!       Abdessamad Qaddouri -  initilal version 2019
!       Abdessamad Qaddouri,   May 2020 (opentop)
!       Rabah Aider,           July 2021 (optimization)
!
      logical, save :: first_time = .true.
      integer j, i, id, k, halox, haloy
      real(kind=REAL64), parameter :: one=1.d0, zero=0.d0, half=0.5d0
      real(kind=REAL64), dimension(l_minx:l_maxx, l_miny:l_maxy,Nk,15) :: A1, B1 ,A2, B2,C1
      real(kind=REAL64), dimension(:,:,:,:), allocatable, save :: stencil
      real, dimension(l_minx:l_maxx, l_miny:l_maxy,Nk+1) :: fdg2
      integer  km, kp,k0,k0t
      integer sol_pil_w_ext, sol_pil_e_ext, sol_pil_s_ext, sol_pil_n_ext
!
!     ---------------------------------------------------------------
!
!!
!(i,j,k)   =>stencil(:,:,:,1)   (i,j-1,k)  =>stencil(:,:,:,10)
!(i-1,j,k) =>stencil(:,:,:,2)   (i,j+1,k)  =>stencil(:,:,:,11)
!(i+1,j,k) =>stencil(:,:,:,3)   (i,j-1,k-1)=>stencil(:,:,:,12)
!(i,j,k-1) =>stencil(:,:,:,4)   (i,j-1,k+1)=>stencil(:,:,:,13)
!(i,j,k+1) =>stencil(:,:,:,5)   (i,j+1,k-1)=>stencil(:,:,:,14)
!(i-1,j,k-1)=>stencil(:,:,:,6)  (i,j+1,k+1)=>stencil(:,:,:,15)
!(i-1,j,k+1)=>stencil(:,:,:,7)
!(i+1,j,k-1)=>stencil(:,:,:,8)
!(i+1,j,k+1)=>stencil(:,:,:,9)

      k0=1+Lam_gbpil_T
      k0t=k0
      if (Schm_opentop_L) k0t=k0-1
      F_Rhs = 0.0d0

      halox=1
      haloy=halox

      fdg2(:,:,:) = 0.

      sol_pil_s_ext=sol_pil_s-1
      sol_pil_n_ext= sol_pil_n
      sol_pil_w_ext=sol_pil_w-1
      sol_pil_e_ext=sol_pil_e

      if(.not.Grd_yinyang_L) then
         if (l_west) sol_pil_w_ext=sol_pil_w
         if (l_east) sol_pil_e_ext=sol_pil_e+1
         if (l_south) sol_pil_s_ext=sol_pil_s
         if (l_north) sol_pil_n_ext= sol_pil_n+1
      endif

      do k = k0, nk
         fdg2(:,:,k) = 0.
         do j=1+sol_pil_s, njl-sol_pil_n
            do i=1+sol_pil_w, nil-sol_pil_e
               fdg2(i,j,k)=F_Sol(i,j,k)
            end do
         end do
      end do

      do j=1+sol_pil_s, njl-sol_pil_n
         do i=1+sol_pil_w, nil-sol_pil_e
               fdg2(i,j,Nk+1) = mc_alfas_H_8(i,j) * F_sol(i,j,NK)   &
                              - mc_betas_H_8(i,j) * F_sol(i,j,NK-1)
         end do
      end do
      if (Schm_opentop_L) then
      do j=1+sol_pil_s, njl-sol_pil_n
         do i=1+sol_pil_w, nil-sol_pil_e
           fdg2(i,j,k0t) =  mc_alfat_8(i,j)* F_sol(i,j,k0)
         end do
      end do
      endif

      if ( Grd_yinyang_L) then
          call yyg_xchng (fdg2, l_minx,l_maxx,l_miny,l_maxy, &
                               l_ni,l_nj, nk+1, .false., 'CUBIC', .true.)
      else
          call rpn_comm_xch_halo(fdg2,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,nk+1, &
                             G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      endif

! Compute stencils
      if(first_time) then
         allocate (stencil(1+sol_pil_w:nil-sol_pil_e,1+sol_pil_s:njl-sol_pil_n,Nk,15))
         do k = 1,Nk
            do i=1,15
               A1(:,:,k,i)=0.d0
               A2(:,:,k,i)=0.d0
               B1(:,:,k,i)=0.d0
               B2(:,:,k,i)=0.d0
               C1(:,:,k,i)=0.d0
               stencil(:,:,k,i)=0.d0
            enddo
         enddo
!
         k=k0
         do j=1+sol_pil_s, njl-sol_pil_n
            do i=1+sol_pil_w, nil-sol_pil_e
               C1(i,j,k,1)=-gama_8*(mc_iJz_8(i,j,k ) &
                            + mu_8*half)*(Ver_idz_8%m(k)+(mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k)) - gg_8
               C1(i,j,k,5)= gama_8*(mc_iJz_8(i,j,k ) &
                            - mu_8*half)*(Ver_idz_8%m(k)+(mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
               if (Schm_opentop_L) then
                   C1(i,j,k,1)=-gama_8*(mc_iJz_8(i,j,k ) + mc_iJz_8(i,j,k-1) )*Ver_idz_8%m(k) &
                              +(mc_Iz_8(i,j,k)-epsi_8)*gama_8*( Ver_wm_8%m(k)*(mc_iJz_8(i,j,k-1) -mu_8*half) &
                                                               -Ver_wp_8%m(k)*(mc_iJz_8(i,j,k )  +mu_8*half) ) - gg_8
                   C1(i,j,k,4)= gama_8*(mc_iJz_8(i,j,k-1) +  mu_8*half)*(Ver_idz_8%m(k) &
                               - (mc_Iz_8(i,j,k)-epsi_8)*Ver_wm_8%m(k))
                   C1(i,j,k,5)= gama_8*(mc_iJz_8(i,j,k ) -mu_8*half)*(Ver_idz_8%m(k) + &
                                        (mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
               endif
            end do
         end do

         do k = k0+1,Nk
            do j=1+sol_pil_s, njl-sol_pil_n
               do i=1+sol_pil_w, nil-sol_pil_e
                   C1(i,j,k,1)=-gama_8*(mc_iJz_8(i,j,k ) + mc_iJz_8(i,j,k-1) )*Ver_idz_8%m(k) &
                              +(mc_Iz_8(i,j,k)-epsi_8)*gama_8*( Ver_wm_8%m(k)*(mc_iJz_8(i,j,k-1) -mu_8*half) &
                                                               -Ver_wp_8%m(k)*(mc_iJz_8(i,j,k )  +mu_8*half) ) - gg_8
                   C1(i,j,k,4)= gama_8*(mc_iJz_8(i,j,k-1) +  mu_8*half)*(Ver_idz_8%m(k) &
                               - (mc_Iz_8(i,j,k)-epsi_8)*Ver_wm_8%m(k))
                   C1(i,j,k,5)= gama_8*(mc_iJz_8(i,j,k ) -mu_8*half)*(Ver_idz_8%m(k) + &
                                        (mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
               end do
            end do
         end do

         do k = k0,Nk
            km=max(k-1,1)
            kp=k+1
            do j=1+sol_pil_s_ext, l_nj-sol_pil_n
               do i=1+sol_pil_w_ext, l_ni-sol_pil_e_ext
                  A1(i,j,k,1)= -geomh_invDX_8(j) + half*mc_Jx_8(i,j,k)*   &
                                (Ver_wp_8%m(k)*mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i,j,km))
                  A2(i,j,k,1)=  geomh_invDX_8(j) + half*mc_Jx_8(i-1,j,k)* &
                                (Ver_wp_8%m(k)*mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i,j,km))
                  A2(i,j,k,2)= -geomh_invDX_8(j) + half*mc_Jx_8(i-1,j,k)* &
                                (Ver_wp_8%m(k)*mc_iJz_8(i-1,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i-1,j,km))
                  A1(i,j,k,3)=  geomh_invDX_8(j) + half*mc_Jx_8(i,j,k)*   &
                                (Ver_wp_8%m(k)*mc_iJz_8(i+1,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i+1,j,km))
                  A1(i,j,k,4)= half*mc_Jx_8(i,j,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j,km)
                  A2(i,j,k,4)= half*mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j,km)
                  A1(i,j,k,5)=-half*mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j,k)
                  A2(i,j,k,5)=-half*mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j,k)
                  A2(i,j,k,6)= half*mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*mc_iJz_8(i-1,j,km)
                  A2(i,j,k,7)=-half*mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*mc_iJz_8(i-1,j,k)
                  A1(i,j,k,8)=half*mc_Jx_8(i,j,k)*Ver_wm_8%m(k)*mc_iJz_8(i+1,j,km)
                  A1(i,j,k,9)=-half*mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*mc_iJz_8(i+1,j,k)
               end do
            end do

            do j=1+sol_pil_s_ext, l_nj-sol_pil_n_ext
               do i=1+sol_pil_w_ext, l_ni-sol_pil_e
                  B1(i,j,k,1)= -geomh_invDYMv_8(j) + half*mc_Jy_8(i,j,k)*     &
                                (Ver_wp_8%m(k)*mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i,j,km))
                  B2(i,j,k,1)=  geomh_invDYMv_8(j-1) + half*mc_Jy_8(i,j-1,k)* &
                                (Ver_wp_8%m(k)*mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i,j,km))
                  B1(i,j,k,4)= half*mc_Jy_8(i,j,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j,km)
                  B2(i,j,k,4)= half*mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j,km)
                  B1(i,j,k,5)= -half*mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j,k)
                  B2(i,j,k,5)= -half*mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j,k)
                  B2(i,j,k,10)=-geomh_invDYMv_8(j-1) + half*mc_Jy_8(i,j-1,k)* &
                               (Ver_wp_8%m(k)*mc_iJz_8(i,j-1,k ) - Ver_wm_8%m(k)*mc_iJz_8(i,j-1,km))
                  B1(i,j,k,11)= geomh_invDYMv_8(j) + half*mc_Jy_8(i,j,k)*  &
                                (Ver_wp_8%m(k)*mc_iJz_8(i,j+1,k ) - Ver_wm_8%m(k)*mc_iJz_8(i,j+1,km))
                  B2(i,j,k,12) =half*mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j-1,km)
                  B2(i,j,k,13)=-half*mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j-1,k)
                  B1(i,j,k,14)= half*mc_Jy_8(i,j,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j+1,km)
                  B1(i,j,k,15)= -half*mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j+1,k)
               end do
            end do

            if(.not.Grd_yinyang_L) then
               do j=1+sol_pil_s, l_nj-sol_pil_n
                  do i=1+sol_pil_w, l_ni-sol_pil_e
                     A2(i,j,k,1)=  geomh_invDX_8(j) + half*mc_Jx_8(i-1,j,k)* &
                                   (Ver_wp_8%m(k)*mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i,j,km))
                     A2(i,j,k,2)= -geomh_invDX_8(j) + half*mc_Jx_8(i-1,j,k)* &
                                   (Ver_wp_8%m(k)*mc_iJz_8(i-1,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i-1,j,km))
                     A2(i,j,k,4)= half*mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j,km)
                     A2(i,j,k,5)= -half*mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j,k)
                     A2(i,j,k,6)= half*mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*mc_iJz_8(i-1,j,km)
                     A2(i,j,k,7)= -half*mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*mc_iJz_8(i-1,j,k)
                     B2(i,j,k,1)=  geomh_invDYMv_8(j-1)  + half*mc_Jy_8(i,j-1,k)* &
                                   (Ver_wp_8%m(k)*mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*mc_iJz_8(i,j,km))
                     B2(i,j,k,4)= half*mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j,km)
                     B2(i,j,k,5)= -half*mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j,k)
                     B2(i,j,k,10)=-geomh_invDYMv_8(j-1) + half*mc_Jy_8(i,j-1,k)* &
                                   (Ver_wp_8%m(k)*mc_iJz_8(i,j-1,k ) - Ver_wm_8%m(k)*mc_iJz_8(i,j-1,km))
                     B2(i,j,k,12) =half*mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*mc_iJz_8(i,j-1,km)
                     B2(i,j,k,13)=-half*mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*mc_iJz_8(i,j-1,k)

                     if(l_west) then
                        A2(1+sol_pil_w,j,k,1) = 0.d0
                        A2(1+sol_pil_w,j,k,2) = 0.d0
                        A2(1+sol_pil_w,j,k,4) = 0.d0
                        A2(1+sol_pil_w,j,k,5) = 0.d0
                        A2(1+sol_pil_w,j,k,6) = 0.d0
                        A2(1+sol_pil_w,j,k,7) = 0.d0
                     endif

                     if(l_south) then
                        B2(i,1+sol_pil_s,k,1) = 0.d0
                        B2(i,1+sol_pil_s,k,4) = 0.d0
                        B2(i,1+sol_pil_s,k,5) = 0.d0
                        B2(i,1+sol_pil_s,k,10)= 0.d0
                        B2(i,1+sol_pil_s,k,12)= 0.d0
                        B2(i,1+sol_pil_s,k,13)= 0.d0
                     endif
                  enddo
               enddo
            endif
         end do

         do k =k0, nk
            do j=1+sol_pil_s, njl-sol_pil_n
               do i=1+sol_pil_w, nil-sol_pil_e
                  do id=1,15
                  stencil (i,j,k,id) =Cstv_hco0_8* ( (A1 (i,j,k,id)-A2 (i,j,k,id))*geomh_invDXM_8(j) &
                            + half * ( mc_Ix_8(i,j,k)*(A1(i,j,k,id)+A2(i,j,k,id)))    &
                            + (B1 (i,j,k,id)*geomh_cyM_8(j)-B2 (i,j,k,id)*&
                              geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                            + half * (mc_Iy_8(i,j,k)*(B1(i,j,k,id)+B2(i,j,k,id)) ) &
                            + C1(i,j,k,id) )
                  enddo
               end do
            end do
         end do

         first_time = .false.

      endif

         do k=k0,NK
            km=max(k-1,1)
            kp=k+1
            do j=1+sol_pil_s, njl-sol_pil_n
               do i=1+sol_pil_w, nil-sol_pil_e
                  F_Rhs(i,j,k)=  stencil (i,j,k,1)*fdg2(i,j,k) &
                              +  stencil (i,j,k,2)*fdg2(i-1,j,k) &
                              +  stencil (i,j,k,3)*fdg2(i+1,j,k)&
                              +  stencil (i,j,k,4)*fdg2(i,j,km) &
                              +  stencil (i,j,k,5)*fdg2(i,j,kp) &
                              +  stencil (i,j,k,6)*fdg2(i-1,j,km) &
                              +  stencil (i,j,k,7)*fdg2(i-1,j,kp) &
                              +  stencil (i,j,k,8)*fdg2(i+1,j,km) &
                              +  stencil (i,j,k,9)*fdg2(i+1,j,kp) &
                              +  stencil (i,j,k,10)*fdg2(i,j-1,k) &
                              +  stencil (i,j,k,11)*fdg2(i,j+1,k) &
                              +  stencil (i,j,k,12)*fdg2(i,j-1,km) &
                              +  stencil (i,j,k,13)*fdg2(i,j-1,kp) &
                              +  stencil (i,j,k,14)*fdg2(i,j+1,km) &
                              +  stencil (i,j,k,15)*fdg2(i,j+1,kp)
                end do
            end do
         end do

      return

      end

