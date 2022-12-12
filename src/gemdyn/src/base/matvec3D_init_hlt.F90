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

!** matvec3D_init compute Sol_stencils for Matrix-vector product subroutines (P & H coordinates)
!
      subroutine matvec3D_init_hlt()
      use cstv
      use geomh
      use gem_options
      use HORgrid_options
      use glb_ld
      use ldnh
      use dynkernel_options
      use opr
      use sol_mem
      use ver
      use metric
      use mem_tstp
      use lam_options
      use dyn_fisl_options
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      integer j, jj, i, ii, id, k
      real(kind=REAL64)  :: di_8
      real(kind=REAL64)  :: xxx, yyy
      real(kind=REAL64), parameter :: one=1.d0, zero=0.d0, half=0.5d0
      integer, parameter :: IDX_POINT=1, IDX_WEST=2, IDX_EAST=3, IDX_NORTH=4, IDX_SOUTH=5, IDX_TOP=6, IDX_BOTTOM=7

      integer  km, kp,k0,k0t
      integer sol_pil_w_ext, sol_pil_e_ext, sol_pil_s_ext, sol_pil_n_ext
!
!     ---------------------------------------------------------------
!
      if (.not. FISLH_LHS_metric_L ) then

!$omp single
         allocate (Sol_stencilp_8(1+sol_pil_w:l_ni-sol_pil_e, 1+sol_pil_s:l_nj-sol_pil_n, 7, l_nk))

         xxx = - Cstv_hco2_8
         yyy = - Cstv_hco1_8

         do k=1, l_nk
            do j=1+sol_pil_s, l_nj-sol_pil_n
               jj=j+l_j0-1
               di_8 = Opr_opsyp0_8(G_nj+jj) * geomh_invcy2_8(j)
               do i=1+sol_pil_w, l_ni-sol_pil_e
                  ii=i+l_i0-1

                  Sol_stencilp_8(i,j,IDX_POINT,k) = Cstv_hco0_8 * (Cstv_hco3_8*Opr_opszp2_8(G_nk+k) + Cstv_hco3_8*Opr_opszpl_8(G_nk+k) &
                              + xxx * Opr_opszpm_8(G_nk+k) + yyy * Opr_opszp0_8(G_nk+k)) &
                              + Opr_opszp0_8(G_nk+k) * (Opr_opsxp2_8(G_ni+ii) * di_8     &
                              + Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(G_nj+jj))           &
                              / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

                  Sol_stencilp_8(i,j,IDX_WEST,k) = Opr_opsxp2_8(ii) * Opr_opszp0_8(G_nk+k) * Opr_opsyp0_8(G_nj+jj) &
                              / (cos( G_yg_8 (jj) )**2) / (Opr_opsxp0_8(G_ni+ii)*Opr_opsyp0_8(G_nj+jj))

                  Sol_stencilp_8(i,j,IDX_EAST,k) = Opr_opsxp2_8(2*G_ni+ii) * Opr_opszp0_8(G_nk+k) * Opr_opsyp0_8(G_nj+jj) &
                             / (cos( G_yg_8 (jj) )**2) / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

                  Sol_stencilp_8(i,j,IDX_SOUTH,k) = Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(jj) * Opr_opszp0_8(G_nk+k) &
                             / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

                  Sol_stencilp_8(i,j,IDX_NORTH,k) = Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(2*G_nj+jj) * Opr_opszp0_8(G_nk+k) &
                              / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

                  Sol_stencilp_8(i,j,IDX_TOP,k) = Cstv_hco0_8 * (Cstv_hco3_8*Opr_opszp2_8(k) + Cstv_hco3_8*Opr_opszpl_8(k) + xxx * Opr_opszpm_8(k))

                  Sol_stencilp_8(i,j,IDX_BOTTOM,k) = Cstv_hco0_8 * (Cstv_hco3_8*Opr_opszp2_8(2*G_nk+k) + Cstv_hco3_8*Opr_opszpl_8(2*G_nk+k) + xxx * Opr_opszpm_8(2*G_nk+k))

               end do
            end do
         end do
!$omp end single

      else

         k0=1+Lam_gbpil_T
         k0t=k0
         if (Schm_opentop_L) k0t=k0-1

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

!$omp do
         do id=1,15
            do k = 1,l_nk
               do j = l_miny,l_maxy
                  do i = l_minx,l_maxx
                     A1(i,j,k,id)=0.d0
                     A2(i,j,k,id)=0.d0
                     B1(i,j,k,id)=0.d0
                     B2(i,j,k,id)=0.d0
                     C1(i,j,k,id)=0.d0
                  enddo
               enddo
            enddo
         enddo
!$omp end do

         k=k0
!$omp do
         do j=1+sol_pil_s, l_nj-sol_pil_n
            do i=1+sol_pil_w, l_ni-sol_pil_e
               C1(i,j,k,1)=-gama_8*(GVM%mc_iJz_8(i,j,k ) &
                            + mu_8*half)*(Ver_idz_8%m(k)+(GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k)) - gg_8
               C1(i,j,k,5)= gama_8*(GVM%mc_iJz_8(i,j,k ) &
                            - mu_8*half)*(Ver_idz_8%m(k)+(GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
               if (Schm_opentop_L) then
                   C1(i,j,k,1)=-gama_8*(GVM%mc_iJz_8(i,j,k ) + GVM%mc_iJz_8(i,j,k-1) )*Ver_idz_8%m(k) &
                              +(GVM%mc_Iz_8(i,j,k)-epsi_8)*gama_8*( Ver_wm_8%m(k)*(GVM%mc_iJz_8(i,j,k-1) -mu_8*half) &
                                                               -Ver_wp_8%m(k)*(GVM%mc_iJz_8(i,j,k )  +mu_8*half) ) - gg_8
                   C1(i,j,k,4)= gama_8*(GVM%mc_iJz_8(i,j,k-1) +  mu_8*half)*(Ver_idz_8%m(k) &
                               - (GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wm_8%m(k))
                   C1(i,j,k,5)= gama_8*(GVM%mc_iJz_8(i,j,k ) -mu_8*half)*(Ver_idz_8%m(k) + &
                                        (GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
               endif
            end do
         end do
!$omp enddo nowait

!$omp do
         do k = k0+1,l_nk
            do j=1+sol_pil_s, l_nj-sol_pil_n
               do i=1+sol_pil_w, l_ni-sol_pil_e
                   C1(i,j,k,1)=-gama_8*(GVM%mc_iJz_8(i,j,k ) + GVM%mc_iJz_8(i,j,k-1) )*Ver_idz_8%m(k) &
                              +(GVM%mc_Iz_8(i,j,k)-epsi_8)*gama_8*( Ver_wm_8%m(k)*(GVM%mc_iJz_8(i,j,k-1) -mu_8*half) &
                                                               -Ver_wp_8%m(k)*(GVM%mc_iJz_8(i,j,k )  +mu_8*half) ) - gg_8
                   C1(i,j,k,4)= gama_8*(GVM%mc_iJz_8(i,j,k-1) +  mu_8*half)*(Ver_idz_8%m(k) &
                               - (GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wm_8%m(k))
                   C1(i,j,k,5)= gama_8*(GVM%mc_iJz_8(i,j,k ) -mu_8*half)*(Ver_idz_8%m(k) + &
                                        (GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
               end do
            end do
         end do
!$omp enddo nowait

!$omp do
         do k = k0,l_nk
            km=max(k-1,1)
            kp=k+1
            do j=1+sol_pil_s_ext, l_nj-sol_pil_n

               do i=1+sol_pil_w_ext, l_ni-sol_pil_e_ext
                  A1(i,j,k,1)= -geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
                                (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
                  A2(i,j,k,1)=  geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
                                (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
                  A2(i,j,k,2)= -geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
                                (Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,km))
                  A1(i,j,k,3)=  geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
                                (Ver_wp_8%m(k)*GVM%mc_iJz_8(i+1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i+1,j,km))
                  A1(i,j,k,4)= half*GVM%mc_Jx_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
                  A2(i,j,k,4)= half*GVM%mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
                  A1(i,j,k,5)=-half*GVM%mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
                  A2(i,j,k,5)=-half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
                  A2(i,j,k,6)= half*GVM%mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,km)
                  A2(i,j,k,7)=-half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k)
                  A1(i,j,k,8)=half*GVM%mc_Jx_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i+1,j,km)
                  A1(i,j,k,9)=-half*GVM%mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i+1,j,k)
               end do
            end do

            do j=1+sol_pil_s_ext, l_nj-sol_pil_n_ext
               do i=1+sol_pil_w_ext, l_ni-sol_pil_e
                  B1(i,j,k,1)= -geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*     &
                                (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
                  B2(i,j,k,1)=  geomh_invDYMv_8(j-1) + half*GVM%mc_Jy_8(i,j-1,k)* &
                                (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
                  B1(i,j,k,4)= half*GVM%mc_Jy_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
                  B2(i,j,k,4)= half*GVM%mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
                  B1(i,j,k,5)= -half*GVM%mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
                  B2(i,j,k,5)= -half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
                  B2(i,j,k,10)=-geomh_invDYMv_8(j-1) + half*GVM%mc_Jy_8(i,j-1,k)* &
                               (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j-1,k ) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j-1,km))
                  B1(i,j,k,11)= geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*  &
                                (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j+1,k ) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j+1,km))
                  B2(i,j,k,12) =half*GVM%mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j-1,km)
                  B2(i,j,k,13)=-half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j-1,k)
                  B1(i,j,k,14)= half*GVM%mc_Jy_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j+1,km)
                  B1(i,j,k,15)= -half*GVM%mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j+1,k)
               end do
            end do

            if(.not.Grd_yinyang_L) then

               do j=1+sol_pil_s, l_nj-sol_pil_n
                  do i=1+sol_pil_w, l_ni-sol_pil_e
                     A2(i,j,k,1)=  geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
                                   (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
                     A2(i,j,k,2)= -geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
                                   (Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,km))
                     A2(i,j,k,4)= half*GVM%mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
                     A2(i,j,k,5)= -half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
                     A2(i,j,k,6)= half*GVM%mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,km)
                     A2(i,j,k,7)= -half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k)
                     B2(i,j,k,1)=  geomh_invDYMv_8(j-1)  + half*GVM%mc_Jy_8(i,j-1,k)* &
                                   (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
                     B2(i,j,k,4)= half*GVM%mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
                     B2(i,j,k,5)= -half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
                     B2(i,j,k,10)=-geomh_invDYMv_8(j-1) + half*GVM%mc_Jy_8(i,j-1,k)* &
                                   (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j-1,k ) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j-1,km))
                     B2(i,j,k,12) =half*GVM%mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j-1,km)
                     B2(i,j,k,13)=-half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j-1,k)

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

!$omp enddo

!$omp do collapse(3)
         do id=1,15
            do k =k0, l_nk
               do j=1+sol_pil_s, l_nj-sol_pil_n
                  do i=1+sol_pil_w, l_ni-sol_pil_e
                     Sol_stencilh_8 (i,j,k,id) =Cstv_hco0_8* ( (A1 (i,j,k,id)-A2 (i,j,k,id))*geomh_invDXM_8(j) &
                                                 + half * ( GVM%mc_Ix_8(i,j,k)*(A1(i,j,k,id)+A2(i,j,k,id)))    &
                                                 + (B1 (i,j,k,id)*geomh_cyM_8(j)-B2 (i,j,k,id)*                &
                                                    geomh_cyM_8(j-1))*geomh_invDYM_8(j)                        &
                                                 + half * (GVM%mc_Iy_8(i,j,k)*(B1(i,j,k,id)+B2(i,j,k,id)) )    &
                                                 + C1(i,j,k,id) )
                  enddo
               end do
            end do
         end do
!$omp enddo

      endif
!
!     ---------------------------------------------------------------
!
      return
      end subroutine matvec3D_init_hlt

