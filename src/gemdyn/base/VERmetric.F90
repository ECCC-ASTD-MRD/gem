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

!**s/r VERmetric - Compute vertical metric coefficients
      
      subroutine VERmetric ()
      use mtn_options
      use HORgrid_options
      use gmm_geof
      use geomh
      use tdpack
      use glb_ld
      use metric
      use cstv
      use ver
      use fislh_sol
      use gem_options
      use lam_options
      use dyn_fisl_options  
      implicit none

      integer :: i,j,k, k0
      real, parameter :: one=1.d0, half=.5d0
!
!     ---------------------------------------------------------------
      k0= 1+Lam_gbpil_T

!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               zmom_8(i,j,k)=ver_z_8%m(k)+Cstv_bar1_8*(Ver_b_8%m(k)*fis0(i,j)+Ver_c_8%m(k)*sls(i,j))/grav_8
               ztht_8(i,j,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*fis0(i,j)+Ver_c_8%t(k)*sls(i,j))/grav_8
            end do
         end do
      end do
!$omp enddo
!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            ztht_8(i,j,0)=ver_z_8%m(0)
            zmom_8(i,j,0)=ver_z_8%m(0)
            zmom_8(i,j,G_nk+1)= Cstv_bar1_8*fis0(i,j)/grav_8
            ztht_8(i,j,G_nk  )= zmom_8(i,j,G_nk+1)
            ztht_8(i,j,G_nk+1)= zmom_8(i,j,G_nk+1)
            lg_pstar_8(i,j,G_nk+1)=log(1.d5)-grav_8*zmom_8(i,j,G_nk+1)/(rgasd_8*Cstv_Tstr_8)
         end do
      end do
!$omp enddo
!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do k=G_nk,1,-1
            do i=1-G_halox,l_ni+G_halox
               lg_pstar_8(i,j,k)=lg_pstar_8(i,j,k+1)+grav_8*(zmom_8(i,j,k+1)-zmom_8(i,j,k))/(rgasd_8*Cstv_Tstr_8)
            end do
         end do
      end do
!$omp enddo
      
!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               mc_Jx_8 (i,j,k)=(zmom_8(i+1,j,k)-zmom_8(i,j,k))*geomh_invDX_8(j)
               mc_Jy_8 (i,j,k)=(zmom_8(i,j+1,k)-zmom_8(i,j,k))*geomh_invDY_8
               mc_iJz_8(i,j,k)=one/(zmom_8(i,j,k+1)-zmom_8(i,j,k))
               mc_Ix_8(i,j,k)=log( (ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(ztht_8(i-1,j,k)-ztht_8(i-1,j,k-1)) )*0.5d0*geomh_invDX_8(j)
               mc_Iy_8(i,j,k)=log( (ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(ztht_8(i,j-1,k)-ztht_8(i,j-1,k-1)) )*0.5d0*geomh_invDY_8
               mc_Iz_8(i,j,k)=log( (zmom_8(i,j,k+1)-zmom_8(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) &
                                  /(zmom_8(i,j,k)-zmom_8(i,j,k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)) )*Ver_idz_8%m(k)
               mc_logJz_8(i,j,k)= 0.0
            end do
         end do
      end do
!$omp enddo
!$omp do
      do j=1-G_haloy+1,l_nj+G_haloy-1
!DIR$ SIMD
         do i=1-G_halox+1,l_ni+G_halox-1
            ztht_8(i,j,G_nk)= ver_z_8%t(G_nk)+Cstv_bar1_8*(Ver_b_8%t(G_nk)*fis0(i,j)+Ver_c_8%t(G_nk)*sls(i,j))/grav_8
            mc_css_H_8(i,j) = one/(gama_8*(isol_i*mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk)-half*mu_8))
            me_full (i,j) = fis0(i,j) / grav_8
            me_large(i,j) = sls (i,j) / grav_8
         end do
      end do
!$omp enddo

      if (Schm_opentop_L) then
!$omp do
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               mc_cst_8(i,j)= one / (-(mu_8* Cstv_tau_nh_8)*(isol_d*Ver_idz_8%t(k0-1) &
                              +isol_i*mc_iJz_8(i,j,k0-1)) & 
                              + half* one/(Cstv_tau_8*cpd_8*Cstv_Tstr_8))
            end do
         end do
!$omp enddo
      endif
!$omp do
      do j=1-G_haloy+1,l_nj+G_haloy-1
         do i=1-G_halox+1,l_ni+G_halox-1
            mc_alfas_H_8(i,j) = mc_css_H_8(i,j)*(gama_8*(isol_i*mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk) &
                                + half*mu_8) + Ver_wmstar_8(G_nk)*gama_8*(isol_i*mc_iJz_8(i,j,G_nk-1) &
                                               + isol_d*Ver_idz_8%t(G_nk-1)-half*mu_8) )

            mc_betas_H_8(i,j) = mc_css_H_8(i,j)* Ver_wmstar_8(G_nk)*gama_8* &
                                (isol_i*mc_iJz_8(i,j,G_nk-1)+isol_d*Ver_idz_8%t(G_nk-1)+half*mu_8)

            mc_cssp_H_8(i,j) = gama_8*(Ver_idz_8%m(G_nk)-epsi_8*Ver_wp_8%m(G_nk))*&
                               (isol_i*mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk)-half*mu_8)*mc_css_H_8(i,j)
         enddo
      enddo
!$omp enddo
      if (Schm_opentop_L) then
!$omp do
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               mc_alfat_8(i,j)  = (-(mu_8* Cstv_tau_nh_8)*(isol_d*Ver_idz_8%t(k0-1)+isol_i*mc_iJz_8(i,j,k0-1))  &
                                 - half* one/(Cstv_tau_8*cpd_8*Cstv_Tstr_8))*mc_cst_8(i,j)
               mc_cstp_8(i,j)   = mc_cst_8(i,j) * gama_8 * &
                                 ((isol_d*Ver_idz_8%t(k0-1)+isol_i*mc_iJz_8(i,j,k0-1))*&
                                 (Ver_idz_8%m(k0)+Ver_wm_8%m(k0)*epsi_8) + &
                                  mu_8*half*(Ver_idz_8%m(k0)+Ver_wm_8%m(k0)*epsi_8) )
            end do
         end do
!$omp enddo
      endif
      
      if (Schm_autobar_L) then
!$omp single
         mc_Jx_8   (:,:,:) = 0.
         mc_Jy_8   (:,:,:) = 0.
         mc_iJz_8  (:,:,:) = 1.
         mc_Ix_8   (:,:,:) = 0.
         mc_Iy_8   (:,:,:) = 0.
         mc_Iz_8   (:,:,:) = 0.
         mc_logJz_8(:,:,:) = 0.

         mc_alfas_H_8(:,:) = 1.0d0
         mc_betas_H_8(:,:) = 0.0d0
         mc_css_H_8  (:,:) = 0.0d0
         mc_cssp_H_8 (:,:) = 0.0d0
         
         if (Schm_opentop_L) then 
            mc_alfat_8  (:,:) = 1.0d0
            mc_cst_8    (:,:) = 0.0d0
            mc_cstp_8   (:,:) = 0.0d0
         end if
!$omp end single
      end if

!      call heights_uv ()
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine VERmetric
