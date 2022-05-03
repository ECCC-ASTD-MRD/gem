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

!**s/r vertical_metric - Compute vertical metric coefficients
      
      subroutine vertical_metric (F_metric, F_topo, F_sls, Minx,Maxx,Miny,Maxy)
      use gem_options
      use lam_options
      use dyn_fisl_options
      use dynkernel_options
      use HORgrid_options
      use geomh
      use tdpack
      use glb_ld
      use mem_nest
      use mem_tstp
      use cstv
      use ver
      use fislh_sol
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy
      real, dimension (Minx:Maxx,Miny:Maxy), intent(IN) :: F_topo, F_sls
      type(Vmetric) , intent(INOUT) :: F_metric
      
      integer :: i,j,k, k0, err, dim
      real, parameter :: one=1.d0, half=.5d0
      real, dimension(:,:,:), pointer :: thick_m, thick_t
!
!     ---------------------------------------------------------------
!
      k0= 1+Lam_gbpil_T
      dim= (l_ni+2*G_halox)*(l_nj+2*G_haloy)*l_nk
      thick_m(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,1:l_nk) => WS1(1:)
      thick_t(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,1:l_nk) => WS1(dim+1:)

      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               F_metric%zmom_8(i,j,k)=ver_z_8%m(k)+Cstv_bar1_8*(Ver_b_8%m(k)*F_topo(i,j)+Ver_c_8%m(k)*F_sls(i,j))/grav_8
               F_metric%ztht_8(i,j,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*F_topo(i,j)+Ver_c_8%t(k)*F_sls(i,j))/grav_8
            end do
         end do
      end do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            F_metric%ztht_8(i,j,0)=ver_z_8%m(0)
            F_metric%zmom_8(i,j,0)=ver_z_8%m(0)
            F_metric%zmom_8(i,j,G_nk+1)= Cstv_bar1_8*F_topo(i,j)/grav_8
            F_metric%ztht_8(i,j,G_nk+1)= F_metric%zmom_8(i,j,G_nk+1)
            F_metric%lg_pstar_8(i,j,G_nk+1)=log(1.d5)-grav_8*F_metric%zmom_8(i,j,G_nk+1)/(rgasd_8*Cstv_Tstr_8)
         end do
      end do
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               thick_m(i,j,k)=F_metric%zmom_8(i,j,k)-F_metric%zmom_8(i,j,k+1)
               thick_t(i,j,k)=F_metric%ztht_8(i,j,k)-F_metric%ztht_8(i,j,k+1)
            end do
         end do
      end do
      err=0
      if (minval(thick_m)<0. .or. minval(thick_t)<0. ) err=-1
      call gem_error (err,'vertical_metric','Heights NOT monotonically decreasing from model top')
      do j=1-G_haloy,l_nj+G_haloy
         do k=G_nk,1,-1
            do i=1-G_halox,l_ni+G_halox
               F_metric%lg_pstar_8(i,j,k)=F_metric%lg_pstar_8(i,j,k+1)+grav_8*(F_metric%zmom_8(i,j,k+1)-F_metric%zmom_8(i,j,k))/(rgasd_8*Cstv_Tstr_8)
            end do
         end do
         do i=1-G_halox,l_ni+G_halox
            F_metric%ztht_8(i,j,G_nk)= F_metric%zmom_8(i,j,G_nk+1) !temporary for mc_Ix_8 and mc_Iy_8 below
         end do
      end do
      do k=1,G_nk
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               F_metric%mc_Jx_8 (i,j,k)=(F_metric%zmom_8(i+1,j,k)-F_metric%zmom_8(i,j,k))*geomh_invDX_8(j)
               F_metric%mc_Jy_8 (i,j,k)=(F_metric%zmom_8(i,j+1,k)-F_metric%zmom_8(i,j,k))*geomh_invDY_8
               F_metric%mc_iJz_8(i,j,k)=one/(F_metric%zmom_8(i,j,k+1)-F_metric%zmom_8(i,j,k))
               F_metric%mc_Ix_8(i,j,k)=log( (F_metric%ztht_8(i+1,j,k)-F_metric%ztht_8(i+1,j,k-1))/(F_metric%ztht_8(i-1,j,k)-F_metric%ztht_8(i-1,j,k-1)) )*0.5d0*geomh_invDX_8(j)
               F_metric%mc_Iy_8(i,j,k)=log( (F_metric%ztht_8(i,j+1,k)-F_metric%ztht_8(i,j+1,k-1))/(F_metric%ztht_8(i,j-1,k)-F_metric%ztht_8(i,j-1,k-1)) )*0.5d0*geomh_invDY_8
               F_metric%mc_Iz_8(i,j,k)=log( (F_metric%zmom_8(i,j,k+1)-F_metric%zmom_8(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) &
                                  /(F_metric%zmom_8(i,j,k)-F_metric%zmom_8(i,j,k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)) )*Ver_idz_8%m(k)
               F_metric%mc_logJz_8(i,j,k)= 0.0
            end do
         end do
      end do
      do j=1-G_haloy+1,l_nj+G_haloy-1
!DIR$ SIMD
         do i=1-G_halox+1,l_ni+G_halox-1
            F_metric%ztht_8(i,j,G_nk)= ver_z_8%t(G_nk)+Cstv_bar1_8*(Ver_b_8%t(G_nk)*F_topo(i,j)+Ver_c_8%t(G_nk)*F_sls(i,j))/grav_8
            F_metric%mc_css_H_8(i,j) = one/(gama_8*(isol_i*F_metric%mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk)-half*mu_8))
         end do
      end do

      if (Schm_opentop_L) then
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               F_metric%mc_cst_8(i,j)= one / (-(mu_8* Cstv_tau_nh_8)*(isol_d*Ver_idz_8%t(k0-1) &
                              +isol_i*F_metric%mc_iJz_8(i,j,k0-1)) & 
                              + half* one/(Cstv_tau_8*cpd_8*Cstv_Tstr_8))
            end do
         end do
      endif
      do j=1-G_haloy+1,l_nj+G_haloy-1
         do i=1-G_halox+1,l_ni+G_halox-1
            F_metric%mc_alfas_H_8(i,j) = F_metric%mc_css_H_8(i,j)*(gama_8*(isol_i*F_metric%mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk) &
                                + half*mu_8) + Ver_wmstar_8(G_nk)*gama_8*(isol_i*F_metric%mc_iJz_8(i,j,G_nk-1) &
                                               + isol_d*Ver_idz_8%t(G_nk-1)-half*mu_8) )

            F_metric%mc_betas_H_8(i,j) = F_metric%mc_css_H_8(i,j)* Ver_wmstar_8(G_nk)*gama_8* &
                                (isol_i*F_metric%mc_iJz_8(i,j,G_nk-1)+isol_d*Ver_idz_8%t(G_nk-1)+half*mu_8)

            F_metric%mc_cssp_H_8(i,j) = gama_8*(Ver_idz_8%m(G_nk)-epsi_8*Ver_wp_8%m(G_nk))*&
                               (isol_i*F_metric%mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk)-half*mu_8)*F_metric%mc_css_H_8(i,j)
         enddo
      enddo
      if (Schm_opentop_L) then
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               F_metric%mc_alfat_8(i,j)  = (-(mu_8* Cstv_tau_nh_8)*(isol_d*Ver_idz_8%t(k0-1)+isol_i*F_metric%mc_iJz_8(i,j,k0-1))  &
                                 - half* one/(Cstv_tau_8*cpd_8*Cstv_Tstr_8))*F_metric%mc_cst_8(i,j)
               F_metric%mc_cstp_8(i,j)   = F_metric%mc_cst_8(i,j) * gama_8 * &
                                 ((isol_d*Ver_idz_8%t(k0-1)+isol_i*F_metric%mc_iJz_8(i,j,k0-1))*&
                                 (Ver_idz_8%m(k0)+Ver_wm_8%m(k0)*epsi_8) + &
                                  mu_8*half*(Ver_idz_8%m(k0)+Ver_wm_8%m(k0)*epsi_8) )
            end do
         end do
      endif
      
      if (Schm_autobar_L) then
         F_metric%mc_Jx_8   (:,:,:) = 0.
         F_metric%mc_Jy_8   (:,:,:) = 0.
         F_metric%mc_iJz_8  (:,:,:) = 1.
         F_metric%mc_Ix_8   (:,:,:) = 0.
         F_metric%mc_Iy_8   (:,:,:) = 0.
         F_metric%mc_Iz_8   (:,:,:) = 0.
         F_metric%mc_logJz_8(:,:,:) = 0.

         F_metric%mc_alfas_H_8(:,:) = 1.0d0
         F_metric%mc_betas_H_8(:,:) = 0.0d0
         F_metric%mc_css_H_8  (:,:) = 0.0d0
         F_metric%mc_cssp_H_8 (:,:) = 0.0d0
         
         if (Schm_opentop_L) then 
            F_metric%mc_alfat_8  (:,:) = 1.0d0
            F_metric%mc_cst_8    (:,:) = 0.0d0
            F_metric%mc_cstp_8   (:,:) = 0.0d0
         end if
      end if

!      call heights_uv ()
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine vertical_metric


      subroutine vertical_metric_omp (F_metric, F_topo, F_sls, Minx,Maxx,Miny,Maxy)
      use gem_options
      use lam_options
      use dyn_fisl_options
      use dynkernel_options
      use HORgrid_options
      use geomh
      use tdpack
      use glb_ld
      use mem_nest
      use mem_tstp
      use cstv
      use ver
      use fislh_sol
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy
      real, dimension (Minx:Maxx,Miny:Maxy), intent(IN) :: F_topo, F_sls
      type(Vmetric) , intent(INOUT) :: F_metric
      
      integer :: i,j,k, k0, err, dim
      real, parameter :: one=1.d0, half=.5d0
      real, dimension(:,:,:), pointer :: thick_m, thick_t
!
!     ---------------------------------------------------------------
!
      k0= 1+Lam_gbpil_T
      dim= (l_ni+2*G_halox)*(l_nj+2*G_haloy)*l_nk
      thick_m(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,1:l_nk) => WS1(1:)
      thick_t(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,1:l_nk) => WS1(dim+1:)

!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               F_metric%zmom_8(i,j,k)=ver_z_8%m(k)+Cstv_bar1_8*(Ver_b_8%m(k)*F_topo(i,j)+Ver_c_8%m(k)*F_sls(i,j))/grav_8
               F_metric%ztht_8(i,j,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*F_topo(i,j)+Ver_c_8%t(k)*F_sls(i,j))/grav_8
            end do
         end do
      end do
!$omp enddo
!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            F_metric%ztht_8(i,j,0)=ver_z_8%m(0)
            F_metric%zmom_8(i,j,0)=ver_z_8%m(0)
            F_metric%zmom_8(i,j,G_nk+1)= Cstv_bar1_8*F_topo(i,j)/grav_8
            F_metric%ztht_8(i,j,G_nk+1)= F_metric%zmom_8(i,j,G_nk+1)
            F_metric%lg_pstar_8(i,j,G_nk+1)=log(1.d5)-grav_8*F_metric%zmom_8(i,j,G_nk+1)/(rgasd_8*Cstv_Tstr_8)
         end do
      end do
!$omp enddo
!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               thick_m(i,j,k)=F_metric%zmom_8(i,j,k)-F_metric%zmom_8(i,j,k+1)
               thick_t(i,j,k)=F_metric%ztht_8(i,j,k)-F_metric%ztht_8(i,j,k+1)
            end do
         end do
      end do
!$omp enddo

!$omp single
      err=0
      if (minval(thick_m)<0. .or. minval(thick_t)<0. ) err=-1
      call gem_error (err,'vertical_metric','Heights NOT monotonically decreasing from model top')
!$omp end single

!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do k=G_nk,1,-1
            do i=1-G_halox,l_ni+G_halox
               F_metric%lg_pstar_8(i,j,k)=F_metric%lg_pstar_8(i,j,k+1)+grav_8*(F_metric%zmom_8(i,j,k+1)-F_metric%zmom_8(i,j,k))/(rgasd_8*Cstv_Tstr_8)
            end do
         end do
         do i=1-G_halox,l_ni+G_halox
            F_metric%ztht_8(i,j,G_nk)= F_metric%zmom_8(i,j,G_nk+1) !temporary for mc_Ix_8 and mc_Iy_8 below
         end do
      end do
!$omp enddo
      
!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               F_metric%mc_Jx_8 (i,j,k)=(F_metric%zmom_8(i+1,j,k)-F_metric%zmom_8(i,j,k))*geomh_invDX_8(j)
               F_metric%mc_Jy_8 (i,j,k)=(F_metric%zmom_8(i,j+1,k)-F_metric%zmom_8(i,j,k))*geomh_invDY_8
               F_metric%mc_iJz_8(i,j,k)=one/(F_metric%zmom_8(i,j,k+1)-F_metric%zmom_8(i,j,k))
               F_metric%mc_Ix_8(i,j,k)=log( (F_metric%ztht_8(i+1,j,k)-F_metric%ztht_8(i+1,j,k-1))/(F_metric%ztht_8(i-1,j,k)-F_metric%ztht_8(i-1,j,k-1)) )*0.5d0*geomh_invDX_8(j)
               F_metric%mc_Iy_8(i,j,k)=log( (F_metric%ztht_8(i,j+1,k)-F_metric%ztht_8(i,j+1,k-1))/(F_metric%ztht_8(i,j-1,k)-F_metric%ztht_8(i,j-1,k-1)) )*0.5d0*geomh_invDY_8
               F_metric%mc_Iz_8(i,j,k)=log( (F_metric%zmom_8(i,j,k+1)-F_metric%zmom_8(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) &
                                  /(F_metric%zmom_8(i,j,k)-F_metric%zmom_8(i,j,k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)) )*Ver_idz_8%m(k)
               F_metric%mc_logJz_8(i,j,k)= 0.0
            end do
         end do
      end do
!$omp enddo
!$omp do
      do j=1-G_haloy+1,l_nj+G_haloy-1
!DIR$ SIMD
         do i=1-G_halox+1,l_ni+G_halox-1
            F_metric%ztht_8(i,j,G_nk)= ver_z_8%t(G_nk)+Cstv_bar1_8*(Ver_b_8%t(G_nk)*F_topo(i,j)+Ver_c_8%t(G_nk)*F_sls(i,j))/grav_8
            F_metric%mc_css_H_8(i,j) = one/(gama_8*(isol_i*F_metric%mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk)-half*mu_8))
         end do
      end do
!$omp enddo

      if (Schm_opentop_L) then
!$omp do
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               F_metric%mc_cst_8(i,j)= one / (-(mu_8* Cstv_tau_nh_8)*(isol_d*Ver_idz_8%t(k0-1) &
                              +isol_i*F_metric%mc_iJz_8(i,j,k0-1)) & 
                              + half* one/(Cstv_tau_8*cpd_8*Cstv_Tstr_8))
            end do
         end do
!$omp enddo
      endif
!$omp do
      do j=1-G_haloy+1,l_nj+G_haloy-1
         do i=1-G_halox+1,l_ni+G_halox-1
            F_metric%mc_alfas_H_8(i,j) = F_metric%mc_css_H_8(i,j)*(gama_8*(isol_i*F_metric%mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk) &
                                + half*mu_8) + Ver_wmstar_8(G_nk)*gama_8*(isol_i*F_metric%mc_iJz_8(i,j,G_nk-1) &
                                               + isol_d*Ver_idz_8%t(G_nk-1)-half*mu_8) )

            F_metric%mc_betas_H_8(i,j) = F_metric%mc_css_H_8(i,j)* Ver_wmstar_8(G_nk)*gama_8* &
                                (isol_i*F_metric%mc_iJz_8(i,j,G_nk-1)+isol_d*Ver_idz_8%t(G_nk-1)+half*mu_8)

            F_metric%mc_cssp_H_8(i,j) = gama_8*(Ver_idz_8%m(G_nk)-epsi_8*Ver_wp_8%m(G_nk))*&
                               (isol_i*F_metric%mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk)-half*mu_8)*F_metric%mc_css_H_8(i,j)
         enddo
      enddo
!$omp enddo
      if (Schm_opentop_L) then
!$omp do
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               F_metric%mc_alfat_8(i,j)  = (-(mu_8* Cstv_tau_nh_8)*(isol_d*Ver_idz_8%t(k0-1)+isol_i*F_metric%mc_iJz_8(i,j,k0-1))  &
                                 - half* one/(Cstv_tau_8*cpd_8*Cstv_Tstr_8))*F_metric%mc_cst_8(i,j)
               F_metric%mc_cstp_8(i,j)   = F_metric%mc_cst_8(i,j) * gama_8 * &
                                 ((isol_d*Ver_idz_8%t(k0-1)+isol_i*F_metric%mc_iJz_8(i,j,k0-1))*&
                                 (Ver_idz_8%m(k0)+Ver_wm_8%m(k0)*epsi_8) + &
                                  mu_8*half*(Ver_idz_8%m(k0)+Ver_wm_8%m(k0)*epsi_8) )
            end do
         end do
!$omp enddo
      endif
      
      if (Schm_autobar_L) then
!$omp single
         F_metric%mc_Jx_8   (:,:,:) = 0.
         F_metric%mc_Jy_8   (:,:,:) = 0.
         F_metric%mc_iJz_8  (:,:,:) = 1.
         F_metric%mc_Ix_8   (:,:,:) = 0.
         F_metric%mc_Iy_8   (:,:,:) = 0.
         F_metric%mc_Iz_8   (:,:,:) = 0.
         F_metric%mc_logJz_8(:,:,:) = 0.

         F_metric%mc_alfas_H_8(:,:) = 1.0d0
         F_metric%mc_betas_H_8(:,:) = 0.0d0
         F_metric%mc_css_H_8  (:,:) = 0.0d0
         F_metric%mc_cssp_H_8 (:,:) = 0.0d0
         
         if (Schm_opentop_L) then 
            F_metric%mc_alfat_8  (:,:) = 1.0d0
            F_metric%mc_cst_8    (:,:) = 0.0d0
            F_metric%mc_cstp_8   (:,:) = 0.0d0
         end if
!$omp end single
      end if

!      call heights_uv ()
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine vertical_metric_omp
