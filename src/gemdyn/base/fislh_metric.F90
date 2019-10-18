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

!**   s/r fislh_metric - calculate metric coefficients for GEM in height-base coordinates
      subroutine fislh_metric()
      use mtn_options
      use HORgrid_options
      use gmm_geof
      use geomh
      use tdpack
      use glb_ld
      use lun
      use metric
      use cstv
      use ver
      use fislh_sol
      implicit none
#include <arch_specific.hf>
!
!Author: Claude Girard, July 2017 (initial version)
!        Syed Husain, June 2019 (revision)

      integer :: i,j,k
      real, parameter :: zero=0.d0, one=1.d0, half=.5d0
      logical,save :: done = .false.
!
!     ---------------------------------------------------------------

      if(.not.done) then
         allocate ( zmom(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    ztht(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    zmom_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    ztht_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                lg_pstar_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) )

         allocate ( mc_Jx_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_Jy_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_iJz_8(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                  mc_logJz_8(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_Ix_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_Iy_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_Iz_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk) )

         allocate ( mc_css_H_8   (l_minx:l_maxx,l_miny:l_maxy), &
                    mc_alfas_H_8 (l_minx:l_maxx,l_miny:l_maxy), &
                    mc_betas_H_8 (l_minx:l_maxx,l_miny:l_maxy), &
                    mc_cssp_H_8  (l_minx:l_maxx,l_miny:l_maxy))
         done=.true.
      end if

      ztht_8(:,:,0)=ver_z_8%m(0)
      zmom_8(:,:,0)=ver_z_8%m(0)

      do k=1,G_nk
         zmom_8(:,:,k)=ver_z_8%m(k)+Cstv_bar1_8*(Ver_b_8%m(k)*fis0(:,:)+Ver_c_8%m(k)*sls(:,:))/grav_8
         ztht_8(:,:,k)=ver_z_8%t(k)+Cstv_bar1_8*(Ver_b_8%t(k)*fis0(:,:)+Ver_c_8%t(k)*sls(:,:))/grav_8
      end do
      zmom_8(:,:,G_nk+1)=Cstv_bar1_8*fis0(:,:)/grav_8
      ztht_8(:,:,G_nk+1)=Cstv_bar1_8*fis0(:,:)/grav_8

      lg_pstar_8(:,:,G_nk+1)=log(1.d5)-grav_8*zmom_8(:,:,G_nk+1)/(rgasd_8*Cstv_Tstr_8)

      do k=G_nk,1,-1
         lg_pstar_8(:,:,k)=lg_pstar_8(:,:,k+1)+grav_8*(zmom_8(:,:,k+1)-zmom_8(:,:,k))/(rgasd_8*Cstv_Tstr_8)
      end do

      do k=1,G_nk
         do j=l_miny+1,l_maxy-1
            do i=l_minx+1,l_maxx-1
               mc_Jx_8 (i,j,k)=(zmom_8(i+1,j,k)-zmom_8(i,j,k))*geomh_invDX_8(j)
               mc_Jy_8 (i,j,k)=(zmom_8(i,j+1,k)-zmom_8(i,j,k))*geomh_invDY_8
               mc_iJz_8(i,j,k)=one/(zmom_8(i,j,k+1)-zmom_8(i,j,k))
            end do
         end do
      end do
      ztht_8(:,:,G_nk)=Cstv_bar1_8*fis0(:,:)/grav_8
      do k=1,G_nk
         do j=l_miny+1,l_maxy-1
            do i=l_minx+1,l_maxx-1
               mc_Ix_8(i,j,k)=log( (ztht_8(i+1,j,k)-ztht_8(i+1,j,k-1))/(ztht_8(i-1,j,k)-ztht_8(i-1,j,k-1)) )*0.5d0*geomh_invDX_8(j)
               mc_Iy_8(i,j,k)=log( (ztht_8(i,j+1,k)-ztht_8(i,j+1,k-1))/(ztht_8(i,j-1,k)-ztht_8(i,j-1,k-1)) )*0.5d0*geomh_invDY_8
               mc_Iz_8(i,j,k)=log( (zmom_8(i,j,k+1)-zmom_8(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) &
                            /(zmom_8(i,j,k)-zmom_8(i,j,k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)) )*Ver_idz_8%m(k)
      !        mc_Ix_8(i,j,k)=0.0
      !        mc_Iy_8(i,j,k)=0.0
      !        mc_Iz_8(i,j,k)=0.0
      !        mc_logJz_8(i,j,k)=log( (ztht_8(i,j,k)-ztht_8(i,j,k-1))/(Ver_z_8%x(k)-Ver_z_8%x(k-1)) )
               mc_logJz_8(i,j,k)=0.0
            end do
         end do
      end do
      ztht_8(:,:,G_nk)=ver_z_8%t(G_nk)+Cstv_bar1_8*(Ver_b_8%t(G_nk)*fis0(:,:)+Ver_c_8%t(G_nk)*sls(:,:))/grav_8

      do j=l_miny+1,l_maxy-1
         do i=l_minx+1,l_maxx-1

            mc_css_H_8(i,j)   = one/(gama_8*(isol_i*mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk)-half*mu_8))

            mc_alfas_H_8(i,j) = mc_css_H_8(i,j)*(gama_8*(isol_i*mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk) &
                                + half*mu_8) + Ver_wmstar_8(G_nk)*gama_8*(isol_i*mc_iJz_8(i,j,G_nk-1) &
                                               + isol_d*Ver_idz_8%t(G_nk-1)-half*mu_8) )

            mc_betas_H_8(i,j) = mc_css_H_8(i,j)* Ver_wmstar_8(G_nk)*gama_8* &
                                (isol_i*mc_iJz_8(i,j,G_nk-1)+isol_d*Ver_idz_8%t(G_nk-1)+half*mu_8)

            mc_cssp_H_8(i,j) = gama_8*(Ver_idz_8%m(G_nk)-epsi_8*Ver_wp_8%m(G_nk))*&
                               (isol_i*mc_iJz_8(i,j,G_nk)+isol_d*Ver_idz_8%t(G_nk)-half*mu_8)*mc_css_H_8(i,j)
         enddo
      enddo

      if (Schm_autobar_L) then
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
      end if

      ! We keep a copy in single precision only for cascade outputs
      do k=1,G_nk+1
         zmom(:,:,k) = real(zmom_8(:,:,k))
         ztht(:,:,k) = real(ztht_8(:,:,k))
      end do

!     ---------------------------------------------------------------
!
      return
      end
