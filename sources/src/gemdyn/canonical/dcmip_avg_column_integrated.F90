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

!**s/r dcmip_avg_column_integrated - Evaluate Average Column Integrated of F_tr

      subroutine dcmip_avg_column_integrated (F_avg_tr,F_tr,Minx,Maxx,Miny,Maxy,F_nk)

      use dynkernel_options
      use glb_ld
      use gmm_geof
      use gmm_itf_mod
      use gmm_pw
      use gmm_vt1
      use lun
      use tdpack, only : rgasd_8, grav_8
      use ver

      implicit none

      !Arguments
      !---------
      integer,                                   intent(in) :: Minx,Maxx,Miny,Maxy  !Dimension H
      integer,                                   intent(in) :: F_nk                 !Number of vertical levels
      real, dimension(Minx:Maxx,Miny:Maxy)     , intent(out):: F_avg_tr             !Average Column Integrated of F_tr
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(in) :: F_tr                 !Tracer

      !object
      !===============================================
      !     Evaluate Average Column Integrated of F_tr
      !===============================================

      !----------------------------------------------------------

      integer :: i,j,k,istat
      real*8  :: avg_8,dp_8,dz_8
      logical :: GEM_P_L
      real, dimension(Minx:Maxx,Miny:Maxy)         :: pr_p0
      real, dimension(Minx:Maxx,Miny:Maxy,1:F_nk+1):: pr_m,pr_t
      real, dimension(Minx:Maxx,Miny:Maxy,0:F_nk+1):: zmom 

      !----------------------------------------------------------

      GEM_P_L = trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P'

      if (GEM_P_L) then

         !Get GMM variables
         !-----------------
         istat = gmm_get(gmmk_st1_s,st1)

         !Evaluate Pressure based on pw_update_GPW
         !----------------------------------------
         call calc_pressure ( pr_m, pr_t, pr_p0, st1, Minx, Maxx, Miny, Maxy, F_nk ) 

         pr_m(1:l_ni,1:l_nj,F_nk+1) = pr_p0(1:l_ni,1:l_nj)

         F_avg_tr(:,:) = 0.

         do j=1,l_nj
         do i=1,l_ni
            avg_8 = 0.0d0
            do k=1,F_nk
               dp_8          = (pr_m(i,j,k+1) - pr_m(i,j,k))
               avg_8         = avg_8         + dp_8
               F_avg_tr(i,j) = F_avg_tr(i,j) + F_tr(i,j,k) * dp_8
            end do
               F_avg_tr(i,j) = F_avg_tr(i,j)/avg_8
         end do
         end do

      else

         !Get GMM variables
         !-----------------
         istat = gmm_get(gmmk_tt1_s,tt1)
         istat = gmm_get(gmmk_pw_pt_plus_s,pw_pt_plus)
         istat = gmm_get (gmmk_fis0_s,fis0)

         zmom(:,:,0) = Ver_z_8%m(0)

         do k=1,F_nk
            zmom(:,:,k) = Ver_z_8%m(k) + Ver_b_8%m(k)*fis0(:,:)/grav_8
         end do

         zmom(:,:,F_nk+1) = fis0(:,:)/grav_8

         F_avg_tr(:,:) = 0.

         do j=1,l_nj
         do i=1,l_ni
            avg_8 = 0.0d0
            do k=1,F_nk
               dz_8          = ( pw_pt_plus(i,j,k) / (rgasd_8 * tt1(i,j,k)) ) * (zmom(i,j,k+1) - zmom(i,j,k))
               avg_8         = avg_8         + dz_8
               F_avg_tr(i,j) = F_avg_tr(i,j) + F_tr(i,j,k) * dz_8
            end do
               F_avg_tr(i,j) = F_avg_tr(i,j)/avg_8
         end do
         end do

      end if

      return
      end
