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

!**s/r get_air_dens - Evaluez la densite reduite de l'air:
!                                                    densite*jacobien

      subroutine get_air_density_hlt (F_time)

      use dyn_fisl_options
      use geomh
      use glb_ld
      use gmm_pw
      use gmm_hzd
      use gmm_vt1
      use tr3d
      use mem_tracers
      use tdpack, only: grav_8
      use ver
      use gem_options
      use lam_options
      use metric
!      
      use, intrinsic :: iso_fortran_env
      implicit none

      !arguments
      !---------
      integer,    intent(in) :: F_time    !Time 1/0
    
      !object
      !==========================================
      !     Evaluate Air Mass at appropriate time
      !==========================================

      integer :: i,j,k,k0,km
      real(kind=REAL64) :: dens,inv_grav
      real, pointer, dimension(:    ) :: tr
      real, pointer, dimension(:,:,:) :: pr_m,hu
!
!---------------------------------------------------------------------
!
      !k0= 1+Lam_gbpil_T
      k0= 1

      !Obtain Momentum pressure levels at appropriate time
      !---------------------------------------------------
      if (F_time == 1) then
         pr_m => pw_pm_plus
         tr   => trt1
         hu   => tracers_P(Tr3d_hu)%pntr
      endif
      if (F_time == 0) then
         pr_m => pw_pm_moins
         tr   => trt0
         hu   => tracers_M(Tr3d_hu)%pntr
      endif

      !Evaluate water tracers at appropriate time if dry mixing ratio
      !--------------------------------------------------------------

      call sumhydro_hlt (sumq_8,l_minx,l_maxx,l_miny,l_maxy,l_nk,Tr3d_ntr,&
                         tr, Schm_dry_mixing_ratio_L.and.Schm_wload_L)
      if (Schm_dry_mixing_ratio_L) then
!$omp do collapse(2)
         do k=k0,l_nk
           do j=1,l_nj
             do i=1,l_ni
               sumq_8(i,j,k) = sumq_8(i,j,k) + hu(i,j,k)
             enddo
           enddo
         enddo
!$omp end do
      end if

      !Evaluate Air Mass
      !-----------------
      inv_grav= 1.d0 / grav_8
      air_dens(:,:,1:k0-1)= 0.
      air_dens_m(:,:,1:k0-1)= 0.
!$omp do
      do k=k0,l_nk
         air_dens(:,l_miny:0     ,k)= 0.
         air_dens_m(:,l_miny:0     ,k)= 0.
         air_dens(:,l_nj+1:l_maxy,k)= 0.
         air_dens_m(:,l_nj+1:l_maxy,k)= 0.
         do j=1,l_nj
         air_dens(l_minx:0     ,j,k)= 0.
         air_dens_m(l_minx:0     ,j,k)= 0.
         air_dens(l_ni+1:l_maxx,j,k)= 0.
         air_dens_m(l_ni+1:l_maxx,j,k)= 0.
         do i=1,l_ni
            air_dens(i,j,k)=  (pr_m(i,j,k) - pr_m(i,j,k+1)) * (1.-sumq_8(i,j,k)) * Ver_idz_8%t(k) * inv_grav
            air_dens(i,j,k)= &
                   air_dens(i,j,k)*(GVM%zmom_8(i,j,k+1)-GVM%zmom_8(i,j,k))/(ver_z_8%m(k+1)-ver_z_8%m(k))
         end do
         end do
      end do
!$omp end do

!$omp do
      do k=1,l_nk
         km=max(k-1,1)
         do j=1,l_nj
            do i=1,l_ni
               air_dens_m(i,j,k) = Ver_wp_8%m(k)*air_dens(i,j,k) &
                                 + Ver_wm_8%m(k)*air_dens(i,j,km)
            end do
         end do
      end do
!$omp end do

!$omp single
         call rpn_comm_xch_halo (air_dens,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo (air_dens_m,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
!$omp end single
!---------------------------------------------------------------------
!
      return
      end subroutine get_air_density_hlt
