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

!**s/r pressure - Compute pressure and log of pressure on momentum
!                 and thermodynamic levels and also surface pressure.
!
      subroutine pressure ( F_pm_4,F_pt_4,F_p0_4,F_log_pm_4,F_log_pt_4,F_pm_8,F_p0_8, &
                            F_minx,F_maxx,F_miny,F_maxy,F_nk,F_time )

      use cstv
      use dyn_fisl_options
      use dynkernel_options
      use gem_options
      use glb_ld
      use gmm_geof
      use gmm_vt0
      use gmm_vt1
      use metric
      use tdpack
      use ver

      implicit none

      !arguments
      !---------
      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_nk,F_time
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk+1), intent(out) :: F_pm_4,F_pt_4,&
                                                                          F_log_pm_4,F_log_pt_4
      real, dimension(F_minx:F_maxx,F_miny:F_maxy),        intent(out) :: F_p0_4
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk+1), intent(out) :: F_pm_8
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy),        intent(out) :: F_p0_8

      integer :: i, j, k, i0,in,j0,jn
      real(kind=REAL64) :: pres_m, pres_t, log_pt
      real, pointer, dimension(:,:,:) :: qt
      real, pointer, dimension(:,:)   :: st
!
!     ________________________________________________________________
!
      if (F_time==1) st => st1
      if (F_time==1) qt => qt1
      if (F_time==0) st => st0
      if (F_time==0) qt => qt0
      i0= 1-G_halox ; in= l_ni+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

         do k=1,l_nk+1
            do j= j0, jn
               do i= i0, in
                  F_pm_8    (i,j,k) = qt(i,j,k)/(rgasd_8*Cstv_Tstr_8)+GVM%lg_pstar_8(i,j,k)
                  F_log_pm_4(i,j,k) = F_pm_8(i,j,k)
               end do
            end do
         end do
         do k=1,l_nk
            do j= j0, jn
               do i= i0, in
                  log_pt= 0.5d0*(F_pm_8(i,j,k+1)+F_pm_8(i,j,k))
                  F_pm_8    (i,j,k) = exp(F_pm_8(i,j,k))
                  F_pm_4    (i,j,k) = F_pm_8(i,j,k)
                  F_log_pt_4(i,j,k) = log_pt
                  F_pt_4    (i,j,k) = exp(log_pt)
               end do
            end do
         end do

         k= l_nk+1
         do j= j0, jn
            do i= i0, in
               F_pm_8    (i,j,k) = exp(F_pm_8(i,j,k))
               F_pm_4    (i,j,k) = F_pm_8    (i,j,k)
               F_pt_4    (i,j,k) = F_pm_4    (i,j,k)
               F_log_pt_4(i,j,k) = F_log_pm_4(i,j,k)
            end do
         end do

         do j= j0, jn
            do i= i0, in
               F_p0_8(i,j) = F_pm_8(i,j,l_nk+1)
               F_p0_4(i,j) = F_pm_4(i,j,l_nk+1)
            end do
         end do

      else

         do j= j0, jn
            do k=1,l_nk+1
            do i= i0, in
               pres_m = Ver_a_8%m(k) + Ver_b_8%m(k) * st(i,j)&
                       +Ver_c_8%m(k) * sls(i,j)
               pres_t = Ver_a_8%t(k) + Ver_b_8%t(k) * st(i,j)&
                       +Ver_c_8%t(k) * sls(i,j)
               F_log_pm_4(i,j,k) = pres_m
               F_log_pt_4(i,j,k) = pres_t
               F_pm_8    (i,j,k) = exp(pres_m)
               F_pm_4    (i,j,k) = F_pm_8(i,j,k)
               F_pt_4    (i,j,k) = exp(pres_t)
            end do
            end do
         end do

         k= l_nk+1
         do j= j0, jn
         do i= i0, in
            F_p0_8(i,j)= F_pm_8(i,j,k)
            F_p0_4(i,j)= F_pm_8(i,j,k)
         end do
         end do

      end if
!
!     ________________________________________________________________
!
      return
      end subroutine pressure
