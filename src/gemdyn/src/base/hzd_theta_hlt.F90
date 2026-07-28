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

!**s/r hzd_theta - applies horizontal diffusion on theta
!
      subroutine hzd_theta_hlt
      use hzd_exp_hlt
      use gmm_pw
      use gmm_vt1
      use gmm_hzd
      use tdpack
      use gem_options
      use glb_ld
      use hvdif_options
      use mem_tstp
      use hzd_mod
      use cstv

      implicit none

      integer i,j,k,dim
      real, parameter :: p_naught=100000., eps=1.0e-5

!-------------------------------------------------------------------

!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy, l_nj+G_haloy
            do i=1-G_halox, l_ni+G_halox
               pres_pt(i,j,k)= (p_naught/pw_pt_plus(i,j,k))**cappa_8
               theta    (i,j,k)= tt1(i,j,k) * pres_pt(i,j,k)
               theta0   (i,j,k)= theta(i,j,k)
            end do
         end do
      end do
!$omp end do
      if (hzd_conserv_th) then
        call hzd_expc_deln ( theta ,air_dens, Hzd_pwr_theta, Hzd_lnR_theta, &
                            l_minx,l_maxx,l_miny,l_maxy, G_nk )
      else
         call hzd_exp_deln ( theta, Hzd_pwr_theta, Hzd_lnR_theta, &
                             l_minx,l_maxx,l_miny,l_maxy, G_nk )
      endif
      if(hzd_apply_th_tend) then 
!$omp do collapse(2)
      	 do k=1,G_nk
            do j=1, l_nj
               do i=1, l_ni
                  hzd_th_tend(i,j,k)= (theta(i,j,k) - theta0 (i,j,k))/Cstv_dt_8 
                  hzd_th_tend(i,j,k)= hzd_th_tend(i,j,k) / theta0(i,j,k) 
               end do
            end do
         end do
!$omp end do
      else
!$omp do collapse(2)
         do k=1,G_nk
            do j=1, l_nj
               do i=1, l_ni
                  tt1(i,j,k)= theta(i,j,k) / pres_pt(i,j,k)
               end do
            end do
         end do
!$omp end do
      endif

!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_theta_hlt
