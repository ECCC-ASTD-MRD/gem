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
      subroutine hzd_theta_z()
      use hzd_exp_hlt
      use gmm_pw
      use gmm_vt1
      use gmm_hzd
      use tdpack
      use gem_options
      use glb_ld
      use hvdif_options
      use mem_tstp
      use cstv
      use dcst
      use ptopo
      implicit none

      integer i,j,k,ik,dim,dim1
      real, parameter :: p_naught=100000.
!
!-------------------------------------------------------------------
!
!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy, l_nj+G_haloy
            do i=1-G_halox, l_ni+G_halox
               pres_pt(i,j,k) = (p_naught/pw_pt_plus(i,j,k))**cappa_8
               theta    (i,j,k) = tt1(i,j,k) * pres_pt(i,j,k)
               theta0   (i,j,k) = theta(i,j,k)
            end do
         end do
      end do
!$omp end do

      !Hybrid diffusion if hzd_hyb_bot >0
      ! Diffusion on constant z for levels hzd_hyb_bot->hzd_hyb_top
      if(hzd_hyb_bot > 0) then 
!$omp do collapse(2)
      do ik=1,hzd_hyb_bot
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
                  wrkt1 (i,j,ik) = theta   (i,j,l_nk+1-ik) 
                  wrkd1 (i,j,ik) = air_dens(i,j,l_nk+1-ik)
            end do
         end do
      end do
!$omp end do
         if (hzd_conserv_th) then
            call hzd_uvwzd_alh(theta,Hzd_lnR_theta_z,Hzd_pwr_theta_z, &
                              l_minx,l_maxx,l_miny,l_maxy,G_nk,4)

            call hzd_expc_deln ( wrkt1 ,wrkd1, Hzd_pwr_theta, Hzd_lnR_theta,&
                                l_minx,l_maxx,l_miny,l_maxy, hzd_hyb_bot )
!$omp do collapse(2)
         do ik=1,hzd_hyb_bot
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  theta(i,j,l_nk+1-ik) = wrkt1 (i,j,ik)
               end do
            end do
         end do
!$omp end do

         if(hzd_hyb_top.ne.1 ) then
!$omp do collapse(2)
            do ik=1,hzd_hyb_top
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox
                     wrkt2(i,j,ik) = theta   (i,j,ik)
                     wrkd2(i,j,ik) = air_dens(i,j,ik)
                  end do
               end do
            end do
!$omp end do
            call hzd_expc_deln ( wrkt2,wrkd2, Hzd_pwr_theta, Hzd_lnR_theta,&
                                l_minx,l_maxx,l_miny,l_maxy, hzd_hyb_top )
!$omp do collapse(2)
            do ik=1,hzd_hyb_top
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox
                     theta(i,j,ik) = wrkt2 (i,j,ik)
                  end do
               end do
            end do
!$omp end do
            endif
         else
            call hzd_uvwzd_alh(theta,Hzd_lnR_theta_z,Hzd_pwr_theta_z,&
                                l_minx,l_maxx,l_miny,l_maxy,G_nk,3)
            call hzd_exp_deln ( wrkt1, Hzd_pwr_theta, Hzd_lnR_theta,&
                                l_minx,l_maxx,l_miny,l_maxy, hzd_hyb_bot )
!$omp do collapse(2)
            do ik=1,hzd_hyb_bot
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox
                     theta(i,j,l_nk+1-ik) = wrkt1 (i,j,ik)
                  end do
               end do
            end do
!$omp end do
         endif
      else
      ! Diffusion on constant z 
         if (hzd_conserv_th) then
            call hzd_uvwzd_alh(theta,Hzd_lnR_theta_z,Hzd_pwr_theta_z,&
                               l_minx,l_maxx,l_miny,l_maxy,G_nk,4)
         else
            call hzd_uvwzd_alh(theta,Hzd_lnR_theta_z,Hzd_pwr_theta_z,&
                                l_minx,l_maxx,l_miny,l_maxy,G_nk,3)
         endif
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
      end subroutine hzd_theta_z
