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
      real, dimension(:,:,:), pointer :: pres_t, th, th0, tmp, tmp1, tmp2, wk
!
!-------------------------------------------------------------------
!
      dim1=(l_maxx-l_minx+1)*(l_maxy-l_miny+1)*(max(hzd_hyb_bot,hzd_hyb_top))
      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*l_nk
      pres_t (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(      1:)
      th     (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(  dim+1:)
      th0    (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(2*dim+1:)
      wk     (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(3*dim+1:)
      tmp    (l_minx:l_maxx,l_miny:l_maxy,1:max(hzd_hyb_bot,hzd_hyb_top)) => WS1(4*dim+1:)
      tmp1   (l_minx:l_maxx,l_miny:l_maxy,1:max(hzd_hyb_bot,hzd_hyb_top)) => WS1(4*dim+dim1+1:)
      tmp2   (l_minx:l_maxx,l_miny:l_maxy,1:max(hzd_hyb_bot,hzd_hyb_top)) => WS1(4*dim+2*dim1+1:)

!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy, l_nj+G_haloy
            do i=1-G_halox, l_ni+G_halox
               pres_t(i,j,k) = (p_naught/pw_pt_plus(i,j,k))**cappa_8
               th    (i,j,k) = tt1(i,j,k) * pres_t(i,j,k)
               th0   (i,j,k) = th(i,j,k)
            end do
         end do
      end do
!$omp end do

      !Hybrid diffusion if hzd_hyb_bot >0
      ! Diffusion on constant z for levels hzd_hyb_bot->hzd_hyb_top
      if(hzd_hyb_bot > 0) then 
!$omp do 
         do ik=1,max(hzd_hyb_bot,hzd_hyb_top )
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  tmp (i,j,ik) = th(i,j,l_nk+1-ik) 
                  tmp1(i,j,ik) = air_dens(i,j,l_nk+1-ik)
                  tmp2(i,j,ik) = th(i,j,ik)
               end do
            end do
         end do
!$omp end do
         if (hzd_conserv_th) then
!$omp single
            call hzd_uvwzd_alh(th,Hzd_lnR_theta_z,Hzd_pwr_theta_z, &
                              l_minx,l_maxx,l_miny,l_maxy,G_nk,4)
!$omp end single
!$omp single
            call hzd_expc_deln ( tmp ,tmp1, Hzd_pwr_theta, Hzd_lnR_theta, wk,&
                                l_minx,l_maxx,l_miny,l_maxy, hzd_hyb_bot )
!$omp end single

!$omp do collapse(2)
         do ik=1,hzd_hyb_bot
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  th(i,j,l_nk+1-ik) = tmp (i,j,ik)
               end do
            end do
         end do
!$omp end do

     if(hzd_hyb_top.ne.1 ) then
         do ik=1,hzd_hyb_top
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  tmp1(i,j,ik) = air_dens(i,j,ik)
               end do
            end do
         end do
!$omp single
            call hzd_expc_deln ( tmp2,tmp1, Hzd_pwr_theta, Hzd_lnR_theta, wk,&
                                l_minx,l_maxx,l_miny,l_maxy, hzd_hyb_top )
!$omp end single

!$omp do collapse(2)
         do ik=1,hzd_hyb_top
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  th(i,j,ik) = tmp2 (i,j,ik)
               end do
            end do
         end do
!$omp end do

      endif

         else
!$omp single
            call hzd_uvwzd_alh(th,Hzd_lnR_theta_z,Hzd_pwr_theta_z,&
                                l_minx,l_maxx,l_miny,l_maxy,G_nk,3)

!$omp end single
            call hzd_exp_deln ( tmp, Hzd_pwr_theta, Hzd_lnR_theta, wk,&
                                l_minx,l_maxx,l_miny,l_maxy, hzd_hyb_bot )
!         endif
!$omp do collapse(2)
         do ik=1,hzd_hyb_bot
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  th(i,j,l_nk+1-ik) = tmp (i,j,ik)
               end do
            end do
         end do
!$omp end do
      endif
      else
      ! Diffusion on constant z 
         if (hzd_conserv_th) then
!$omp single
            call hzd_uvwzd_alh(th,Hzd_lnR_theta_z,Hzd_pwr_theta_z,&
                               l_minx,l_maxx,l_miny,l_maxy,G_nk,4)
!$omp end single
         else
!$omp single
            call hzd_uvwzd_alh(th,Hzd_lnR_theta_z,Hzd_pwr_theta_z,&
                                l_minx,l_maxx,l_miny,l_maxy,G_nk,3)
!$omp end single
         endif
      endif

      if(hzd_apply_th_tend) then 
!$omp do collapse(2)
      	 do k=1,G_nk
            do j=1, l_nj
               do i=1, l_ni
                  hzd_th_tend(i,j,k)= (th(i,j,k) - th0 (i,j,k))/Cstv_dt_8 
                  hzd_th_tend(i,j,k)= hzd_th_tend(i,j,k) / th0(i,j,k) 
               end do
            end do
         end do
!$omp end do
      else
!$omp do collapse(2)
         do k=1,G_nk
            do j=1, l_nj
               do i=1, l_ni
                  tt1(i,j,k)= th(i,j,k) / pres_t(i,j,k)
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
