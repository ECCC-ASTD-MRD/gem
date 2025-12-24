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
      use tdpack
      use gem_options
      use glb_ld
      use hvdif_options
      use mem_tstp
!
      use cstv
      use dcst
!
      use ptopo

      implicit none



      integer i,j,k,ik,dim
      real, parameter :: p_naught=100000., eps=1.0e-5
      real, dimension(:,:,:), pointer :: pres_t, th, wk , tmp
!
!-------------------------------------------------------------------
!
      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*l_nk
      pres_t (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(      1:)
      th     (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(  dim+1:)
      wk     (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(2*dim+1:)
      tmp    (l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_nk) => WS1(3*dim+1:)

      Hzd_lnr_theta_z= min(max(0.,Hzd_lnr_theta_z),0.9999999)

!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy, l_nj+G_haloy
            do i=1-G_halox, l_ni+G_halox
               pres_t(i,j,k)= (p_naught/pw_pt_plus(i,j,k))**cappa_8
               th    (i,j,k)= tt1(i,j,k) * pres_t(i,j,k)
            end do
         end do
      end do
!$omp end do

      !Hybrid diffusion if hzd_hyb_nk >0
      if(hzd_hyb_nk > 0) then 
         do ik=1,hzd_hyb_nk
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  tmp (i,j,ik) = th(i,j,l_nk+1-ik) 
               end do
            end do
         end do
         call hzd_theta_alh (th,Hzd_lnr_theta_z,l_minx,l_maxx,l_miny,l_maxy, &
                            G_nk,hzd_Theta_ALH_it)
         call hzd_exp_deln ( tmp, Hzd_pwr_theta, Hzd_lnR_theta, wk,&
                          l_minx,l_maxx,l_miny,l_maxy, hzd_hyb_nk )
         do ik=1,hzd_hyb_nk
            do j=1-G_haloy,l_nj+G_haloy
               do i=1-G_halox,l_ni+G_halox
                  th(i,j,l_nk+1-ik) = tmp (i,j,ik)
               end do
            end do
         end do
      else
      !	Diffusion on constant z 
         call hzd_theta_alh (th,Hzd_lnr_theta_z,l_minx,l_maxx,l_miny,l_maxy, &
                               G_nk,hzd_Theta_ALH_it)
      endif

!$omp do collapse(2)
      do k=1,G_nk
         do j=1, l_nj
            do i=1, l_ni
               tt1(i,j,k)= th(i,j,k) / pres_t(i,j,k)
            end do
         end do
      end do
!$omp end do
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_theta_z
