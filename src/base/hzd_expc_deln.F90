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

!
!**s/r hzd_exp_visco - applies explicit 9pt del N horizontal filtering operator
!
      subroutine hzd_expc_deln( F_f2hzd, F_rho, F_pwr, F_lnr,&
                                   Minx,Maxx,Miny,Maxy, NK )
      use hzd_mod
      use hvdif_options
      use gem_options
      use HORgrid_options
      use glb_ld
      use glb_pil
      use dcst
      use tdpack
      use mem_tstp
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_pwr,Minx,Maxx,Miny,Maxy,Nk
      real,    intent(in) :: F_lnr
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_f2hzd
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) ::  F_rho 
      real, dimension(:,:,:), pointer :: F_wk

      integer :: i,j,k, nn,mm, i0,in,j0,jn
      real(kind=REAL64) :: nu_dif,visco
      real(kind=REAL64), parameter :: epsilon=1.0d-12, pt25=0.25d0
!
!-------------------------------------------------------------------
!
      F_wk(minx:maxx,miny:maxy,1:nk) => WS1(1:)

      if (Grd_yinyang_L) then
         i0 = 1    + 2*west
         j0 = 1    + 2*south
         in = l_ni - 2*east
         jn = l_nj - 2*north
      else
         i0 = 1    + pil_w
         j0 = 1    + pil_s
         in = l_ni - pil_e
         jn = l_nj - pil_n
      end if


      nu_dif = 0.0d0
      if (F_pwr ==2) then
         nu_dif = pt25*dble(F_lnr)
      else
         nu_dif= pt25*sqrt(F_lnr)
     endif

      nn = F_pwr/2

!$omp single
      call rpn_comm_xch_halo ( F_f2hzd, l_minx,l_maxx,l_miny,l_maxy,&
           l_ni,l_nj, Nk, G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
!$omp end single

!$omp do collapse(2)
      do k=1,nk
         do j=1-G_haloy, l_nj+G_haloy
            do i=1-G_halox, l_ni+G_halox
               F_wk(i,j,k)= F_f2hzd(i,j,k)
            end do
         end do
      end do
!$omp end do

      do mm=1, nn
         if (mm == 2) then
!$omp do collapse(2)
            do k=1, nk
               do j=j0-1, jn+1
                  do i=i0-1, in+1
                     F_wk(i,j,k)= F_f2hzd(i,j,k) - F_wk(i,j,k)
                  end do
               end do
            end do
!$omp end do
         else if (mm > 2) then
!$omp do collapse(2)
            do k=1, nk
               do j=j0-1+south, jn+1-north
                  do i=i0-1+west, in+1-east
                     F_wk(i,j,k)= F_f2hzd(i,j,k) - F_wk(i,j,k)
                  end do
               end do
            end do
!$omp end do
         end if

         call hzd_CvDel2_flt9pt (F_f2hzd,F_wk,F_rho,l_minx,l_maxx,l_miny,l_maxy,nk,&
                                   nu_dif,mm, nn,i0,in,j0,jn)
         if (mm /= nn) then
!$omp single
              call rpn_comm_xch_halo( F_wk, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj, Nk, G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
!$omp end single
         end if
      end do
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_expc_deln

