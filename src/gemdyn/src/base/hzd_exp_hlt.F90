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

module hzd_exp_hlt
  use hzd_mod
  use hvdif_options
  use gem_options
  use HORgrid_options
  use glb_ld
  use glb_pil
  use, intrinsic :: iso_fortran_env
  implicit none

contains

!**s/r hzd_exp_del2 - 5 points explicit del 'n' horizontal diffusion
!                     for LAM configuration

      subroutine hzd_exp_del2 ( F_c1, F_hgrid_S, Minx,Maxx,Miny,Maxy,Nk,&
                                F_geom, F_vv)
      use, intrinsic :: iso_fortran_env
      implicit none

      character(len=*), intent(in) :: F_hgrid_S
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,NK),           intent (inout) :: F_c1
      real, dimension(Minx:Maxx,Miny:Maxy,NK), optional, intent (inout) :: F_vv
      real(kind=REAL64), dimension(Miny:Maxy,*), intent(IN) :: F_geom
!author
!    Abdessamad Qaddouri - summer 2015
!
      integer :: iter, i0,in,j0,jn,inv,jnv,lni, available
!
!     ---------------------------------------------------------------
!
      lni = l_ni ; if ( F_hgrid_S == 'U' ) lni = l_niu
      i0  = 2         - G_halox * (1 - west )
      j0  = 2         - G_haloy * (1 - south)
      in  = lni   - 1 + G_halox * (1 - east )
      jn  = l_nj  - 1 + G_haloy * (1 - north)
      inv = l_ni  - 1 + G_halox * (1 - east )
      jnv = l_njv - 1 + G_haloy * (1 - north)
      available= 0

      do iter= 1, Vspng_niter
         if (available<1)   then
!$omp single
         if (Grd_yinyang_L) then
         if (present(F_vv)) then
            call yyg_xchng_vec_uv2uv (F_c1, F_vv,l_minx,l_maxx,l_miny,l_maxy,G_nk)
            if (Glb_pilotcirc_L) then
                call rpn_comm_propagate_pilot_circular(F_c1, &
                l_minx,l_maxx,l_miny,l_maxy, &
                l_niu,l_nj,NK,Glb_pil_e,Glb_pil_s,G_halox,G_haloy)
                call rpn_comm_propagate_pilot_circular(F_vv, &
                l_minx,l_maxx,l_miny,l_maxy, &
                l_ni,l_njv,NK,Glb_pil_e,Glb_pil_s,G_halox,G_haloy)
            else
                call rpn_comm_xch_halo(F_c1,l_minx,l_maxx,l_miny,l_maxy,&
               l_niu,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
                call rpn_comm_xch_halo(F_vv,l_minx,l_maxx,l_miny,l_maxy,&
               l_ni,l_njv,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
            endif
         else
            call yyg_xchng (F_c1, l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj,&
                            Nk, .false., 'CUBIC', .true. )
         end if
         else
         if (present(F_vv)) then
            call rpn_comm_xch_halo(F_c1,l_minx,l_maxx,l_miny,l_maxy,&
               l_niu,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
            call rpn_comm_xch_halo(F_vv,l_minx,l_maxx,l_miny,l_maxy,&
               l_ni,l_njv,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         else
            call rpn_comm_xch_halo(F_c1,l_minx,l_maxx,l_miny,l_maxy,&
               l_ni,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         end if
         end if
!$omp end single
         available=G_halox
         endif

         call hzd_flt5pt ( F_c1, F_geom, l_minx,l_maxx,l_miny,l_maxy,&
                           Nk, Vspng_coef_8, i0,in,j0,jn )
         if (present(F_vv)) then
           call hzd_flt5pt (F_vv,Hzd_geom_v,l_minx,l_maxx,l_miny,l_maxy,&
                            Nk, Vspng_coef_8, i0,inv,j0,jnv)
         end if
         available=available-1
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine hzd_exp_del2
!
!**s/r hzd_exp_visco - applies explicit 9pt del N horizontal filtering operator
!
      subroutine hzd_exp_deln( F_f2hzd, F_pwr, F_lnr, F_wk,&
                                   Minx,Maxx,Miny,Maxy, NK )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_pwr,Minx,Maxx,Miny,Maxy,Nk
      real,    intent(in) :: F_lnr
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_f2hzd,F_wk

      integer :: i,j,k, nn,mm, i0,in,j0,jn
      real(kind=REAL64) :: nu_dif,visco
      real(kind=REAL64), parameter :: epsilon=1.0d-12, pt25=0.25d0
!
!-------------------------------------------------------------------
!
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
      if (F_pwr > 0) nu_dif = pt25*dble(F_lnr)**(2.d0/dble(F_pwr))
      nu_dif = min ( nu_dif, pt25-epsilon )
      visco  = min ( nu_dif, pt25 )
      if (nu_dif < 1.0e-10) return

      nn = F_pwr/2

!$omp single
      call rpn_comm_xch_halo ( F_f2hzd, l_minx,l_maxx,l_miny,l_maxy,&
           l_ni,l_nj, Nk, G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
!$omp end single

!$omp do collapse(2)
      do k=1, l_nk
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
            do k=1, l_nk
               do j=j0-1, jn+1
                  do i=i0-1, in+1
                     F_wk(i,j,k)= F_f2hzd(i,j,k) - F_wk(i,j,k)
                  end do
               end do
            end do
!$omp end do
         else if (mm > 2) then
!$omp do collapse(2)
            do k=1, l_nk
               do j=j0-1+south, jn+1-north
                  do i=i0-1+west, in+1-east
                     F_wk(i,j,k)= F_f2hzd(i,j,k) - F_wk(i,j,k)
                  end do
               end do
            end do
!$omp end do
         end if

         call hzd_flt9pt (F_f2hzd, F_wk, l_minx,l_maxx,l_miny,l_maxy,&
                                    l_nk, visco, mm, nn, i0,in,j0,jn)

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
      end subroutine hzd_exp_deln

end module hzd_exp_hlt
