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

module hzd_exp
  use hzd_mod
  use hvdif_options
  use gem_options
  use HORgrid_options
  use glb_ld
  use glb_pil
  use, intrinsic :: iso_fortran_env
  implicit none
#include <arch_specific.hf>
  private
  public :: hzd_exp_deln, hzd_exp_visco

contains

!**s/r hzd_exp_deln - 5 points explicit del 'n' horizontal diffusion
!                     for LAM configuration

      subroutine hzd_exp_deln ( F_c1, F_hgrid_S, Minx,Maxx,Miny,Maxy,Nk,&
                                F_vv, F_type_S )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=*)          , intent(in) :: F_hgrid_S
      character(len=*), optional, intent(in) :: F_type_S
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,NK),           intent (inout) :: F_c1
      real, dimension(Minx:Maxx,Miny:Maxy,NK), optional, intent (inout) :: F_vv
!author
!    Abdessamad Qaddouri - summer 2015
!
      integer iter1,mm,dpwr,itercnt,Niter
      real c1(Minx:Maxx,Miny:Maxy,Nk), c2(Minx:Maxx,Miny:Maxy,Nk)
      real(kind=REAL64) coef_8(Nk)
!
!     ---------------------------------------------------------------
!
      dpwr = Hzd_pwr
      niter= Hzd_niter
      if (niter > 0) coef_8(1:NK) = Hzd_coef_8(1:Nk)

      if (present(F_type_S)) then
         if (F_type_S == 'S_THETA') then
            dpwr = Hzd_pwr_theta
            niter= Hzd_niter_theta
            if (niter > 0) coef_8(1:NK) = Hzd_coef_8_theta(1:Nk)
         end if
         if (F_type_S == 'S_TR') then
            dpwr = Hzd_pwr_tr
            niter= Hzd_niter_tr
            if (niter > 0) coef_8(1:NK) = Hzd_coef_8_tr(1:Nk)
         end if
         if (F_type_S == 'VSPNG') then
            dpwr = 2
            niter= Vspng_niter
            if (niter > 0)coef_8(1:NK) = Vspng_coef_8(1:Nk)
         end if
      end if

      if (niter <= 0) return
      dpwr=dpwr/2

!     Fill all halo regions

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
           c2 = F_vv
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
               c2 = F_vv
         else
            call rpn_comm_xch_halo(F_c1,l_minx,l_maxx,l_miny,l_maxy,&
               l_ni,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         end if
      end if

      c1 = F_c1
      itercnt=0

      do iter1= 1, niter

         do mm=1,dpwr

            call hzd_exp5p ( F_c1, c1, l_minx,l_maxx,l_miny,l_maxy,&
                             Nk, coef_8, F_hgrid_S, mm,dpwr )
            if (present(F_vv)) then
               call hzd_exp5p ( F_vv, c2, l_minx,l_maxx,l_miny,l_maxy,&
                                Nk, coef_8, 'V' , mm,dpwr )
            end if
            itercnt = itercnt + 1

            if (itercnt == G_halox) then
               if (Grd_yinyang_L) then
                   if (present(F_vv)) then !2 vecteurs
                      call yyg_xchng_vec_uv2uv (c1,c2,l_minx,l_maxx,l_miny,l_maxy,NK)
                      if (Glb_pilotcirc_L) then
                          call rpn_comm_propagate_pilot_circular(c1, &
                          l_minx,l_maxx,l_miny,l_maxy, &
                          l_niu,l_nj,NK,Glb_pil_e,Glb_pil_s,G_halox,G_haloy)
                          call rpn_comm_propagate_pilot_circular(c2, &
                          l_minx,l_maxx,l_miny,l_maxy, &
                          l_ni,l_njv,NK,Glb_pil_e,Glb_pil_s,G_halox,G_haloy)
                      else
                          call rpn_comm_xch_halo(c1,l_minx,l_maxx,l_miny,l_maxy,&
                          l_niu,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         
                          call rpn_comm_xch_halo(c2,l_minx,l_maxx,l_miny,l_maxy,&
                          l_ni,l_njv,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
                      endif
                   else ! 1 champ scalaire
                       call yyg_xchng ( c1, l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj,&
                                        Nk, .false., 'CUBIC', .true. )
                   end if
               else
                   if (present(F_vv)) then !2vecteurs
                       call rpn_comm_xch_halo(c1,l_minx,l_maxx,l_miny,l_maxy,&
                               l_niu,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
                       call rpn_comm_xch_halo(c2,l_minx,l_maxx,l_miny,l_maxy,&
                               l_ni,l_njv,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
                   else ! 1 champ scalaire
                       call rpn_comm_xch_halo(c1,l_minx,l_maxx,l_miny,l_maxy,&
                               l_ni,l_nj,Nk,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
                   end if
               end if
               itercnt=0
            end if

         end do

      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine hzd_exp_deln
!
!**s/r hzd_exp_visco - applies horizontal explicit diffusion 9pt
!
      subroutine hzd_exp_visco(F_f2hzd, F_grd_S, Minx,Maxx,Miny,Maxy, NK)
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_grd_S
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_f2hzd

      integer :: nn, mm
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,Nk) :: wk1
      real :: rnr,pwr
      real(kind=REAL64) :: nu_dif,lnr,visco
      real(kind=REAL64), parameter :: epsilon = 1.0d-12
      real(kind=REAL64), parameter :: pt25=0.25d0

!     __________________________________________________________________
!
      rnr = log(1.- Hzd_lnR)
      pwr = Hzd_pwr

      if (F_grd_S == 'S_THETA') then
         rnr = log(1.- Hzd_lnR_theta)
         pwr = Hzd_pwr_theta
      end if
      if (F_grd_S == 'S_TR') then
         rnr = log(1.- Hzd_lnR_tr)
         pwr = Hzd_pwr_tr
      end if

      if ((F_grd_S == 'S_THETA').or.(F_grd_S == 'S_TR')) then

         lnr    = 1.0d0 - exp(rnr)
         nu_dif = 0.0d0
         if (pwr > 0) nu_dif = pt25*lnr**(2.d0/pwr)
         nu_dif = min ( nu_dif, pt25-epsilon )
         visco  = min ( nu_dif, pt25 )
         if (nu_dif < 1.0e-10) return

      else
!from 4.8: pwr=Hzd_del, lnr=Hzd_visco
         lnr    = 1.0d0 - exp(rnr)
         nu_dif = 0.0d0
         if (pwr > 0) nu_dif = pt25*lnr**(2.d0/pwr)
         nu_dif  = min ( nu_dif, pt25-epsilon )
         visco= min ( nu_dif, pt25 )
         if (nu_dif < 1.0e-10) return

      end if

      nn = pwr/2

      call rpn_comm_xch_halo ( F_f2hzd, l_minx,l_maxx,l_miny,l_maxy,&
           l_ni,l_nj, Nk, G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      do mm=1,nn

         call hzd_nudeln(F_f2hzd, wk1, l_minx,l_maxx,l_miny,l_maxy,&
                                                   Nk, visco, mm,nn )

         if (mm /= nn) then
              call rpn_comm_xch_halo( wk1, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj, Nk, G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         end if

      end do
!     __________________________________________________________________
!
      return
      end subroutine hzd_exp_visco

end module hzd_exp
