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
!---------------------------------- LICENCE END ----------------------

!*s/r set_smago - computing vertically varying coefficients for smag. diffusion

      subroutine set_smago

      use dynkernel_options
      use hzd_mod
      use HORgrid_options
      use hvdif_options
      use geomh
      use tdpack
      use glb_ld
      use cstv
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

!
!author
!     Syed Husain
!

      integer k
      real(kind=REAL64) top_m, bot_m, top_t, bot_t, pi2
      real(kind=REAL64) Hzd_smago_bot_lev, Hzd_smago_top_lev
      real(kind=REAL64) Hzd_smago_bot_lnr, Hzd_smago_mid_lnr, Hzd_smago_top_lnr
!

      allocate( Hzd_smago_lnrM_8(G_nk), Hzd_smago_lnrT_8(G_nk))

      if (trim(Dynamics_Kernel_S)=='DYNAMICS_FISL_H') then
         call set_smago_fislh()
         return
      end if

      if ((Hzd_smago_lev(1) < 0.) .and. (Hzd_smago_lev(2) < 0.) &
          .and. (Hzd_smago_lnr(1) < 0.)) then
         Hzd_smago_lnrM_8(1:G_nk)= 0.d0
         Hzd_smago_lnrT_8(1:G_nk)= 0.d0
         return
      end if

      pi2=pi_8/2.d0

      Hzd_smago_bot_lev=Hzd_smago_lev(1)
      Hzd_smago_top_lev=Hzd_smago_lev(2)

      Hzd_smago_bot_lnr=Hzd_smago_lnr(1)
      Hzd_smago_mid_lnr=Hzd_smago_lnr(2)
      Hzd_smago_top_lnr=Hzd_smago_lnr(3)

      top_m= max( Hzd_smago_top_lev, dble(Ver_hyb%m(  1 ) ))
      bot_m= min( Hzd_smago_bot_lev, dble(Ver_hyb%m(G_nk) ))

      top_t= max( Hzd_smago_top_lev, dble(Ver_hyb%t(  1 ) ))
      bot_t= min( Hzd_smago_bot_lev, dble(Ver_hyb%t(G_nk) ))

      if ((Hzd_smago_bot_lev < 0. .or. Hzd_smago_top_lev < 0.) &
          .and. (Hzd_smago_bot_lnr > 0.)) then
         do k=1,G_nk
            Hzd_smago_lnrM_8(k)= Hzd_smago_bot_lnr
            Hzd_smago_lnrT_8(k)= Hzd_smago_bot_lnr
         end do

      else
         do k=1,G_nk
            ! For momentum levels
            if (Ver_hyb%m(k) <= bot_m .and. Ver_hyb%m(k) >= top_m) then
               Hzd_smago_lnrM_8(k) = cos(pi2-pi2*(bot_m - Ver_hyb%m(k))/(bot_m-top_m))
               Hzd_smago_lnrM_8(k) = Hzd_smago_bot_lnr + &
                                     Hzd_smago_lnrM_8(k)*Hzd_smago_lnrM_8(k)* &
                                     (Hzd_smago_mid_lnr - Hzd_smago_bot_lnr)
            elseif (Ver_hyb%m(k) < top_m) then
               Hzd_smago_lnrM_8(k) = cos(pi2-pi2*(top_m - Ver_hyb%m(k))/(top_m- Ver_hyb%m(1)))
               Hzd_smago_lnrM_8(k) = Hzd_smago_mid_lnr + &
                                     Hzd_smago_lnrM_8(k)*Hzd_smago_lnrM_8(k)* &
                                     (Hzd_smago_top_lnr - Hzd_smago_mid_lnr)
            else
               Hzd_smago_lnrM_8(k)= Hzd_smago_bot_lnr
            end if

            ! For thermodynamic levels
            if (Ver_hyb%t(k) <= bot_t .and. Ver_hyb%t(k) >= top_t) then
               Hzd_smago_lnrT_8(k) = cos(pi2-pi2*(bot_t - Ver_hyb%t(k))/(bot_t-top_t))
               Hzd_smago_lnrT_8(k) = Hzd_smago_bot_lnr + &
                                     Hzd_smago_lnrT_8(k)*Hzd_smago_lnrT_8(k)* &
                                     (Hzd_smago_mid_lnr - Hzd_smago_bot_lnr)
            elseif (Ver_hyb%t(k) < top_t) then
               Hzd_smago_lnrT_8(k) = cos(pi2-pi2*(top_t - Ver_hyb%t(k))/(top_t-Ver_hyb%t(1)))
               Hzd_smago_lnrT_8(k) = Hzd_smago_mid_lnr + &
                                     Hzd_smago_lnrT_8(k)*Hzd_smago_lnrT_8(k)* &
                                     (Hzd_smago_top_lnr - Hzd_smago_mid_lnr)
            else
               Hzd_smago_lnrT_8(k)= Hzd_smago_bot_lnr
            end if

            Hzd_smago_lnrM_8(k)= max(0.d0, Hzd_smago_lnrM_8(k))
            Hzd_smago_lnrT_8(k)= max(0.d0, Hzd_smago_lnrT_8(k))

         end do
      end if
!
      return
      end subroutine set_smago
