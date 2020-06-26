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

!**s/r hzd_exp_set

      subroutine hzd_exp_set
      use hzd_mod
      use hvdif_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      use ver
      implicit none
#include <arch_specific.hf>

!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,1002)

      Hzd_lnr = min(max(0.,Hzd_lnr),0.9999999)
      Hzd_pwr = Hzd_pwr / 2
      Hzd_pwr = min(max(2,Hzd_pwr*2),8)

      Hzd_lnr_theta= min(max(0.,Hzd_lnr_theta),0.9999999)
      Hzd_pwr_theta= Hzd_pwr_theta / 2
      Hzd_pwr_theta= min(max(2,Hzd_pwr_theta*2),8)

      if (Hzd_lnr_tr < 0.) Hzd_lnr_tr = Hzd_lnr
      if (Hzd_pwr_tr < 0 ) Hzd_pwr_tr = Hzd_pwr
      Hzd_lnr_tr = min(max(0.,Hzd_lnr_tr),0.9999999)
      Hzd_pwr_tr = Hzd_pwr_tr / 2
      Hzd_pwr_tr = min(max(2,Hzd_pwr_tr*2),8)

      if ((Hzd_lnr <= 0.).and.(Hzd_lnr_theta <= 0.)  &
                         .and.(Hzd_lnr_tr <= 0.)) then
         if((Hzd_smago_param <= 0.).and.(Hzd_smago_lnr(2) == 0.)) then
            if (Lun_out > 0) write(Lun_out,1003)
         else
            if (Lun_out > 0) then
               write(Lun_out,1004) Hzd_smago_param,100*Hzd_smago_lnr(2)
            end if
         end if
      end if

      call hzd_exp_geom ()

      call hzd_exp5p_set ()

 1002 format(/,'INITIALIZATING HIGH ORDER HORIZONTAL DIFFUSION ',  &
               '(S/R HZD_SET)',/,60('='))
 1003 format(/,'NO HORIZONTAL DIFFUSION REQUESTED',/,33('='))
 1004 format(/,'  HORIZONTAL DIFFUSION A LA SMAGORINSKY',/,2x,37('=')// &
              ,'  PARAMETER =',f5.2,'  BACKGROUND =',f4.1,' %/TIMESTEP')
!
!     ---------------------------------------------------------------
!
      return
      end
