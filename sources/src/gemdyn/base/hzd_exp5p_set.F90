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

!**s/r hzd_exp5p_set - Horizontal diffusion delN setup for LAMs

      subroutine hzd_exp5p_set
      use dcst
      use hzd_mod
      use HORgrid_options
      use hvdif_options
      use gem_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      implicit none
#include <arch_specific.hf>

!author
!    Abdessamad Qaddouri - summer 2015
!
!revision
! v4_80 - Qaddouri A.      - initial version
! v4_80 - Lee   - optimization
!


      real*8 coef_8,coef_theta_8,coef_tr_8,nutop_8,c_8,deg2rad_8
!
!     ---------------------------------------------------------------
!
      deg2rad_8 = pi_8 / 180.0d0
      c_8= min(Grd_dx,Grd_dy)
      c_8= c_8 * deg2rad_8

      Hzd_Niter= 0 ; Hzd_Niter_theta= 0 ; Hzd_Niter_tr= 0
      coef_8= 0. ; coef_theta_8= 0. ; coef_tr_8= 0.
      !for U,V,W,Zd
      if ( (Hzd_lnR > 0.) .and. (Hzd_pwr > 0) ) then
         nutop_8 = 1./4. * Hzd_lnR**(2./Hzd_pwr)
         Hzd_Niter= max(int(8.d0*nutop_8+0.9999999),1)
         coef_8= nutop_8/max(1.,float(HZD_niter))* &
                                  ((Dcst_rayt_8*c_8)**2)/Cstv_dt_8
         allocate( Hzd_coef_8(G_nk))
         Hzd_coef_8(1:G_nk) = coef_8*(Dcst_inv_rayt_8**2)*Cstv_dt_8
      end if

      !for Theta
      if ( (Hzd_lnR_theta > 0) .and. (Hzd_pwr_theta > 0) )  then
         nutop_8 = 1./4. * Hzd_lnR_theta**(2./Hzd_pwr_theta)
         Hzd_Niter_theta = max(int(8.d0*nutop_8+0.9999999),1)
         coef_theta_8=nutop_8/max(1.,float(hzd_niter_theta))* &
                                  ((Dcst_rayt_8*c_8)**2)/Cstv_dt_8
         allocate( Hzd_coef_8_theta(G_nk))
         Hzd_coef_8_theta(1:G_nk) = coef_theta_8*(Dcst_inv_rayt_8**2)*Cstv_dt_8
      end if

      !for Tracers
      if ( (Hzd_lnR_tr > 0) .and. (Hzd_pwr_tr > 0) )  then
         nutop_8 = 1./4. * Hzd_lnR_tr**(2./Hzd_pwr_tr)
         Hzd_Niter_tr = max(int(8.d0*nutop_8+0.9999999),1)
         coef_tr_8=nutop_8/max(1.,float(hzd_niter_tr))* &
                                  ((Dcst_rayt_8*c_8)**2)/Cstv_dt_8
         allocate( Hzd_coef_8_tr(G_nk))
         Hzd_coef_8_tr(1:G_nk) = coef_tr_8*(Dcst_inv_rayt_8**2)*Cstv_dt_8
      end if

      if (Lun_out > 0) then
         write(Lun_out,1010) coef_8      ,Hzd_pwr/2      ,'Winds ',Hzd_Niter
         write(Lun_out,1010) coef_theta_8,Hzd_pwr_theta/2,'Theta ',Hzd_Niter_theta
         write(Lun_out,1010) coef_tr_8   ,Hzd_pwr_tr/2   ,'Tracer',Hzd_Niter_tr, ' if specified'
      end if

1010 format (3X,'Diffusion Coefficient =  (',e17.10,' m**2)**',i1,'/sec ',a,' Niter=',i2,a )
!
!     ---------------------------------------------------------------
!
      return
      end
