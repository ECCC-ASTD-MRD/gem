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

!**s/r init_bar - prepare data for autobarotropic runs (Williamson cases)

      subroutine init_bar ()
      use dynkernel_options
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use inp_mod
      use dyn_fisl_options
      use ctrl
      use gem_options
      use wil_options
      use glb_ld
      use tr3d
      use step_options
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>

      !object
      !============================================================
      !Prepare data for autobarotropic runs (Williamson cases)
      !------------------------------------------------------------
      !NOTE: U,V output on Staggered grids
      !============================================================

      integer :: istat, k
      real, dimension (:,:,:), pointer :: hu
      real, dimension (l_minx:l_maxx,l_miny:l_maxy,G_nk) :: gz_t
!
!-------------------------------------------------------------------
!
      if (Vtopo_L)      call gem_error (-1,'INIT_BAR','Vtopo_L not available YET')

      if (Schm_sleve_L) call gem_error (-1,'INIT_BAR','  SLEVE not available YET')

      istat = gmm_get (gmmk_pw_uu_plus_s, pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s, pw_vv_plus)
      istat = gmm_get (gmmk_pw_tt_plus_s, pw_tt_plus)
      istat = gmm_get (gmmk_ut1_s ,ut1 )
      istat = gmm_get (gmmk_vt1_s ,vt1 )
      istat = gmm_get (gmmk_wt1_s ,wt1 )
      istat = gmm_get (gmmk_tt1_s ,tt1 )
      istat = gmm_get (gmmk_zdt1_s,zdt1)
      istat = gmm_get (gmmk_st1_s ,st1 )
      istat = gmm_get (gmmk_sls_s ,sls )
      istat = gmm_get (gmmk_fis0_s,fis0)
      istat = gmm_get (gmmk_qt1_s ,qt1 )

      !Setup Williamson Case 7: The 21 December 1978 Initial conditions are read
      !-------------------------------------------------------------------------
      if (Williamson_case==7) &
          call inp_data ( ut1, vt1, wt1, tt1  ,&
                          zdt1,st1,fis0,l_minx,l_maxx,l_miny,l_maxy,&
                          G_nk,.true. ,'TR/',':P',Step_runstrt_S )

      !Initialize T/ZD/W/Q
      !-------------------
      tt1 = Cstv_Tstr_8 ; zdt1 = 0. ; wt1 = 0. ; qt1 = 0.

      !Initialize HU
      !-------------
      istat = gmm_get ('TR/HU:P',hu)
      hu = 0.

      !Initialize d(Zeta)dot and dz/dt
      !-------------------------------
      Inp_zd_L = .true.
      Inp_w_L  = .true.

      !Prepare initial conditions (staggered u-v,gz,s,topo) for Williamson cases
      !-------------------------------------------------------------------------
      call wil_init (ut1,vt1,gz_t,st1,fis0,'TR/',':P',l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

      !Required for CASE5 LAM version
      !------------------------------
      istat = gmm_get (gmmk_topo_low_s , topo_low )
      istat = gmm_get (gmmk_topo_high_s, topo_high)

      topo_high(1:l_ni,1:l_nj) =      fis0(1:l_ni,1:l_nj)
      topo_low (1:l_ni,1:l_nj) = topo_high(1:l_ni,1:l_nj)

      !Estimate U-V and T on scalar grids
      !----------------------------------
      call hwnd_stag ( pw_uu_plus,pw_vv_plus,ut1,vt1, &
                       l_minx,l_maxx,l_miny,l_maxy,G_nk,.false. )

      pw_tt_plus = tt1

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H' .and. .not.Ctrl_testcases_adv_L) then
         do k=1,G_nk
            qt1(1:l_ni ,1:l_nj, k) = max(gz_t(1:l_ni,1:l_nj,1), 0.)
         end do
      end if
!
!-------------------------------------------------------------------
!
      return
      end subroutine init_bar
