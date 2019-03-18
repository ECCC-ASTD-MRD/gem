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

      subroutine init_bar ( F_u, F_v, F_w, F_t, F_zd, F_s, F_q, F_topo,&
                            Mminx,Mmaxx,Mminy,Mmaxy, Nk               ,&
                            F_stag_L,F_trprefix_S, F_trsuffix_S, F_datev )

      use dynkernel_options
      use gmm_geof
      use inp_mod
      use gmm_pw
      use dyn_fisl_options
      use ctrl
      use gem_options
      use wil_options
      use glb_ld
      use cstv
      use tr3d
      use gmm_itf_mod

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      character(len=*) F_trprefix_S, F_trsuffix_S, F_datev
      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk
      real F_u (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_v (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_w (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_t (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_zd(Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_q (Mminx:Mmaxx,Mminy:Mmaxy,Nk+1), &
           F_topo(Mminx:Mmaxx,Mminy:Mmaxy)

      logical F_stag_L

      !object
      !============================================================
      !     prepare data for autobarotropic runs (Williamson cases)
      !============================================================

      !---------------------------------------------------------------

      integer :: istat, k
      real, dimension (:,:,:), pointer :: hu
      real, dimension (Mminx:Mmaxx,Mminy:Mmaxy,Nk) :: gz_t

      !---------------------------------------------------------------

      if (Vtopo_L)      call handle_error (-1,'INIT_BAR','Vtopo_L not available YET')

      if (Schm_sleve_L) call handle_error (-1,'INIT_BAR','  SLEVE not available YET')

      !Setup Williamson Case 7: The 21 December 1978 Initial conditions are read
      !-------------------------------------------------------------------------
      if (Williamson_case==7) &
      call inp_data ( F_u, F_v, F_w, F_t, F_zd, F_s, F_q, F_topo,&
                      Mminx,Mmaxx,Mminy,Mmaxy, Nk,F_stag_L      ,&
                      F_trprefix_S, F_trsuffix_S, F_datev )

      !Initialize T/ZD/W/Q
      !-------------------
      F_t = Cstv_Tstr_8 ; F_zd = 0. ; F_w = 0. ; F_q = 0.

      !Initialize HU
      !-------------
      istat = gmm_get (trim(F_trprefix_S)//"HU"//trim(F_trsuffix_S),hu)

      hu = 0.

      !Initialize d(Zeta)dot and dz/dt
      !-------------------------------
      Inp_zd_L = .true.
      Inp_w_L  = .true.

      !Prepare initial conditions (staggered u-v,gz,s,topo) for Williamson cases
      !-------------------------------------------------------------------------
      call wil_init (F_u,F_v,gz_t,F_s,F_topo,F_trprefix_S,F_trsuffix_S,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

      !Required for CASE5 LAM version
      !------------------------------
      istat = gmm_get (gmmk_topo_low_s , topo_low )
      istat = gmm_get (gmmk_topo_high_s, topo_high)

      topo_high(1:l_ni,1:l_nj) =    F_topo(1:l_ni,1:l_nj)
      topo_low (1:l_ni,1:l_nj) = topo_high(1:l_ni,1:l_nj)

      !Estimate U-V and T on scalar grids
      !----------------------------------
      if (F_stag_L) then

         istat = gmm_get (gmmk_pw_uu_plus_s, pw_uu_plus)
         istat = gmm_get (gmmk_pw_vv_plus_s, pw_vv_plus)
         istat = gmm_get (gmmk_pw_tt_plus_s, pw_tt_plus)

         call hwnd_stag ( pw_uu_plus,pw_vv_plus,F_u,F_v, &
                          Mminx,Mmaxx,Mminy,Mmaxy,Nk,.false. )

         pw_tt_plus = F_t

      end if

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H' .and. .not.Ctrl_testcases_adv_L) then
         do k=1,G_nk
            F_q(1:l_ni ,1:l_nj, k) = max(gz_t(1:l_ni,1:l_nj,1), 0.)
         end do
      end if

      return

      end
