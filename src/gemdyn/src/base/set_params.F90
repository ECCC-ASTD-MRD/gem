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

!**   s/r set_params - initialize some constant parameters

      subroutine set_params (F_check_and_stop_L, F_errcode)
      use dynkernel_options
      use dyn_fisl_options
      use tdpack
      use cstv
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
      
      logical, intent(in) :: F_check_and_stop_L
      integer, intent(inout) :: F_errcode

      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0
!
!     ---------------------------------------------------------------

      Cstv_tau_8   = Cstv_dt_8 * Cstv_bA_8
      Cstv_invT_8  = one/Cstv_tau_8
      Cstv_Beta_8  = (one-Cstv_bA_8)/Cstv_bA_8

      Cstv_tau_m_8   = Cstv_dt_8 * Cstv_bA_m_8
      Cstv_invT_m_8  = one/Cstv_tau_m_8
      Cstv_Beta_m_8  = (one-Cstv_bA_m_8)/Cstv_bA_m_8

!     Parameters for the nonhydrostatic case
      Cstv_tau_nh_8   = Cstv_dt_8 * Cstv_bA_nh_8
      Cstv_invT_nh_8  = one/Cstv_tau_nh_8
      Cstv_Beta_nh_8  = (one-Cstv_bA_nh_8)/Cstv_bA_nh_8

      if (Schm_advec == 1) then ! traditional advection
         Cstv_dtA_8  = Cstv_dt_8 * 0.5d0
         Cstv_dtzA_8 = Cstv_dt_8 * 0.5d0
      end if
      if (Schm_advec == 2) then ! consistent advection
         Cstv_dtA_8  = Cstv_tau_m_8
         Cstv_dtzA_8 = Cstv_tau_8
      end if
      if (Schm_advec == 3) then ! reversed advection
         Cstv_dtA_8  = (one-Cstv_bA_m_8)*Cstv_dt_8
         Cstv_dtzA_8 = (one-Cstv_bA_8)*Cstv_dt_8
      end if

      Cstv_dtD_8  = Cstv_dt_8 - Cstv_dtA_8
      Cstv_dtzD_8 = Cstv_dt_8 - Cstv_dtzA_8

      if (Schm_advec == 0) then ! no advection
         Cstv_dtA_8  = 0.d0
         Cstv_dtD_8  = 0.d0
         Cstv_dtzA_8 = 0.d0
         Cstv_dtzD_8 = 0.d0
      end if

      Ver_igt_8    = Cstv_invT_8/grav_8
      Ver_ikt_8    = Cstv_invT_m_8/cappa_8
      if(Dynamics_hydro_L) Ver_igt_8=zero

      Ver_igt2_8   = Ver_igt_8*(Cstv_invT_nh_8/grav_8)

      call set_dync ( F_check_and_stop_L, F_errcode )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine set_params
