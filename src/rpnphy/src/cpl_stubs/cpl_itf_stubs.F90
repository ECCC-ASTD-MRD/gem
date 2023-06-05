!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer, 
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms 
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer 
!version 3 or (at your option) any later version that should be found at: 
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html 
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software; 
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec), 
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

module cpl_itf

  ! Subprogram elements of the API
!  use cpl_nml_mod
!  use cpl_init_mod
!  use cpl_step_mod
!  use cpl_update_mod
!  use cpl_restart_mod
!  use cpl_terminate_mod

  implicit none

  !@authors Desgagne, Roy -- Spring 2015
  !@revision

  !@objective 
  ! Fill the external interface (API) for the coupler package

  ! Allow access to the full API
  public

  ! Externally-accessible constants and variables
  integer, parameter :: CPL_COMPATIBILITY_LVL = 1

contains
      integer function cpl_nml (F_namelistf_S, F_unout)
      implicit none
!!!#include <arch_specific.hf>
      character(len=*), intent(in) :: F_namelistf_S
      integer, intent(in)          :: F_unout
      cpl_nml= 0
      if (F_namelistf_S == '' .and. F_unout >=0) print *,'cpl_nml (stub)',F_namelistf_S
      return
      end function cpl_nml

      integer function cpl_init (F_path_S, F_print_L, F_unout, F_dateo, F_dt, &
                                 F_glb_ni, F_glb_nj, F_lcl_ni, F_lcl_nj     , &
                                 F_i0, F_j0, F_in, F_jn, F_nk               , &
                                 F_z0mtype, F_z0ttype, F_z0lat)
      implicit none
!!!#include <arch_specific.hf>
      character(len=*), intent(in) :: F_path_S,F_z0mtype,F_z0ttype
      logical, intent(in)          :: F_print_L
      integer, intent(in)          :: F_dateo, F_unout                , &
                                   F_glb_ni,F_glb_nj,F_lcl_ni,F_lcl_nj, &
                                   F_i0, F_j0, F_in, F_jn, F_nk
      real   , intent(in)          :: F_dt,F_z0LAT(:)
      cpl_init= 0
      if (F_path_S == '' .and. F_print_L .and. F_unout > 0) print *,'cpl_init (stub)', &
           F_path_S, F_print_L, F_unout, F_dateo, F_dt, &
           F_glb_ni, F_glb_nj, F_lcl_ni, F_lcl_nj     , &
           F_i0, F_j0, F_in, F_jn, F_nk               , &
           F_z0mtype, F_z0ttype, F_z0lat
      return
      end function cpl_init

      subroutine cpl_step (F_stepcount, F_stepdriver)
      implicit none
!!!#include <arch_specific.hf>
      integer, intent(in) :: F_stepcount, F_stepdriver
      if (F_stepcount < 0 .and. F_stepdriver < 0) print *,'cpl_step (stub)', F_stepcount, F_stepdriver
      return
      end subroutine cpl_step

      subroutine cpl_update (F_f2u, F_name_S, F_pos, F_ni, rho, u, v, vmod, cmu, cplu)
      implicit none
!!!#include <arch_specific.hf>

      character(len=*), intent(in) :: F_name_S
      integer, intent(in)          :: F_ni
      integer, dimension(2,F_ni),           intent(in)    :: F_pos
      real   , dimension  (F_ni),           intent(inout) :: F_f2u
      real   , dimension  (F_ni), optional, target,&
                                            intent(in)    :: rho, u, v, vmod
      real   , dimension  (F_ni), optional, target,&
                                            intent(inout) :: cmu
      logical, optional, intent(out) :: cplu

      cplu = .false.
      if (CPL_COMPATIBILITY_LVL == 0) print *,'cpl_update (stub)', F_f2u, F_name_S, F_pos, F_ni, rho, u, v, vmod, cmu
      return
      end subroutine cpl_update

      subroutine cpl_restart (F_WorR_S)
      implicit none
!!!#include <arch_specific.hf>
      character(len=*),intent(in) :: F_WorR_S
      if (CPL_COMPATIBILITY_LVL == 0) print *,'cpl_restart (stub)', F_WorR_S
      return
      end subroutine cpl_restart

      subroutine cpl_terminate (F_pcomm_L)
      implicit none
!!!#include <arch_specific.hf>
      logical, intent(in) :: F_pcomm_L
      if (CPL_COMPATIBILITY_LVL == 0) print *,'cpl_terminate (stub)', F_pcomm_L
      return
      end subroutine cpl_terminate

      subroutine cpl_snapshot (F_mode)
      implicit none
!!!#include <arch_specific.hf>
      character(len=*), intent(in) :: F_mode
      if (CPL_COMPATIBILITY_LVL == 0) print *,'cpl_snapshot (stub)', F_mode
      return
      end subroutine cpl_snapshot

end module cpl_itf
