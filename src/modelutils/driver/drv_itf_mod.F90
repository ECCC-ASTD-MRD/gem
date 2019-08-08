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

!/@*
module drv_itf_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use clib_itf_mod, only: clib_tolower, clib_getenv
   use wb_itf_mod
   use env_utils
   use str_mod
   use config_mod, only: config_init
   use drv_ptopo_mod, only: drv_ptopo_init,drv_ptopo_terminate
   use drv_path_mod, only: drv_path_config_dir0_S
   use drv_time_mod, only: drv_time_config,drv_time_init,drv_time_info,drv_time_increment
   use drv_grid_mod, only: drv_grid_config,drv_grid_init
   use drv_levels_mod, only: drv_levels_config,drv_levels_init
   implicit none
   private
   !@objective Main Driver API
   !@author  Stephane Chamberland, 2012-03
   ! Public functions
   public :: drv_config,drv_init,drv_verbosity, &
        drv_ptopo_init,drv_ptopo_terminate, &
        drv_time_info,drv_time_increment
   ! 
   ! Public constants
   ! 
   ! Public var
   !@Description
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <mu_gmm.hf>
#include <msg.h>

contains

   !/@*
   function drv_config(F_cfg_basename_S) result(F_istat)
      implicit none
      !@arguments
      character(len=*),intent(out) :: F_cfg_basename_S
      !@return
      integer :: F_istat
      !*@/
      !---------------------------------------------------------------------
      F_istat = clib_getenv('UM_EXEC_CONFIG_BASENAME',F_cfg_basename_S)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(drv_config) UM_EXEC_CONFIG_BASENAME not defiened, trying with default settings file, model_settings')
         F_cfg_basename_S = './model_settings'
      endif
      call config_init(drv_path_config_dir0_S)
      F_istat = drv_grid_config(F_cfg_basename_S)
      F_istat = min(drv_levels_config(F_cfg_basename_S),F_istat)
      F_istat = min(drv_time_config(F_cfg_basename_S),F_istat)
      !---------------------------------------------------------------------
      return
   end function drv_config


   !/@*
   function drv_init(F_dateo_S,F_dt_8,F_stepno,F_ntimelevels) result(F_istat)
      implicit none
      !@arguments
      character(len=*),intent(out) :: F_dateo_S
      real(REAL64),intent(out) :: F_dt_8
      integer,intent(out) :: F_stepno
      integer,intent(in) :: F_ntimelevels
      !@return
      integer :: F_istat
      !*@/
      logical :: is_chkpt,is_last
      !---------------------------------------------------------------------
      F_istat = drv_grid_init()
      F_istat = min(drv_levels_init(),F_istat)
      F_istat = min(drv_time_init(F_dateo_S,F_dt_8,F_stepno,is_chkpt,is_last,F_ntimelevels),F_istat)
      !---------------------------------------------------------------------
      return
   end function drv_init


   !/@*
   subroutine drv_verbosity(F_prefix_S)
      implicit none
      !@objective Set verbosity
      !@argument
      character(len=*),intent(in) :: F_prefix_S
      !@author  Stephane Chamberland, 2011-09
      !*@/
      integer :: istat,myproc,mycol,myrow
      character(len=256) :: tmp_S,prefix_S
      logical :: debug_L
      !---------------------------------------------------------------------
      prefix_S = F_prefix_S
      if (len_trim(prefix_S) > 0) prefix_S = trim(prefix_S)//'_'

      istat = clib_getenv(trim(prefix_S)//'VERBOSITY_PROC',tmp_S)
      if (RMN_IS_OK(istat)) then
         call str_tab2space(tmp_S)
         tmp_S = adjustl(tmp_S)
         istat = clib_tolower(tmp_S)
         select case(tmp_S(1:1))
         case('0') !- print from pe0 only
            istat = rpn_comm_mype(myproc,mycol,myrow)
         case('a') !- print from all proc
            myproc = 0
         case('b') !- print from bloc master only
            call rpn_comm_rank(RPN_COMM_BLOC_COMM,myproc,istat)
         end select
      else
         istat = rpn_comm_mype(myproc,mycol,myrow)
      endif
      call msg_set_p0only(max(0,myproc))

      istat = env_get('UM_EXEC_DEBUG',debug_L,.false.)
      if (debug_L) then
         tmp_S = 'debug'
      else
         istat = clib_getenv(trim(prefix_S)//'VERBOSITY',tmp_S)
      endif
      if (.not.RMN_IS_OK(istat)) tmp_S='info'
      call str_tab2space(tmp_S)
      tmp_S = adjustl(tmp_S)
      istat = clib_tolower(tmp_S)
      select case(tmp_S(1:1))
      case('d')
         tmp_S='debug'
         istat = fstopc('MSGLVL','DEBUG',RMN_OPT_SET)
         istat = wb_verbosity(WB_MSG_DEBUG)
         istat = gmm_verbosity(GMM_MSG_DEBUG)
         call msg_set_minMessageLevel(MSG_DEBUG)
         call handle_error_setdebug(debug_L)
      case('p')
         tmp_S='plus'
         istat = fstopc('MSGLVL','INFORM',RMN_OPT_SET)
         istat = wb_verbosity(WB_MSG_WARN)
         istat = gmm_verbosity(GMM_MSG_ERROR)
         call msg_set_minMessageLevel(MSG_INFOPLUS)
      case('w')
         tmp_S='warn'
         istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
         istat = wb_verbosity(WB_MSG_WARN)
         istat = gmm_verbosity(GMM_MSG_ERROR)
         call msg_set_minMessageLevel(MSG_WARNING)
      case('e')
         tmp_S='error'
         istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
         istat = wb_verbosity(WB_MSG_ERROR)
         istat = gmm_verbosity(GMM_MSG_ERROR)
         call msg_set_minMessageLevel(MSG_ERROR)
      case('c')
         tmp_S='critical'
         istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
         istat = wb_verbosity(WB_MSG_ERROR)
         istat = gmm_verbosity(GMM_MSG_FATAL)
         call msg_set_minMessageLevel(MSG_CRITICAL)
      case default
         tmp_S='info'
         istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
         istat = wb_verbosity(WB_MSG_WARN)
         istat = gmm_verbosity(GMM_MSG_ERROR)
         call msg_set_minMessageLevel(MSG_INFO)
      end select
      call msg(MSG_INFO,'(drv) Set Verbosity Level='//trim(tmp_S))
      !---------------------------------------------------------------------
      return
   end subroutine drv_verbosity


end module drv_itf_mod
