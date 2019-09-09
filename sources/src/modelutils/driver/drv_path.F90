!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

!/@
module drv_path_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use clib_itf_mod, only: clib_getenv, clib_mkdir, clib_isdir, clib_getcwd
   use wb_itf_mod
   implicit none
   private
   !@objective Define/create path task_setup style (following GEM 4.4)
   !@author  Stephane Chamberland, 2012-02
   ! Public functions
   public :: drv_path_set_basedir,drv_path_set
   ! Public constants
   integer,parameter :: DRV_PATH_LEN = 512
   ! Public var
   character(len=DRV_PATH_LEN),public,save :: &
        drv_path_basedir_S = '.', &
        drv_path_config_dir0_S = '.', &
        drv_path_config_dir_S = '.', &
        drv_path_work_domaine_S = '.', &
        drv_path_work_grid_S = '.', &
        drv_path_work_S = '.', &
        drv_path_input_domaine_S = '.', &
        drv_path_input_S = '.', &
        drv_path_output_S = '.'
   !@Description
   !@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

contains

   !/@
   function drv_path_set_basedir(F_basedir_S) result(F_istat)
      implicit none
      !@arguments
      character(len=*),intent(in),optional :: F_basedir_S
      !@return
      integer :: F_istat
      !@/
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (present(F_basedir_S)) then
         F_istat = clib_isdir(trim(adjustl(F_basedir_S)))
         if (RMN_IS_OK(F_istat)) drv_path_basedir_S = adjustl(F_basedir_S)
      else if (drv_path_basedir_S == '.') then
         F_istat = clib_getenv('TASK_BASEDIR'  ,drv_path_basedir_S)
         if (.not.RMN_IS_OK(F_istat)) then
            F_istat = clib_getcwd(drv_path_basedir_S)
         else
            F_istat = clib_isdir(drv_path_basedir_S)
            if (.not.RMN_IS_OK(F_istat)) drv_path_basedir_S = '.'
         endif
      else
         F_istat = RMN_OK
         return
      endif
      if (RMN_IS_OK(F_istat)) then
         F_istat = wb_put('path/basedir',drv_path_basedir_S,WB_REWRITE_MANY)
      endif
      !---------------------------------------------------------------------
      return
   end function drv_path_set_basedir


   !/@
   function drv_path_set(F_mydomain,F_ngrids,F_igrid,F_ipex,F_ipey) result(F_istat)
      implicit none
      !@arguments
      integer,intent(in) :: F_mydomain
      integer,intent(in),optional :: F_ngrids,F_igrid,F_ipex,F_ipey
      !@return
      integer :: F_istat
      !@/
      character(len=512) :: tmp_S
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] drv_path_set')
      F_istat = clib_getenv('TASK_WORK',drv_path_work_S)
      F_istat = min(clib_getenv('TASK_INPUT' ,drv_path_input_S),F_istat)
      F_istat = min(clib_getenv('TASK_OUTPUT',drv_path_output_S),F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(drv_path_set) Undefined: TASK_WORK, TASK_INPUT, TASK_OUTPUT')
         return
      endif

      tmp_S = ' '
      call priv_add_path_element('input',drv_path_input_S,tmp_S)
      call priv_add_path_element('output',drv_path_output_S,tmp_S)
      call priv_add_path_element('work',drv_path_work_S,tmp_S)


      write(tmp_S,'(a,i4.4)') 'cfg_',F_mydomain
      call priv_add_path_element('input',drv_path_input_S,tmp_S)
      call priv_add_path_element('output',drv_path_output_S,tmp_S)
      call priv_add_path_element('work',drv_path_work_S,tmp_S)
      drv_path_work_domaine_S = drv_path_work_S
      call priv_add_path_element('work_domaine',drv_path_work_domaine_S,' ')
      drv_path_input_domaine_S = drv_path_input_S
      call priv_add_path_element('input_domaine',drv_path_input_domaine_S,' ')
      drv_path_config_dir0_S = drv_path_input_domaine_S
      call priv_add_path_element('config_dir0',drv_path_config_dir0_S,' ')

      if (present(F_ngrids) .and. present(F_igrid)) then
         if (F_ngrids > 1) then
            write(tmp_S,'(i3.3)') F_igrid
            call priv_add_path_element('input',drv_path_input_S,tmp_S)
            call priv_add_path_element('output',drv_path_output_S,tmp_S)
            call priv_add_path_element('work',drv_path_work_S,tmp_S)
         endif
      endif
      drv_path_work_grid_S = drv_path_work_S
      call priv_add_path_element('work_grid',drv_path_work_grid_S,' ')

      if (present(F_ipex) .and. present(F_ipey)) then
         write(tmp_S,'(i3.3,a1,i3.3)') F_ipex,'-',F_ipey
         call priv_add_path_element('output',drv_path_output_S,tmp_S)
         call priv_add_path_element('work',drv_path_work_S,tmp_S)
      endif
      drv_path_config_dir_S = drv_path_work_S
      call priv_add_path_element('config_dir',drv_path_config_dir_S,' ')

!!$      path_ind_S = trim(path_input_S)//'/MODEL_INPUT' !TODO-later:
!!$      path_ind_S = trim(path_input_S)//'/MODEL_INPUT/'//trim(Grd_yinyang_S)
!!$      path_phy_S = trim(path_input_S)//'/'
!!$      err = clib_getcwd(pwd_S)
!!$      Path_nml_S    = trim(pwd_S)//'/model_settings'
!!$      Path_outcfg_S = trim(pwd_S)//'/output_settings'
!!$      Path_phyincfg_S = trim(pwd_S)//'/physics_input_table'
      call msg(MSG_DEBUG,'[END] drv_path_set')
      !---------------------------------------------------------------------
      return
   end function drv_path_set


   !/@
   subroutine priv_add_path_element(F_name_S,F_path_S,F_part_S)
      implicit none
      !@arguments
      character(len=*),intent(inout) ::F_path_S
      character(len=*),intent(in) ::  F_name_S,F_part_S
      !@/
      integer :: istat
      character(len=512) :: part_S
      !---------------------------------------------------------------------
      if (F_part_S /= ' ') then
         part_S = adjustl(F_part_S)
         F_path_S  = adjustl(trim(F_path_S) // '/' // trim(part_S))
      endif
      istat = clib_mkdir(trim(F_path_S))
      istat = clib_isdir(trim(F_path_S))
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(drv_path_set) Dir Not Found: '//trim(F_path_S))
      endif
      istat = wb_put('path/'//trim(F_name_S),F_path_S,WB_REWRITE_MANY)
      !---------------------------------------------------------------------
      return
   end subroutine priv_add_path_element

end module drv_path_mod
