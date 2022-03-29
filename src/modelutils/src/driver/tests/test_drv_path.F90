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

!/@
subroutine test_drv_path()
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod
   use wb_itf_mod
   use testutils
   use drv_path_mod
   implicit none
   !@objective
   !@author Stephane Chamberland, 2012-02
!@/
#include <rmnlib_basics.hf>
#include <msg.h>
   character(len=RMN_PATH_LEN) :: basedir_S,work_S,input_S,output_S,work_domaine_S,work_grid_S,tmp_S,input_domaine_S

   integer :: istat,npx,npy
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_drv_path')

   istat = clib_getcwd(basedir_S)
   basedir_S = trim(basedir_S)//'/__test_drv_path-to-delete__'
   istat = clib_mkdir(trim(basedir_S))
   istat = clib_chdir(trim(basedir_S))

   istat = drv_path_set_basedir()
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'basedir status')
   call testutils_assert_eq(drv_path_basedir_S,basedir_S,'basedir value')

   !-
   istat = drv_path_set(5)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'set no var')
 
   !-
   work_S = trim(basedir_S)//'/work'
   input_S = trim(basedir_S)//'/input'
   output_S = trim(basedir_S)//'/output'
   istat = clib_putenv('TASK_WORK='//trim(work_S))
   istat = clib_putenv('TASK_INPUT='//trim(input_S))
   istat = clib_putenv('TASK_OUTPUT='//trim(output_S))

   tmp_S = '/cfg_0005'
   work_S = trim(work_S)//trim(tmp_S)
   input_S = trim(input_S)//trim(tmp_S)
   output_S = trim(output_S)//trim(tmp_S)
   work_domaine_S = work_S
   work_grid_S = work_S
   input_domaine_S = input_S

   istat = drv_path_set(5)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'set dom')
   call testutils_assert_eq(drv_path_work_domaine_S,work_domaine_S,'set dom wd')
   call testutils_assert_eq(drv_path_work_grid_S,work_grid_S,'set dom wg')
   call testutils_assert_eq(drv_path_work_S,work_S,'set dom w')
   call testutils_assert_eq(drv_path_input_domaine_S,input_domaine_S,'set dom i')
   call testutils_assert_eq(drv_path_input_S,input_S,'set dom id')
   call testutils_assert_eq(drv_path_output_S,output_S,'set dom o')
   call testutils_assert_eq(drv_path_config_dir0_S,work_domaine_S,'set dom cfg0')
   call testutils_assert_eq(drv_path_config_dir_S,work_S,'set dom cfg')
   istat = clib_isdir(trim(drv_path_work_domaine_S))
   istat = min(clib_isdir(trim(drv_path_work_grid_S)),istat)
   istat = min(clib_isdir(trim(drv_path_work_S)),istat)
   istat = min(clib_isdir(trim(drv_path_input_S)),istat)
   istat = min(clib_isdir(trim(drv_path_output_S)),istat)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'set dom dirs')

   !-
   work_S = trim(basedir_S)//'/work'
   input_S = trim(basedir_S)//'/input'
   output_S = trim(basedir_S)//'/output'
   tmp_S = '/cfg_0005'
   input_S = trim(input_S)//trim(tmp_S)
   tmp_S = '/cfg_0005/003-007'
   work_S = trim(work_S)//trim(tmp_S)
   output_S = trim(output_S)//trim(tmp_S)
   input_domaine_S = input_S

   istat = drv_path_set(5,1,1,3,7)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'set 1 grid status')
   call testutils_assert_eq(drv_path_work_domaine_S,work_domaine_S,'set 1 grid wd')
   call testutils_assert_eq(drv_path_work_grid_S,work_grid_S,'set 1 grid wg')
   call testutils_assert_eq(drv_path_work_S,work_S,'set 1 grid w')
   call testutils_assert_eq(drv_path_input_domaine_S,input_domaine_S,'set 1 grid id')
   call testutils_assert_eq(drv_path_input_S,input_S,'set 1 grid i')
   call testutils_assert_eq(drv_path_output_S,output_S,'set 1 grid o')
   call testutils_assert_eq(drv_path_config_dir0_S,work_domaine_S,'set 1 grid cfg0')
   call testutils_assert_eq(drv_path_config_dir_S,work_S,'set 1 grid cfg')
   istat = clib_isdir(trim(drv_path_work_domaine_S))
   istat = min(clib_isdir(trim(drv_path_work_grid_S)),istat)
   istat = min(clib_isdir(trim(drv_path_work_S)),istat)
   istat = min(clib_isdir(trim(drv_path_input_S)),istat)
   istat = min(clib_isdir(trim(drv_path_output_S)),istat)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'set 1 grid dirs')

   !-
   work_S = trim(basedir_S)//'/work'
   input_S = trim(basedir_S)//'/input'
   output_S = trim(basedir_S)//'/output'
   tmp_S = '/cfg_0005/001'
   work_grid_S = trim(work_S)//trim(tmp_S)
   input_S = trim(input_S)//trim(tmp_S)
   tmp_S = '/cfg_0005/001/003-007'
   work_S = trim(work_S)//trim(tmp_S)
   output_S = trim(output_S)//trim(tmp_S)

   istat = drv_path_set(5,2,1,3,7)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'set 2 grid status')
   call testutils_assert_eq(drv_path_work_domaine_S,work_domaine_S,'set 2 grid wd')
   call testutils_assert_eq(drv_path_work_grid_S,work_grid_S,'set 2 grid wg')
   call testutils_assert_eq(drv_path_work_S,work_S,'set 2 grid w')
   call testutils_assert_eq(drv_path_input_domaine_S,input_domaine_S,'set 2 grid id')
   call testutils_assert_eq(drv_path_input_S,input_S,'set 2 grid i')
   call testutils_assert_eq(drv_path_output_S,output_S,'set 2 grid o')
   call testutils_assert_eq(drv_path_config_dir0_S,work_domaine_S,'set 2 grid cfg0')
   call testutils_assert_eq(drv_path_config_dir_S,work_S,'set 2 grid cfg')
   istat = clib_isdir(trim(drv_path_work_domaine_S))
   istat = min(clib_isdir(trim(drv_path_work_grid_S)),istat)
   istat = min(clib_isdir(trim(drv_path_work_S)),istat)
   istat = min(clib_isdir(trim(drv_path_input_S)),istat)
   istat = min(clib_isdir(trim(drv_path_output_S)),istat)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'set 2 grid dirs')
   ! ---------------------------------------------------------------------
   return
end subroutine test_drv_path
