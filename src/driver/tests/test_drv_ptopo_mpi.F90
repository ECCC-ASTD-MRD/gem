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
subroutine test_drv_ptopo_mpi()
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod
   use testutils
   use drv_ptopo_mod
   use env_utils, only: env_get
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-02
!@/
#include <rmnlib_basics.hf>
#include <msg.h>

   integer,parameter :: NGRIDS_DEFAULT = 1
   integer,parameter :: NGRIDS_MIN = 1
   integer,parameter :: NGRIDS_MAX = 10

!!$   character(len=RMN_PATH_LEN) :: basedir0_S,basedir_S,work_S,input_S,output_S,work_domaine_S,work_grid_S,tmp_S,input_domaine_S
   integer :: istat,ngrids !,npx,npy,ndoms
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_drv_ptopo')

!!$   istat = testutils_getenv_int('TEST_NDOMS',ndoms)
!!$
!!$   istat = clib_getcwd(basedir0_S)
!!$   basedir_S = trim(basedir0_S)//'/__test_drv_ptopo-to-delete__'
!!$   istat = clib_mkdir(trim(basedir_S))
!!$   istat = clib_chdir(trim(basedir_S))
!!$
!!$   work_S = trim(basedir_S)//'/work'
!!$   input_S = trim(basedir_S)//'/input'
!!$   output_S = trim(basedir_S)//'/output'
!!$   istat = clib_putenv('TASK_WORK='//trim(work_S))
!!$   istat = clib_putenv('TASK_INPUT='//trim(input_S))
!!$   istat = clib_putenv('TASK_OUTPUT='//trim(output_S))
!!$   if (ndoms == 1) then
!!$      istat = clib_putenv('UM_EXEC_NDOMAINS=1')
!!$   else
!!$      istat = clib_putenv('UM_EXEC_NDOMAINS=2:3')
!!$   endif
!!$   istat = clib_putenv('UM_EXEC_CONFIG_DICT='//trim(basedir0_S)//'test_drv_ptopo_mpi.dict')
!!$   istat = clib_putenv('UM_EXEC_CONFIG_USER='//trim(basedir0_S)//'test_drv_ptopo_mpi.cfg')

   istat = env_get('UM_EXEC_NGRIDS',ngrids,NGRIDS_DEFAULT,NGRIDS_MIN,NGRIDS_MAX)
   istat = drv_ptopo_init(ngrids)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'init')

   !TODO: test drv_ptopo_init result

   call drv_ptopo_terminate()
   ! ---------------------------------------------------------------------
   return
end subroutine test_drv_ptopo_mpi
