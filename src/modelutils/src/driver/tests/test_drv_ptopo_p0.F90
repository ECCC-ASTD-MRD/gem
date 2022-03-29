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
subroutine test_drv_ptopo_p0()
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod
   use wb_itf_mod
   use testutils
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-02
!@/
#include <rmnlib_basics.hf>
#include <msg.h>
   character(len=*),parameter :: section_S = 'ptopo_cfgs'
   character(len=*),parameter :: statusfile_S = './status_mod.dot'

   integer :: istat,npx,npy
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_drv_ptopo_p0')

   istat = wb_put('path/output','.',WB_REWRITE_MANY)
   istat = wb_put('path/config_dir0','.',WB_REWRITE_MANY)

   istat = clib_putenv('UM_EXEC_CONFIG_BASENAME=__no-such-file__')
   call drv_ptopo_p0(npx,npy)
   call testutils_assert_eq((/npx,npy/),(/0,0/),'no dict')
   call close_status_file3()
   istat = clib_isfile(statusfile_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'no dict status')
   istat = clib_unlink(statusfile_S)

!!$   istat = clib_putenv('UM_EXEC_CONFIG_BASENAME=test_drv_ptopo_p0')
!!$   call drv_ptopo_p0(npx,npy)
!!$   call testutils_assert_eq((/npx,npy/),(/1,1/),'no cfg')

   istat = clib_putenv('UM_EXEC_CONFIG_BASENAME=test_drv_ptopo_p0')
   call drv_ptopo_p0(npx,npy)
   call testutils_assert_eq((/npx,npy/),(/2,3/),'full cfg')
   call close_status_file3()
   istat = clib_isfile(statusfile_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'full cfg status')
   istat = clib_unlink(statusfile_S)

   istat = wb_put(section_S//'/npx',6,WB_REWRITE_MANY)
   istat = wb_put(section_S//'/npy',9,WB_REWRITE_MANY)
   call drv_ptopo_p0(npx,npy)
   call testutils_assert_eq((/npx,npy/),(/6,9/),'valid wb values')
   call close_status_file3()
   istat = clib_isfile(statusfile_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'valid wb status')
   istat = clib_unlink(statusfile_S)

   istat = wb_put(section_S//'/npx',-1,WB_REWRITE_MANY)
   istat = wb_put(section_S//'/npy',0,WB_REWRITE_MANY)
   call drv_ptopo_p0(npx,npy)
   call testutils_assert_eq((/npx,npy/),(/0,0/),'invalid wb values')
   call close_status_file3()
   istat = clib_isfile(statusfile_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'invalid wb values status')
   istat = clib_unlink(statusfile_S)
   ! ---------------------------------------------------------------------
   return
end subroutine test_drv_ptopo_p0
