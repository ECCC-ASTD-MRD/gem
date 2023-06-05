!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 
!/@*
subroutine test_config_mpi2
   use testutils
   use config_mod
   use wb_itf_mod
   use clib_itf_mod
   implicit none
   !@objective Test config_mod module functions
   !@author Stephane Chamberland, July 2010
   !@revisions
   !  2012-02, Stephane Chamberland: to new config_mod and new testutils
!*@/
#include <rmnlib_basics.hf>

   integer :: istat,myproc,npx
   character(len=1024) :: basename_S,secname_S,wkdir_S,cwd_S
   !---------------------------------------------------------------------
   myproc = testutils_initmpi()
   call testutils_verbosity()
   call testutils_set_name('test_config_mpi2')

   basename_S = 'test_config'
   secname_S = 'ptopo_cfgs'

   istat = clib_getcwd(cwd_S)
   call config_init(cwd_S)

   write(wkdir_S,'(a,i3.3)') '__'//trim(basename_S)//'-to-remove__',myproc
   istat = clib_mkdir(trim(wkdir_S))
   istat = clib_chdir(trim(wkdir_S))

   istat = config_read(basename_S,secname_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_read status')
   istat = wb_get(trim(secname_S)//'/npx',npx)
   call testutils_assert_eq(npx,2,'config_read npx')

   istat = clib_chdir(trim(cwd_S))
   istat = clib_unlink(trim(wkdir_S)//'/'//trim(basename_S)//'.dict')
   istat = clib_unlink(trim(wkdir_S)//'/'//trim(basename_S)//'.cfg')
   istat = clib_rmdir(trim(wkdir_S))

   call rpn_comm_finalize(istat)
   !---------------------------------------------------------------------
end subroutine test_config_mpi2
