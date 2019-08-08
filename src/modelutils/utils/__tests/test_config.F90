!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 
!/@*
subroutine test_config
   use testutils
   use config_mod
   use wb_itf_mod
   use clib_itf_mod
   implicit none
   !@objective Test config_mod module functions
   !@author Stephane Chamberland, Feb 2012
!*@/
#include <rmnlib_basics.hf>
   character(len=RMN_PATH_LEN) :: basename_S, dictname_S, secname_S,cwd_S,wkdir_S,path_S
   integer :: istat,npx
   !---------------------------------------------------------------------
   call testutils_verbosity()
   basename_S = 'test_config'
   call testutils_set_name(basename_S)

   dictname_S = trim(basename_S)//'.dict'
   secname_S = 'ptopo_cfgs'

   istat = config_read('__no-such_file__',secname_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'config_read no-such_file')


   istat = config_read(basename_S,'__no-such-section__')
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'config_read no-such-section')

   istat = config_read(basename_S,secname_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_read status')
   istat = wb_get(trim(secname_S)//'/npx',npx)
   call testutils_assert_eq(npx,2,'config_read npx')
   call config_print(secname_S)

   istat = clib_getcwd(cwd_S)

   istat = config_cp2localdir(dictname_S, cwd_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_cp cwd status')

   wkdir_S = '__'//trim(basename_S)//'-to-remove__'
   istat = clib_mkdir(trim(wkdir_S))
   istat = clib_chdir(trim(wkdir_S))

   istat = clib_unlink(trim(wkdir_S)//'/'//trim(dictname_S))
   istat = clib_isfile(trim(dictname_S))
   call testutils_assert_eq(RMN_IS_OK(istat),.false.,'config_cp pre isfile')

   istat = config_path(path_S,basename_S,'dict',cwd_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_path')
   istat = clib_isfile(trim(dictname_S))
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_path isfile')

   istat = clib_unlink(trim(wkdir_S)//'/'//trim(dictname_S))

   istat = config_cp2localdir(dictname_S, cwd_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_cp status')
   istat = clib_isfile(dictname_S)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_cp isfile')

!!$   istat = clib_unlink(trim(wkdir_S)//'/'//trim(dictname_S))
!!$
!!$   istat = config_read(basename_S,secname_S,F_basedir_S=cwd_S)
!!$   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_read 2 status')
!!$   istat = clib_isfile(dictname_S)
!!$   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'config_cp isfile')
!!$   istat = wb_get('ptopo_cfgs/npx',npx)
!!$   call testutils_assert_eq(npx,2,'config_read 2 npx')

   istat = clib_chdir(trim(cwd_S))
   istat = clib_unlink(trim(wkdir_S)//'/'//trim(dictname_S))
   istat = clib_rmdir(trim(wkdir_S))
   !---------------------------------------------------------------------
end subroutine test_config
