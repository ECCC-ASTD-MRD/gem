!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 

!/@*
subroutine drv_ptopo_p0(F_npx,F_npy)
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_isdir, clib_isfile, clib_mkdir, clib_getenv
   use wb_itf_mod
   use config_mod
   implicit none
   !@objective Initialize processor topology
   !@arguments
   integer :: F_npx !O, number of procesor along X
   integer :: F_npy !O, number of procesor along Y
   !@author
   !  Michel Desgagne, Feb 2008
   !  Ron McTaggart-Cowan, Feb 2008
   !  Stephane Chamberland, Feb 2008
   !@revision
   !  2008-10, Stephane Chamberland: replace config_read by wb_read
   !  2012-02, Stephane Chamberland: update for rpnphy offline
   !@description
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

   logical,parameter :: IS_BEGIN_L = .true.
   character(len=*),parameter :: section_S = 'ptopo_cfgs'
   character(len=*),parameter :: statusfile_S = 'status_mod.dot'

   integer ::istat, nblocx, nblocy,ninblocx,ninblocy
   character(len=RMN_PATH_LEN) :: dict_S,user_S,tmp_S,tmp2_S,path_output_S
   logical :: is_ok
   !---------------------------------------------------------------------
!   if (F_npx>0 .and. F_npy>0) return

   !- Create status file
   istat = wb_get('path/output',path_output_S)
   if (.not.RMN_IS_OK(istat)) path_output_S = '.'
   istat = clib_isdir(path_output_S)
   if (.not.RMN_IS_OK(istat)) istat = clib_mkdir(path_output_S)
   call open_status_file3(trim(path_output_S)//'/'//trim(statusfile_S))
   call write_status_file3('_status=ABORT')

   call drv_print_banner(IS_BEGIN_L)

   !- Get processor topology
   F_npx  = 0
   F_npy  = 0
   nblocx = 0
   nblocy = 0
   is_ok = .true.

   istat = wb_get(section_S//'/npx',F_npx)
   istat = min(istat,wb_get(section_S//'/npy',F_npy))
   READ_NML: if (RMN_IS_OK(istat)) then
      call msg(MSG_INFO,'(drv_ptopo_p0) Got topology from Whiteboard')
   else
      istat = clib_getenv('UM_EXEC_CONFIG_BASENAME',tmp_S)
      if (.not.RMN_IS_OK(istat)) &
           istat = wb_get('UM_EXEC_CONFIG_BASENAME',tmp_S)
      if (.not.RMN_IS_OK(istat)) tmp_S = 'drv_settings'
      istat = clib_getenv('UM_EXEC_CONFIG_DIR',tmp2_S)
      if (.not.RMN_IS_OK(istat)) &
           istat = wb_get('path/config_dir0',tmp2_S)
      if (.not.RMN_IS_OK(istat)) tmp2_S = '.'
      tmp_S = trim(adjustl(tmp2_S))//'/'//trim(adjustl(tmp_S))
      user_S = trim(tmp_S)//'.cfg'
      dict_S = trim(tmp_S)//'.dict'

      istat = clib_isfile(dict_S)
      if (RMN_IS_OK(istat)) then
         call msg(MSG_INFO,'(drv_ptopo_p0) Reading topology from: '//trim(user_S)//'::'//trim(section_S))
         istat = wb_read(section_S//'/',dict_S,section_S,WB_DICT_FILE)
      endif
      READ_DICT: if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(drv_ptopo_p0) Unable to read: '//trim(dict_S)//'::'//trim(section_S))
         is_ok = .false.
      else
         istat = clib_isfile(user_S)
         if (RMN_IS_OK(istat)) &
              istat = wb_read(section_S//'/',user_S,section_S,WB_CONF_FILE)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING,'(drv_ptopo_p0) Using default npx,npy - Unable to read: '//trim(user_S)//'::'//trim(section_S))
         endif

         istat = wb_get(section_S//'/npx',F_npx)
         istat = min(istat,wb_get(section_S//'/npy',F_npy))
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,'(drv_ptopo_p0) Unable to get processor topology from: '//trim(user_S)//'::'//trim(section_S))
            is_ok = .false.
         endif
      endif READ_DICT
   endif READ_NML

   if (any((/F_npx,F_npy/) < 1)) is_ok = .false.

   if (.not.is_ok) then
      call msg(MSG_CRITICAL,'(drv_ptopo_p0) Topology not defined or not valid\n==== ABORT ====')
      F_npx = 0
      F_npy = 0
      return
   endif

   istat = wb_get(section_S//'/nblocx',nblocx)
   if (nblocx < 1 .or. .not.RMN_IS_OK(istat)) nblocx = F_npx
   nblocx = min(max(1,nblocx),F_npx)
   istat = wb_put(section_S//'/nblocx',nblocx,WB_REWRITE_MANY)

   istat = wb_get(section_S//'/nblocy',nblocy)
   if (nblocy < 1 .or. .not.RMN_IS_OK(istat)) nblocy = F_npy
   nblocy = min(max(1,nblocy),F_npy)
   istat = wb_put(section_S//'/nblocy',nblocy,WB_REWRITE_MANY)

   istat = wb_get(section_S//'/ninblocx',ninblocx)
   if (ninblocx < 1 .or. .not.RMN_IS_OK(istat)) ninblocx = nblocx
   ninblocx = min(max(1,ninblocx),F_npx)
   istat = wb_put(section_S//'/ninblocx',ninblocx,WB_REWRITE_MANY)

   istat = wb_get(section_S//'/ninblocy',ninblocy)
   if (ninblocy < 1 .or. .not.RMN_IS_OK(istat)) ninblocy = nblocy
   ninblocy = min(max(1,ninblocy),F_npy)
   istat = wb_put(section_S//'/ninblocy',ninblocy,WB_REWRITE_MANY)

   call config_print(section_S)

   write(tmp_S,'(a,2i4,a,2i4,a,2i4)') &
        "(drv_ptopo_p0) Topology defined: npex,y=",F_npx,F_npy,'; nblocx,y=:', &
        nblocx, nblocy,'; ninblocx,y:',ninblocx, ninblocy
   call msg(MSG_INFO,tmp_S)
   !---------------------------------------------------------------------
   return
end subroutine drv_ptopo_p0


