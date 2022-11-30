!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!--------------------------------------------------------------------------

!/@*
module config_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use clib_itf_mod, only: clib_realpath, clib_realpath, clib_isfile, clib_getcwd
   use wb_itf_mod
   implicit none
   private
   !@objective Read options config [dictionary + users config] for MPI jobs
   !@author
   !  Michel Desgagne, Feb 2008
   !  Ron McTaggart-Cowan, Feb 2008
   !  Stephane Chamberland, Feb 2008
   !@revisions
   !  2012-02, Stephane Chamberland: config_print, config_path, update for Phyoffline
   !@public_parameters
   integer, parameter,public :: CFG_MAX_FILE_SIZE = 100000
   !@public_functions
   public :: config_init, config_path, config_read, config_cp2localdir, config_print
   !@public_vars
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   character(len=RMN_PATH_LEN),save :: m_basedir_S = '.'
   integer,save :: m_me_pe     = -1
   logical,save :: m_is_init_L = .false.

contains

   !/@*
   subroutine config_init(F_basedir_S)
      implicit none
      !@objective Initialize config module
      !@arguments
      character(len=*),intent(in),optional :: F_basedir_S !-
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
   !/@*
      integer :: istat
      !---------------------------------------------------------------------
      m_is_init_L = .true.
      if (present(F_basedir_S)) then
         istat = clib_realpath(F_basedir_S,m_basedir_S) !- make sure F_basedir_S is an absolute path
      else
         istat = clib_realpath('.',m_basedir_S)
      endif
      call msg(MSG_DEBUG,'(config_init) basedir set to: '//trim(m_basedir_S))
      call rpn_comm_rank(RPN_COMM_GRID,m_me_pe,istat)
      !---------------------------------------------------------------------
      return
   end subroutine config_init


   !/@*
   function config_path(F_path_S,F_cfg_S,F_ext_S,F_basedir_S,F_force_cp_L) result(F_istat)
      implicit none
      !@objective Get config path/name.ext, cp file to local dir if needed
      !@arguments
      character(len=*),intent(out) :: F_path_S !- config file basename [dict/user], no path, no .ext
      character(len=*),intent(in) :: F_cfg_S !- config file basename [dict/user], no path, no .ext
      character(len=*),intent(in),optional :: F_ext_S !- config filename ext [dict/cfg]
      character(len=*),intent(in),optional :: F_basedir_S  !- dir of orig. file
      logical,intent(in),optional :: F_force_cp_L !- copy even if local copy exists
      !@returns
      integer :: F_istat !- RMN_OK on success, RMN_ERR otherwise
      !@author
      !  Stephane Chamberland, 2012-02
   !/@*
      character(len=RMN_PATH_LEN) :: basedir_S,path_S,ext_S
      integer :: istat
      logical :: force_cp_L
      !---------------------------------------------------------------------
      basedir_S = m_basedir_S
      if (present(F_basedir_S)) istat = clib_realpath(F_basedir_S,basedir_S)

      if (.not.m_is_init_L) then
         call config_init(basedir_S)
      endif

      force_cp_L = .false.
      if (present(F_force_cp_L)) force_cp_L = F_force_cp_L

      ext_S = 'cfg'
      if (present(F_ext_S)) ext_S = F_ext_S
      path_S = trim(adjustl(F_cfg_S))//'.'//trim(ext_S)

      F_istat = config_cp2localdir(path_S,basedir_S,force_cp_L)
      if (RMN_IS_OK(F_istat)) then
         F_istat = clib_realpath('./'//path_S,F_path_S)
      endif
      !---------------------------------------------------------------------
      return
   end function config_path


   !/@*
   function config_read(F_cfg_S,F_sec_S,F_pkg_S,F_nocfg_istat,F_basedir_S,F_force_cp_L) result(F_istat)
      implicit none
      !@objective Read dict/config into global whiteboard
      !@arguments
      character(len=*),intent(in) :: F_cfg_S !- config file basename [dict/user], no path, no .ext
      character(len=*),intent(in) :: F_sec_S !- config file section name [as for namelist], relate to dictionary file 
      character(len=*),intent(in),optional :: F_pkg_S !- wb name prefix
      integer,intent(in),optional :: F_nocfg_istat    !- return F_nocfg_istat if user cfg file does not exists
      character(len=*),intent(in),optional :: F_basedir_S  !- dir of orig. file
      logical,intent(in),optional :: F_force_cp_L !- copy even if local copy exists
      !@returns
      integer :: F_istat !- RMN_OK on success, RMN_ERR otherwise
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
      !@revisions
      !  2012-02, Stephane Chamberland: F_nocfg_istat,F_basedir_S,F_force_cp
   !/@*
      character(len=RMN_PATH_LEN) :: pkg_S,dict_S,user_S,basedir_S
      integer :: istat,nocfg_istat,ii
      logical :: force_cp_L
      !---------------------------------------------------------------------
      basedir_S = m_basedir_S
      if (present(F_basedir_S)) istat = clib_realpath(F_basedir_S,basedir_S)

      if (.not.m_is_init_L) then
         call config_init(basedir_S)
      endif

      force_cp_L = .false.
      if (present(F_force_cp_L)) force_cp_L = F_force_cp_L

      pkg_S = adjustl(F_sec_S)
      if (present(F_pkg_S)) pkg_S = adjustl(F_pkg_S)
      if (pkg_S /= '') then
         ii = len_trim(pkg_S)
         if (pkg_S(ii:ii) /= '/') pkg_S = trim(pkg_S)//'/'
      endif
      nocfg_istat = RMN_OK
      if (present(F_nocfg_istat)) nocfg_istat = F_nocfg_istat
      dict_S = trim(adjustl(F_cfg_S))//'.dict'
      user_S = trim(adjustl(F_cfg_S))//'.cfg'

      F_istat = config_cp2localdir(dict_S,basedir_S,force_cp_L)
      if (RMN_IS_OK(F_istat)) &
           F_istat = wb_read(pkg_S,dict_S,F_sec_S,WB_DICT_FILE)
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(config_read) Problem reading: '//trim(dict_S)//'::'//trim(F_sec_S))
         return
      endif
      call msg(MSG_INFO,'(config_read) Reading: '//trim(user_S)//'::'//trim(F_sec_S))
 
      F_istat = config_cp2localdir(user_S,basedir_S,force_cp_L)
      if (.not.RMN_IS_OK(F_istat)) then
         F_istat = nocfg_istat
      else
         F_istat = wb_read(pkg_S,user_S,F_sec_S,WB_CONF_FILE)
      endif
      if (RMN_IS_OK(F_istat)) call config_print(pkg_S)
      !---------------------------------------------------------------------
      return
   end function config_read


   !/@*
   function config_cp2localdir(F_filename_S, F_basedir_S, F_force_cp_L) result(F_istat)
      implicit none
      !@objective Copy a file from the experiement dir to the local process dir
      !@arguments
      character(len=*),intent(in)          :: F_filename_S !- fileName to copy
      character(len=*),intent(in),optional :: F_basedir_S  !- dir of orig. file
      logical,intent(in),optional :: F_force_cp_L !- copy even if local copy exists
      !@returns
      integer :: F_istat !- Error Status
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
      !@revisions
      !  2012-02, Stephane Chamberland: F_basedir_S now override m_basedir_S; F_force_cp
   !*@/
      character(len=RMN_PATH_LEN) :: basedir_S, cwd_S,fullpath0_S,fullpath_S
      integer :: bufnml(CFG_MAX_FILE_SIZE), istat
      logical :: cp_flag_L,force_cp_L
      !---------------------------------------------------------------------
      F_istat = RMN_OK
      basedir_S = m_basedir_S
      if (present(F_basedir_S)) istat = clib_realpath(F_basedir_S,basedir_S)

      if (.not.m_is_init_L) then
         call config_init(basedir_S)
      endif

      force_cp_L = .false.
      if (present(F_force_cp_L)) force_cp_L = F_force_cp_L
      istat = clib_getcwd(cwd_S)

      F_istat = RMN_ERR
      if (.not.force_cp_L) then
         fullpath0_S = trim(cwd_S)//'/'//trim(F_filename_S)
         F_istat = clib_realpath(fullpath0_S,fullpath_S)
         F_istat = min(clib_isfile(fullpath_S),F_istat)
      endif
      call collect_error(F_istat)
      if (RMN_IS_OK(F_istat)) return

      fullpath0_S = trim(basedir_S)//'/'//trim(F_filename_S)

      F_istat = clib_realpath(fullpath0_S,fullpath_S)
      F_istat = min(clib_isfile(fullpath_S),F_istat)
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(config_cp2localdir) No such File: '//trim(fullpath0_S))
         return
      endif

      !TODO-later: implement this as soon a clib_filemtime exists [int stat (const char *filename, struct stat *buf)]
      !alt-implementation 1: store [fullpath,pwd] ;
      !   cp_flag_L = all(list.ne.'fullpath,pwd') [?Safer since not depedent
      !   on risk of user changing fullpath runtime]
      !alt-implementation 2: store start date [when doing the init above];
      !   cp_flag_L = filemtime(F_filename) < start_time [?Safer since not
      !   depedent on risk of user changing fullpath runtime]
      cp_flag_L = .true.
!!$      cp_flag_L = .false.
!!$      if (clib_isfile(F_filename_S).ne.CLIB_OK) then
!!$         cp_flag_L = .true.
!!$      else
!!$         mtime_src = clib_filemtime(F_filename_S)
!!$         mtime_dst = clib_filemtime(fullpath_S)
!!$         if (mtime_src > mtime_dst) cp_flag_L = .true.
!!$      endif

      if (cp_flag_L) then
         if (m_me_pe == RPN_COMM_MASTER) &
              call array_from_file(bufnml,size(bufnml),fullpath_S)
         call rpn_comm_bcast(bufnml,size(bufnml),RPN_COMM_INTEGER,RPN_COMM_MASTER,RPN_COMM_MULTIGRID,istat)
         if (basedir_S == cwd_S) return !src == dst, nothing to do (should be done after rpn_comm_bcast to avoid mpi hang
         if (RMN_IS_OK(istat)) then
            call array_to_file(bufnml,size(bufnml),F_filename_S)
         else
            call msg(MSG_ERROR,'(config_cp2localdir) rpn_comm_bcast error')
            F_istat = RMN_ERR
         endif
         fullpath_S = trim(cwd_S)//'/'//trim(F_filename_S)
         F_istat = clib_isfile(fullpath_S)
         if (.not.RMN_IS_OK(F_istat)) then
            call msg(MSG_ERROR,'(config_cp2localdir) Copy not succesfull - No such Config File: '//trim(fullpath_S))
         else
            call msg(MSG_INFOPLUS,'(config_cp2localdir) '//trim(fullpath0_S)//' ==> '//trim(fullpath_S))
         endif
      endif
      !---------------------------------------------------------------------
      return
   end function config_cp2localdir


   !/@*
   subroutine config_print(F_pkg_S)
      implicit none
      !@objective Print wb sub set values
      !@arguments
      character(len=*),intent(in) :: F_pkg_S
      !@author
      !  Stephane Chamberland, 2012-02
   !/@*
      integer,parameter :: MAX_KEYS = 1024
      integer,parameter :: MAX_VALUES = 2048
      integer,parameter :: MAX_PRINT = 64
      integer :: istat,nkeys,ikey,element_type,typelen,array_size,options, &
           nvals,ii,msgLevelMin,msgUnit
      character(len=WB_MAXNAMELENGTH) :: keys(MAX_KEYS)
      character(len=WB_MAXSTRINGLENGTH) :: cval,caval(MAX_VALUES)
      character(len=MSG_MAXLEN) :: msg_S
      integer :: ival,iaval(MAX_VALUES)
      integer(INT64) :: i8val,i8aval(MAX_VALUES)
      real :: rval,raval(MAX_VALUES)
      real(REAL64) :: r8val,r8aval(MAX_VALUES)
      logical :: bval,baval(MAX_VALUES),canWrite_L
      !---------------------------------------------------------------------
      call msg_getInfo(canWrite_L,msgLevelMin,msgUnit,msg_S)
      if (MSG_INFO < msgLevelMin .or. .not.canWrite_L) return

      if (F_pkg_S /= '') then
         ii = len_trim(F_pkg_S)
         if (F_pkg_S(ii:ii) == '/') ii = ii-1
         write(msgUnit, '(a)') '@'//trim(F_pkg_S(1:ii))
      else
         write(msgUnit,'(a)') '@'
      endif
      istat = wb_keys(keys,nkeys,F_pkg_S)
      do ikey = 1, nkeys
         istat = wb_get_meta(keys(ikey),element_type,typelen,array_size,options)
         select case(element_type)
         case(WB_FORTRAN_REAL)
            if (typelen == 4) then
               if (array_size == 0) then
                  istat = wb_get(keys(ikey),rval)
                  write(msgUnit,'(a,g0)') trim(keys(ikey))//' = ',rval
               else
                  istat = wb_get(keys(ikey),raval,nvals)
                  if (nvals <= MAX_PRINT) then
                     write(msgUnit,*) trim(keys(ikey))//' = ',raval(1:nvals)
                  else
                     write(msgUnit,*) trim(keys(ikey))//' = ', &
                          raval(1:MAX_PRINT/2),', ..., ', &
                          raval(nvals-MAX_PRINT/2:nvals)
                  endif
               endif
            else
               if (array_size == 0) then
                  istat = wb_get(keys(ikey),r8val)
                  write(msgUnit,*) trim(keys(ikey))//' = ',r8val
               else
                  istat = wb_get(keys(ikey),r8aval,nvals)
                  if (nvals <= MAX_PRINT) then
                     write(msgUnit,*) trim(keys(ikey))//' = ', &
                          r8aval(1:nvals)
                  else
                     write(msgUnit,*) trim(keys(ikey))//' = ', &
                          r8aval(1:MAX_PRINT/2),', ..., ',r8aval(nvals-MAX_PRINT/2:nvals)
                  endif
               endif
            endif
         case(WB_FORTRAN_INT)
            if (typelen == 4) then
               if (array_size == 0) then
                  istat = wb_get(keys(ikey),ival)
                  write(msgUnit,*) trim(keys(ikey))//' = ',ival
               else
                  istat = wb_get(keys(ikey),iaval,nvals)
                  if (nvals <= MAX_PRINT) then
                     write(msgUnit,*) trim(keys(ikey))//' = ',iaval(1:nvals)
                  else
                     write(msgUnit,*) trim(keys(ikey))//' = ', &
                          iaval(1:MAX_PRINT/2),', ..., ', &
                          iaval(nvals-MAX_PRINT/2:nvals)
                  endif
               endif
            else
               if (array_size == 0) then
                  istat = wb_get(keys(ikey),i8val)
                  write(msgUnit,*) trim(keys(ikey))//' = ',i8val
               else
                  istat = wb_get(keys(ikey),i8aval,nvals)
                  if (nvals <= MAX_PRINT) then
                     write(msgUnit,*) trim(keys(ikey))//' = ',i8aval(1:nvals)
                  else
                     write(msgUnit,*) trim(keys(ikey))//' = ', &
                          i8aval(1:MAX_PRINT/2),', ..., ', &
                          i8aval(nvals-MAX_PRINT/2:nvals)
                  endif
               endif
            endif
         case(WB_FORTRAN_CHAR)
            if (array_size == 0) then
               istat = wb_get(keys(ikey),cval)
               write(msgUnit,*) trim(keys(ikey))//' = ',trim(cval)
            else
               istat = wb_get(keys(ikey),caval,nvals)
               if (nvals <= MAX_PRINT) then
                  write(msgUnit,*) trim(keys(ikey))//' = ', &
                       ('"'//trim(caval(ii))//'", ',ii=1,nvals)
               else
                  write(msgUnit,*) trim(keys(ikey))//' = ', &
                       ('"'//trim(caval(ii))//'", ',ii=1,MAX_PRINT/2),', ..., ',&
                       ('"'//trim(caval(ii))//'", ',ii=nvals-MAX_PRINT/2,nvals)
               endif
            endif
         case(WB_FORTRAN_BOOL)
            if (array_size == 0) then
               istat = wb_get(keys(ikey),bval)
               write(msgUnit,*) trim(keys(ikey))//' = ',bval
            else
               write(msgUnit,*) trim(keys(ikey))//' = ',baval(1:nvals)
               if (nvals <= MAX_PRINT) then
                  write(msgUnit,*) trim(keys(ikey))//' = ',baval(1:nvals)
               else
                  write(msgUnit,*) trim(keys(ikey))//' = ',baval(1:MAX_PRINT/2),', ..., ',baval(nvals-MAX_PRINT/2:nvals)
               endif
            endif
!!$               case defautl
         end select
      enddo
      write(msgUnit,*) '@'
      !---------------------------------------------------------------------
      return
   end subroutine config_print


end module config_mod
