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
module inputio_files_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_toupper, clib_isdir, clib_isfile, clib_realpath
   use fstmpio_mod
   use ptopo_utils, only: PTOPO_IODIST !#TODO: , PTOPO_BLOC
   implicit none
   private
   !@objective
   !@author Stephane Chamberland,2013-03
   !@description
   ! Public functions
   public :: inputio_files_new, inputio_files_set_name, inputio_files_set_basedir, &
        inputio_files_open, inputio_files_close, inputio_files_get_idx, &
        inputio_files_get_type, inputio_files_set_iotype
   ! Public constants
   integer, parameter, public :: INPUT_FILES_NMAX = 16
   integer, parameter, public :: INPUT_FILES_ANAL = 1
   integer, parameter, public :: INPUT_FILES_CLIM = 2
   integer, parameter, public :: INPUT_FILES_GEOP = 3
   integer, parameter, public :: INPUT_FILES_NTYPES = 3
   ! Public Types
   public :: INPUTIO_FILES_T
   !@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   integer, parameter :: NMAX_NLIST = 16
   integer, parameter :: FILE_AVAIL = -1
   integer, parameter :: FILE_RESERV = 0
   integer, parameter :: DEFAULT_LIST_ID = 0

   type :: INFILE_T
      integer :: unit, ftype
      logical :: isdir_L
      character(len=8) :: tag_S
      character(len=1024) :: name_S
   end type INFILE_T

   type :: INPUTIO_FILES_T
      logical :: init_L = .false.
      character(len=1024) :: basedir_S
      integer :: iotype = PTOPO_IODIST
      integer :: nfiles = 0
      type(INFILE_T) :: files(INPUT_FILES_NMAX)
   end type INPUTIO_FILES_T

contains

   !/@*
   function inputio_files_new(F_inputobj, F_basedir_S, F_iotype) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INPUTIO_FILES_T), intent(inout) :: F_inputobj
      character(len=*), intent(in), optional :: F_basedir_S
      integer, intent(in), optional :: F_iotype
      !@return
      integer :: F_istat
      !*@/
      !----------------------------------------------------------------------
      F_istat = RMN_OK
      F_inputobj%init_L = .false.
      call priv_init(F_inputobj)
      if (present(F_basedir_S)) then
         F_istat = inputio_files_set_basedir(F_inputobj, F_basedir_S)
      else
         F_istat = inputio_files_set_basedir(F_inputobj, '.')
      endif
      if (present(F_iotype)) &
           F_inputobj%iotype = F_iotype !#TODO: check validity
      !----------------------------------------------------------------------
   end function inputio_files_new


   !/@*
   function inputio_files_set_name(F_inputobj, F_key_S, F_filename_S, &
        F_isdir_L, F_type, F_check_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INPUTIO_FILES_T), intent(inout) :: F_inputobj
      character(len=*), intent(in) :: F_key_S, F_filename_S
      logical, optional, intent(in) :: F_isdir_L
      integer, optional, intent(in) :: F_type
      logical, optional, intent(in) :: F_check_L
      !@return
      integer :: F_istat
      !*@/
      character(len=8) :: key_S
      character(len=64) :: msg_S
      integer :: istat, ifile, ftype
      logical :: check_L
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio_files) set_name [BEGIN]')
      call priv_init(F_inputobj)
      F_istat = RMN_ERR
      ftype = INPUT_FILES_ANAL
      check_L = .true.
      if (present(F_type)) ftype = min(max(1, F_type), INPUT_FILES_NTYPES)
      if (present(F_check_L)) check_L = F_check_L

      ifile = inputio_files_get_idx(F_inputobj, F_key_S)
      if (ifile < 1) then
!!$         call msg(MSG_WARNING, &
!!$              '(inputio_files) set_name, unknown file key name: '//trim(F_key_S))
!!$         return
         if (F_inputobj%nfiles >= INPUT_FILES_NMAX) then
            call msg(MSG_WARNING, &
                 '(inputio_files) set_name, too many file key name, canot set: ' &
                 //trim(F_key_S))
            return
         endif
         key_S = F_key_S
         istat = clib_toupper(key_S)
         ifile = F_inputobj%nfiles + 1
         F_inputobj%files(ifile)%tag_S = key_S
         F_inputobj%files(ifile)%name_S = F_filename_S
         F_inputobj%files(ifile)%isdir_L = .false.
      endif
      if (present(F_isdir_L)) F_inputobj%files(ifile)%isdir_L = F_isdir_L

      if (check_L) then
         if (F_inputobj%files(ifile)%isdir_L) then
            F_istat = clib_isdir(trim(F_inputobj%basedir_S)//'/'//trim(F_filename_S))
            if (.not.RMN_IS_OK(F_istat)) then
               F_istat = clib_isfile(trim(F_inputobj%basedir_S)//'/'//trim(F_filename_S))
            endif
         else
            F_istat = clib_isfile(trim(F_inputobj%basedir_S)//'/'//trim(F_filename_S))
         endif
      else
         F_istat = RMN_OK
      endif
      if (RMN_IS_OK(F_istat)) then
         F_istat = ifile
         F_inputobj%files(ifile)%name_S = trim(F_filename_S)
         F_inputobj%files(ifile)%ftype = ftype
         call msg(MSG_INFO, '(inputio_files) set_name: '//trim(F_key_S) &
              //' => '//trim(F_filename_S))
      else
         call msg(MSG_WARNING, '(inputio_files) set_name, not such dir/file: ' &
              //trim(F_inputobj%basedir_S)//'/'//trim(F_filename_S))
         return
      endif

      F_inputobj%nfiles = max(F_inputobj%nfiles, ifile)
      write(msg_S, '(i0,1x,i0)') ifile, F_inputobj%nfiles
      call msg(MSG_DEBUG, '(inputio_files) set_name [END] '//trim(F_key_S)//':'//trim(msg_S))
      !----------------------------------------------------------------------
      return
   end function inputio_files_set_name


   !/@*
   function inputio_files_set_basedir(F_inputobj, F_basedir_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(INPUTIO_FILES_T), intent(inout) :: F_inputobj
      character(len=*), intent(in) :: F_basedir_S
      !@return
      integer :: F_istat
      !*@/
      character(len=1024) :: basedir_S
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio_files) set_basedir [BEGIN]')
      call priv_init(F_inputobj)
      F_istat = RMN_ERR
      F_istat = clib_isdir(trim(F_basedir_S))
      F_istat = min(clib_realpath(trim(F_basedir_S), basedir_S), F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR, &
              '(inputio_files) set_basedir, not such dir: '//trim(F_basedir_S))
      else
         call msg(MSG_INFO, &
              '(inputio_files) set_basedir: '//trim(F_basedir_S))
         F_inputobj%basedir_S = trim(basedir_S)
         F_istat = RMN_OK
      endif
      call msg(MSG_DEBUG, '(inputio_files) set_basedir [END]')
      !----------------------------------------------------------------------
      return
   end function inputio_files_set_basedir


   !/@*
   function inputio_files_get_idx(F_inputobj, F_name_S) result(F_fileidx)
      implicit none
      !@objective
      !@arguments
      type(INPUTIO_FILES_T), intent(in) :: F_inputobj
      character(len=*), intent(in) :: F_name_S
      !@return
      integer :: F_fileidx
      !*@/
      character(len=8) :: name_S
      integer :: istat, ifile
      !------------------------------------------------------------------
!!$      call priv_init(F_inputobj)
      F_fileidx = RMN_ERR
      name_S = F_name_S
      istat = clib_toupper(name_S)
      do ifile=1, F_inputobj%nfiles
         if (name_S == F_inputobj%files(ifile)%tag_S) then
            F_fileidx = ifile
            exit
         endif
      enddo
      !# Truncated names search kept for backward compatibility
      if (.not.RMN_IS_OK(F_fileidx)) then
         do ifile=1, F_inputobj%nfiles
            if (name_S(1:4) == F_inputobj%files(ifile)%tag_S(1:4)) then
               F_fileidx = ifile
               exit
            endif
         enddo
      endif
      !------------------------------------------------------------------
      return
   end function inputio_files_get_idx


   !/@*
   function inputio_files_get_type(F_inputobj, F_file_idx) result(F_type)
      implicit none
      !@objective
      !@arguments
      type(INPUTIO_FILES_T), intent(in) :: F_inputobj
      integer, intent(in) :: F_file_idx
      !@return
      integer :: F_type
      !*@/
      integer :: fileidx
      !------------------------------------------------------------------
!!$      call priv_init(F_inputobj)
      F_type = RMN_ERR
      fileidx = min(max(1, F_file_idx), F_inputobj%nfiles)
      F_type = F_inputobj%files(fileidx)%ftype
      !------------------------------------------------------------------
      return
   end function inputio_files_get_type


  !/@*
   subroutine inputio_files_set_iotype(F_inputobj, F_iotype)
      implicit none
      type(INPUTIO_FILES_T), intent(inout) :: F_inputobj
      integer, intent(in), optional :: F_iotype
      !*@/
      !----------------------------------------------------------------------
      if (present(F_iotype)) then
         F_inputobj%iotype = F_iotype !#TODO: check validity
      endif
      call fstmpio_set_iotype(F_inputobj%iotype)
      !----------------------------------------------------------------------
      return
   end subroutine inputio_files_set_iotype


   !/@*
   function inputio_files_open(F_inputobj, F_file_idx) result(F_funit)
      implicit none
      type(INPUTIO_FILES_T), intent(inout) :: F_inputobj
      integer, intent(in) :: F_file_idx
      integer :: F_funit
      !*@/
      character(len=1024) :: filename_S, filename0_S, unit_S
      integer :: istat, fileidx
      !------------------------------------------------------------------
      call priv_init(F_inputobj)
      call inputio_files_set_iotype(F_inputobj)
      call msg(MSG_DEBUG, '(inputio_files) open [BEGIN]')
      F_funit = RMN_ERR
      fileidx = min(max(1, F_file_idx), F_inputobj%nfiles)
      if (F_inputobj%files(fileidx)%unit < 1) then
         if (F_inputobj%basedir_S == ' ') F_inputobj%basedir_S = '.'
         filename0_S = trim(F_inputobj%basedir_S)//'/'// &
              F_inputobj%files(fileidx)%name_S
         istat = clib_realpath(filename0_S, filename_S)
         if (.not.RMN_IS_OK(istat)) filename_S = filename0_S
         F_inputobj%files(fileidx)%unit = fstmpio_open(trim(filename_S), &
              FST_READONLY, F_inputobj%files(fileidx)%isdir_L)
         if (F_inputobj%files(fileidx)%unit > 1) then
            write(unit_S, '(i0)') F_inputobj%files(fileidx)%unit
            call msg(MSG_INFOPLUS, '(inputio_files) opened file '//trim(unit_S)//': ' &
                 //trim(F_inputobj%basedir_S)//'/'//F_inputobj%files(fileidx)%name_S)
         endif
      endif
      if (F_inputobj%files(fileidx)%unit < 1) then
         F_inputobj%files(fileidx)%unit = FILE_AVAIL
         call msg(MSG_ERROR, '(inputio_files) Unable to open: ' &
              //trim(F_inputobj%basedir_S)//'/'//F_inputobj%files(fileidx)%name_S)
      else
         F_funit = F_inputobj%files(fileidx)%unit
      endif
      write(unit_S, '(i0)') F_funit
      call msg(MSG_DEBUG, '(inputio_files) open [END]'//trim(unit_S))
      !------------------------------------------------------------------
      return
   end function inputio_files_open


   !/@*
   function inputio_files_close(F_inputobj) result(F_istat)
      implicit none
      !@objective
      !@return
      type(INPUTIO_FILES_T), intent(inout) :: F_inputobj
      integer :: F_istat
      !*@/
      character(len=1024) :: unit_S
      integer :: ii, istat
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG, '(inputio_files) close [BEGIN]')
      call priv_init(F_inputobj)
      call inputio_files_set_iotype(F_inputobj)
      F_istat = RMN_OK
      do ii = 1, F_inputobj%nfiles
         if (F_inputobj%files(ii)%unit > 0) then
            write(unit_S, '(i0)') F_inputobj%files(ii)%unit
            call msg(MSG_INFOPLUS, '(inputio_files) closing file '//trim(unit_S)//': ' &
                 //trim(F_inputobj%basedir_S)//'/'//F_inputobj%files(ii)%name_S)
            istat = fstmpio_close(F_inputobj%files(ii)%unit)
            if (RMN_IS_OK(istat)) F_inputobj%files(ii)%unit = FILE_AVAIL
            F_istat = min(F_istat, istat)
         endif
      enddo
      call msg(MSG_DEBUG, '(inputio_files) close [END]')
      !----------------------------------------------------------------------
      return
   end function inputio_files_close


   !==== Private Functions =================================================


   !/@*
   subroutine priv_init(F_inputobj)
      implicit none
      type(INPUTIO_FILES_T), intent(inout) :: F_inputobj
      !*@/
      logical, parameter :: NO_CHECK_L = .false.
      integer, parameter :: NBASE_FILES = 7
      character(len=8), parameter :: BASEFILE_TAG_S(NBASE_FILES) = &
           (/'GEOPHY  ', 'CLIMATO ', 'ANALYSIS', &
             'INREP   ', 'MANAL   ', 'MINREP  ', &
             'MINPUT  '/)
      character(len=16), parameter :: BASEFILE_NAME_S(NBASE_FILES) = &
           (/'GEOPHY          ', 'CLIMATO         ', 'ANALYSIS        ', &
             'INREP           ', 'MODEL_ANALYSIS  ', 'MODEL_INREP     ', &
             'MODEL_INPUT     '/)
      logical, parameter :: BASEFILE_ISDIR_L(NBASE_FILES) = &
           (/.true., .true., .true., &
             .true., .true., .true., &
             .true./)
      integer, parameter :: BASEFILE_TYPE(NBASE_FILES) = &
           (/INPUT_FILES_GEOP, INPUT_FILES_CLIM, INPUT_FILES_ANAL, &
             INPUT_FILES_ANAL, INPUT_FILES_ANAL, INPUT_FILES_ANAL, &
             INPUT_FILES_ANAL/)
      integer :: ifile, istat
      !------------------------------------------------------------------
      if (F_inputobj%init_L) return
      call msg(MSG_DEBUG, '(inputio_files) init [BEGIN]')
      F_inputobj%init_L = .true.
      F_inputobj%iotype = PTOPO_IODIST
      F_inputobj%nfiles = 0
      F_inputobj%basedir_S = '.'
      do ifile=1, INPUT_FILES_NMAX
         F_inputobj%files(ifile)%tag_S = ' '
         F_inputobj%files(ifile)%name_S = ' '
         F_inputobj%files(ifile)%unit = FILE_AVAIL
         F_inputobj%files(ifile)%ftype = INPUT_FILES_ANAL
         F_inputobj%files(ifile)%isdir_L = .false.
      enddo
      F_inputobj%nfiles = 0
      do ifile=1, NBASE_FILES
         istat = inputio_files_set_name(F_inputobj, BASEFILE_TAG_S(ifile), &
              BASEFILE_NAME_S(ifile), BASEFILE_ISDIR_L(ifile), &
              BASEFILE_TYPE(ifile), NO_CHECK_L)
      enddo
      call msg(MSG_DEBUG, '(inputio_files) init [END]')
      !------------------------------------------------------------------
      return
   end subroutine priv_init

end module inputio_files_mod

