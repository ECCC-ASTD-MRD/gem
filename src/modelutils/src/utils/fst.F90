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
#include <rmn/msg.h>

!/@
module fst_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_isfile, clib_isdir, clib_isreadok
   use mu_fglob_mod
   use fst_read_mod
   use fst_write_mod
   implicit none
   private
   !@objective 
   !@author  Stephane Chamberland, 2011-04
   !@description
   ! Public functions
   public :: fst_open, fst_close, fst_find, fst_read, fst_rdhint, fst_write, &
        fst_get_gridid, fst_getmeta, fst_get_hgridid, fst_get_vgrid, &
        fst_write_opt, fst_checkalloc, fst_checkalloc2
   public :: fst_find_vect, fst_find_0, fst_find_3d_vect, fst_find_3d_0
   ! Public constants
   logical, parameter, public :: FST_READONLY   = .true.
   public :: FST_FIND_LT, FST_FIND_LE, FST_FIND_NEAR, FST_FIND_GE,  &
        FST_FIND_GT, FST_NPAK_DEFAULT, FST_NPAK_FULL32, FST_OPT_OUTGRID_TYPE, &
        FST_ZGRID, FST_DIEZEGRID, FST_FIND_DIAG_T, FST_FIND_DIAG_M
!@/

#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

contains

   !/@
   function fst_open(F_filename_S, F_readonly_L, F_dir_ok_L) result(F_fileid)
      implicit none
      !@objective Open rpn std file
      !@arguments
      character(len=*), intent(in) :: F_filename_S
      logical, intent(in),optional :: F_readonly_L
      logical, intent(in),optional :: F_dir_ok_L
      !@author
      !@return
      integer :: F_fileid
      !@/
      integer,parameter :: NMAXFILES = 2048, FSTD89 = 1, FSTD98 = 33
 
      logical :: readonly_L,dir_ok_L,is_file_L,is_dir_L,can_read_L
      integer :: istat,nfiles,nfiles2,ifile,unitlist(NMAXFILES),ftype
      character(len=32) :: type_S
      character(len=512) :: msg_S
      character(len=RMN_PATH_LEN) :: filelist_S(NMAXFILES)

      integer, external :: wkoffit
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fst) open [BEGIN]')
      F_fileid = 0
      readonly_L = .not.FST_READONLY
      dir_ok_L = .false.
      if (present(F_readonly_L)) readonly_L = F_readonly_L
      if (present(F_dir_ok_L)) dir_ok_L = F_dir_ok_L

      is_file_L  = RMN_IS_OK(clib_isfile(trim(F_filename_S)))
      is_dir_L   = RMN_IS_OK(clib_isdir(trim(F_filename_S)))
      can_read_L = RMN_IS_OK(clib_isreadok(trim(F_filename_S)))
      if (is_dir_L) readonly_L = .true.

      if ((readonly_L.or.is_file_L.or.is_dir_L).and..not.can_read_L) then
         F_fileid = RMN_ERR
         call msg(MSG_WARNING,'(fst_open) not readable: '//trim(F_filename_S))
         return
      endif
 
      if (is_dir_L.and..not.dir_ok_L) then
         F_fileid = RMN_ERR
         call msg(MSG_WARNING,'(fst_open) not a file: '//trim(F_filename_S))
         return
      endif

      type_S = 'RND+R/W'
      if (readonly_L) type_S = 'RND+OLD+R/O'

!!$      if (.not.F_readonly_L) then
!!$         istat = min(clib_iswriteok(trim(F_filename_S)),istat)
!!$      endif

      !#TODO: accept also filename pattern to be resolved by glob
      if (.not.is_dir_L) then
         istat = fnom(F_fileid,trim(F_filename_S),type_S,0)
         if (RMN_IS_OK(istat) .and. F_fileid > 0) istat = fstouv(F_fileid,'RND')
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING,'(fst_open) probleme opening file: '//trim(F_filename_S))
            istat = fstfrm(F_fileid)
            istat = fclos(F_fileid)
            F_fileid = RMN_ERR
            return
         endif
         write(msg_S,'(a,i6,a)') '(fst_open) id=',F_fileid,' name='//trim(F_filename_S)
         call msg(MSG_INFO,msg_S)
      else
         !#TODO: make recurse an option/arg
         nfiles = mu_fglob(filelist_S, trim(F_filename_S), &
                           F_recurse=1, F_filter=MU_FGLOB_FILE_ONLY)
         nfiles2 = 0
         do ifile = 1,nfiles
            is_file_L  = RMN_IS_OK(clib_isfile(trim(filelist_S(ifile))))
            can_read_L = RMN_IS_OK(clib_isreadok(trim(filelist_S(ifile))))
            ftype = wkoffit(trim(filelist_S(ifile)))
            if (.not.(is_file_L.and.can_read_L.and.(ftype==FSTD89.or.ftype==FSTD98))) then
               call msg(MSG_WARNING,'(fst_open) ignoring: '//trim(filelist_S(ifile)))
               cycle
            endif
            F_fileid = 0
            istat = fnom(F_fileid,trim(filelist_S(ifile)),type_S,0)
            if (RMN_IS_OK(istat) .and. F_fileid > 0) istat = fstouv(F_fileid,'RND')
            if (.not.RMN_IS_OK(istat)) then
               call msg(MSG_WARNING,'(fst_open) probleme opening file, ignoring: '//trim(filelist_S(ifile)))
               istat = fstfrm(F_fileid)
               istat = fclos(F_fileid)
               cycle
            endif
            nfiles2 = nfiles2 + 1
            unitlist(nfiles2) = F_fileid
            write(msg_S,'(a,i6,a)') '(fst_open) id=',F_fileid,' name='//trim(filelist_S(ifile))
            call msg(MSG_INFO,msg_S)

         enddo
         if (nfiles2 > 0) then
            F_fileid = unitlist(1)
            if (nfiles2 > 1) then
               istat = fstlnk(unitlist(1:nfiles2),nfiles2)
               if (.not.RMN_IS_OK(istat)) then
                  call msg(MSG_ERROR,'(fst_open) Problem linking openedfiles from dir: '//trim(F_filename_S))
                  !TODO-later: fstfrm/fclos all files
                  F_fileid = RMN_ERR
                  return
               endif
            endif
         else
            call msg(MSG_ERROR,'(fst_open) Was not able to open any file from dir: '//trim(F_filename_S))
            F_fileid = RMN_ERR
            return
         endif
      endif
      call msg(MSG_DEBUG,'(fst) open [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_open


   !/@
   function fst_close(F_fileid) result(F_istat)
      implicit none
      !@objective Close rpn std file
      !@arguments
      integer, intent(in) :: F_fileid
      !@author
      !@return
      integer :: F_istat
      !@/
      character(len=512) :: msg_S
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(fst) close [BEGIN]')
      F_istat = RMN_ERR
      if (F_fileid > 0) then
         F_istat = fstfrm(F_fileid)
         F_istat = min(fclos(F_fileid),F_istat)
         write(msg_S,'(a,i6)') '(fst_close) id=',F_fileid
         call msg(MSG_INFO,msg_S)
      endif
      !TODO-later: anything special to close a fstlnk list of files?
      call msg(MSG_DEBUG,'(fst) end [BEGIN]')
      ! ---------------------------------------------------------------------
      return
   end function fst_close


end module fst_mod
