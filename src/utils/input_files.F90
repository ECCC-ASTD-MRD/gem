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
module input_files_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower, clib_isdir, clib_isfile, clib_realpath
   use mu_jdate_mod
   use fstmpi_mod
   use ezgrid_mod
   use vGrid_Descriptors
   use vgrid_wb
   use vgrid_from_file_mod
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland,2013-03
   !@description
   ! Public functions
   public :: input_files_get_idx, input_files_get_type, input_files_read, &
        input_files_vgrid,input_files_sfcref
   public :: input_files_set_name, input_files_set_basedir, input_files_open, input_files_close
   public :: input_files_set_name0, input_files_set_basedir0, input_files_open0, input_files_close0
   public :: input_files_set_name1, input_files_set_basedir1, input_files_open1, input_files_close1
   ! Public constants
   integer,parameter,public :: INPUT_FILES_NMAX = 16
   integer,parameter,public :: INPUT_FILES_ANAL = 1
   integer,parameter,public :: INPUT_FILES_CLIM = 2
   integer,parameter,public :: INPUT_FILES_GEOP = 3
   integer,parameter,public :: INPUT_FILES_NTYPES = 3
   !@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface input_files_set_name
      module procedure input_files_set_name0
      module procedure input_files_set_name1
   end interface

   interface input_files_set_basedir
      module procedure input_files_set_basedir0
      module procedure input_files_set_basedir1
   end interface

   interface input_files_open
      module procedure input_files_open0
      module procedure input_files_open1
   end interface

   interface input_files_close
      module procedure input_files_close0
      module procedure input_files_close1
   end interface

   interface input_files_vgrid
      module procedure input_files_vgrid_4
      module procedure input_files_vgrid_8
   end interface

   interface input_files_sfcref
      module procedure input_files_sfcref_4
      module procedure input_files_sfcref_8
   end interface

   interface input_files_read
      module procedure input_files_read_4
      module procedure input_files_read_8
   end interface


   integer,parameter :: NMAX_NLIST = 16
   integer,parameter :: FILE_AVAIL = -1
   integer,parameter :: FILE_RESERV = 0
   integer,parameter :: DEFAULT_LIST_ID = 0
   real,parameter :: MB2PA = 100.

   type :: input_file_T
!!$      sequence
      integer :: id,type
      logical :: isdir_L
      character(len=8) :: tag_S
      character(len=1024) :: name_S
!!$      type(vgrid_descriptor),pointer :: vgrid
   end type input_file_T

   type(input_file_T),save :: m_files(0:NMAX_NLIST,INPUT_FILES_NMAX)
   integer,save :: m_nfiles(0:NMAX_NLIST)
   character(len=1024),save :: m_basedir_S(0:NMAX_NLIST)

contains

   !/@*
   function input_files_set_name0(F_key_S,F_filename_S,F_isdir_L,F_type) result(F_istat)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_key_S,F_filename_S
      logical,optional,intent(in) :: F_isdir_L
      integer,optional,intent(in) :: F_type
      !@return
      integer :: F_istat
      !*@/
      integer :: type
      logical :: isdir_L
      !----------------------------------------------------------------------
      isdir_L = .false.
      type = INPUT_FILES_ANAL
      if (present(F_isdir_L)) isdir_L = F_isdir_L
      if (present(F_type)) type = F_type
      F_istat = input_files_set_name1(DEFAULT_LIST_ID,F_key_S,F_filename_S,isdir_L,type)
      !----------------------------------------------------------------------
      return
   end function input_files_set_name0


   !/@*
   function input_files_set_name1(F_id,F_key_S,F_filename_S,F_isdir_L,F_type) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_key_S,F_filename_S
      logical,optional,intent(in) :: F_isdir_L
      integer,optional,intent(in) :: F_type
      !@return
      integer :: F_istat
      !*@/
      character(len=8) :: key_S
      integer :: istat,ifile,type
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input_files) set_name [BEGIN]')
      call priv_init()
      F_istat = RMN_ERR
      if (F_id < 0 .or. F_id > NMAX_NLIST) then
         call msg(MSG_WARNING,'(input_files) set_name, index out of range, ignoring')
         return
      endif
      type = INPUT_FILES_ANAL
      if (present(F_type)) type = min(max(1,F_type),INPUT_FILES_NTYPES)

      ifile = input_files_get_idx(F_id,F_key_S)
      if (ifile < 1) then
!!$         call msg(MSG_WARNING,'(input_files) set_name, unknown file key name: '//trim(F_key_S))
!!$         return
         if (m_nfiles(F_id) >= INPUT_FILES_NMAX) then
            call msg(MSG_WARNING,'(input_files) set_name, too many file key name, canot set: '//trim(F_key_S))
            return
         endif
         key_S = F_key_S
         istat = clib_tolower(key_S)
         ifile = m_nfiles(F_id) + 1
         m_files(F_id,ifile)%tag_S = key_S
         m_files(F_id,ifile)%name_S = F_filename_S
         m_files(F_id,ifile)%isdir_L = .false.
      endif
      if (present(F_isdir_L)) m_files(F_id,ifile)%isdir_L = F_isdir_L

      if (m_files(F_id,ifile)%isdir_L) then
         F_istat = clib_isdir(trim(m_basedir_S(F_id))//'/'//trim(F_filename_S))
         if (.not.RMN_IS_OK(F_istat)) then
            F_istat = clib_isfile(trim(m_basedir_S(F_id))//'/'//trim(F_filename_S))
         endif
      else
         F_istat = clib_isfile(trim(m_basedir_S(F_id))//'/'//trim(F_filename_S))
      endif
      if (RMN_IS_OK(F_istat)) then
         m_files(F_id,ifile)%name_S = trim(F_filename_S)
         m_files(F_id,ifile)%type = type
         call msg(MSG_INFOPLUS,'(input_files) set_name: '//trim(F_key_S)//' => '//trim(F_filename_S))
      else
         call msg(MSG_WARNING,'(input_files) set_name, not such dir/file: '// &
              trim(m_basedir_S(F_id))//'/'//trim(F_filename_S))
         return
      endif

      m_nfiles(F_id) = max(m_nfiles(F_id),ifile)
      call msg(MSG_DEBUG,'(input_files) set_name [END]')
      !----------------------------------------------------------------------
      return
   end function input_files_set_name1


   !/@*
   function input_files_set_basedir0(F_basedir_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_basedir_S
      !@return
      integer :: F_istat
      !*@/
      !----------------------------------------------------------------------
      F_istat = input_files_set_basedir1(DEFAULT_LIST_ID,F_basedir_S)
      !----------------------------------------------------------------------
      return
   end function input_files_set_basedir0


   !/@*
   function input_files_set_basedir1(F_id,F_basedir_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_basedir_S
      !@return
      integer :: F_istat
      !*@/
      character(len=1024) :: basedir_S
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input_files) set_basedir [BEGIN]')
      call priv_init()
      F_istat = RMN_ERR
      if (F_id < 0 .or. F_id > NMAX_NLIST) then
         call msg(MSG_WARNING,'(input_files) set_basedir, index out of range, ignoring')
         return
      endif
      F_istat = clib_isdir(trim(F_basedir_S))
      F_istat = min(clib_realpath(trim(F_basedir_S),basedir_S),F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(input_files) set_basedir, not such dir: '//trim(F_basedir_S))
         return
      endif
      m_basedir_S(F_id) = trim(basedir_S)
      F_istat = RMN_OK
      call msg(MSG_DEBUG,'(input_files) set_basedir [END]')
      !----------------------------------------------------------------------
      return
   end function input_files_set_basedir1


   !/@*
   function input_files_get_idx(F_id,F_name_S) result(F_fileidx)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id
      character(len=*),intent(in) :: F_name_S
      !@return
      integer :: F_fileidx
      !*@/
      character(len=8) :: name_S
      integer :: istat,ifile
      !------------------------------------------------------------------
      call priv_init()
      F_fileidx = RMN_ERR
      if (F_id < 0 .or. F_id > NMAX_NLIST) then
         call msg(MSG_WARNING,'(input_files) get_idx, index out of range')
         return
      endif
      name_S = F_name_S
      istat = clib_tolower(name_S)
      do ifile=1,m_nfiles(F_id)
         if (name_S == m_files(F_id,ifile)%tag_S) then
            F_fileidx = ifile
            exit
         endif
      enddo
      !# Truncated names search kept for backward compatibility
      if (.not.RMN_IS_OK(F_fileidx)) then
         do ifile=1,m_nfiles(F_id)
            if (name_S(1:4) == m_files(F_id,ifile)%tag_S(1:4)) then
               F_fileidx = ifile
               exit
            endif
         enddo
      endif
      !------------------------------------------------------------------
      return
   end function input_files_get_idx


   !/@*
   function input_files_get_type(F_id,F_file_idx) result(F_type)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_file_idx
      !@return
      integer :: F_type
      !*@/
      integer :: fileidx
      !------------------------------------------------------------------
      call priv_init()
      F_type = RMN_ERR
      if (F_id < 0 .or. F_id > NMAX_NLIST) then
         call msg(MSG_WARNING,'(input_files) get_type, index out of range')
         return
      endif
      fileidx = min(max(1,F_file_idx),m_nfiles(F_id))
      F_type = m_files(F_id,fileidx)%type
      !------------------------------------------------------------------
      return
   end function input_files_get_type


   !/@*
   function input_files_open0(F_file_idx) result(F_istat)
      implicit none
      integer,intent(in) :: F_file_idx
      integer :: F_istat
      !*@/
      !------------------------------------------------------------------
      F_istat = input_files_open1(DEFAULT_LIST_ID,F_file_idx)
      !------------------------------------------------------------------
      return
   end function input_files_open0


   !/@*
   function input_files_open1(F_id,F_file_idx) result(F_istat)
      implicit none
      integer,intent(in) :: F_id,F_file_idx
      integer :: F_istat
      !*@/
      character(len=1024) :: filename_S,filename0_S
      integer :: istat,fileidx
      !------------------------------------------------------------------
      call priv_init()
      call msg(MSG_DEBUG,'(input_files) open [BEGIN]')
      F_istat = RMN_ERR
      if (F_id < 0 .or. F_id > NMAX_NLIST) then
         call msg(MSG_WARNING,'(input_files) open, index out of range, ignoring')
         return
      endif
      fileidx = min(max(1,F_file_idx),m_nfiles(F_id))
      if (m_files(F_id,fileidx)%id < 1) then
         if (m_basedir_S(F_id) == ' ') m_basedir_S(F_id) = '.'
         filename0_S = trim(m_basedir_S(F_id))//'/'//m_files(F_id,fileidx)%name_S
         istat = clib_realpath(filename0_S,filename_S)
         if (.not.RMN_IS_OK(istat)) filename_S = filename0_S
         m_files(F_id,fileidx)%id = fstmpi_open(trim(filename_S),FST_READONLY,m_files(F_id,fileidx)%isdir_L)
      endif
      if (m_files(F_id,fileidx)%id < 1) then
         m_files(F_id,fileidx)%id = FILE_AVAIL
         call msg(MSG_ERROR,'(input_files) Unable to open: '//trim(m_basedir_S(F_id))// &
              '/'//m_files(F_id,fileidx)%name_S)
      else
         F_istat = m_files(F_id,fileidx)%id
!!$         call msg(MSG_INFOPLUS,'(input_files) opened file: '//trim(m_basedir_S(F_id))// &
!!$         '/'//m_files(F_id,fileidx)%name_S)
      endif
      call msg(MSG_DEBUG,'(input_files) open [END]')
      !------------------------------------------------------------------
      return
   end function input_files_open1


   !/@*
   function input_files_close0() result(F_istat)
      implicit none
      !@objective
      !@return
      integer :: F_istat
      !*@/
      !----------------------------------------------------------------------
      F_istat = input_files_close1(DEFAULT_LIST_ID)
      !----------------------------------------------------------------------
      return
   end function input_files_close0


   !/@*
   function input_files_close1(F_id) result(F_istat)
      implicit none
      !@objective
      !@return
      integer,intent(in) :: F_id
      integer :: F_istat
      !*@/
      integer :: ii,istat
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input_files) close [BEGIN]')
      call priv_init()
      F_istat = RMN_ERR
      if (F_id < 0 .or. F_id > NMAX_NLIST) then
         call msg(MSG_WARNING,'(input_files) close, index out of range, ignoring')
         return
      endif
      F_istat = RMN_OK
      do ii = 1,m_nfiles(F_id)
         if (m_files(F_id,ii)%id > 0) then
            istat = fstmpi_close(m_files(F_id,ii)%id)
            if (RMN_IS_OK(istat)) m_files(F_id,ii)%id = FILE_AVAIL
            F_istat = min(F_istat,istat)
         endif
      enddo
      call msg(MSG_DEBUG,'(input_files) close [END]')
      !----------------------------------------------------------------------
      return
   end function input_files_close1


   !/@*
   function input_files_vgrid_4(F_id, F_file_idx, F_varname_S, F_ip2, &
        F_datev, F_vgrid_S, F_sfcfld, F_sfcgridid, &
        F_sfcfld2, F_sfcgridid2) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_file_idx,F_ip2
      integer,intent(in) :: F_datev
      character(len=*),intent(in) :: F_varname_S
      character(len=*),intent(out) :: F_vgrid_S
      real,pointer,optional :: F_sfcfld(:,:,:)
      integer,intent(out),optional :: F_sfcgridid
      real,pointer,optional :: F_sfcfld2(:,:,:)
      integer,intent(out),optional :: F_sfcgridid2
      !@return
      integer :: F_istat
      !*@/
      integer(INT64) :: jdatev
      !------------------------------------------------------------------
      jdatev = jdate_from_cmc(F_datev)
      if (present(F_sfcfld).and.present(F_sfcgridid)) then
         if (present(F_sfcfld2).and.present(F_sfcgridid2)) then
            F_istat = input_files_vgrid_8(F_id, F_file_idx, F_varname_S, &
                 F_ip2, jdatev, F_vgrid_S, F_sfcfld, F_sfcgridid, &
                 F_sfcfld2, F_sfcgridid2)
         else
            F_istat = input_files_vgrid_8(F_id, F_file_idx, F_varname_S, &
                 F_ip2, jdatev, F_vgrid_S, F_sfcfld, F_sfcgridid)
         endif
      else
         F_istat = input_files_vgrid_8(F_id,F_file_idx,F_varname_S,F_ip2,jdatev,F_vgrid_S)
      endif
      !------------------------------------------------------------------
      return
   end function input_files_vgrid_4


   !/@*
   function input_files_vgrid_8(F_id, F_file_idx, F_varname_S, F_ip2, &
        F_jdatev, F_vgrid_S, F_sfcfld, F_sfcgridid, &
        F_sfcfld2, F_sfcgridid2) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_file_idx,F_ip2
      integer(INT64),intent(in) :: F_jdatev
      character(len=*),intent(in) :: F_varname_S
      character(len=*),intent(out) :: F_vgrid_S
      real,pointer,optional :: F_sfcfld(:,:,:)
      integer,intent(out),optional :: F_sfcgridid
      real,pointer,optional :: F_sfcfld2(:,:,:)
      integer,intent(out),optional :: F_sfcgridid2
      !@return
      integer :: F_istat
      !*@/
      type(vgrid_descriptor) :: vgrid, vgrid0
      integer,pointer :: ip1list(:)
      integer :: fileid, istat, sfcRefKey, sfcRefKey2
      character(len=4) :: levtype_S, sfcfld_S, sfcfld2_S
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input_files) vgrid [BEGIN]')
      call priv_init()
      F_istat = RMN_ERR
      F_vgrid_S = ' '
      fileid = input_files_open1(F_id,F_file_idx)
      if (.not.RMN_IS_OK(fileid)) return
      !TODO: for climatic field, should use F_ip2 instead of F_datev
      if (F_ip2 /= RMN_ANY_I) then
         call msg(MSG_WARNING,'(input_files) vgrid, provided ip2 ignored')
      endif
      nullify(ip1list)
      if (present(F_sfcfld).and.present(F_sfcgridid)) then
         nullify(F_sfcfld)
         F_sfcgridid = RMN_ERR
         istat = vgrid_from_file_mpi(fileid, F_varname_S, F_jdatev, vgrid, &
              ip1list, levtype_S, F_sfcRefKey=sfcRefKey, &
              F_sfcRefKey2=sfcRefKey2)
         if (RMN_IS_OK(sfcRefKey)) then
            istat = fstmpi_read(sfcRefKey, F_sfcfld, fileid, F_sfcgridid)
            if (.not.RMN_IS_OK(istat)) &
                 call msg(MSG_WARNING,'(input_files) Problem reading vgrid Sfc Ref field')
         endif
         if (present(F_sfcfld2) .and. present(F_sfcgridid2 ).and. &
              RMN_IS_OK(sfcRefKey2)) then
            nullify(F_sfcfld)
            F_sfcgridid2 = RMN_ERR
            istat = fstmpi_read(sfcRefKey2, F_sfcfld2, fileid, F_sfcgridid2)
            if (.not.RMN_IS_OK(istat)) &
                 call msg(MSG_WARNING,'(input_files) Problem reading vgrid Sfc Ref field')
         endif
      else
         istat = vgrid_from_file_mpi(fileid, F_varname_S, F_jdatev, vgrid, &
              ip1list, levtype_S)
      endif
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(input_files) vgrid, problem getting vgrid from file for '//trim(F_varname_S))
         return
      endif
      istat = vgd_get(vgrid, 'RFLD', sfcfld_S)
      if (istat /= VGD_OK) sfcfld_S = ' '
      istat = vgd_get(vgrid, 'RFLS', sfcfld2_S)
      if (istat /= VGD_OK) sfcfld2_S = ' '
      F_vgrid_S = 'in/'//trim(F_varname_S)
      istat = vgrid_wb_get(F_vgrid_S, vgrid0)
      if (.not.RMN_IS_OK(istat)) then
         istat = vgrid_wb_put(F_vgrid_S, vgrid, ip1list, sfcfld_S, sfcfld2_S)
      elseif (.not.(vgrid0 == vgrid)) then
         call msg(MSG_WARNING,'(input_files) vgrid, ignoring an inconsistant vgrid for '//trim(F_varname_S))
         return
      endif
      istat = vgd_free(vgrid0)
      if (present(F_sfcfld).and.present(F_sfcgridid)) then
         if (associated(F_sfcfld) .and. any(sfcfld_S == (/'P0','p0'/))) then
            F_sfcfld = F_sfcfld * MB2PA
         endif
      endif
      if (present(F_sfcfld2).and.present(F_sfcgridid2)) then
         if (associated(F_sfcfld2) .and. any(sfcfld2_S == (/'P0','p0'/))) then
            F_sfcfld2 = F_sfcfld2 * MB2PA
         endif
      endif
      if (associated(ip1list)) deallocate(ip1list,stat=istat)
      F_istat = RMN_OK
      call msg(MSG_DEBUG,'(input_files) vgrid [END]')
      !------------------------------------------------------------------
      return
   end function input_files_vgrid_8


   !/@*
   function input_files_sfcref_4(F_id, F_file_idx, F_datev, F_vgrid, &
        F_sfcfld, F_sfcgridid, F_rfls_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_file_idx
      integer,intent(in) :: F_datev
      type(vgrid_descriptor),intent(in) :: F_vgrid
      real,pointer :: F_sfcfld(:,:,:)
      integer,intent(out) :: F_sfcgridid
      logical,intent(in), optional :: F_rfls_L
      !@return
      integer :: F_istat
      !*@/
      integer(INT64) :: jdatev
      logical :: rfls_L
      !------------------------------------------------------------------
      jdatev = jdate_from_cmc(F_datev)
      rfls_L = .false.
      if (present(F_rfls_L)) rfls_L = F_rfls_L
      F_istat = input_files_sfcref_8(F_id, F_file_idx, jdatev, F_vgrid, &
           F_sfcfld, F_sfcgridid, rfls_L)
      !------------------------------------------------------------------
      return
   end function input_files_sfcref_4


   !/@*
   function input_files_sfcref_8(F_id, F_file_idx, F_jdatev, F_vgrid, &
        F_sfcfld, F_sfcgridid, F_rfls_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_file_idx
      integer(INT64),intent(in) :: F_jdatev
      type(vgrid_descriptor),intent(in) :: F_vgrid
      real,pointer :: F_sfcfld(:,:,:)
      integer,intent(out) :: F_sfcgridid
      logical,intent(in), optional :: F_rfls_L
      !@return
      integer :: F_istat
      !*@/
      integer :: fileid,sfcRefKey,istat
      character(len=4) :: sfcfld_S
      logical :: rfls_L
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input_files) sfcref  [BEGIN]')
      call priv_init()
      F_istat = RMN_ERR
      F_sfcgridid = RMN_ERR
      rfls_L = .false.
      if (present(F_rfls_L)) rfls_L = F_rfls_L
      fileid = input_files_open1(F_id, F_file_idx)
      nullify(F_sfcfld)
      sfcRefKey = vgrid_from_file_rfld_key_mpi(fileid, F_jdatev, F_vgrid, &
           rfls_L)
      if (sfcRefKey == VGRID_FROM_FILE_NORFLD) then
         F_istat = RMN_OK
         return
      endif
      F_istat = fstmpi_read(sfcRefKey, F_sfcfld, fileid, F_sfcgridid)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(input_files) Problem reading vgrid Sfc Ref field')
      endif
      if (rfls_L) then
         istat = vgd_get(F_vgrid, key='RFLS', value=sfcfld_S)
      else
         istat = vgd_get(F_vgrid, key='RFLD', value=sfcfld_S)
      endif
      if (associated(F_sfcfld) .and. any(sfcfld_S == (/'P0','p0'/))) then
         F_sfcfld = F_sfcfld * MB2PA
      endif
      call msg(MSG_DEBUG,'(input_files) sfcref [END]')
      !------------------------------------------------------------------
      return
   end function input_files_sfcref_8


!!$   !/@*
!!$   function input_files_vgrid(F_fld,F_id,F_file_idx) result(F_istat)
!!$      implicit none
!!$      type(input_var_T),intent(inout) :: F_fld
!!$      integer,intent(in) :: F_file_idx
!!$      integer :: F_istat
!!$      !*@/
!!$      integer :: istat,istat2,fileid
!!$      character(len=4) :: levtyp_S
!!$      !------------------------------------------------------------------
!!$      call msg(MSG_DEBUG,'(input_files) vgrid_from_file [BEGIN]')
!!$      F_istat = RMN_ERR
!!$      fileid = priv_openfile(F_file_idx)
!!$      nullify(F_fld%ip1list)
!!!      istat = vgrid_from_file_mpi(fileid,F_fld%vn1_S,RMN_ANY_DATE,m_files(F_id,F_file_idx)%vgrid,F_fld%ip1list,levtyp_S)
!!$      if (.not.RMN_IS_OK(istat)) then
!!$         call msg(MSG_WARNING,'(input) Cannot get vgrid in file for: '//trim(F_fld%vn1_S))
!!$         return
!!$      endif
!!$      if (.not.any(levtyp_S(1:1) == (/'M','T','m','t'/))) then
!!$         call msg(MSG_WARNING,'(input_files) Cannot vertically interpolate arbitrary levels, for: '//trim(F_fld%vn1_S))
!!$         return
!!$      endif
!!$      F_fld%vgrid => m_files(F_id,F_file_idx)%vgrid
!!$      F_fld%nk = size(F_fld%ip1list)
!!$      nullify(F_fld%d1,F_fld%d2)
!!$      allocate(F_fld%d1(F_fld%ni,F_fld%nj,F_fld%nk),stat=istat)
!!$      istat2 = 0
!!$      if (F_fld%vn2_S /= ' ') &
!!$           allocate(F_fld%d2(F_fld%ni,F_fld%nj,F_fld%nk),stat=istat2)
!!$      if (istat /= 0 .or. istat2 /= 0 .or. .not.associated(F_fld%d1)) then
!!$         call msg(MSG_WARNING,'(input_files) Problem allocating memory')
!!$         return
!!$      endif
!!$      F_istat = RMN_OK
!!$      call msg(MSG_DEBUG,'(input_files) vgrid_from_file [END]')
!!$      !------------------------------------------------------------------
!!$      return
!!$   end function input_files_vgrid


   !/@*
   function input_files_read_4(F_id,F_file_idx,F_varname_S,F_ip1,F_ip2,F_datev, &
        F_datevfuzz,F_fuzztype,F_ingridid,F_data,F_typvar_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_id,F_file_idx,F_datevfuzz,F_fuzztype
      integer,intent(inout) :: F_ip1,F_ip2,F_ingridid
      integer,intent(inout) :: F_datev
      character(len=*),intent(in) :: F_varname_S
      character(len=*),intent(in),optional :: F_typvar_S
      !@return
      integer :: F_istat
      !*@/
      integer(INT64) :: jdatev
      !------------------------------------------------------------------
      jdatev = jdate_from_cmc(F_datev)
      if (present(F_typvar_S)) then
         F_istat = input_files_read_8(F_id,F_file_idx,F_varname_S,F_ip1,F_ip2, &
              jdatev,F_datevfuzz,F_fuzztype,F_ingridid,F_data,F_typvar_S)
      else
         F_istat = input_files_read_8(F_id,F_file_idx,F_varname_S,F_ip1,F_ip2, &
              jdatev,F_datevfuzz,F_fuzztype,F_ingridid,F_data)
      endif
      if (RMN_IS_OK(F_istat)) F_datev = jdate_to_cmc(jdatev)
      !------------------------------------------------------------------
      return
   end function input_files_read_4


   !/@*
   function input_files_read_8(F_id,F_file_idx,F_varname_S,F_ip1,F_ip2,F_jdatev, &
        F_datevfuzz,F_fuzztype,F_ingridid,F_data,F_typvar_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: F_id,F_file_idx,F_datevfuzz,F_fuzztype
      integer,intent(inout) :: F_ip1,F_ip2,F_ingridid
      integer(INT64),intent(inout) :: F_jdatev
      character(len=*),intent(in) :: F_varname_S
      character(len=*),intent(in),optional :: F_typvar_S
      !@return
      integer :: F_istat
      !*@/
      integer :: key,istat,dateo,datev,deet,npas,ip3,ingridid,k,msgLevelMin,msgUnit,fileid
      character(len=512) :: dummy_S,dummy2_S,nomvar_S,etiket_S,typvar0_S,typvar_S,datev_S,msgFormat_S
      real :: nhours
      logical :: canWrite_L
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input_files) read [BEGIN]')
      call msg_getInfo(canWrite_L,msgLevelMin,msgUnit,msgFormat_S)
      F_istat = RMN_ERR
      typvar0_S = RMN_ANY_TYP
      if (present(F_typvar_S)) typvar0_S = F_typvar_S
      fileid = input_files_open1(F_id,F_file_idx)
      if (.not.RMN_IS_OK(fileid)) return
      datev = jdate_to_cmc(F_jdatev)
      key = fstmpi_find(fileid,F_varname_S,datev,F_ip1,F_ip2,RMN_ANY_I,F_datevfuzz, &
           F_fuzztype,F_typvar_S=typvar0_S)
      if (.not.RMN_IS_OK(key)) then
         datev_S = '-1'
         if (F_jdatev /= MU_JDATE_ANY) datev_S = jdate_to_print(F_jdatev)
         dummy_S = ' '
         write(dummy_S,'(a4,a,3I9,a,i10,i3,a)') trim(F_varname_S), &
              ' (datev='//trim(adjustl(datev_S))//') [ip123=',F_ip1, F_ip2, RMN_ANY_I, &
              ', tv='//typvar0_S(1:1)//'] [datefuzz=',F_datevfuzz,F_fuzztype,']'
         call msg(MSG_INFO,'(input_files) Could not find: '//trim(dummy_S))
         return
      endif

      F_jdatev = jdate_from_cmc(datev)
      datev_S = '-1'
      if (F_jdatev /= MU_JDATE_ANY) datev_S = jdate_to_print(F_jdatev)
      istat = fstmpi_getmeta(key,nomvar_S,dateo,deet,npas,F_ip1, F_ip2, ip3, etiket_S,typvar_S)
      dummy_S = ' '
      if (RMN_IS_OK(istat)) then
         if (any(typvar_S(1:1)==(/'C','c'/))) then
            write(dummy_S,'(a4,a,3I9,a)') trim(nomvar_S),' (datev=-1) [ip123=', &
                 F_ip1, F_ip2, ip3,', tv='//typvar0_S(1:1)//']'
         else
            nhours = real((dble(deet)*dble(npas))/3600.)
            write(dummy_S,'(a4,a,F8.1,a,3I9,a)') trim(nomvar_S),' (datev='// &
                 trim(adjustl(datev_S))//'; nhours=',nhours,') [ip123=', &
                 F_ip1, F_ip2, ip3,', tv='//typvar0_S(1:1)//']'
         endif
      else
         write(dummy_S,'(a4,a,3I9,a,i10,i3,a)') trim(F_varname_S),' (datev='// &
              trim(adjustl(datev_S))//') [ip123=',F_ip1, F_ip2, RMN_ANY_I,', tv='// &
              typvar0_S(1:1)//'] [datefuzz=',F_datevfuzz,F_fuzztype,']'
      endif

      istat = fstmpi_read(key,F_data)
      if (RMN_IS_OK(istat)) then
         if (MSG_INFOPLUS >= msgLevelMin) then
            do k=lbound(F_data,3),ubound(F_data,3)
               write(dummy2_S,'(a,a,i3,a,1pe14.7,a,1pe14.7,a,1pe14.7,a)') &
                    '('//'input'//') Read:',trim(dummy_S),k,&
                    ' [Mean: ',sum(F_data(:,:,k))/max(1.,real(size(F_data(:,:,k)))),&
                    ' ; Min: ',minval(F_data(:,:,k)),&
                    ' ; Max: ',maxval(F_data(:,:,k)),']'
               call msg(MSG_INFO,dummy2_S)
            enddo
         else
            call msg(MSG_INFO,'(input_files) Read:'//trim(dummy_S))
         endif
      else
         call msg(MSG_WARNING,'(input_files) Problem reading field: '//trim(dummy_S))
         return
      endif

      ingridid = fstmpi_get_gridid(fileid,key)
      if (.not.RMN_IS_OK(ingridid)) then
         call msg(MSG_WARNING,'(input_files) Could not get grid for field: '//trim(dummy_S))
         return
      endif
      if (F_ingridid >= 0 .and. F_ingridid /= ingridid) then
         if (.not.ezgrid_samegrid(F_ingridid,ingridid)) then
            call msg(MSG_WARNING,'(input_files) Inconsitent/wrong grid for field: '//trim(dummy_S))
            return
         else
            istat = gdrls(ingridid)
            ingridid = F_ingridid
         endif
      endif
      F_ingridid = ingridid
      F_istat = RMN_OK
      write(dummy2_S, '(i0,1x,i0,1x,i0,1x,i0,1x,i0)') F_istat, F_ip1, F_ip2, F_ingridid, F_jdatev
      call msg(MSG_DEBUG,'(input_files) read [END]'//dummy2_S)
      !------------------------------------------------------------------
      return
   end function input_files_read_8


   !==== Private Functions =================================================


   !/@*
   subroutine priv_init()
      implicit none
      !*@/
      integer,parameter :: NBASE_FILES = 7
      character(len=8),parameter :: BASEFILE_TAG_S(NBASE_FILES) = &
           (/'GEOPHY  ', 'CLIMATO ', 'ANALYSIS', &
             'INREP   ', 'MANAL   ', 'MINREP  ', &
             'MINPUT  '/)
      character(len=16),parameter :: BASEFILE_NAME_S(NBASE_FILES) = &
           (/'GEOPHY          ', 'CLIMATO         ', 'ANALYSIS        ',&
             'INREP           ', 'MODEL_ANALYSIS  ', 'MODEL_INREP     ',&
             'MODEL_INPUT     '/)
      logical,parameter :: BASEFILE_ISDIR_L(NBASE_FILES) = &
           (/.true., .true., .true., &
             .true., .true., .true., &
             .true./)
      integer,parameter :: BASEFILE_TYPE(NBASE_FILES) = &
           (/INPUT_FILES_GEOP, INPUT_FILES_CLIM, INPUT_FILES_ANAL, &
             INPUT_FILES_ANAL, INPUT_FILES_ANAL, INPUT_FILES_ANAL, &
             INPUT_FILES_ANAL/)
      logical,save :: init_L = .false.
      integer :: ifile,istat
      character(len=8) :: name_S
      !------------------------------------------------------------------
      if (init_L) return
      m_nfiles = 0
      m_basedir_S = '.'
      do ifile=1,INPUT_FILES_NMAX
         m_files(:,ifile)%tag_S = ' '
         m_files(:,ifile)%name_S = ' '
         m_files(:,ifile)%id = FILE_AVAIL
         m_files(:,ifile)%type = INPUT_FILES_ANAL
         m_files(:,ifile)%isdir_L = .false.
      enddo
      m_nfiles = NBASE_FILES
      do ifile=1,NBASE_FILES
         name_S = BASEFILE_TAG_S(ifile)
         istat = clib_tolower(name_S)
         m_files(:,ifile)%tag_S = name_S
         m_files(:,ifile)%name_S = BASEFILE_NAME_S(ifile)
         m_files(:,ifile)%type = BASEFILE_TYPE(ifile)
         m_files(:,ifile)%isdir_L = BASEFILE_ISDIR_L(ifile)
      enddo
      init_L = .true.
      !------------------------------------------------------------------
      return
   end subroutine priv_init

end module input_files_mod

