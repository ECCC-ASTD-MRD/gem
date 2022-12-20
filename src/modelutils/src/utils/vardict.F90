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

!/@*
module vardict_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_isfile, clib_isreadok
   use wb_itf_mod
   use str_mod
   implicit none
   !@author Stephane Chamberland, 2012-05
   !@Objective Var/Fields dictionary: Name, Context/Pkg, Units, Description
   private
   public :: vardict_add_string,vardict_add,vardict_get,vardict_fromfile,vardict_getlist
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>


   character(len=*),parameter :: PREFIX = 'VD#/'
   character(len=*),parameter :: PREFIX2 = 'VD%/'
   character(len=*),parameter :: WBFMT0 = '(a)'
   character(len=*),parameter :: WBFMT1 = '(a,a,"/")'
   character(len=*),parameter :: WBFMT2 = '(a,a,"/",a,"/")'
   character(len=*),parameter :: PKG_DEFAULT = '-'
   integer,parameter :: IDX_NAME = 1
   integer,parameter :: IDX_PKG = 2
   integer,parameter :: IDX_UNITS = 3
   integer,parameter :: IDX_DESC = 4
   integer,parameter :: IDX_NAMEIN = 5
   integer,parameter :: IDX_NAMEOUT = 6
   integer,parameter :: NIDX = 6
   integer,parameter :: NIDX_MIN = 1
   integer,parameter :: NIDX_EXTRA = 32 - NIDX
   integer,parameter :: NIDX_MAX = NIDX + NIDX_EXTRA
   integer,parameter :: KEY = 1
   integer,parameter :: VAL = 2
   logical,parameter :: NAMEFIRST = .true.
   logical,parameter :: PKGFIRST = .false.

   character(len=WB_MAXNAMELENGTH),parameter :: KNOWN_KEYS(NIDX) = (/&
        'vn     ', & !# mandatory
        'vb     ', &
        'vu     ', &
        'vd     ', &
        'in     ', & 
        'on     ' &
        /)

contains

   !/@*
   function vardict_add_string(F_meta_S,F_print_L) result(F_istat)
      implicit none
      !@objective Register a new var with meta in the system
      !@arguments
      character(len=*),intent(in) :: F_meta_S
      logical,intent(in),optional :: F_print_L
      !@return
      integer :: F_istat
      !*@/
      integer,external :: str_split2keyval
      integer :: nkeys,n,nextra
      logical :: print_L
      character(len=WB_MAXSTRINGLENGTH) :: kv_S(2,NIDX_MAX),extra_S(NIDX_EXTRA)
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_meta_S == ' ') then
         call msg(MSG_WARNING,'(vardict) New Var from string, no Meta')
         return
      endif
      print_L = .true.
      if (present(F_print_L)) print_L = F_print_L
      kv_S = ' '
      kv_S(KEY,1:NIDX) = KNOWN_KEYS(1:NIDX)
      nkeys = str_split2keyval(kv_S,F_meta_S,NIDX_MAX)
      if (nkeys < NIDX_MIN .or. any(kv_S(VAL,1:NIDX_MIN) == ' ')) then
         call msg(MSG_ERROR,'(vardict) newvar - Mandatory params not all provided for: '//trim(F_meta_S))
         return
      endif
      extra_S = ' '
      nextra = 0
      do n=NIDX+1,nkeys
         nextra = nextra+1
         extra_S(nextra) = trim(kv_S(KEY,n))//'='//trim(kv_S(VAL,n))
      enddo
      nextra = max(1,nextra)
      F_istat = vardict_add(kv_S(VAL,IDX_NAME),kv_S(VAL,IDX_PKG),kv_S(VAL,IDX_NAMEIN), &
           kv_S(VAL,IDX_NAMEOUT),kv_S(VAL,IDX_UNITS),kv_S(VAL,IDX_DESC),extra_S(1:nextra),print_L)
      !---------------------------------------------------------------
      return
   end function vardict_add_string


   !/@*
   function vardict_add(F_varname_S,F_pkg_S,F_namein_S,F_nameout_S, &
        F_units_S,F_desc_S,F_extra_S,F_print_L) result(F_istat)
      implicit none
      !@objective Register a new var with meta in the system
      !@arguments
      character(len=*),intent(in) :: F_varname_S
      character(len=*),intent(in),optional :: F_pkg_S,F_namein_S,F_nameout_S, &
           F_units_S,F_desc_S,F_extra_S(:)
      logical,intent(in),optional :: F_print_L
      !@return
      integer :: F_istat
      !*@/
      integer :: n,n1,n2,istat,nextra,iverb
      logical :: print_L
      character(len=WB_MAXNAMELENGTH) :: name_S
      character(len=WB_MAXSTRINGLENGTH) :: meta_S(NIDX_MAX), meta0_S(NIDX_MAX),msg_S
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_varname_S == ' ') then
         call msg(MSG_WARNING,'(vardict) New Var, invalid varname: '//trim(F_varname_S))
         return
      endif
      print_L = .true.
      if (present(F_print_L)) print_L = F_print_L
      meta_S = ' '
      meta_S(IDX_NAME) = F_varname_S
      meta_S(IDX_PKG) = PKG_DEFAULT
      if (present(F_pkg_S)) meta_S(IDX_PKG) = F_pkg_S
      if (present(F_namein_S)) meta_S(IDX_NAMEIN) = F_namein_S
      if (present(F_nameout_S)) meta_S(IDX_NAMEOUT) = F_nameout_S
      if (present(F_units_S)) meta_S(IDX_UNITS) = F_units_S
      if (present(F_desc_S)) meta_S(IDX_DESC) = F_desc_S
      nextra = 0
      if (present(F_extra_S)) then
         nextra = min(size(F_extra_S),NIDX_EXTRA)
         meta_S(NIDX+1:NIDX+nextra) = F_extra_S(1:nextra)
      endif
      n1 = max(8,len_trim(meta_S(IDX_NAME)))
      n2 = max(8,len_trim(meta_S(IDX_UNITS)))
      do n=1,NIDX_MAX
         call str_normalize(meta_S(n))
      enddo
      call priv_name(name_S,meta_S(IDX_PKG),meta_S(IDX_NAME),PKGFIRST)
      iverb = wb_verbosity(WB_MSG_FATAL)
      istat = wb_get(name_S,meta0_S,n)
      iverb = wb_verbosity(iverb)
      if (RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(vardict) Cannot redefine'//&
              ': VN='//meta_S(IDX_NAME)(1:n1)//&
              '; VB='//meta_S(IDX_PKG)(1:4))
         return
      endif
      F_istat = wb_put(name_S,meta_S)
      call priv_name(name_S,meta_S(IDX_PKG),meta_S(IDX_NAME),NAMEFIRST)
      F_istat = min(wb_put(name_S,meta_S),F_istat)
      if (RMN_IS_OK(F_istat)) then
         if (print_L) then
            msg_S = '(vardict) VN='//meta_S(IDX_NAME)(1:n1)//&
                 '; VB='//meta_S(IDX_PKG)(1:4)//&
                 '; IN='//meta_S(IDX_NAMEIN)(1:4)//&
                 '; ON='//meta_S(IDX_NAMEOUT)(1:4)//&
                 '; VU='//meta_S(IDX_UNITS)(1:n2)//&
                 '; VD="'//trim(meta_S(IDX_DESC))//'"'
            do n=1,nextra
               msg_S = trim(msg_S)//'; '//trim(meta_S(NIDX+n))
            enddo
            call msg(MSG_INFO,msg_S)
         endif
      else
         call msg(MSG_WARNING,'(vardict) Probleme registering New Var'//&
              ': VN='//meta_S(IDX_NAME)(1:n1)//&
              '; VB='//meta_S(IDX_PKG)(1:4))
      endif
      !---------------------------------------------------------------
      return
   end function vardict_add


   !/@*
   function vardict_get(F_varname_S,F_pkg_S,F_namein_S,F_nameout_S,F_units_S, &
        F_desc_S,F_extra_S,F_idx0) result(F_istat)
      implicit none
      !@objective Return Metadata for matching elements (wildcard = ' ')
      !@arguments
      character(len=*),intent(inout) :: F_varname_S,F_pkg_S
      character(len=*),intent(inout),optional :: F_namein_S,F_nameout_S
      character(len=*),intent(out),optional :: F_units_S,F_desc_S,F_extra_S(:)
      integer,intent(inout),optional :: F_idx0 !- Search from provided idx, return idx of item found
      !@return
      integer :: F_istat
      !*@/
      integer,parameter :: NMAX_VARS = 4096
      integer :: n,nval,i,nvar,idx0,istat,nextra,iverb
      character(len=WB_MAXNAMELENGTH) :: name_S,varlist_S(NMAX_VARS)
      character(len=WB_MAXSTRINGLENGTH) :: meta0_S(NIDX_MAX),meta_S(NIDX_MAX),meta1_S(NIDX_MAX)
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      if (present(F_units_S)) F_units_S = ' '
      if (present(F_desc_S)) F_desc_S = ' '
      if (present(F_extra_S)) F_extra_S = ' '

      idx0 = 1
      if (present(F_idx0)) then
         idx0 = max(1,F_idx0)
         F_idx0 = -1
      endif
      meta0_S = ' '
      meta0_S(IDX_NAME) = F_varname_S
      meta0_S(IDX_PKG) = F_pkg_S
      if (present(F_namein_S)) meta0_S(IDX_NAMEIN) = F_namein_S
      if (present(F_nameout_S)) meta0_S(IDX_NAMEOUT) = F_nameout_S
      do n=1,NIDX
         call str_normalize(meta0_S(n))
      enddo

      IF_NAME: if (F_varname_S /= ' '.and. F_pkg_S /= ' ') then

         call priv_name(name_S,meta0_S(IDX_PKG),meta0_S(IDX_NAME),PKGFIRST)
         meta_S = ' '
         iverb = wb_verbosity(WB_MSG_FATAL)
         F_istat = wb_get(name_S,meta_S,nval)
         iverb = wb_verbosity(iverb)
         if (RMN_IS_OK(F_istat)) then
            meta1_S = ' '
            do i=1,NIDX
               if (meta0_S(i) /= ' ') meta1_S(i) = meta_S(i)
            enddo
            if (any(meta0_S /= meta1_S)) F_istat = RMN_ERR
         endif
         idx0 = -1

      else !IF_NAME

         if (F_varname_S /= ' ') then
            call priv_name(name_S,meta0_S(IDX_PKG),meta0_S(IDX_NAME),NAMEFIRST)
         else
            call priv_name(name_S,meta0_S(IDX_PKG),meta0_S(IDX_NAME),PKGFIRST)
         endif
         F_istat = wb_keys(varlist_S,nvar,name_S)
         if (.not.RMN_IS_OK(F_istat)) return
         F_istat = RMN_ERR
         DOVAR: do n=idx0,nvar
            meta_S = ' '
            istat = wb_get(trim(varlist_S(n)),meta_S,nval)
            meta1_S = ' '
            do i=1,NIDX
               if (meta0_S(i) /= ' ') meta1_S(i) = meta_S(i)
            enddo
            if (all(meta0_S == meta1_S)) then
               idx0 = n
               F_istat = RMN_OK
               exit
            endif
         enddo DOVAR

      endif IF_NAME

      if (RMN_IS_OK(F_istat)) then
         F_varname_S = meta_S(IDX_NAME)
         F_pkg_S = meta_S(IDX_PKG)
         if (present(F_namein_S)) F_namein_S = meta_S(IDX_NAMEIN)
         if (present(F_nameout_S)) F_nameout_S = meta_S(IDX_NAMEOUT)
         if (present(F_units_S)) F_units_S = meta_S(IDX_UNITS)
         if (present(F_desc_S)) F_desc_S = meta_S(IDX_DESC)
         if (present(F_idx0)) F_idx0 = idx0
         if (present(F_extra_S)) then
            nextra = min(size(F_extra_S),NIDX_EXTRA)
            F_extra_S(1:nextra) = meta_S(NIDX+1:NIDX+nextra)
         endif
      endif
      !---------------------------------------------------------------
      return
   end function vardict_get


   !/@*
   function vardict_getlist(F_varlist_S,F_varname_S,F_pkg_S,F_namein_S,F_nameout_S) result(F_nvar)
      implicit none
      !@objective Return list of var for matching elements (wildcard = ' ')
      !@arguments
      character(len=*) :: F_varlist_S(:)
      character(len=*),intent(in),optional :: F_varname_S,F_pkg_S
      character(len=*),intent(in),optional :: F_namein_S,F_nameout_S
      !@return
      integer :: F_nvar
      !*@/
      integer,parameter :: NMAX_VARS = 4096
      character(len=WB_MAXSTRINGLENGTH) :: meta0_S(NIDX_MAX),meta_S(NIDX_MAX),meta1_S(NIDX_MAX)
      character(len=WB_MAXNAMELENGTH) :: name_S,varlist_S(NMAX_VARS)
      integer :: n,nval,nvar,nvarmax,i,istat
      !---------------------------------------------------------------
      F_nvar = 0
      F_varlist_S = ' '
      nvarmax = size(F_varlist_S)

      meta0_S = ' '
      if (present(F_varname_S)) meta0_S(IDX_NAME) = F_varname_S
      if (present(F_pkg_S)) meta0_S(IDX_PKG) = F_pkg_S
      if (present(F_namein_S)) meta0_S(IDX_NAMEIN) = F_namein_S
      if (present(F_nameout_S)) meta0_S(IDX_NAMEOUT) = F_nameout_S
      do n=1,NIDX
         call str_normalize(meta0_S(n))
      enddo

      if (F_varname_S /= ' ') then
         call priv_name(name_S,meta0_S(IDX_PKG),meta0_S(IDX_NAME),NAMEFIRST)
      else
         call priv_name(name_S,meta0_S(IDX_PKG),meta0_S(IDX_NAME),PKGFIRST)
      endif

      istat = wb_keys(varlist_S,nvar,name_S)
      if (.not.RMN_IS_OK(istat)) return

      DOVAR: do n=1,nvar
         meta_S = ' '
         istat = wb_get(trim(varlist_S(n)),meta_S,nval)
         meta1_S = ' '
         do i=1,NIDX
            if (meta0_S(i) /= ' ') meta1_S(i) = meta_S(i)
         enddo
         if (all(meta0_S == meta1_S)) then
            F_nvar = min(F_nvar + 1,nvarmax)
            F_varlist_S(F_nvar) = meta_S(IDX_NAME)
         endif
      enddo DOVAR
      !---------------------------------------------------------------
      return
   end function vardict_getlist


!!$   !/@*
!!$   subroutine vardict_printlist(F_varname_S,F_pkg_S,F_namein_S,F_nameout_S) result(F_nvar)
!!$      implicit none
!!$      !@objective Print list of var for matching elements (wildcard = ' ')
!!$      !@arguments
!!$      character(len=*),intent(in),optional :: F_varname_S,F_pkg_S
!!$      character(len=*),intent(in),optional :: F_namein_S,F_nameout_S
!!$      !*@/
!!$      integer,parameter :: NMAX_VARS = 4096
!!$      character(len=WB_MAXSTRINGLENGTH) :: meta0_S(NIDX),meta_S(NIDX)
!!$      character(len=WB_MAXNAMELENGTH) :: name_S,varlist_S(NMAX_VARS)
!!$      integer :: n,nval,nvar,istat
!!$      !---------------------------------------------------------------
!!$      meta0_S = ' '
!!$      if (present(F_varname_S)) meta0_S(IDX_NAME) = F_varname_S
!!$      if (present(F_pkg_S)) meta0_S(IDX_PKG) = F_pkg_S
!!$      if (present(F_namein_S)) meta0_S(IDX_NAMEIN) = F_namein_S
!!$      if (present(F_nameout_S)) meta0_S(IDX_NAMEOUT) = F_nameout_S
!!$
!!$      nvar = vardict_getlist(varlist_S,meta0_S(IDX_NAME),meta0_S(IDX_PKG),meta0_S(IDX_NAMEIN),meta0_S(IDX_NAMEOUT))
!!$
!!$      if (F_pkg_S /= ' ') then
!!$         DOVAR: do n=1,nvar
!!$            call priv_name(name_S,meta0_S(IDX_PKG),varlist_S(n),PKGFIRST))
!!$            istat = wb_get(name_S,meta_S,nval)
!!$
!!$         enddo DOVAR
!!$      else
!!$
!!$      endif
!!$      !---------------------------------------------------------------
!!$      return
!!$   end subroutine vardict_printlist


   !/@*
   function vardict_fromfile(F_filename_S,F_print_L) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      logical,intent(in),optional :: F_print_L
      integer :: F_istat
      !*@/
      logical :: print_L
      integer :: istat,istat2,fileid
      character(len=WB_MAXSTRINGLENGTH) :: string_S
      !---------------------------------------------------------------
      F_istat = RMN_ERR

      istat = clib_isfile(trim(F_filename_S))
      istat = min(clib_isreadok(trim(F_filename_S)),istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(vardict) read dict fromFile, File not found or not readable: ' &
              //trim(F_filename_S))
         return
      endif
 
      fileid = 0
      istat = fnom(fileid,F_filename_S,'SEQ/FMT+R/O+OLD',0)
      if (.not.RMN_IS_OK(istat) .or. fileid <= 0) then
         call msg(MSG_WARNING,'(vardict) read dict fromFile, Problem opening file: '// &
              trim(F_filename_S))
         return
      endif

      print_L = .false.
      if (present(F_print_L)) print_L = F_print_L

      F_istat = RMN_OK
      DOFILE: do
         read(fileid,'(a)',iostat=istat2) string_S
         if (istat2 /= 0) exit DOFILE
         if (string_S == ' ') cycle DOFILE
         F_istat = min(vardict_add_string(string_S,print_L),F_istat)
      enddo DOFILE
      !---------------------------------------------------------------
      return
   end function vardict_fromfile


   !/@*
   function vardict_from_rdict(F_filename_S,F_pkg_S) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      character(len=*),intent(in),optional :: F_pkg_S
      integer :: F_istat
      !*@/
      logical,parameter :: NOPRINT = .false.
      integer,parameter :: UNIT_I0 = 67
      integer :: istat,istat2,fileid,n
      character(len=WB_MAXNAMELENGTH) :: varname_S,units_S,pkg_S
      character(len=WB_MAXSTRINGLENGTH) :: string_S,desc_S
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      pkg_S = ' '
      if (present(F_pkg_S)) pkg_S = F_pkg_S

      istat = clib_isfile(trim(F_filename_S))
      istat = min(clib_isreadok(trim(F_filename_S)),istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(vardict) read dict fromFile, File not found or not readable: ' &
              //trim(F_filename_S))
         return
      endif
 
      fileid = 0
      istat = fnom(fileid,F_filename_S,'SEQ/FMT+R/O+OLD',0)
      if (.not.RMN_IS_OK(istat) .or. fileid <= 0) then
         call msg(MSG_WARNING,'(vardict) read dict fromFile, Problem opening file: '// &
              trim(F_filename_S))
         return
      endif

      F_istat = RMN_OK
      DOFILE: do
         read(fileid,'(a)',iostat=istat2) string_S
         if (istat2 /= 0) exit DOFILE
         n = len_trim(string_S)
         if (n < UNIT_I0) cycle DOFILE
         varname_S = string_S(1:4)
         desc_S = string_S(5:UNIT_I0)
         n = min(max(UNIT_I0,len_trim(string_S)),UNIT_I0+len(units_S)-1)
         units_S = string_S(UNIT_I0:n)
         F_istat = min(vardict_add(varname_S,pkg_S,' ',' ',units_S,desc_S,F_print_L=NOPRINT),F_istat)
      enddo DOFILE
      !---------------------------------------------------------------
      return
   end function vardict_from_rdict


   !==== Private functions =================================================


   !/@*
   subroutine priv_name(F_name_S,F_pkg_S,F_varname_S,F_namefirst_L)
      implicit none
      character(len=*),intent(out) :: F_name_S
      character(len=*),intent(in) :: F_pkg_S,F_varname_S
      logical,intent(in) :: F_namefirst_L
      !*@/
      character(len=WB_MAXNAMELENGTH) :: prefix_S,part1_S,part2_S
      !---------------------------------------------------------------
      prefix_S = PREFIX
      part1_S = F_pkg_S
      part2_S = F_varname_S
      if (F_namefirst_L) then
         prefix_S = PREFIX2
         part1_S = F_varname_S
         part2_S = F_pkg_S
      endif
      if (part2_S == ' ') then
         if (part1_S == ' ') then
            write(F_name_S,WBFMT0) trim(prefix_S)
         else
            write(F_name_S,WBFMT1) trim(prefix_S),trim(part1_S)
         endif
      else
         write(F_name_S,WBFMT2) trim(prefix_S),trim(part1_S),trim(part2_S)
      endif
      !---------------------------------------------------------------
      return
   end subroutine priv_name


end module vardict_mod
