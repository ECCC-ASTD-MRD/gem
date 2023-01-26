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
module input_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower, clib_toupper
   use fstmpi_mod
   use incfg_mod
   use input_files_mod
   use ezgrid_mod
   use vGrid_Descriptors
   use vgrid_wb
   use hinterp4yy_mod
   use vinterp_mod
   use time_interp_mod
   use mu_jdate_mod, only: jdate_year, jdate_month, jdate_midmonth, MU_JDATE_ANY
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland,2011-04
   !@description
   ! Public functions
   public :: input_new, input_add, input_nbvar, input_meta, input_isvarstep, &
        input_varindex,input_setgridid,input_getgridid,input_close_files, &
        input_set_basedir, input_set_filename, INPUT_FILES_ANAL, &
        INPUT_FILES_CLIM, INPUT_FILES_GEOP

   public :: input_get
   ! Public constants
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface input_new
      module procedure incfg_new_4
      module procedure incfg_new_8
   end interface

   interface input_add
      module procedure incfg_add_string
      module procedure incfg_add_kv
   end interface

   interface input_nbvar
      module procedure incfg_nbvar
   end interface

   interface input_meta
      module procedure incfg_meta
   end interface

   interface input_isvarstep
      module procedure incfg_isvarstep
   end interface

   interface input_varindex
      module procedure incfg_varindex
   end interface

   interface input_getgridid
      module procedure incfg_getgridid
   end interface

   interface input_setgridid
      module procedure incfg_setgridid
   end interface

   interface input_set_filename
      module procedure input_files_set_name0
      module procedure input_files_set_name1
   end interface

   interface input_set_basedir
      module procedure input_files_set_basedir0
      module procedure input_files_set_basedir1
   end interface

   interface input_close_files
      module procedure input_files_close0
      module procedure input_files_close1
   end interface

   interface input_get
      module procedure input_get_scalar
      module procedure input_get_vect
      module procedure input_get3d_scalar
      module procedure input_get3d_vect
   end interface

   integer,parameter :: NMAX_LEVELS = 9999
   integer,parameter :: VTYPE_TDIAG = -2
   integer,parameter :: VTYPE_MDIAG = -10

   type :: input_ptr_T
      real,pointer,dimension(:,:,:) :: p
   end type input_ptr_T

   type :: input_var_T
!!$      sequence
      logical :: needonce_L
      integer :: ni,nj,nk,hgridid,hgrididcore,vtype,dt,hstat
      integer(INT64) :: jdatev,jdateo
      integer,pointer :: ip1list(:)
      type(input_ptr_T) :: d(2)
      character(len=INCFG_STRLEN) :: vn(2),hint_S,vint_S,tint_S,sfc_S,vgrid_S,typvar_S
      type(vgrid_descriptor),pointer :: vgrid
   end type input_var_T

   logical,save :: m_allow_per_level_grid_L = .true.

contains

   !TODO-later: support for other array type and rank

   !TODO-later: input_free...incfg_free

   !/@*
   function input_get_scalar(F_id,F_index,F_step,F_hgridid,F_data,F_varname_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_index,F_step,F_hgridid
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(out),optional :: F_varname_S
      !@return
      integer :: F_istat
      !*@/
      character(len=32) :: varname_S
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input) get_scalar [BEGIN]')
      F_istat = input_get3d_scalar(F_id,F_index,F_step,F_hgridid,' ',F_data,varname_S)
      if (present(F_varname_S)) F_varname_S = varname_S
      call msg(MSG_DEBUG,'(input) get_scalar [END]')
      !----------------------------------------------------------------------
      return
   end function input_get_scalar


   !/@*
   function input_get3d_scalar(F_id,F_index,F_step,F_hgridid,F_vgridid_S,F_data,F_varname_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_index,F_step,F_hgridid
      character(len=*),intent(in) :: F_vgridid_S
      real,pointer :: F_data(:,:,:)
      character(len=*),intent(out),optional :: F_varname_S
      !@return
      integer :: F_istat
      !*@/
      character(len=32) :: varname_S,varname2_S
      real,pointer :: data2(:,:,:)
      integer :: istat
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input) get3d_scalar [BEGIN]')
      nullify(data2)
      F_istat = input_get3d_vect(F_id,F_index,F_step,F_hgridid,F_vgridid_S,F_data,data2,varname_S,varname2_S)
      if (associated(data2)) then
         call msg(MSG_WARNING,'(input) requesting vectorial field as a scalar, 2nd component ignored: ' &
              //trim(varname_S)//' : '//trim(varname2_S))
         deallocate(data2,stat=istat)
         nullify(data2)
      endif
      if (present(F_varname_S)) F_varname_S = varname_S
      call msg(MSG_DEBUG,'(input) get3d_scalar [END]')
      !----------------------------------------------------------------------
      return
   end function input_get3d_scalar


   !/@*
   function input_get_vect(F_id,F_index,F_step,F_hgridid,F_data,F_data2,F_varname_S,F_varname2_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_index,F_step,F_hgridid
      real,pointer :: F_data(:,:,:),F_data2(:,:,:)
      character(len=*),intent(out),optional :: F_varname_S,F_varname2_S
      !@return
      integer :: F_istat
      !*@/
      character(len=32) :: varname_S,varname2_S
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input) get_vect [BEGIN]')
      F_istat = input_get3d_vect(F_id,F_index,F_step,F_hgridid,' ',F_data,F_data2,varname_S,varname2_S)
      if (present(F_varname_S)) F_varname_S = varname_S
      if (present(F_varname2_S)) F_varname2_S = varname2_S
      call msg(MSG_DEBUG,'(input) get_vect [END]')
      !----------------------------------------------------------------------
      return
   end function input_get_vect


   !/@*
   function input_get3d_vect(F_id,F_index,F_step,F_hgridid,F_vgridid_S,F_data, &
        F_data2,F_varname_S,F_varname2_S,F_ovname1_S, F_ovname2_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_id,F_index,F_step,F_hgridid
      character(len=*),intent(in) :: F_vgridid_S
      real,pointer :: F_data(:,:,:),F_data2(:,:,:)
      character(len=*),intent(out),optional :: F_varname_S,F_varname2_S
      character(len=*),intent(in),optional :: F_ovname1_S, F_ovname2_S  !#input name overrides
      !@return
      integer :: F_istat
      !*@/
      character(len=INCFG_STRLEN) :: files_S(INPUT_FILES_NMAX),grtyp_S,lvl_type_S
      character(len=256) :: tmp_S,vn1_S,fn1_S,fn2_S
      integer :: istat,ii,ig1,ig2,ig3,ig4,fileidx,filetype
      integer,target :: ip1list(NMAX_LEVELS)
      logical :: ip1list_alloc_L
      type(vgrid_descriptor),target :: vgrid
      type(input_var_T) :: fld
      !----------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input) get_vect [BEGIN]')
      F_istat = RMN_ERR
      if (F_hgridid < 0) then
         call msg(MSG_WARNING,'(input) invalid gridid')
         return
      endif

      istat = incfg_meta(F_id,F_index,fld%vn(1),fld%vn(2),files_S, &
           fld%hint_S,fld%vint_S,fld%tint_S,F_lvl_type_S=lvl_type_S, &
           F_ip1list=ip1list,F_nip1=fld%nk,F_needonce_L=fld%needonce_L, &
           F_typvar_S=fld%typvar_S)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(input) Problem getting input var metadata')
         return
      endif
      if (present(F_ovname1_S)) then
         if (F_ovname1_S /= '') then
            fld%vn(1) = F_ovname1_S
            istat = clib_tolower(fld%vn(1))
            if (present(F_ovname2_S)) then
               if (F_ovname2_S /= '') then
                  fld%vn(2) = F_ovname2_S
                  istat = clib_tolower(fld%vn(2))
               endif
            endif
         endif
      endif
      if (present(F_varname_S)) F_varname_S = fld%vn(1)
      if (present(F_varname2_S)) F_varname2_S = fld%vn(2)

      ip1list_alloc_L = .false.
      fld%ip1list => ip1list
      fld%vtype = VGRID_GROUND_TYPE
      if (lvl_type_S(1:5) == 'tdiag') then
         fld%vtype = VTYPE_TDIAG
      else if (lvl_type_S(1:5) == 'mdiag') then
         fld%vtype = VTYPE_MDIAG
      endif
      IF_V_INT: if (fld%vint_S /= 'none') then
         IF_VGRID: if (F_vgridid_S == ' ') then
            call msg(MSG_WARNING,'(input) Cannot interpolate, no dest vgrid provided for: '// &
                 trim(fld%vn(1))//' [v_int='//trim(fld%vint_S)//']')
            return
         else
            fld%vgrid_S = F_vgridid_S
            nullify(fld%ip1list)
            ip1list_alloc_L = .true.
            istat = vgrid_wb_get(F_vgridid_S,vgrid,fld%ip1list,fld%vtype,fld%sfc_S)
            if (.not.RMN_IS_OK(istat)) then
               call msg(MSG_WARNING,'(input) Problem getting vgrid meta, cannot get: '// &
                    trim(fld%vn(1))//' [vgrid='//trim(F_vgridid_S)//']')
               return
            endif
            fld%nk = 0
            if (associated(fld%ip1list)) fld%nk = size(fld%ip1list)
            fld%vgrid => vgrid
         endif IF_VGRID
      else
         if (fld%nk == 1 .and. ip1list(1) == -1) then
            if (F_vgridid_S /= ' ') then
               nullify(fld%ip1list)
               ip1list_alloc_L = .true.
               istat = vgrid_wb_get(F_vgridid_S,vgrid,fld%ip1list)
               if (.not.RMN_IS_OK(istat)) then
                  call msg(MSG_WARNING,'(input) Problem getting vgrid meta, cannot get: '//&
                       trim(fld%vn(1))//' [vgrid='//trim(F_vgridid_S)//']')
                  return
               endif
               fld%nk = 0
               if (associated(fld%ip1list)) fld%nk = size(fld%ip1list)
            else
               call msg(MSG_WARNING,'(input) Cannot get data, level list unknown, no dest '// &
                    'vgrid provided for: '//trim(fld%vn(1)))
               return
            endif
         endif
      endif IF_V_INT
      if (fld%nk < 1) then
         call msg(MSG_WARNING,'(input) Problem getting vgrid meta, nk < 1, cannot get: '// &
              trim(fld%vn(1))//' [vgrid='//trim(F_vgridid_S)//']')
         return
      endif

      istat = ezgprm(F_hgridid,grtyp_S,fld%ni,fld%nj,ig1,ig2,ig3,ig4)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(input) Problem getting grid info')
         return
      endif
      fld%hgridid = F_hgridid
      fld%hgrididcore = incfg_getgridid(F_id)
      if (fld%hgrididcore < 0) fld%hgrididcore = fld%hgridid

      istat = incfg_time(F_id,F_step,fld%jdatev,fld%jdateo,fld%dt)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(input) Problem getting time info')
         return
      endif

      nullify(fld%d(1)%p,fld%d(2)%p)
      if (.not.associated(F_data)) then
         allocate(F_data(fld%ni,fld%nj,fld%nk),stat=istat)
         if (istat /= 0 .or. .not.associated(F_data)) then
            call msg(MSG_WARNING,'(input) Problem allocating memory')
            return
         endif
      endif
      if (any(shape(F_data) /= (/fld%ni,fld%nj,fld%nk/))) then
         call msg(MSG_WARNING,'(input) Wrong shape of provided pointer')
         return
      endif
      fld%d(1)%p => F_data
      IF_VN2: if (fld%vn(2) /= ' ') then
         if (.not.associated(F_data2)) then
            allocate(F_data2(fld%ni,fld%nj,fld%nk),stat=istat)
            if (istat /= 0 .or. .not.associated(F_data2)) then
               call msg(MSG_WARNING,'(input) Problem allocating memory')
               return
            endif
         endif
         if (any(shape(F_data2) /= (/fld%ni,fld%nj,fld%nk/))) then
            call msg(MSG_WARNING,'(input) Wrong shape of 2nd provided pointer')
            return
         endif
         fld%d(2)%p => F_data2
      endif IF_VN2

      ii = 0
      DOFILES: do
         ii = ii + 1
         if (ii > size(files_S) .or. RMN_IS_OK(F_istat)) exit DOFILES
         if (files_S(ii) == ' ') cycle DOFILES
         fileidx = input_files_get_idx(F_id,files_S(ii))
         if (.not.RMN_IS_OK(fileidx)) then
            call msg(MSG_INFO,'(input) Ignoring unknown file: '//trim(files_S(ii)))
            cycle
         endif
         vn1_S = fld%vn(1); istat = clib_toupper(vn1_S)
         fn1_S = files_S(ii); istat = clib_toupper(fn1_S)
         if (ii == 1) call msg(MSG_INFO,'(input) Looking for '//trim(vn1_S)//' in file "'// &
              trim(fn1_S)//'" (interpolation h/v/t: '//trim(fld%hint_S)//'/'// &
              trim(fld%vint_S)//'/'//trim(fld%tint_S)//')')
         if (ii > 1) then
            fn2_S = files_S(ii-1); istat = clib_toupper(fn2_S)
            call msg(MSG_INFO,'(input) '//trim(vn1_S)//' not found in "'//trim(fn2_S)// &
                 '"; looking in "'//trim(fn1_S)//'"')
         endif
         filetype = input_files_get_type(F_id,fileidx)
         F_istat = priv_input_data(fld,F_id,fileidx,filetype)
      enddo DOFILES
      if (associated(fld%ip1list) .and. ip1list_alloc_L) then
         if (fld%vint_S /= 'none') istat = vgd_free(fld%vgrid)
         deallocate(fld%ip1list,stat=istat)
      endif
      write(tmp_S, '(i0)') F_istat
      call msg(MSG_DEBUG,'(input) get_vect [END] '//trim(tmp_S))
      !----------------------------------------------------------------------
      return
   end function input_get3d_vect


   !==== Private Functions =================================================


   !/@*
   function priv_input_data(F_fld,F_id,F_file_idx,F_filetype) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(input_var_T),intent(inout) :: F_fld
      integer,intent(in) :: F_id,F_file_idx,F_filetype
      !@return
      integer :: F_istat
      !*@/
      character(len=INCFG_STRLEN) :: tint_S
      character(len=256) :: tmp_S
      integer :: istat,datevfuzz,t_int_type,ip2, &
           find_type,trials, &
           ip2m1,ip2p1,hstat_p,hstat_n,hstat_ps,hstat_ns,&
           lijk(3),uijk(3)
      integer(INT64) :: findjdatev,myjdatev,jdatevm1,jdatevp1,jdatevp,jdatevn
      logical :: inrestart_L
      type(input_var_T) :: fld2
      type(input_ptr_T) :: hdata(2),hdatap(2),hdatan(2)
      real,pointer,dimension(:,:,:) :: sfcfld,sfcfldp,sfcfldn,ptr11,ptr12,ptr21,ptr22
      character(len=32) :: vgrid_S,sfc_S,vgridp_S,vgridn_S
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input) input_data [BEGIN] '//trim(F_fld%vn(1)))
      F_istat = RMN_ERR

      tint_S = F_fld%tint_S
      if (F_filetype == INPUT_FILES_GEOP) tint_S = 'any '
      call priv_clim_date(F_filetype,F_fld%jdatev,findjdatev,jdatevm1,jdatevp1,&
           ip2,ip2m1,ip2p1)
      datevfuzz = 0
      vgrid_S = ' '
      sfc_S = ' '

      nullify(hdata(1)%p,hdata(2)%p,sfcfld)
      IF_TINT: select case(tint_S(1:4))
      case('none')

         F_istat = priv_read_hinterp(F_fld,F_id,F_file_idx,ip2,findjdatev,&
              datevfuzz,FST_FIND_NEAR,hdata,sfcfld,vgrid_S,sfc_S)
         if (F_fld%vint_S /= 'none') &
              F_istat = min(priv_vinterp(F_fld,hdata,sfcfld,vgrid_S),F_istat)

      case('any ')

         ip2 = RMN_ANY_I
         F_istat = priv_read_hinterp(F_fld,F_id,F_file_idx,ip2,MU_JDATE_ANY,&
              datevfuzz,FST_FIND_NEAR,hdata,sfcfld,vgrid_S,sfc_S)
         if (F_fld%vint_S /= 'none') &
              F_istat = min(priv_vinterp(F_fld,hdata,sfcfld,vgrid_S),F_istat)

      case default !# line, near, step, next

         inrestart_L = .not.F_fld%needonce_L
         t_int_type = time_interp_typecode(tint_S)
         datevfuzz = huge(datevfuzz)
         find_type = FST_FIND_LE
         trials = 0
         DO_TRIALS: do

            F_istat = time_interp_status(F_fld%vn(1),F_fld%jdatev,t_int_type)
            if (F_istat == TIME_INTERP_NOT_FOUND &
                 .or. F_istat == TIME_INTERP_NEED_PREV) then
               find_type = FST_FIND_LE
               if (tint_S(1:4) == 'next') find_type = FST_FIND_GE
            elseif (F_istat == TIME_INTERP_NEED_NEXT) then
               find_type = FST_FIND_GE
!!$               if (any(tint_S(1:4) == (/'next','step'/))) find_type = FST_FIND_GT
               if (trials >1) find_type = FST_FIND_GT
            else
               exit DO_TRIALS
            endif

            nullify(hdata(1)%p,hdata(2)%p,sfcfld)
            if (F_filetype == INPUT_FILES_CLIM) then
               istat = priv_read_hinterp_clim(F_fld,F_id,F_file_idx, &
                    find_type,hdata,sfcfld,vgrid_S,sfc_S,myjdatev)
            else
               istat = priv_read_hinterp(F_fld,F_id,F_file_idx,ip2,findjdatev,&
                    datevfuzz,find_type,hdata,sfcfld,vgrid_S,sfc_S,myjdatev)
            endif

            if (istat /= RMN_ERR .and. associated(hdata(1)%p)) then
               if (associated(sfcfld) .and. sfc_S /= ' ') then
                  sfc_S = 'rfld/'//sfc_S
                  istat = time_interp_set(sfcfld,sfc_S,myjdatev,inrestart_L,&
                       F_flag=F_fld%hstat)
               endif
               istat = time_interp_set(hdata(1)%p,F_fld%vn(1),myjdatev, &
                    vgrid_S,sfc_S,inrestart_L,F_flag=F_fld%hstat)
               if (associated(hdata(2)%p) .and. F_fld%vn(2) /= ' ') &
                    istat = time_interp_set(hdata(2)%p,F_fld%vn(2),myjdatev, &
                         vgrid_S,sfc_S,inrestart_L,F_flag=F_fld%hstat)
            endif

            trials = trials + 1
            if (trials > 2) exit

         enddo DO_TRIALS

         IF_VINT: if (F_fld%vint_S == 'none') then

            F_istat = time_interp_get(F_fld%d(1)%p,F_fld%vn(1), &
                 F_fld%jdatev,t_int_type,F_fld%dt)
            if (F_fld%vn(2) /= ' ') &
                 F_istat = min(time_interp_get(F_fld%d(2)%p,F_fld%vn(2), &
                     F_fld%jdatev,t_int_type,F_fld%dt),F_istat)

         else !IF_VINT

            F_istat = RMN_ERR
            nullify(hdatap(1)%p,hdatap(2)%p,hdatan(1)%p,hdatan(2)%p,sfcfldp,sfcfldn)
            istat = time_interp_retrieve(F_fld%vn(1),TIME_INTERP_PREV, &
                 jdatevp,F_data=hdatap(1)%p,F_vgrid_S=vgridp_S, &
                 F_sfcfld_S=sfc_S,F_flag=hstat_p)
            if (sfc_S /= '') &
                 istat = time_interp_retrieve(sfc_S,TIME_INTERP_PREV, &
                 myjdatev,F_data=sfcfldp,F_flag=hstat_ps)
            istat = time_interp_retrieve(F_fld%vn(1),TIME_INTERP_NEXT, &
                 jdatevn,F_data=hdatan(1)%p,F_vgrid_S=vgridn_S, &
                 F_sfcfld_S=sfc_S,F_flag=hstat_n)
            if (sfc_S /= '') &
                 istat = time_interp_retrieve(sfc_S,TIME_INTERP_NEXT, &
                 myjdatev,F_data=sfcfldn,F_flag=hstat_ns)
            if (F_fld%vn(2) /= ' ') then
               istat = time_interp_retrieve(F_fld%vn(2),TIME_INTERP_PREV, &
                    myjdatev,F_data=hdatap(2)%p)
               istat = time_interp_retrieve(F_fld%vn(2),TIME_INTERP_NEXT, &
                    myjdatev,F_data=hdatan(2)%p)
            endif

            istat = time_interp_status(F_fld%vn(1),F_fld%jdatev,t_int_type)
            TI_STAT: if (istat == 0) then
               F_istat = priv_vinterp(F_fld,hdatap,sfcfldp,vgridp_S)
            else if (istat == nint(TIME_INTERP_WEIGHT_FACT)) then
               F_istat = priv_vinterp(F_fld,hdatan,sfcfldn,vgridn_S)
            else if (istat > 0 .and. istat < nint(TIME_INTERP_WEIGHT_FACT)) then
               lijk(:) = lbound(F_fld%d(1)%p) ; uijk = ubound(F_fld%d(1)%p)
               fld2 = F_fld
               nullify(fld2%d(1)%p,fld2%d(2)%p,fld2%ip1list)
               !TODO: what about vgrid, ip1list?
               allocate(fld2%d(1)%p(lijk(1):uijk(1),lijk(2):uijk(2),lijk(3):uijk(3)))
               if (F_fld%vn(2) /= ' ') &
                    allocate(fld2%d(2)%p(lijk(1):uijk(1),lijk(2):uijk(2),lijk(3):uijk(3)))
               !# vint: dp,vp,sp => dpv,v0,sp
               !# vint: dn,vn,sn => dnv,v0,sn
               !# tint: dpv,dnv => d0v

               nullify(ptr11,ptr12,ptr21,ptr22)
               F_istat = priv_vinterp(F_fld,hdatap,sfcfldp,vgridp_S)!,.true.)
               if (RMN_IS_OK(F_istat)) then
                  ptr11 => F_fld%d(1)%p
                  ptr12 => F_fld%d(2)%p
               endif
               F_istat = priv_vinterp(fld2,hdatan,sfcfldn,vgridn_S)!,.true.)
               if (RMN_IS_OK(F_istat)) then
                  ptr21 => fld2%d(1)%p
                  ptr22 => fld2%d(2)%p
               endif
               F_istat = time_interp_get(F_fld%d(1)%p, &
                    ptr11,ptr21,F_fld%jdatev,jdatevp,jdatevn,&
                    t_int_type,F_varname_S=F_fld%vn(1))
               if (F_fld%vn(2) /= ' ' .and. RMN_IS_OK(F_istat)) &
                    F_istat = time_interp_get(F_fld%d(2)%p,&
                         ptr12,ptr22,F_fld%jdatev,jdatevp,jdatevn,&
                         t_int_type,F_varname_S=F_fld%vn(2))
               deallocate(fld2%d(1)%p,fld2%d(2)%p,stat=istat)
            else !TI_STAT
               call msg(MSG_ERROR,'(input) Problem getting time-int data; for '// &
                    trim(F_fld%vn(1))//' (interpolation h/v/t: '//trim(F_fld%hint_S)// &
                    '/'//trim(F_fld%vint_S)//'/'//trim(tint_S)//')')
               F_istat = RMN_ERR
            endif TI_STAT

         endif IF_VINT

         if (RMN_IS_OK(F_istat)) then
            istat = time_interp_status(F_fld%vn(1),F_fld%jdatev,t_int_type)
            if (istat > 0 .and. istat < nint(TIME_INTERP_WEIGHT_FACT)) then
               write(tmp_S,'(i6)') (istat*100)/nint(TIME_INTERP_WEIGHT_FACT)
               call msg(MSG_INFO,'(input) Got time interpolated value for '// &
                    trim(F_fld%vn(1))//' (interpolation h/v/t: '//trim(F_fld%hint_S)// &
                    '/'//trim(F_fld%vint_S)//'/'//trim(tint_S)//') (t_int_weight='//trim(tmp_S)//'%)')
            endif
         else
            call msg(MSG_WARNING,'(input) Problem getting: '//trim(F_fld%vn(1))// &
                 ' (interpolation h/v/t: '//trim(F_fld%hint_S)//'/'//trim(F_fld%vint_S)// &
                 '/'//trim(tint_S)//')')
         endif

      end select IF_TINT

      write(tmp_S, '(i0)') F_istat
      call msg(MSG_DEBUG,'(input) input_data [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function priv_input_data


    !/@*
   subroutine priv_clim_date(F_filetype,F_jdatev0,F_jdatev,F_jdatevm1,F_jdatevp1,F_ip2,F_ip2m1,F_ip2p1)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_filetype
      integer(INT64),intent(in) :: F_jdatev0
      integer(INT64),intent(out) :: F_jdatev,F_jdatevm1,F_jdatevp1
      integer,intent(out) :: F_ip2,F_ip2m1,F_ip2p1
      !*@/
      integer :: year
      !------------------------------------------------------------------
      F_jdatev = F_jdatev0
      F_jdatevm1 = F_jdatev0
      F_jdatevp1 = F_jdatev0
      F_ip2 = RMN_ANY_I ; F_ip2m1 = RMN_ANY_I ; F_ip2p1 = RMN_ANY_I
      if (F_filetype /= INPUT_FILES_CLIM) return

      F_jdatev = RMN_ANY_DATE
      F_ip2 = jdate_month(F_jdatev0)
      F_ip2m1 = F_ip2
      year = jdate_year(F_jdatev0)
      F_jdatevm1 = jdate_midmonth(year,F_ip2m1)
      if (F_jdatevm1 > F_jdatev0) then
         F_ip2m1 = F_ip2 - 1
         if (F_ip2m1 < 1) then
            F_ip2m1 = F_ip2m1 + 12
            year = year - 1
         endif
         F_jdatevm1 = jdate_midmonth(year,F_ip2m1)
      endif
      F_ip2p1 = F_ip2m1 + 1
      if (F_ip2p1 > 12) then
         F_ip2p1 = F_ip2p1 - 12
         year = year + 1
      endif
      F_jdatevp1 = jdate_midmonth(year,F_ip2p1)
      !------------------------------------------------------------------
      return
   end subroutine priv_clim_date


   !/@*
   function priv_read_hinterp_clim(F_fld,F_id,F_file_idx,F_fuzztype,F_hdata, &
        F_sfcfld,F_vgrid_S,F_sfc_S,F_foundjdatev) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(input_var_T),intent(inout) :: F_fld
      integer,intent(in) :: F_id,F_file_idx,F_fuzztype
      type(input_ptr_T) :: F_hdata(2)
      real,pointer :: F_sfcfld(:,:,:)
      character(len=*),intent(out) :: F_vgrid_S,F_sfc_S
      integer(INT64),intent(out) :: F_foundjdatev
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: tmp_S
      integer(INT64) :: findjdatev,jdatevm1,jdatevp1
      integer :: ip2,ip2m1,ip2p1
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input) read_interp_clim [BEGIN]')
      call priv_clim_date(INPUT_FILES_CLIM,F_fld%jdatev,findjdatev,jdatevm1,jdatevp1,ip2,ip2m1,ip2p1)
      if (F_fuzztype == FST_FIND_LE) then
         findjdatev = jdatevm1
         ip2 = ip2m1
      else
         findjdatev = jdatevp1
         ip2 = ip2p1
      endif
      F_istat = priv_read_hinterp(F_fld,F_id,F_file_idx,ip2,MU_JDATE_ANY,0, &
           FST_FIND_NEAR,F_hdata,F_sfcfld,F_vgrid_S,F_sfc_S,F_foundjdatev)
      F_foundjdatev = 0
      if (F_istat /= RMN_ERR) F_foundjdatev = findjdatev
      write(tmp_S, '(i0,i0)') F_istat,F_foundjdatev
      call msg(MSG_DEBUG,'(input) read_interp_clim [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function priv_read_hinterp_clim


   !/@*
   function priv_read_hinterp(F_fld,F_id,F_file_idx,F_ip2,F_jdatev,F_datevfuzz, &
        F_fuzztype,F_hdata,F_sfcfld,F_vgrid_S,F_sfc_S,F_foundjdatev) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(input_var_T),intent(inout) :: F_fld
      integer,intent(in) :: F_id,F_file_idx,F_ip2,F_datevfuzz,F_fuzztype
      integer(INT64),intent(in) :: F_jdatev
      type(input_ptr_T) :: F_hdata(2)
      real,pointer :: F_sfcfld(:,:,:)
      character(len=*),intent(out) :: F_vgrid_S,F_sfc_S
      integer(INT64),intent(out),optional :: F_foundjdatev
      !@return
      integer :: F_istat
      !*@/
      integer,parameter :: MYMSG_QUIET = 99
      character(len=5),parameter :: HINT_TYPE_S = 'CUBIC'
      character(len=INCFG_STRLEN) :: vgrid_S
      character(len=256) :: tmp_S, dummy_S
      integer :: datevfuzz,fuzztype,ingridid,ip1,ip2,istat,istat2, &
           k,nk,vtype,lij(3),uij(3),vb0, ikind
      real :: zp1
      integer(INT64) :: jdatev2
      integer,target :: ip1list0(1)
      integer,pointer :: ip1list(:)
      logical :: ip1list_alloc_L
      type(input_ptr_T) :: indata(2)
      type(vgrid_descriptor) :: vgrid
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input) read_interp [BEGIN]')
      F_istat = RMN_ERR
      if (present(F_foundjdatev)) F_foundjdatev = RMN_ERR
      F_fld%hstat = HINTERP4YY_NONE
      if (F_fld%nk > size(F_fld%d(1)%p,3)) then
         call msg(MSG_WARNING,'(input) Too many requested levels for allocated space')
         return
      endif

      jdatev2 = F_jdatev ; datevfuzz = F_datevfuzz ; fuzztype = F_fuzztype
      ingridid = -1
      nk = F_fld%nk ; ip1list => F_fld%ip1list
      F_vgrid_S = ' ' ; F_sfc_S = ' '
      F_hdata(1)%p => F_fld%d(1)%p
      if (associated(F_fld%d(2)%p)) F_hdata(2)%p => F_fld%d(2)%p
      nullify(F_sfcfld,indata(1)%p,indata(2)%p)

      ip1list_alloc_L = .false.
      IF_VINT: if (F_fld%vint_S /= 'none') then
         istat = input_files_vgrid(F_id,F_file_idx,F_fld%vn(1),RMN_ANY_I, &
              RMN_ANY_DATE,F_vgrid_S)
         nullify(ip1list)
         ip1list_alloc_L = .true.
         istat = min(vgrid_wb_get(F_vgrid_S,vgrid,ip1list,vtype,F_sfc_S),istat)
         call collect_error(istat)
         if (F_vgrid_S == '' .or. .not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING,'(input) Problem getting vgrid for '//trim(F_fld%vn(1))//' '//trim(F_fld%vn(2)))
            return
         endif
         nk = size(ip1list) ; lij(:) = lbound(F_fld%d(1)%p) ; uij = ubound(F_fld%d(1)%p)
         allocate(F_hdata(1)%p(lij(1):uij(1),lij(2):uij(2),nk))!,stat=istat2)
         if (associated(F_fld%d(2)%p)) &
              allocate(F_hdata(2)%p(lij(1):uij(1),lij(2):uij(2),nk))!,stat=istat2)
      endif IF_VINT

      IF_DIAG: if (any(F_fld%vtype == (/VTYPE_TDIAG,VTYPE_MDIAG/))) then
         istat = input_files_vgrid(F_id,F_file_idx,F_fld%vn(1),RMN_ANY_I, &
              RMN_ANY_DATE,vgrid_S)
         istat = min(vgrid_wb_get(vgrid_S,vgrid),istat)
         ip1list0(1) = -1
         call msg_verbosity_get(vb0)
         call msg_verbosity(MYMSG_QUIET)
         if (F_fld%vtype == VTYPE_TDIAG) then
            istat = min(vgd_get(vgrid,'DIPT',ip1list0(1)),istat)
         else
            istat = min(vgd_get(vgrid,'DIPM',ip1list0(1)),istat)
         endif
         call msg_verbosity(vb0)
         if (RMN_IS_OK(istat) .and. RMN_IS_OK(ip1list0(1))) then
            nk = 1
            if (associated(ip1list) .and. ip1list_alloc_L) &
                 deallocate(ip1list,stat=istat)
            ip1list => ip1list0
         else
            call msg(MSG_INFOPLUS,'(input) Problem getting vgrid diag for '// &
                 trim(F_fld%vn(1))//' '//trim(F_fld%vn(2)))
            nk = 1
            zp1 = 1.
            ikind = RMN_CONV_HY
            dummy_S = ' '
            call convip_plus(ip1, zp1, ikind, RMN_CONV_P2IPNEW, dummy_S, .not.RMN_CONV_USEFORMAT_L)
            ip1list0(1) = ip1
            if (associated(ip1list) .and. ip1list_alloc_L) &
                 deallocate(ip1list,stat=istat)
            ip1list => ip1list0
         endif
      endif IF_DIAG

      DO_LEVELS: do k=1,nk
         F_istat = RMN_ERR
         ip1 = ip1list(k) ; ip2 = F_ip2
         istat = input_files_read(F_id,F_file_idx,F_fld%vn(1),ip1,ip2, &
              jdatev2,datevfuzz,fuzztype,ingridid,indata(1)%p, &
              F_typvar_S=F_fld%typvar_S)
         datevfuzz = 0 ; fuzztype = FST_FIND_NEAR
         if (F_fld%vn(2) /= ' ') &
              istat = min(input_files_read(F_id,F_file_idx,F_fld%vn(2),ip1,ip2,&
                   jdatev2,datevfuzz,fuzztype,ingridid,indata(2)%p, &
                   F_typvar_S=F_fld%typvar_S),istat)
         call collect_error(istat)
         if (.not.RMN_IS_OK(istat)) exit DO_LEVELS
         istat = hinterp4yy2d(F_hdata(1)%p,F_hdata(2)%p, &
              indata(1)%p,indata(2)%p,k,ingridid, &
              F_fld%hgridid,F_fld%hgrididcore,F_fld%hint_S, &
              F_fld%vn(1),F_fld%vn(2))
         call collect_error(istat)
         F_fld%hstat = min(istat,F_fld%hstat)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_WARNING,'(input) Problem interpolating '//trim(F_fld%vn(1))//' '//trim(F_fld%vn(2)))
            exit DO_LEVELS
         endif
         if (m_allow_per_level_grid_L .or. k==nk) then
            ingridid = -1
            if (associated(indata(1)%p)) deallocate(indata(1)%p,stat=istat2)
            if (associated(indata(2)%p)) deallocate(indata(2)%p,stat=istat2)
            nullify(indata(1)%p,indata(2)%p)
         endif
         F_istat = RMN_OK  !#max(0,datev2)
         if (present(F_foundjdatev)) F_foundjdatev = jdatev2
      enddo DO_LEVELS

      IF_RFLD: if (F_fld%vint_S /= 'none' .and. F_sfc_S /= ' ') then
         istat = input_files_sfcref(F_id,F_file_idx,jdatev2,vgrid, &
              indata(1)%p,ingridid)
         if (associated(indata(1)%p)) then
            allocate(F_sfcfld(lij(1):uij(1),lij(2):uij(2),1))!,stat=istat2)
            istat = hinterp4yy2d(F_sfcfld,indata(2)%p,indata(1)%p,indata(2)%p, &
                 1,ingridid,F_fld%hgridid,F_fld%hgrididcore,'cubic','vRFLD',' ')
            if (.not.RMN_IS_OK(istat)) then
               F_istat = RMN_ERR
               call msg(MSG_WARNING,'(input) Problem interpolating sfc ref field for '// &
                    trim(F_fld%vn(1))//' '//trim(F_fld%vn(2)))
               if (associated(F_sfcfld)) deallocate(F_sfcfld,stat=istat2)
               nullify(F_sfcfld)
            endif
            deallocate(indata(1)%p,stat=istat)
            nullify(indata(1)%p)
         endif
         if (.not.RMN_IS_OK(istat)) F_istat = RMN_ERR
      endif IF_RFLD

      if (associated(ip1list) .and. ip1list_alloc_L) &
           deallocate(ip1list,stat=istat)

      write(tmp_S, '(i0)') F_istat
      call msg(MSG_DEBUG,'(input) read_interp [END] '//trim(tmp_S))
      !----------------------------------------------------------------------
      return
   end function priv_read_hinterp


   !/@*
   function priv_vinterp(F_fld,F_hdata,F_sfcfld,F_vgrid_S,F_same_sfc_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      type(input_var_T),intent(inout) :: F_fld
      type(input_ptr_T) :: F_hdata(2)
      real,pointer :: F_sfcfld(:,:,:)
      character(len=*),intent(in) :: F_vgrid_S
      logical,intent(in),optional :: F_same_sfc_L
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: tmp_S
      logical :: same_sfc_L
      integer :: nlinbot, ivar, istat
      real,pointer :: sfcfld(:,:),sfcfldout(:,:),sfcfldls(:,:),sfcfldoutls(:,:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(input) vinterp [BEGIN]')
      F_istat = RMN_ERR
      if (.not.(associated(F_fld%d(1)%p).and.associated(F_hdata(1)%p))) return
      if (F_fld%vn(2) /= ' '.and. &
           .not.(associated(F_fld%d(2)%p).and.associated(F_hdata(2)%p))) return

      F_istat = RMN_OK
      IF_VINT: if (F_fld%vint_S == 'none') then

         do ivar = 1,2
            if (F_fld%vn(ivar) /= ' ' .and. associated(F_fld%d(ivar)%p) &
                 .and. associated(F_hdata(ivar)%p)) then
               if (all(shape(F_fld%d(ivar)%p) == shape(F_hdata(ivar)%p))) then
                  F_fld%d(ivar)%p = F_hdata(ivar)%p
               else
                  F_istat = RMN_ERR
               endif
            endif
         enddo

      else !IF_VINT

         nlinbot = 0
         same_sfc_L = .false.
         if (any(F_fld%vint_S(1:4) == (/'l-co','line'/))) &
              nlinbot = size(F_fld%d(1)%p,3)
         if (any(F_fld%vint_S(1:4) == (/'l-co','c-co'/)) .and. &
              F_fld%hstat == HINTERP4YY_NONE) same_sfc_L = .true.
         nullify(sfcfld)
         if (associated(F_sfcfld)) sfcfld => F_sfcfld(:,:,1)

         nullify(sfcfldout, sfcfldls, sfcfldoutls)
         if (present(F_same_sfc_L)) then
            if (F_same_sfc_L) then
               if (associated(F_sfcfld)) sfcfldout => F_sfcfld(:,:,1)
               !#TODO: if (associated(F_sfcfldls)) sfcfldoutls => F_sfcfldls(:,:,1)
            endif
         endif

         do ivar = 1,2
            if (F_fld%vn(ivar) /= ' ' .and. associated(F_hdata(ivar)%p)) then
               if (associated(sfcfldout)) then
                  istat = vinterp( &
                       F_fld%d(ivar)%p, F_fld%vgrid_S,  &
                       F_hdata(ivar)%p, F_vgrid_S, &
                       sfcfldout, sfcfld,  &
                       nlinbot, F_fld%vn(ivar), same_sfc_L, &
                       sfcfldoutls, sfcfldls)
                  F_istat = min(istat, F_istat)
               else
                  istat = vinterp( &
                       F_fld%d(ivar)%p, F_fld%vgrid_S,  &
                       F_hdata(ivar)%p, F_vgrid_S,  &
                       F_sfcfldin=sfcfld,  &
                       F_nlinbot=nlinbot, F_msg_S=F_fld%vn(ivar),  &
                       F_use_same_sfcfld_L=same_sfc_L, &
                       F_sfcfldin2=sfcfldls)
                  F_istat = min(istat, F_istat)
               endif
            endif
         enddo

      endif IF_VINT

      write(tmp_S, '(i0)') F_istat
      call msg(MSG_DEBUG,'(input) vinterp [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function priv_vinterp


end module input_mod

