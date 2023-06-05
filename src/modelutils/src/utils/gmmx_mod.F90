!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer, 
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms 
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer 
!version 3 or (at your option) any later version that should be found at: 
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html 
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software; 
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec), 
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------


!/@
module gmmx_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower
   use wb_itf_mod
   !use rmn_gmm
   use vardict_mod
   use str_mod
   use hgrid_wb
   use vgrid_wb
   use vGrid_Descriptors
   use rmn_gmm
   implicit none
   private
   !@objective GMM/Buses data/metadata managment tool
   !@author Stephane Chamberland, 2012-01
   !@description
   ! Public functions
   public :: gmmx_new, gmmx_meta, gmmx_data, gmmx_list
!!$, gmmx_print_list
   ! Public constants
   !
!@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface gmmx_new
      module procedure gmmx_new0
      module procedure gmmx_new1
   end interface

   interface gmmx_data
      module procedure gmmx_get_data_r4_2d
      module procedure gmmx_get_data_r4_3d
      module procedure gmmx_get_data_r4_bus
   end interface

   interface gmmx_list
      module procedure vardict_getlist
   end interface

   integer,parameter :: IDX_NAME = 1
   integer,parameter :: IDX_HGRID = 2
   integer,parameter :: IDX_VGRID = 3
   integer,parameter :: IDX_VB = 4
   integer,parameter :: NIDX = 4
   integer,parameter :: NIDX_MIN = 3
   integer,parameter :: NIDX_EXTRA = 32 - NIDX
   integer,parameter :: NIDX_MAX = NIDX + NIDX_EXTRA
   integer,parameter :: KEY = 1
   integer,parameter :: VAL = 2
   character(len=WB_MAXNAMELENGTH),parameter :: KNOWN_KEYS(NIDX) = (/&
        'vn     ', & !# mandatory
        'hgrid  ', & !# mandatory
        'vgrid  ', & !# mandatory
        'vb     ' &
        /)

contains


   !/@*
   function gmmx_new1(F_fullmeta_S,F_flags) result(F_istat)
      implicit none
      !@objective Create a new var
      !@arguments
      character(len=*),intent(in) :: F_fullmeta_S !# metadata
      integer,intent(in),optional :: F_flags      !# gmm flags
      !@return
      integer :: F_istat !- exit status or bus index
      !@author  S. Chamberland, 2012-01
      !*@/
      integer,external :: str_split2keyval
      character(len=WB_MAXSTRINGLENGTH) :: meta_S,kv_S(2,NIDX_MAX)
      integer :: n,flags,nkeys
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_fullmeta_S == ' ') then
         call msg(MSG_WARNING,'(gmmx) New Var from string, no Meta')
         return
      endif
      flags = GMM_NULL_FLAGS
      if (present(F_flags)) flags = F_flags
      kv_S = ' '
      kv_S(KEY,1:NIDX) = KNOWN_KEYS(1:NIDX)
      nkeys = str_split2keyval(kv_S,F_fullmeta_S,NIDX_MAX)
      if (nkeys < NIDX_MIN .or. any(kv_S(VAL,1:NIDX_MIN) == ' ')) then
         call msg(MSG_ERROR,'(gmmx) newvar - Mandatory params not all provided for: '// &
              trim(F_fullmeta_S))
         return
      endif
      meta_S = ' '
      do n=1,nkeys
         call str_normalize(kv_S(KEY,n))
         if (any(kv_S(KEY,n) == KNOWN_KEYS)) cycle
         meta_S = trim(meta_S)//trim(kv_S(KEY,n))//'="'//trim(kv_S(VAL,n))//'";'
      enddo
      F_istat = gmmx_new(kv_S(VAL,IDX_NAME),kv_S(VAL,IDX_HGRID),kv_S(VAL,IDX_VGRID), &
           kv_S(VAL,IDX_VB),meta_S,flags)
      !---------------------------------------------------------------------
      return
   end function gmmx_new1


   !/@*
   function gmmx_new0(F_name_S,F_hgrid_S,F_vgrid_S,F_vb_S,F_meta_S,F_flags) result(F_istat)
      implicit none
      !@objective Create a new var
      !@arguments
      character(len=*),intent(in) :: F_name_S !# Key (internal var name)
      character(len=*),intent(in) :: F_hgrid_S,F_vgrid_S !# horiz/verti grid desc
      character(len=*),intent(in),optional :: F_vb_S   !# gmm or bus name
      character(len=*),intent(in),optional :: F_meta_S !# vardict type of metadata
      integer,intent(in),optional :: F_flags  !# gmm flags
      !@return
      integer :: F_istat !- exit status or bus index
      !@author  S. Chamberland, 2012-01
      !*@/
      character(len=WB_MAXNAMELENGTH) :: name_S,vb_S
      character(len=WB_MAXSTRINGLENGTH) :: meta_S
      integer :: istat,flags
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (len_trim(F_name_S) == 0) then
         call msg(MSG_ERROR,'(gmmx) new - Need to provide an internal name')
         return
      endif
      if (len_trim(F_hgrid_S) == 0) then
         call msg(MSG_ERROR,'(gmmx) new - Need to provide an horizontal description '// &
              'hgrid, ignoring: '//trim(F_name_S))
         return
      endif
      if (len_trim(F_vgrid_S) == 0) then
         call msg(MSG_ERROR,'(gmmx) new - Need to provide a vertical description '// &
              'vgrid, ignoring: '//trim(F_name_S))
         return
      endif
      vb_S = 'g'
      if (present(F_vb_S)) then
         if (F_vb_S /= ' ') vb_S = F_vb_S
      endif
      name_S = F_name_S
      istat = vardict_get(name_S,vb_S)
      if (RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(gmmx) newvar - Var already exists: '//trim(F_name_S)// &
              ' (VB='//trim(vb_S)//')')
         return
      endif

      call str_normalize(vb_S)
      flags = GMM_NULL_FLAGS
      if (present(F_flags)) flags = F_flags
      if (vb_S == 'g') then
         F_istat = priv_newvar_gmm(F_name_S,F_hgrid_S,F_vgrid_S,flags)
      else
         call msg(MSG_ERROR,'(gmmx) newvar - Bus type not yet supported, ignoring: '//trim(F_name_S))
         F_istat = RMN_ERR
!!$         F_istat = priv_newvar_bus(F_name_S,F_hgrid_S,F_vgrid_S,flags,vb_S)
      endif
      if (.not.RMN_IS_OK(F_istat)) return

      meta_S = 'VN='//trim(F_name_S)//';'//&
           'VB='//trim(vb_S)//';'//&
           'hgrid='//trim(F_hgrid_S)//';'//&
           'vgrid='//trim(F_vgrid_S)//';'
      if (present(F_meta_S)) meta_S = trim(meta_S)//trim(F_meta_S)
      F_istat = vardict_add_string(meta_S)

!!$      mylen = min(GMM_MAXNAMELENGTH,WB_MAXNAMELENGTH) - len_trim(PREFIX_S)
!!$      call msg(MSG_INFO,'(gmmx) VN='//F_name_S(1:mylen)//'; IN='// &
!!$           kv_S(VAL,KEY_IN)(1:MAX_IONAME_LEN)//'; ON='//kv_S(VAL,KEY_ON)(1:MAX_IONAME_LEN)// &
!!$           '; HG='//trim(kv_S(VAL,KEY_HGRID))//'; VG='//trim(kv_S(VAL,KEY_VGRID))// &
!!$           '; UNITS='//trim(kv_S(VAL,KEY_UNITS)))
      !---------------------------------------------------------------------
      return
   end function gmmx_new0


   !/@*
   function gmmx_meta(F_name_S,F_namein_S,F_nameout_S,F_vb_S,F_hgrid_S, &
        F_vgrid_S,F_units_S,F_desc_S,F_idx0) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      character(len=*),intent(inout) :: F_name_S !- Key (internal var name)
      character(len=*),intent(inout),optional :: F_namein_S,F_nameout_S,F_vb_S
      character(len=*),intent(out),optional :: F_hgrid_S,F_vgrid_S,F_units_S,F_desc_S
      integer,intent(inout),optional :: F_idx0 !- Search from provided idx, return idx of item found
      !@return
      integer :: F_istat !- exit status or bus index
      !@author  S. Chamberland, 2012-01
      !*@/
      character(len=WB_MAXNAMELENGTH) :: namein_S,nameout_S,vb_S,extra_S(NIDX_EXTRA),units_S,name_S
      character(len=WB_MAXSTRINGLENGTH) :: desc_S
      integer :: idx0,n,n1,n2
      !---------------------------------------------------------------------
      namein_S = ' '
      nameout_S = ' '
      vb_S = ' '
      idx0 = 1
      if (present(F_namein_S)) namein_S = F_namein_S
      if (present(F_nameout_S)) nameout_S = F_nameout_S
      if (present(F_vb_S)) vb_S = F_vb_S
      if (present(F_idx0)) idx0 = max(1,F_idx0)

      F_istat = vardict_get(F_name_S,vb_S,namein_S,nameout_S,units_S,desc_S,extra_S,idx0)

      if (RMN_IS_OK(F_istat)) then
         if (present(F_vb_S)) F_vb_S = vb_S
         if (present(F_namein_S)) F_namein_S = namein_S
         if (present(F_nameout_S)) F_nameout_S = nameout_S
         if (present(F_units_S)) F_units_S = units_S
         if (present(F_desc_S)) F_desc_S = desc_S
         if (present(F_idx0)) F_idx0 = idx0
         if (present(F_hgrid_S)) then
            name_S = 'hgrid='
            n1 = len_trim(name_S)
            do n=1,size(extra_S)
               if (extra_S(n)(1:n1) == name_S) then
                  n2 = max(n1+1,min(len_trim(extra_S(n)),len(F_hgrid_S)))
                  F_hgrid_S = extra_S(n)(n1+1:n2)
                  exit
               endif
            enddo
         endif
         if (present(F_vgrid_S)) then
            name_S = 'vgrid='
            n1 = len_trim(name_S)
            do n=1,size(extra_S)
               if (extra_S(n)(1:n1) == name_S) then
                  n2 = max(n1+1,min(len_trim(extra_S(n)),len(F_vgrid_S)))
                  F_vgrid_S = extra_S(n)(n1+1:n2)
                  exit
               endif
            enddo
         endif
      endif

      !TODO-later: return GMM meta/flags (or passed in flags)

!!$      if (.not.RMN_IS_OK(F_istat) .and. present(F_idx0)) F_idx0 = m_nbvar
      !---------------------------------------------------------------------
      return
   end function gmmx_meta


!!$   !/@*
!!$   subroutine gmmx_print_list()
!!$      implicit none
!!$      !*@/
!!$      integer,parameter :: NMAX_VARS = 4096
!!$      character(len=WB_MAXNAMELENGTH) :: name_S,varlist_S(NMAX_VARS)
!!$      character(len=128) :: msg_S
!!$      integer :: istat,ivar,nvar
!!$      real,pointer :: ptr2d(:,:),ptr3d(:,:,:)
!!$      type(gmm_metadata) :: mymeta
!!$      !---------------------------------------------------------------------
!!$      nvar = gmmx_get_list(varlist_S)
!!$      do ivar = 1,nvar
!!$         nullify(ptr2d,ptr3d)
!!$         name_S = varlist_S(ivar)
!!$         istat = gmm_getmeta(name_S,mymeta)
!!$         if (istat >= 0) then
!!$            if (mymeta%l(3)%n == 0) then
!!$               istat = gmm_get(name_S,ptr2d,mymeta)
!!$               if (RMN_IS_OK(istat).and.associated(ptr2d)) then
!!$                  write(msg_S,'(a,2i7,a,2i7)') '(gmmx) '//trim(name_S)//' size 2d:',lbound(ptr2d),';',ubound(ptr2d)
!!$               endif
!!$            else
!!$               istat = gmm_get(name_S,ptr3d,mymeta)
!!$               if (RMN_IS_OK(istat).and.associated(ptr3d)) then
!!$                  write(msg_S,'(a,3i7,a,3i7)') '(gmmx) '//trim(name_S)//' size 3d:',lbound(ptr3d),';',ubound(ptr3d)
!!$               endif
!!$            endif
!!$         else
!!$            write(msg_S,'(a)') '(gmmx) '//trim(name_S)//' : GMM VAR Not Found'
!!$         endif
!!$         call msg(MSG_INFO,msg_S)
!!$      enddo
!!$      !---------------------------------------------------------------------
!!$      return
!!$   end subroutine gmmx_print_list


   !TODO: gmmx_get_data should require vb?

   !/@*
   function gmmx_get_data_r4_3d(F_name_S,F_ptr3d,F_print_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      character(len=*),intent(in) :: F_name_S !- Key (internal var name)
      real,pointer :: F_ptr3d(:,:,:)
      logical,intent(in),optional :: F_print_L
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-01
      !*@/
      integer :: istat
      character(len=GMM_MAXNAMELENGTH) :: name_S,vb_S
      type(gmm_metadata) :: mymeta
      logical :: print_L
      !---------------------------------------------------------------------
      nullify(F_ptr3d)
      print_L = .true.
      if (present(F_print_L)) print_L = F_print_L
      name_S = F_name_S
      istat = clib_tolower(name_S) !# Note: bug in GMM, NOT case independent
      vb_S = ' '
      F_istat = gmmx_meta(name_S,F_vb_S=vb_S)
      if (.not.RMN_IS_OK(F_istat)) then
         if (print_L) call msg(MSG_WARNING,'(gmmx) get3d - not such var: '//trim(name_S))
         return
      endif
      if (vb_S(1:1) /= 'g') then
         if (print_L) call msg(MSG_WARNING,'(gmmx) get3d - var not right type/rank: '//trim(name_S))
         F_istat = RMN_ERR
         return
      endif
      F_istat = gmm_getmeta(name_S,mymeta)
      if (.not.RMN_IS_OK(F_istat)) then
         if (print_L) call msg(MSG_WARNING,'(gmmx) get3d - GMM VAR not found: '//trim(name_S))
         F_istat = RMN_ERR
         return
      endif
      if (mymeta%l(3)%n < 1) then
         if (print_L) call msg(MSG_WARNING,'(gmmx) get3d - Rank mismatch: '//trim(name_S))
         F_istat = RMN_ERR
         return
      endif
      F_istat = gmm_get(name_S,F_ptr3d,mymeta)
      if (.not.(RMN_IS_OK(F_istat).and.associated(F_ptr3d))) then
         F_istat = RMN_ERR
         if (print_L) call msg(MSG_WARNING,'(gmmx) get3d - problem in gmm_get for: '//trim(name_S))
         return
      endif
!!$      print *,associated(F_ptr3d),F_istat,' (gmmx) Found ptr3d for: '//trim(name_S)
      !---------------------------------------------------------------------
      return
   end function gmmx_get_data_r4_3d


   !/@*
   function gmmx_get_data_r4_2d(F_name_S,F_ptr2d,F_print_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      character(len=*),intent(in) :: F_name_S !- Key (internal var name)
      real,pointer :: F_ptr2d(:,:)
      logical,intent(in),optional :: F_print_L
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-01
      !*@/
      integer :: istat
      character(len=GMM_MAXNAMELENGTH) :: name_S,vb_S
      type(gmm_metadata) :: mymeta
      logical :: print_L
      !---------------------------------------------------------------------
      nullify(F_ptr2d)
      print_L = .true.
      if (present(F_print_L)) print_L = F_print_L
      name_S = F_name_S
      istat = clib_tolower(name_S) !# Note: bug in GMM, NOT case independent
      vb_S = ' '
      F_istat = gmmx_meta(name_S,F_vb_S=vb_S)
      if (.not.RMN_IS_OK(F_istat)) then
         if (print_L) call msg(MSG_WARNING,'(gmmx) get2d - not such var: '//trim(name_S))
         return
      endif
      if (vb_S(1:1) /= 'g') then
         if (print_L) call msg(MSG_WARNING,'(gmmx) get2d - var not right type/rank: '//trim(name_S))
         F_istat = RMN_ERR
         return
      endif
      F_istat = gmm_getmeta(name_S,mymeta)
      if (.not.RMN_IS_OK(F_istat)) then
         if (print_L) call msg(MSG_WARNING,'(gmmx) get2d - GMM var not found: '//trim(name_S))
         F_istat = RMN_ERR
         return
      endif
      if (mymeta%l(3)%n /= 0) then
         if (print_L) call msg(MSG_WARNING,'(gmmx) get2d - Rank mismatch: '//trim(name_S))
         F_istat = RMN_ERR
         return
      endif
      F_istat = gmm_get(name_S,F_ptr2d,mymeta)
      if (.not.(RMN_IS_OK(F_istat).and.associated(F_ptr2d))) then
         F_istat = RMN_ERR
         if (print_L) call msg(MSG_WARNING,'(gmmx) get2d - problem in gmm_get for: '//trim(name_S))
         return
      endif
!!$      print *,associated(F_ptr2d),F_istat,' (gmmx) Found ptr2d for: '//trim(name_S)
      !---------------------------------------------------------------------
      return
   end function gmmx_get_data_r4_2d


   !/@*
   function gmmx_get_data_r4_bus(F_name_S,F_ptr) result(F_index)
      implicit none
      !@objective 
      !@arguments
      character(len=*),intent(in) :: F_name_S !- Key (internal var name)
      real,pointer :: F_ptr(:,:,:,:)
      !@return
      integer :: F_index !- exit status or bus index
      !@author  S. Chamberland, 2012-01
      !*@/
      character(len=GMM_MAXNAMELENGTH) :: name_S,vb_S
      !---------------------------------------------------------------------
      name_S = F_name_S
      F_index = gmmx_meta(name_S,F_vb_S=vb_S)
      if (.not.RMN_IS_OK(F_index)) then
         call msg(MSG_WARNING,'(gmmx_get_data) not such var: '//trim(F_name_S))
         return
      endif
      if (vb_S(1:1) == 'g') then
         call msg(MSG_WARNING,'(gmmx_get_data) var not right type/dim: '//trim(F_name_S))
         F_index = RMN_ERR
         return
      endif
      nullify(F_ptr)
      F_index = RMN_ERR
      call msg(MSG_ERROR,'(gmmx) Bus type not yet supported')

      !TODO-later: allocate bus if not done yet
      !TODO-later: F_ptr => bus_ptr(vb_S)
      !---------------------------------------------------------------------
      return
   end function gmmx_get_data_r4_bus


   !==== Private functions =================================================

   !/@*
   function priv_newvar_gmm(F_name_S,F_hgrid_S,F_vgrid_S,F_flags) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      character(len=*),intent(in) :: F_name_S  !- Key (internal var name)
      character(len=*),intent(in) :: F_hgrid_S,F_vgrid_S
      integer,intent(in) :: F_flags
      !@return
      integer :: F_istat !- exit status or bus index
      !@author  S. Chamberland, 2012-01
      !*@/
      type(gmm_metadata) :: mymeta,mymeta2
      type(vgrid_descriptor) :: myvgrid
      integer :: mytype,istat,kn,k0
      integer,pointer :: ip1list(:)
      real,pointer :: data2d(:,:),data3d(:,:,:)
      logical :: is_2d_L
      character(len=GMM_MAXNAMELENGTH) :: name_S
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
      nullify(ip1list,data2d,data3d)
      F_istat = vgrid_wb_get(F_vgrid_S,myvgrid,ip1list,mytype)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(gmmx) new - Probleme getting vgrid info for: '// &
              trim(F_name_S)//' VGRID='//trim(F_vgrid_S))
         return
      endif
      istat = vgd_free(myvgrid)
      k0 = lbound(ip1list,1)
      kn = ubound(ip1list,1)
      is_2d_L = .false.
      !TODO-later: 2d option if kn==k0 .and. mytype == ?
      deallocate(ip1list,stat=istat)
      F_istat = min(hgrid_wb_gmmmeta(F_hgrid_S,mymeta,kn,k0),F_istat)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(gmmx) new - Probleme getting hgrid info for: '// &
              trim(F_name_S)//' HGRID='//trim(F_hgrid_S))
         return
      endif
      name_S = F_name_S
      istat = clib_tolower(name_S) !# Note: GMM is NOT case independent

      F_istat = gmm_getmeta(name_S,mymeta2)
      if (RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(gmmx) new - Registering previously created GMM var for: ' &
              //trim(F_name_S))
         return
      endif

      if (is_2d_L) then
         F_istat = gmm_create(name_S,data2d,mymeta,F_flags)
         if (RMN_IS_OK(F_istat)) &
              write(msg_s,'(a,2i5,a,2i5,a)') ' [',lbound(data2d),';',ubound(data2d),'] 2d'
      else
         F_istat = gmm_create(name_S,data3d,mymeta,F_flags)
         if (RMN_IS_OK(F_istat)) &
              write(msg_s,'(a,3i5,a,3i5,a)') ' [',lbound(data3d),';',ubound(data3d),'] 3d'
      endif
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_ERROR,'(gmmx) newvar - Probleme creating GMM var for: '//trim(F_name_S))
         return
      endif
      call msg(MSG_INFOPLUS,'(gmmx) New GMM var:'//trim(name_S)//trim(msg_s))
      !---------------------------------------------------------------------
      return
   end function priv_newvar_gmm


!!$   !/@*
!!$   function priv_newvar_bus(F_name_S,F_hgrid_S,F_vgrid_S,F_flags,F_vb_S) result(F_istat)
!!$      implicit none
!!$      !@objective 
!!$      !@arguments
!!$      character(len=*),intent(in) :: F_name_S  !- Key (internal var name)
!!$      character(len=*),intent(in) :: F_hgrid_S,F_vgrid_S !# horiz/verti grid desc
!!$      character(len=*),intent(in) :: F_vb_S   !# bus name
!!$      integer,intent(in) :: F_flags
!!$      !@return
!!$      integer :: F_istat !- exit status or bus index
!!$      !@author  S. Chamberland, 2012-01
!!$      !*@/
!!$      !---------------------------------------------------------------------
!!$      F_istat = RMN_ERR
!!$      call msg(MSG_ERROR,'(gmmx) newvar - Bus type not yet supported, ignoring: '//trim(F_name_S))
!!$      !TODO-later: implement gmmx for bus
!!$      !TODO-later: return if F_kv_S(VAL,KEY_VB) is already allocated
!!$      !---------------------------------------------------------------------
!!$      return
!!$   end function priv_newvar_bus

end module gmmx_mod
