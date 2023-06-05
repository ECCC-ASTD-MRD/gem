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
module vgrid_wb
   use, intrinsic :: iso_fortran_env, only: REAL64, INt64
   use iso_c_binding
   use rpn_comm_itf_mod
   use vGrid_Descriptors
   use vgrid_ov, only: vgrid_nullify
   use wb_itf_mod
   use rmn_gmm
   implicit none
   private
   !@objective Whiteboard (data store) for vgrid + ip1 <=> index association
   !@author Stephane Chamberland, 2012-01
   !@description
   ! Public functions
   public :: vgrid_wb_is_press_kind
   public :: vgrid_wb_exists, vgrid_wb_put, vgrid_wb_get, vgrid_wb_bcast
   public :: vgrid_wb_bcast_s, vgrid_wb_bcast_v, vgrid_wb_bcast_v2
   ! Public constants
   !TODO-FR: (feature request) include these in the vGrid_Descriptors module
   integer, parameter, public :: VGRID_SIGM_KIND = 1 !Sigma
   integer, parameter, public :: VGRID_SIGM_VER  = 1

   integer, parameter, public :: VGRID_ETA_KIND = 1 !Eta
   integer, parameter, public :: VGRID_ETA_VER  = 2

   integer, parameter, public :: VGRID_HYBN_KIND = 1 !Hybrid Normalized
   integer, parameter, public :: VGRID_HYBN_VER  = 3

!!$   integer,parameter,public :: VGRID_????_KIND = 1 !???
!!$   integer,parameter,public :: VGRID_????_VER  = 4

   integer, parameter, public :: VGRID_PRES_KIND = 2 !pressure
   integer, parameter, public :: VGRID_PRES_VER  = 1

   integer, parameter, public :: VGRID_HYB_KIND = 5 !Hybrid Un-staggered
   integer, parameter, public :: VGRID_HYB_VER  = 1

   integer, parameter, public :: VGRID_HYBS_KIND = 5 !Hybrid staggered
   integer, parameter, public :: VGRID_HYBS_VER  = 2

   integer, parameter, public :: VGRID_HYBT_KIND = 5 !Hybrid staggered with unstaggered last Thermo level
   integer, parameter, public :: VGRID_HYBT_VER  = 3

   integer, parameter, public :: VGRID_HYBM_KIND = 5 !Hybrid staggered, first level is a momentum level,
   !        same number of thermo and momentum levels
   integer, parameter, public :: VGRID_HYBM_VER  = 4

   integer, parameter, public :: VGRID_HYBMD_KIND = 5 !Hybrid staggered, first level is a momentum level,
   !        same number of thermo and momentum levels, Diag level heights (m AGL) encoded
   integer, parameter, public :: VGRID_HYBMD_VER  = 5


   integer,parameter,public :: VGRID_HYBSLEVE_KIND = 5  !Hybrid staggered SLEVE like coordinate
   integer,parameter,public :: VGRID_HYBSLEVE_VER  = 100

   integer,parameter,public :: VGRID_HYBO_KIND = 5      !Hybrid unstaggered coordinate of unknown/other origin (like ECMWF)
   integer,parameter,public :: VGRID_HYBO_VER  = 999

   integer,parameter,public :: VGRID_GC_KIND = 21     !Hybrid heights levels on Charney-Philips grid, may be non SLEVE (reference field: ME) or SLEVE (reference fields: ME, MELS) 
   integer,parameter,public :: VGRID_GC_VER  = 1

   integer,parameter,public :: VGRID_LORZ_KIND = 21    !Hybrid heights levels on Lorenz grid, may be non SLEVE (reference field: ME) or SLEVE (reference fields: ME, MELS) 
   integer,parameter,public :: VGRID_LORZ_VER  = 2

   !#TODO: add missing VGRID *_KIND and *_VER

   integer, parameter, public :: VGRID_GROUND_TYPE = -1
   integer, parameter, public :: VGRID_SURF_TYPE  = 0
   integer, parameter, public :: VGRID_UPAIR_TYPE = 1
   integer, parameter, public :: VGRID_UPAIR_M_TYPE = 1
   integer, parameter, public :: VGRID_UPAIR_T_TYPE = 2
   integer, parameter, public :: VGRID_DIAG_TYPE  = 3
   integer, parameter, public :: VGRID_DIAG_M_TYPE  = 3
   integer, parameter, public :: VGRID_DIAG_T_TYPE  = 4

   character(len=6), parameter, public :: VGRID_TYPES_S(-1:4) = &
        (/'Ground', 'Surf  ', 'UpAir ', 'UpAirT', 'Diag  ', 'DiagT '/)

!@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   character(len=*),parameter :: PREFIXI_S = 'VGWBI/'
   character(len=*),parameter :: PREFIXV_S = 'VGWB/'

   interface vgrid_wb_put
      module procedure vgrid_wb_put_i
      module procedure vgrid_wb_put_v
      module procedure vgrid_wb_put_v2
   end interface

   interface vgrid_wb_get
      module procedure vgrid_wb_get_i
      module procedure vgrid_wb_get_s
   end interface

   interface vgrid_wb_bcast
      module procedure vgrid_wb_bcast_v
      module procedure vgrid_wb_bcast_v2
      module procedure vgrid_wb_bcast_s
   end interface vgrid_wb_bcast

   interface vgrid_wb_is_press_kind
      module procedure vgrid_wb_is_press_kind_i
      module procedure vgrid_wb_is_press_kind_s
      module procedure vgrid_wb_is_press_kind_v
   end interface vgrid_wb_is_press_kind

contains

   !/@*
   function vgrid_wb_is_press_kind_i(F_kind) result(F_pkind_L)
      implicit none
      integer, intent(in) :: F_kind  !# vgrid coor kind
      logical :: F_pkind_L           !# True if pressure kind
      F_pkind_L = all(F_kind /= (/VGRID_GC_KIND,  VGRID_LORZ_KIND/))
      return
   end function vgrid_wb_is_press_kind_i

   !/@*
   function vgrid_wb_is_press_kind_v(F_vgrid) result(F_pkind_L)
      implicit none
      type(vgrid_descriptor), intent(in) :: F_vgrid
      logical :: F_pkind_L           !# True if pressure kind
      integer :: istat, ikind
      istat = vgd_get(F_vgrid, key='KIND', value=ikind)
      F_pkind_L = vgrid_wb_is_press_kind_i(ikind)
      return
   end function vgrid_wb_is_press_kind_v
   
   !/@*
   function vgrid_wb_is_press_kind_S(F_vgrid_S) result(F_pkind_L)
      implicit none
      character(len=*), intent(in) :: F_vgrid_S
      logical :: F_pkind_L           !# True if pressure kind
      integer :: istat
      type(vgrid_descriptor) :: vgrid
      istat = vgrid_wb_get(F_vgrid_S, vgrid)
      F_pkind_L = vgrid_wb_is_press_kind_v(vgrid)
      return
   end function vgrid_wb_is_press_kind_S
  
   !/@*
   function vgrid_wb_exists(F_name_S, F_id, F_id_S, F_ip1list, F_itype, F_restart_L) &
        result(F_istat)
      implicit none
      !@objective Check if F_name_S is already in use
      !@arguments
      character(len=*), intent(in)  :: F_name_S
      integer, intent(out), optional :: F_id
      character(len=*), intent(out), optional :: F_id_S
      integer, pointer, optional :: F_ip1list(:)
      integer, intent(out), optional :: F_itype
      logical, intent(out), optional :: F_restart_L
      !@return
      integer :: F_istat
      !@author  S. Chamberland, 2012-05
      !*@/
      integer :: id, istat
      integer, pointer :: i4ptr1d(:)
      type(gmm_metadata) :: i4meta1d
      character(len=GMM_MAXNAMELENGTH) :: id_S
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      nullify(i4ptr1d)
      id_S = priv_name(F_name_S, id)
      istat = gmm_get(trim(PREFIXI_S)//trim(id_S), i4ptr1d, i4meta1d)
      if (associated(i4ptr1d)) F_istat = RMN_OK
      if (present(F_id))      F_id = id
      if (present(F_id_S))    F_id_S = id_S
      if (present(F_ip1list)) F_ip1list => i4ptr1d
      if (present(F_itype))   F_itype = i4meta1d%a%uuid1
      if (present(F_restart_L)) F_restart_L = (iand(i4meta1d%a%flags, GMM_FLAG_READ) > 0)
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_exists


   !/@*
   function vgrid_wb_put_i(F_name_S, F_type, F_ip1list, F_overwrite_L, F_restart_L) &
        result(F_id)
      implicit none
      !@objective Store a new vgrid
      !@arguments
      character(len=*), intent(in) :: F_name_S !- Key (internal var name)
      integer, intent(in) :: F_type
      integer, pointer :: F_ip1list(:)         !- list of ip1
      logical, intent(in), optional :: F_overwrite_L
      logical, intent(in), optional :: F_restart_L
      !@return
      integer :: F_id
      !@author  S. Chamberland, 2012-01
      !*@/
      character(len=256) :: msg_S, name_S, id_S, id2_S
      integer :: istat, lip1, uip1, itype, id, flags
      integer, pointer :: i4ptr1d(:)
      logical :: overwrite_L, same_L, same2_L, exists_L, restart_L
      type(gmm_metadata) :: i4meta1d
      !---------------------------------------------------------------------
      F_id = RMN_ERR
      name_S = F_name_S
      if (len_trim(F_name_S) == 0) then
         call msg(MSG_ERROR, '(vgrid_wb_put) need to provide a internal name')
         return
      endif
      overwrite_L = .false.
      if (present(F_overwrite_L)) overwrite_L = F_overwrite_L
      restart_L = .true.
      if (present(F_restart_L)) restart_L = F_restart_L
      flags = 0
      if (restart_L) flags = GMM_FLAG_RSTR

      if (.not.any(F_type == (/VGRID_GROUND_TYPE, VGRID_SURF_TYPE, VGRID_UPAIR_TYPE, VGRID_UPAIR_M_TYPE, VGRID_UPAIR_T_TYPE/))) then
         call msg(MSG_ERROR, '(vgrid_wb_put) invalid vgrid_wb type: '//trim(F_name_S))
         return
      endif

      nullify(i4ptr1d)
      istat    = vgrid_wb_exists(F_name_S, id, id_S, i4ptr1d, itype)
      exists_L = (RMN_IS_OK(istat) .and. associated(i4ptr1d))
      id2_S    = trim(PREFIXI_S)//trim(id_S)

      IF_EXISTS: if (exists_L) then
         same_L  = (itype == F_type .and. size(i4ptr1d) == size(F_ip1list))
         same2_L = same_L
         if (same_L) same2_L = all(i4ptr1d == F_ip1list)
         if (overwrite_L) then
            if (.not.same_L) then
               istat    = gmm_delete(trim(id2_S))
               exists_L = .false.
            else
               call msg(MSG_INFO, '(vgrid_wb) Put, updating: '//trim(F_name_S))
            endif
         else
            if (.not.(same_L.and.same2_L)) then
               call msg(MSG_ERROR, '(vgrid_wb_put) vgrid already exists with different params1: '//trim(F_name_S))
               return
            else
               F_id = id
               return
            endif
         endif
      endif IF_EXISTS

      lip1 = lbound(F_ip1list, 1)
      uip1 = ubound(F_ip1list, 1)
      itype = min(max(VGRID_GROUND_TYPE, F_type), VGRID_UPAIR_T_TYPE)

      IF_EXISTS2: if (.not.exists_L) then
         i4meta1d = GMM_NULL_METADATA
         i4meta1d%l(1) = gmm_layout(lip1, uip1, 0, 0, uip1-lip1+1)
         i4meta1d%a%uuid1 = itype
         nullify(i4ptr1d)
         istat = gmm_create(trim(id2_S), i4ptr1d, i4meta1d, flags)
         if (.not.associated(i4ptr1d)) then
            call msg(MSG_ERROR, '(vgrid_wb_put) cannot allocate mem for ip1list: '//trim(F_name_S))
            return
         endif
      endif IF_EXISTS2

      i4ptr1d(:) = F_ip1list(:)

      write(msg_S, '(a, " [type=", a, "] ip1[", i4, ":", i4, "] = (", i12, ", ..., ", i12, ") id=", i16.16)') &
           name_S(1:16), VGRID_TYPES_S(itype), lip1, uip1, F_ip1list(lip1), F_ip1list(uip1), id
      call msg(MSG_INFO, '(vgrid_wb) Put: '//trim(msg_S))
      F_id = id
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_put_i


   !/@*
   function vgrid_wb_put_v(F_name_S, F_vgrid, F_ip1list, F_sfcfld_S, &
        F_overwrite_L, F_altfld_S, F_restart_L) result(F_id)
      implicit none
      !@objective Store a new vgrid
      !@arguments
      character(len=*), intent(in) :: F_name_S       !- Key (internal var name)
      type(vgrid_descriptor), intent(in) :: F_vgrid  !- 
      integer, pointer :: F_ip1list(:)               !- list of ip1
      character(len=*), intent(in), optional :: F_sfcfld_S !- Name of ref sfc fields for levels computations
      logical, intent(in), optional :: F_overwrite_L
      character(len=*), intent(in), optional :: F_altfld_S !- Name of alternate 3d field descibing the vertical - in M for Press vertcoor and vice-versa
      logical, intent(in), optional :: F_restart_L
      !@return
      integer :: F_id
      !@author  S. Chamberland, 2012-01
      !*@/
      character(len=GMM_MAXNAMELENGTH) :: sfcfld_S, sfcfld2_S, altfld_S
      logical :: overwrite_L, restart_L
      !---------------------------------------------------------------------
      overwrite_L = .false.
      if (present(F_overwrite_L)) overwrite_L = F_overwrite_L
      restart_L = .true.
      if (present(F_restart_L)) restart_L = F_restart_L
      sfcfld_S = ' '
      if (present(F_sfcfld_S)) sfcfld_S = F_sfcfld_S
      altfld_S = ' '
      if (present(F_altfld_S)) altfld_S = F_altfld_S
      sfcfld2_S = ' '
      F_id = vgrid_wb_put_v2(F_name_S, F_vgrid, F_ip1list, sfcfld_S, &
           sfcfld2_S, overwrite_L, F_altfld_S=altfld_S, F_restart_L=restart_L)
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_put_v


   !/@*
   function vgrid_wb_put_v2(F_name_S, F_vgrid, F_ip1list, F_sfcfld_S, &
        F_sfcfld2_S, F_overwrite_L, F_altfld_S, F_restart_L) result(F_id)
      implicit none
      !@objective Store a new vgrid
      !@arguments
      character(len=*),intent(in) :: F_name_S       !- Key (internal var name)
      type(vgrid_descriptor),intent(in) :: F_vgrid  !- 
      integer,pointer :: F_ip1list(:)               !- list of ip1
      character(len=*),intent(in) :: F_sfcfld_S !- Name of ref sfc fields for levels computations
      character(len=*),intent(in) :: F_sfcfld2_S !- Name of ref sfc fields for levels computations
      logical,intent(in),optional :: F_overwrite_L
      character(len=*), intent(in), optional :: F_altfld_S !- Name of alternate 3d field descibing the vertical - in M for Press vertcoor and vice-versa
      logical, intent(in), optional :: F_restart_L
      !@return
      integer :: F_id
      !@author  S. Chamberland, 2012-01
      !*@/
      character(len=256) :: msg_S
      character(len=WB_MAXSTRINGLENGTH) :: sfcfld_S, sfcfld2_S, id_S, id2_S, altfld_S
      integer :: istat, lijk(3), uijk(3), n, lip1, uip1, flags
      logical :: overwrite_L, exists_L, same_L, same2_L, restart_L
      real(REAL64), pointer :: vtbl(:, :, :), vtbl2(:, :, :)
      type(gmm_metadata) :: r8meta3d, r8meta3d2
      !---------------------------------------------------------------------
      overwrite_L = .false.
      if (present(F_overwrite_L)) overwrite_L = F_overwrite_L
      restart_L = .true.
      if (present(F_restart_L)) restart_L = F_restart_L
      altfld_S = ' ' 
      if (present(F_altfld_S)) altfld_S = F_altfld_S
      flags = 0
      if (restart_L) flags = GMM_FLAG_RSTR

      if (associated(F_ip1list)) then
         lip1 = lbound(F_ip1list,1)
         uip1 = ubound(F_ip1list,1)
         write(msg_S, &
              '(a, " ip1[", i4, ":", i4, "] = (", i12, ", ..., ", i12, ") sfcref=", a)') &
              trim(F_name_S), lip1, uip1, F_ip1list(lip1), F_ip1list(uip1),trim(F_sfcfld_S)//' '//trim(F_sfcfld2_S)//'; alt='//trim(altfld_S)
      else
          write(msg_S, '(a)') trim(F_name_S)//" sfcref="//trim(F_sfcfld_S)//' '//trim(F_sfcfld2_S)
      endif
      call msg(MSG_INFO, '(vgrid_wb) Put: '//trim(msg_S))

      F_id = vgrid_wb_put_i(F_name_S, VGRID_UPAIR_TYPE, F_ip1list, &
           overwrite_L, F_restart_L=restart_L)
      if (.not.RMN_IS_OK(F_id)) return

      nullify(vtbl)
      istat = vgd_get(F_vgrid, key='VTBL', value=vtbl)
      if (istat /= VGD_OK .or. .not.associated(vtbl)) then
         call msg(MSG_WARNING, '(vgrid_wb_put) problem cloning vgrid for: '//trim(F_name_S))
         F_id = RMN_ERR
         return
      endif

      sfcfld_S = ' '
      if (F_sfcfld_S /= ' ') then
         sfcfld_S = F_sfcfld_S
      else
         istat = vgd_get(F_vgrid, key='RFLD', value=sfcfld_S)
         if (istat /= VGD_OK) sfcfld_S = ' '
      endif

      sfcfld2_S = ' '
      if (F_sfcfld2_S /= ' ') then
         sfcfld2_S = F_sfcfld2_S
      else
         istat = vgd_get(F_vgrid, key='RFLS', value=sfcfld2_S)
         if (istat /= VGD_OK) sfcfld2_S = ' '
      endif

      id_S  = priv_name(F_name_S)
      id2_S = trim(PREFIXV_S)//trim(id_S)
      lijk = lbound(vtbl)
      uijk = ubound(vtbl)
      r8meta3d = GMM_NULL_METADATA
      do n=1, 3
         r8meta3d%l(n) = gmm_layout(lijk(n), uijk(n), 0, 0, uijk(n)-lijk(n)+1)
      enddo
      r8meta3d%a%uuid1 = transfer(sfcfld_S(1:8),  r8meta3d%a%uuid1)
      r8meta3d%a%uuid2 = transfer(sfcfld_S(9:16), r8meta3d%a%uuid2)

      nullify(vtbl2)
      istat = gmm_get(trim(id2_S), vtbl2, r8meta3d2)
      exists_L = associated(vtbl2)

      IF_EXISTS3: if (exists_L) then
         same_L = (size(vtbl2) == size(vtbl) .and. &
              r8meta3d2%a%uuid1 == r8meta3d%a%uuid1 .and. &
              r8meta3d2%a%uuid2 == r8meta3d%a%uuid2)
         same2_L = same_L
         if (same_L) same2_L = all(vtbl2 == vtbl)

         if (overwrite_L) then
            if (.not.same_L) then
               istat    = gmm_delete(trim(id2_S))
               exists_L = .false.
            else
               call msg(MSG_INFO, '(vgrid_wb) Put, updating: '//trim(F_name_S))
            endif
         else
            if (.not.(same_L.and.same2_L)) then
               call msg(MSG_ERROR, '(vgrid_wb_put) vgrid already exists with different params2: '//trim(F_name_S))
               F_id = RMN_ERR
               deallocate(vtbl, stat=istat)
               return
!!$            else
!!$               F_id = id
!!$               return
            endif
         endif
      endif IF_EXISTS3

      IF_EXISTS4: if (exists_L) then
         vtbl2(:, :, :) = vtbl(:, :, :)
      else
         istat = gmm_create(trim(id2_S), vtbl, r8meta3d, flags)
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR, '(vgrid_wb_put) problem storing vtbl: '//trim(F_name_S))
            F_id = RMN_ERR
            return
         endif
      endif IF_EXISTS4

      istat = wb_put(trim(id2_S)//'/RFLS', sfcfld2_S, WB_REWRITE_MANY)
      istat = wb_put(trim(id2_S)//'/ALTF', altfld_S, WB_REWRITE_MANY)
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_put_v2


   !/@*
   function vgrid_wb_get_s(F_name_S, F_vgrid, F_ip1list, F_type, &
        F_sfcfld_S, F_sfcfld2_S, F_altfld_S) result(F_istat)
      implicit none
      !@objective Retreive stored vgrid
      !@arguments
      character(len=*), intent(in) :: F_name_S !- Key (internal var name)
      type(vgrid_descriptor), intent(out) :: F_vgrid    !- vgrid struct
      integer, pointer, optional :: F_ip1list(:)         !- list of ip1
      integer, intent(out), optional :: F_type
      character(len=*), intent(out), optional :: F_sfcfld_S  !- Name of ref sfc fields for levels computations
      character(len=*), intent(out), optional :: F_sfcfld2_S !- Name of ref sfc fields for levels computations
      character(len=*), intent(out), optional :: F_altfld_S  !- Name of alternate 3d field descibing the vertical - in M for Press vertcoor and vice-versa
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-01
      !*@/
      integer :: vgrid_idx, itype, lip1, uip1
      character(len=GMM_MAXNAMELENGTH) :: sfcfld_S, sfcfld2_S, altfld_S
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(vgrid_wb) get [BEGIN] '//trim(F_name_S))
      call vgrid_nullify(F_vgrid)
      vgrid_idx = priv_id(F_name_S)

      if (present(F_ip1list)) then
         F_istat = vgrid_wb_get_i(vgrid_idx, F_vgrid, F_ip1list, F_type=itype, &
              F_sfcfld_S=sfcfld_S, F_name_S=F_name_S, F_sfcfld2_S=sfcfld2_S, &
              F_altfld_S=altfld_S)
      else
         F_istat = vgrid_wb_get_i(vgrid_idx, F_vgrid, F_type=itype, &
              F_sfcfld_S=sfcfld_S, F_name_S=F_name_S, F_sfcfld2_S=sfcfld2_S, &
              F_altfld_S=altfld_S)
      endif
      if (.not.RMN_IS_OK(F_istat)) then
         itype = -1
         sfcfld_S = ' '
         sfcfld2_S = ' '
         altfld_S = ' '
      endif
      if (present(F_type)) F_type = itype
      if (present(F_sfcfld_S))  F_sfcfld_S  = sfcfld_S
      if (present(F_sfcfld2_S)) F_sfcfld2_S = sfcfld2_S
      if (present(F_altfld_S))  F_altfld_S  = altfld_S

      if (RMN_IS_OK(F_istat)) then
         if (present(F_ip1list)) then
            lip1 = lbound(F_ip1list,1)
            uip1 = ubound(F_ip1list,1)
            write(msg_S, &
                 '(a, " [type=", i4, "] ip1[", i4, ":", i4, "] = (", i12, ", ..., ", i12, ") sfcref=", a)') &
                 trim(F_name_S), itype, lip1, uip1, F_ip1list(lip1), F_ip1list(uip1),trim(sfcfld_S)//' '//trim(sfcfld2_S)
         else
            write(msg_S, '(a, " [type=", i4, "] sfcref=", a)') trim(F_name_S), itype, trim(sfcfld_S)//' '//trim(sfcfld2_S)
         endif
         msg_S = trim(msg_S)//' ; alt='//trim(altfld_S)
      else
         msg_S = trim(F_name_S) // 'Not Found'
      endif
      call msg(MSG_DEBUG, '(vgrid_wb) get [END] '//trim(msg_S))
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_get_s


   !/@*
   function vgrid_wb_get_i(F_id, F_vgrid, F_ip1list, F_type, F_sfcfld_S, &
        F_name_S, F_sfcfld2_S, F_altfld_S) result(F_istat)
      implicit none
      !@objective Retreive stored vgrid
      !@arguments
      integer, intent(in) :: F_id !- vgrid id returned by vgrid_wb_put
      type(vgrid_descriptor), intent(out) :: F_vgrid    !- vgrid struct
      integer, pointer, optional :: F_ip1list(:)         !- list of ip1
      integer, intent(out), optional :: F_type
      character(len=*), intent(out), optional :: F_sfcfld_S !- Name of ref sfc fields for levels computations
      character(len=*), intent(in), optional :: F_name_S !- Key (internal var name)
      character(len=*),intent(out), optional :: F_sfcfld2_S !- Name of ref sfc fields for levels computations
      character(len=*), intent(out), optional :: F_altfld_S  !- Name of alternate 3d field descibing the vertical - in M for Press vertcoor and vice-versa
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-01
      !*@/
      character(len=GMM_MAXNAMELENGTH) :: name_S, id_S, id2_S, tmp_S, tmp2_S
      character(len=WB_MAXSTRINGLENGTH) :: sfcfld2_S, altfld_S
      integer :: istat, nip1, lip1, uip1
      integer, pointer :: i4ptr1d(:)
      real(REAL64), pointer :: vtbl(:, :, :)
      type(gmm_metadata) :: i4meta1d, r8meta3d
      !---------------------------------------------------------------------
      F_istat = RMN_ERR
      call vgrid_nullify(F_vgrid) 
      write(id_S, '(i16.16)') F_id
      name_S = id_S
      if (present(F_name_S)) name_S = trim(F_name_S)//':'//trim(name_S)
      nullify(i4ptr1d, vtbl)
      istat = gmm_get(trim(PREFIXI_S)//trim(id_S), i4ptr1d, i4meta1d)
      if (.not.associated(i4ptr1d)) then
         call msg(MSG_INFOPLUS, '(vgrid_wb_get) vgrid not found: '//trim(name_S))
         return
      endif
      if (present(F_type)) F_type = i4meta1d%a%uuid1
      nip1 = size(i4ptr1d)
      if (present(F_ip1list)) then
         lip1 = lbound(i4ptr1d, 1)
         uip1 = ubound(i4ptr1d, 1)
         if (.not.associated(F_ip1list)) then
            allocate(F_ip1list(lip1:uip1), stat=istat)
            if (istat /= 0) then
               call msg(MSG_ERROR, '(vgrid_wb_get) Cannot allocate memory for ip1list: '//trim(name_S))
               return
            endif
         endif
         if (lbound(F_ip1list, 1) /= lip1 .or. ubound(F_ip1list, 1) < uip1) then
            if (lbound(F_ip1list, 1) /= lip1) then
               call msg(MSG_ERROR, '(vgrid_wb_get) provided ip1list lbound mismatch: '//trim(name_S))
            else
               call msg(MSG_ERROR, '(vgrid_wb_get) provided ip1list size mismatch: '//trim(name_S))
            endif
            return
         endif
         F_ip1list = -1
         F_ip1list(lip1:uip1) = i4ptr1d(lip1:uip1)
      endif

      id2_S = trim(PREFIXV_S)//trim(id_S)
      istat = gmm_get(trim(id2_S), vtbl, r8meta3d)
      istat = VGD_OK
      if (associated(vtbl)) then
         istat = vgd_new(F_vgrid, vtbl)
         if (istat /= VGD_OK) then
            call msg(MSG_ERROR, '(vgrid_wb_get) problem cloning vgrid for: '//trim(F_name_S))
            return
         endif
         if (present(F_sfcfld_S)) then
            tmp_S = ' '
            tmp2_S = ' '
            tmp_S  = transfer(r8meta3d%a%uuid1, tmp_S)
            tmp2_S = transfer(r8meta3d%a%uuid2, tmp2_S)
            F_sfcfld_S = tmp_S(1:8)//tmp2_S(1:8)
         endif
         if (present(F_sfcfld2_S)) then
            F_sfcfld2_S = ' '
            istat = wb_get(trim(id2_S)//'/RFLS', sfcfld2_S)
            if (RMN_IS_OK(istat)) F_sfcfld2_S = sfcfld2_S
         endif
         if (present(F_altfld_S)) then
            F_altfld_S = ' '
            istat = wb_get(trim(id2_S)//'/ALTF', altfld_S)
            if (RMN_IS_OK(istat)) F_altfld_S = altfld_S
         endif
      endif
      F_istat = nip1
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_get_i


   !/@*
   function vgrid_wb_bcast_s(F_name_S, F_comm_S, F_ipe_master, F_ipe) result(F_istat)
      implicit none
      !@objective  MPI bcast stored vgrid
      !@arguments
      character(len=*), intent(in) :: F_name_S  !- Key (internal var name)
      character(len=*), intent(in) :: F_comm_S  !- RPN_COMM communicator name
      integer, intent(in), optional :: F_ipe_master !- Sending PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      integer, intent(in), optional :: F_ipe    !- PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-08
      !*@/
      character(len=GMM_MAXNAMELENGTH) :: sfcfld_S, sfcfld2_S, altfld_S
      integer :: itype, me, istat, ipe_master
      logical :: ismaster_L
      type(vgrid_descriptor) :: vgrid
      integer, pointer :: ip1list(:)
      !---------------------------------------------------------------------
      F_istat = RMN_OK
      if (present(F_ipe)) then
         me = F_ipe
      else
         call rpn_comm_rank(F_comm_S, me, istat)
      endif
      ipe_master = RPN_COMM_MASTER
      if (present(F_ipe_master)) ipe_master = F_ipe_master
      ismaster_L = (me == ipe_master)
      nullify(ip1list)
      sfcfld_S = ' '
      sfcfld2_S = ' '
      altfld_S = ' '
      if (ismaster_L) then
         F_istat = vgrid_wb_get(F_name_S,vgrid,ip1list,itype,sfcfld_S,sfcfld2_S, F_altfld_S=altfld_S)
      endif
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return
      F_istat = vgrid_wb_bcast(vgrid, ip1list, itype ,sfcfld_S, sfcfld2_S, F_comm_S, ipe_master, me, F_altfld_S=altfld_S)
      if (RMN_IS_OK(F_istat) .and. .not.ismaster_L) then
         if (itype < VGRID_UPAIR_TYPE) then
            F_istat = vgrid_wb_put(F_name_S, itype, ip1list)
         else
            F_istat = vgrid_wb_put(F_name_S, vgrid, ip1list, sfcfld_S, sfcfld2_S, F_altfld_S=altfld_S)
         endif
      endif
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_bcast_s


   !/@*
   function vgrid_wb_bcast_v(F_vgrid, F_ip1list, F_itype, F_sfcfld_S, &
        F_comm_S, F_ipe_master, F_ipe) result(F_istat)
      implicit none
      !@objective  MPI bcast stored vgrid
      !@arguments
      type(vgrid_descriptor),intent(inout) :: F_vgrid 
      integer,pointer :: F_ip1list(:)            !- list of ip1
      integer,intent(inout) :: F_itype
      character(len=*),intent(inout) :: F_sfcfld_S
      character(len=*),intent(in) :: F_comm_S    !- RPN_COMM communicator name
      integer, intent(in), optional :: F_ipe_master !- Sending PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      integer, intent(in), optional :: F_ipe      !- PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-08
      !*@/
      character(len=GMM_MAXNAMELENGTH) :: sfcfld2_S
      integer :: me, ipe_master, istat
      !---------------------------------------------------------------------
      if (present(F_ipe)) then
         me = F_ipe
      else
         call rpn_comm_rank(F_comm_S, me, istat)
      endif
      ipe_master = RPN_COMM_MASTER
      if (present(F_ipe_master)) ipe_master = F_ipe_master
      sfcfld2_S = ' '
      F_istat = vgrid_wb_bcast_v2(F_vgrid, F_ip1list, F_itype, F_sfcfld_S, sfcfld2_S, F_comm_S, ipe_master, me)
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_bcast_v


   !/@*
   function vgrid_wb_bcast_v2(F_vgrid, F_ip1list, F_itype, F_sfcfld_S, &
        F_sfcfld2_S, F_comm_S, F_ipe_master, F_ipe, F_altfld_S) result(F_istat)
      implicit none
      !@objective  MPI bcast stored vgrid
      !@arguments
      type(vgrid_descriptor), intent(inout) :: F_vgrid
      integer, pointer :: F_ip1list(:)            !- list of ip1
      integer, intent(inout) :: F_itype
      character(len=*), intent(inout) :: F_sfcfld_S
      character(len=*), intent(inout) :: F_sfcfld2_S
      character(len=*), intent(in) :: F_comm_S    !- RPN_COMM communicator name
      integer, intent(in), optional :: F_ipe_master !- Sending PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      integer, intent(in), optional :: F_ipe      !- PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      character(len=*), intent(inout), optional :: F_altfld_S !- Name of alternate 3d field descibing the vertical - in M for Press vertcoor and vice-versa
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-08
      !*@/
      integer, parameter :: NMAXIP1 = 2048
      integer, parameter :: CHARPERBYTE = 4
      integer, parameter :: STRLEN = 32
      integer, parameter :: STRSIZE = STRLEN/CHARPERBYTE
      integer, parameter :: ADDINT = 3
      integer, parameter :: IBUFSIZE = ADDINT + NMAXIP1 + 3*STRSIZE
      integer :: me, istat, istat2, nip1, i0, in, ipe_master, ibuf4(STRSIZE)
      logical :: ismaster_L, ok_L
      integer :: ibuf(IBUFSIZE), n123(3)
      character(len=STRLEN) :: sfcfld_S, sfcfld2_S, altfld_S
      real(REAL64), pointer :: vtbl_8(:,:,:)
      !---------------------------------------------------------------------
      F_istat = RMN_OK
      if (present(F_ipe)) then
         me = F_ipe
      else
         call rpn_comm_rank(F_comm_S, me, istat)
      endif
      ipe_master = RPN_COMM_MASTER
      if (present(F_ipe_master)) ipe_master = F_ipe_master
      altfld_S = ' '
      if (present(F_altfld_S)) altfld_S = F_altfld_S

      ismaster_L = (me == ipe_master)
      nullify(vtbl_8)
      if (ismaster_L) then
         F_istat = vgd_get(F_vgrid, 'VTBL', vtbl_8)
         n123 = ubound(vtbl_8)
         sfcfld_S  = F_sfcfld_S
         sfcfld2_S = F_sfcfld2_S
         ibuf = 0
         ibuf(1:STRSIZE) = transfer(sfcfld_S,ibuf4)
         ibuf(1+STRSIZE:2*STRSIZE) = transfer(sfcfld2_S,istat)
         ibuf(1+2*STRSIZE:3*STRSIZE) = transfer(altfld_S,istat)
         ibuf(3*STRSIZE+1) = lbound(F_ip1list,1)
         ibuf(3*STRSIZE+2) = ubound(F_ip1list,1)
         ibuf(3*STRSIZE+3) = F_itype
         nip1 = size(F_ip1list)
         ibuf(3*STRSIZE+(ADDINT+1):3*STRSIZE+(ADDINT+1)+(nip1-1)) = F_ip1list(:)
      endif
      call rpn_comm_bcast(n123, size(n123), RPN_COMM_INTEGER, ipe_master, F_comm_S, istat)
      if (.not.ismaster_L) then
         allocate(vtbl_8(n123(1), n123(2), n123(3)), stat=istat)
      endif
      call rpn_comm_bcast(vtbl_8, size(vtbl_8), RPN_COMM_REAL8, ipe_master, F_comm_S, istat)
      call rpn_comm_bcast(ibuf, size(ibuf), RPN_COMM_INTEGER, ipe_master, F_comm_S, istat2)
      F_istat = min(F_istat, istat, istat2)
      if (.not.ismaster_L) then
         F_istat = vgd_new(F_vgrid, vtbl_8)
         sfcfld_S = ' '
         sfcfld_S = transfer(ibuf(1:STRSIZE), sfcfld_S)
         F_sfcfld_S = sfcfld_S
         sfcfld2_S = ' '
         sfcfld2_S = transfer(ibuf(1+STRSIZE:2*STRSIZE),sfcfld2_S)
         F_sfcfld2_S = sfcfld2_S
         if (present(F_altfld_S)) then
            altfld_S = ' '
            altfld_S = transfer(ibuf(1+2*STRSIZE:3*STRSIZE),altfld_S)
            F_altfld_S = altfld_S
         endif
         i0 = ibuf(3*STRSIZE+1)
         in = ibuf(3*STRSIZE+2)
         F_itype = ibuf(3*STRSIZE+3)
         ok_L = associated(F_ip1list)
         if (ok_L) ok_L = (lbound(F_ip1list,1) == i0 .and. ubound(F_ip1list,1) >= in)
         !#TODO: Memory leak if F_ip1list already allocated but size too small
         if (.not.ok_L)then
            nullify(F_ip1list)
            allocate(F_ip1list(i0:in))
         endif
         F_ip1list(i0:in) = ibuf(3*STRSIZE+(ADDINT+1):3*STRSIZE+(ADDINT+1)+in-i0)
      endif
      if (associated(vtbl_8)) deallocate(vtbl_8, stat=istat)
!!$      call collect_error(F_istat)
!!$      if (.not.RMN_IS_OK(F_istat)) return
      !---------------------------------------------------------------------
      return
   end function vgrid_wb_bcast_v2

   !==== Private Functions =================================================

   function priv_id(F_name_S) result(F_id)
      implicit none
      character(len=*), intent(in) :: F_name_S
      integer :: F_id

      integer, parameter :: NWORD=8
      integer, parameter :: NCHARPERWORD=4
      integer, parameter :: NCHAR=NWORD*NCHARPERWORD
      integer, external :: f_crc32
      character(len=NCHAR) :: name_S
      integer :: crc, lbuf, n, i0, i1
      integer :: buf(NCHAR)  !#deliberatly have bigger buffer to prevent crc32 outof bound
      !---------------------------------------------------------------------
      if (len_trim(F_name_S) == 0) then
         F_id = 0
         return
      endif
      crc = 0
      lbuf = size(buf)
      name_S = F_name_S
      buf = 0
      do n=1, NWORD
         i0 = 1+(n-1)*NCHARPERWORD
         i1 = n*NCHARPERWORD
         buf(n) = transfer(name_S(i0:i1), buf(n))
      enddo
      F_id = abs(f_crc32(crc, buf, lbuf))
      !---------------------------------------------------------------------
      return
   end function priv_id


   function priv_name(F_name_S, F_id) result(F_id_S)
      implicit none
      character(len=*), intent(in) :: F_name_S
      integer, intent(out), optional :: F_id
      character(len=GMM_MAXNAMELENGTH) :: F_id_S
      integer :: id
      !---------------------------------------------------------------------
      id = priv_id(F_name_S)
      if (present(F_id)) F_id = id
      write(F_id_S, '(i16.16)') id
      return
      !---------------------------------------------------------------------
      return
   end function priv_name


end module vgrid_wb
