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
module fstmpio_rdhint_mod
   use iso_c_binding
   use rpn_comm_itf_mod
   use vGrid_Descriptors
   use vgrid_ov, only: vgrid_nullify
   use cmcdate_mod
   use ezgrid_mod
   use fst_mod
   use hinterp4yy_mod
   use ptopo_utils
   use vgrid_wb
   implicit none
   private
   !@objective 
   !@author  Stephane Chamberland, 2017-04
   !@description
   ! Public functions
   public :: fstmpio_find, fstmpio_getmeta, fstmpio_get_vgrid, fstmpio_rdhint
   public :: fstmpio_find_0, fstmpio_find_vect, fstmpio_find_3d_0, fstmpio_find_3d_vect
   public :: fstmpio_rdhint_3d_r4, fstmpio_rdhint_3d_r4_vect, fstmpio_get_hgridid
 
   ! Public constants
   public :: FST_FIND_LT, FST_FIND_LE, FST_FIND_NEAR, FST_FIND_GE, FST_FIND_GT
   public :: FST_FIND_DIAG_T, FST_FIND_DIAG_M
!@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

   interface fstmpio_find
      module procedure fstmpio_find_0
      module procedure fstmpio_find_vect
      module procedure fstmpio_find_3d_0
      module procedure fstmpio_find_3d_vect
   end interface

   interface fstmpio_rdhint
      module procedure fstmpio_rdhint_3d_r4
      module procedure fstmpio_rdhint_3d_r4_vect
   end interface

   integer, parameter :: CHARPERBYTE = 4

contains

   !# Find is done by all ptopo_io_ipe >= 0, give list of keys
   !# Pool N keys to read + hinterp in a 3d buffer... needs list of keys, varname, hinterp
   !# Special case for vector fields, needs vectorial hinterp

   !/@
   function fstmpio_find_0(F_key, F_fileid, F_nomvar, &
        F_datev, F_ip1, F_ip2, F_ip3, F_datevfuzz, F_fuzzopr, F_typvar_S) &
        result(F_istat)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      integer, intent(out) :: F_key
      character(len=*), intent(in) :: F_nomvar
      integer, intent(inout) :: F_datev
      integer, intent(in) :: F_fileid, F_ip1, F_ip2, F_ip3
      integer, intent(in), optional :: F_datevfuzz, F_fuzzopr
      character(len=*), intent(in), optional :: F_typvar_S
      !@return
      integer :: F_istat
      !@author
      !@revision
      !@/
      logical :: isiomaster_L, isiope_L
      integer :: datevfuzz, fuzzopr, istat, mydata(2), mysize, &
           comm_ipe_io_master
      character(len=128) :: msg_S, communicator_S
      character(len=12) :: nomvar_S
      character(len=2) :: typvar_s
      !---------------------------------------------------------------------
      nomvar_S = F_nomvar
      write(msg_S, '(a, 3i10, a)') '(fstmpio) Find, looking for: '// nomvar_S(1:4)// &
           ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=', F_ip1, F_ip2, F_ip3, ']'
      call msg(MSG_DEBUG, msg_S)

      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      F_key = RMN_OK
      mydata = (/F_key, F_datev/)

      if (isiomaster_L) then
         datevfuzz = 0
         fuzzopr = FST_FIND_NEAR
         typvar_S = RMN_ANY_TYP
         if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
         if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
         if (present(F_typvar_S)) typvar_S = F_typvar_S
         F_key = fst_find(F_fileid, F_nomvar, F_datev, F_ip1, F_ip2, F_ip3, &
              datevfuzz, fuzzopr, typvar_S)
         mydata = (/F_key, F_datev/)
      endif

      mysize = size(mydata)
      call rpn_comm_bcast(mydata, mysize, RPN_COMM_INTEGER, comm_ipe_io_master, &
           trim(communicator_S), istat)

      if (.not.isiomaster_L) then
         F_key = mydata(1)
         F_datev = mydata(2)
      endif

      if (RMN_IS_OK(F_key)) then
         write(msg_S, '(a, 3i10, a, i12, a)') '(fstmpio) Found: '//nomvar_S(1:4)// &
              ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=', F_ip1, F_ip2, F_ip3, &
              '] [key=', F_key, ']'
      else
         write(msg_S, '(a, 3i10, a, i12, a)') '(fstmpio) Not Found: '//nomvar_S(1:4)// &
              ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=', F_ip1, F_ip2, F_ip3, &
              '] [key=', F_key, ']'
      endif
      F_istat = F_key
      call msg(MSG_DEBUG, msg_S)
      !---------------------------------------------------------------------
      return
   end function fstmpio_find_0


   !/@
   function fstmpio_find_vect(F_key1, F_key2, F_fileid, F_nomvar1, F_nomvar2, &
        F_datev, F_ip1, F_ip2, F_ip3, F_datevfuzz, F_fuzzopr, F_typvar_S) &
        result(F_istat)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      integer, intent(out) :: F_key1, F_key2
      character(len=*), intent(in) :: F_nomvar1, F_nomvar2
      integer, intent(inout) :: F_datev
      integer, intent(in) :: F_fileid, F_ip1, F_ip2, F_ip3
      integer, intent(in), optional :: F_datevfuzz, F_fuzzopr
      character(len=*), intent(in), optional :: F_typvar_S
      !@return
      integer :: F_istat
      !@author
      !@revision
      !@/
      logical :: isiomaster_L, isiope_L
      integer :: datevfuzz, fuzzopr, istat, mydata(4), mysize, &
           comm_ipe_io_master
      character(len=128) :: msg_S, communicator_S
      character(len=12) :: nomvar1_S, nomvar2_S
      character(len=2) :: typvar_s
      !---------------------------------------------------------------------
      nomvar1_S = F_nomvar1
      nomvar2_S = F_nomvar2
      write(msg_S, '(a,3i10,a)') '(fstmpio) Find, looking for: '//nomvar1_S(1:4)// &
           ' + '//nomvar2_S(1:4)//' [datev='//cmcdate_toprint(F_datev)// &
           '] [ip123=', F_ip1, F_ip2, F_ip3, ']'
      call msg(MSG_DEBUG, msg_S)

      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      F_istat = RMN_OK
      F_key1 = RMN_OK
      F_key2 = RMN_OK
      mydata = (/F_istat, F_key1, F_key2, F_datev/)

      if (isiomaster_L) then
         datevfuzz = 0
         fuzzopr = FST_FIND_NEAR
         typvar_S = RMN_ANY_TYP
         if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
         if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
         if (present(F_typvar_S)) typvar_S = F_typvar_S
         F_istat = fst_find(F_key1, F_key2, F_fileid, F_nomvar1, F_nomvar2, &
              F_datev, F_ip1, F_ip2, F_ip3, datevfuzz, fuzzopr, typvar_S)
         mydata = (/F_istat, F_key1, F_key2, F_datev/)
      endif

      mysize = size(mydata)
      call rpn_comm_bcast(mydata, mysize, RPN_COMM_INTEGER, comm_ipe_io_master, &
           trim(communicator_S), istat)

      if (.not.isiomaster_L) then
         F_istat = min(mydata(1), RMN_OK)
         F_key1 = min(mydata(2), RMN_OK)
         F_key2 = min(mydata(3), RMN_OK)
         F_datev = mydata(4)
      endif

      if (RMN_IS_OK(F_istat)) then
         write(msg_S, '(a,3i10,a,i12,i12,a)') '(fstmpio) Found: '//nomvar1_S(1:4)// &
              ' + '//nomvar2_S(1:4)//' [datev='//cmcdate_toprint(F_datev)// &
              '] [ip123=', F_ip1, F_ip2, F_ip3, '] [keys=', F_key1, F_key2, ']'
      else
         write(msg_S, '(a,3i10,a,i12,i12,a)') '(fstmpio) Not Found: '//nomvar1_S(1:4)// &
              ' + '//nomvar2_S(1:4)//' [datev='//cmcdate_toprint(F_datev)//'] [ip123=', &
              F_ip1, F_ip2, F_ip3, '] [keys=', F_key1, F_key2, ']'
      endif
      call msg(MSG_DEBUG, msg_S)
      !---------------------------------------------------------------------
      return
   end function fstmpio_find_vect


   !/@
   function fstmpio_find_3d_0(F_keys1, F_fileid, F_nomvar1, &
        F_datev, F_ip1s, F_ip2, F_ip3, F_datevfuzz, F_fuzzopr, F_typvar_S, &
        F_vgrid, F_lvltyp_S) result(F_nkeys)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      integer, pointer :: F_keys1(:)  !# out
      integer, intent(in) :: F_fileid
      character(len=*), intent(in) :: F_nomvar1
      integer, intent(inout) :: F_datev
      integer, pointer :: F_ip1s(:)  !# inout
      integer, intent(in) :: F_ip2, F_ip3
      integer, intent(in), optional :: F_datevfuzz, F_fuzzopr
      character(len=*), intent(in), optional :: F_typvar_S  !#TODO:inout
      type(vgrid_descriptor), intent(out), optional :: F_vgrid
      character(len=*), intent(out), optional :: F_lvltyp_S
      !@return
      integer :: F_nkeys
      !@author
      !@revision
      !@/
      integer, parameter :: NMAX = 2000
      logical :: isiomaster_L, isiope_L
      integer :: datevfuzz, fuzzopr, istat, mydata(NMAX), mysize, lvltyp, nip1s, &
           itype, comm_ipe_io_master
      character(len=128) :: msg_S, ip1_S, communicator_S
      character(len=12) :: nomvar_S, lvltyp_S
      character(len=2) :: typvar_s
      type(vgrid_descriptor) :: vgrid
      !---------------------------------------------------------------------
      nomvar_S = F_nomvar1
      ip1_S = '-1'
      if (associated(F_ip1s)) write(ip1_S, '(i10,a)') F_ip1s(1), ', ...'
      write(msg_S, '(a, 2i10, a)') '(fstmpio) Find, looking for: '//nomvar_S(1:4)// &
           ' [datev='//cmcdate_toprint(F_datev)//'] [ip123='//trim(ip1_S), F_ip2, F_ip3, ']'
      call msg(MSG_DEBUG, msg_S)

      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      F_nkeys = RMN_ERR
      if (present(F_vgrid)) call vgrid_nullify(F_vgrid) 

      if (associated(F_keys1)) F_keys1 = RMN_OK
      mydata  = RMN_OK
      mydata(1:4) = (/F_datev, 0, 0, 0/)
      IF_MASTER: if (isiomaster_L) then
         datevfuzz = 0
         fuzzopr = FST_FIND_NEAR
         typvar_S = RMN_ANY_TYP
         lvltyp_S = ' '
         if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
         if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
         if (present(F_typvar_S)) typvar_S = F_typvar_S
         if (present(F_lvltyp_S) .or. present(F_vgrid)) then
            F_nkeys = fst_find_3d_0(F_keys1, F_fileid, nomvar_S, F_datev,  &
                 F_ip1s, F_ip2, F_ip3, datevfuzz, fuzzopr, typvar_S, &
                 vgrid, lvltyp_S)
         else
            F_nkeys = fst_find_3d_0(F_keys1, F_fileid, nomvar_S, F_datev,  &
                 F_ip1s, F_ip2, F_ip3, datevfuzz, fuzzopr, typvar_S)
         endif
         lvltyp = 0
         select case(lvltyp_S)
         case('SFC')
            lvltyp = 1
         case('M')
            lvltyp = 2
         case('T')
            lvltyp = 3
         end select
!!$         if (F_nkeys > size(mydata)/2-4) error
         if (F_nkeys > 0 .and. associated(F_keys1) .and. associated(F_ip1s)) then
            mydata(1:4) = (/F_datev, F_nkeys, size(F_ip1s), lvltyp/)
            mydata(5:5+F_nkeys-1) = F_keys1(1:F_nkeys)
            if (size(F_ip1s) > 0) then
               mydata(NMAX/2:NMAX/2+size(F_ip1s)-1) = F_ip1s(:)
            endif
         else
            mydata(1:4) = (/F_datev, F_nkeys, 0, 0/)
         endif
      endif IF_MASTER

      mysize = size(mydata)
      call rpn_comm_bcast(mydata, mysize, RPN_COMM_INTEGER, comm_ipe_io_master, &
           trim(communicator_S), istat)

      IF_SLAVE: if (.not.isiomaster_L) then
         F_datev = mydata(1)
         F_nkeys = mydata(2)
         nip1s   = mydata(3)
         lvltyp_S = ' '
         select case(mydata(4))
         case(1)
            lvltyp_S = 'SFC'
         case(2)
            lvltyp_S = 'M'
         case(3)
            lvltyp_S = 'T'
         end select
 !!$         if (nkeys > size(mydata)/2-4) error
         if (F_nkeys > 0) then
            if (associated(F_keys1)) then
               if (size(F_keys1) < F_nkeys) then
                  call msg(MSG_WARNING, '(fstmpio) find - provided key list too small')
                  F_nkeys = RMN_ERR
               endif
            else
               allocate(F_keys1(F_nkeys), stat=istat)
            endif
            F_keys1(1:F_nkeys) = mydata(5:5+F_nkeys-1)
         endif
         if (nip1s > 0) then
            if (.not.associated(F_ip1s)) then
               allocate(F_ip1s(nip1s), stat=istat)
            endif
            if (size(F_ip1s) < nip1s) then
               call msg(MSG_WARNING, '(fstmpio) find - provided ip1 list too small')
!!$               F_nkeys = RMN_ERR
            else
               F_ip1s(1:nip1s) = mydata(NMAX/2:NMAX/2+nip1s-1)
            endif
         endif
      endif IF_SLAVE

      if (present(F_lvltyp_S)) F_lvltyp_S = lvltyp_S
      if (present(F_vgrid) .and. any(lvltyp_S == (/'M', 'T'/)) .and. &
           associated(F_ip1s)) then
         itype = 0
         istat = vgrid_wb_bcast(vgrid, F_ip1s, itype, lvltyp_S, communicator_S, &
              comm_ipe_io_master)
         F_vgrid = vgrid
      endif

      ip1_S = '-1'
      if (associated(F_ip1s)) write(ip1_S, '(i10,a)') F_ip1s(1), ', ...'
      if (F_nkeys > 0) then
         write(msg_S, '(a, 2i10, a, i12, a)') '(fstmpio) Found: '//nomvar_S(1:4)// &
              ' [datev='//cmcdate_toprint(F_datev)//'] [ip123='//trim(ip1_S), &
              F_ip2, F_ip3, '] [key=', F_keys1(1), ']'
      else
         write(msg_S, '(a, 2i10, a, i12, a)') '(fstmpio) Not Found: '//nomvar_S(1:4)// &
              ' [datev='//cmcdate_toprint(F_datev)//'] [ip123='//trim(ip1_S), &
              F_ip2, F_ip3, '] [key=', -1, ']'
      endif
      call msg(MSG_DEBUG, msg_S)
      !---------------------------------------------------------------------
      return
   end function fstmpio_find_3d_0


   !/@
   function fstmpio_find_3d_vect(F_keys1, F_keys2, F_fileid, F_nomvar1, F_nomvar2, &
        F_datev, F_ip1s, F_ip2, F_ip3, F_datevfuzz, F_fuzzopr, F_typvar_S, &
        F_vgrid, F_lvltyp_S) result(F_nkeys)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      integer, pointer :: F_keys1(:), F_keys2(:)  !# out
      integer, intent(in) :: F_fileid
      character(len=*), intent(in) :: F_nomvar1, F_nomvar2
      integer, intent(inout) :: F_datev
      integer, pointer :: F_ip1s(:)  !# inout
      integer, intent(in) :: F_ip2, F_ip3
      integer, intent(in), optional :: F_datevfuzz, F_fuzzopr
      character(len=*), intent(in), optional :: F_typvar_S
      type(vgrid_descriptor), intent(out), optional :: F_vgrid
      character(len=*), intent(out), optional :: F_lvltyp_S
      !@return
      integer :: F_nkeys
      !@author
      !@revision
      !@/
      integer, parameter :: NMAX = 3000
      logical :: isiomaster_L, isiope_L
      integer :: datevfuzz, fuzzopr, istat, mydata(NMAX), mysize, lvltyp, nip1s, &
           itype, comm_ipe_io_master
      character(len=128) :: msg_S, ip1_S, communicator_S
      character(len=12) :: nomvar1_S, nomvar2_S, lvltyp_S
      character(len=2) :: typvar_s
      type(vgrid_descriptor) :: vgrid
      !---------------------------------------------------------------------
      nomvar1_S = F_nomvar1
      nomvar2_S = F_nomvar2
      ip1_S = '-1'
      if (associated(F_ip1s)) write(ip1_S, '(i10,a)') F_ip1s(1), ', ...'
      write(msg_S, '(a, 2i10, a)') '(fstmpio) Find, looking for: '//nomvar1_S(1:4)// &
           ' + '//nomvar2_S(1:4)//' [datev='//cmcdate_toprint(F_datev)// &
           '] [ip123='//trim(ip1_S), F_ip2, F_ip3, ']'
      call msg(MSG_DEBUG, msg_S)

      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      F_nkeys = RMN_ERR
      if (present(F_vgrid)) call vgrid_nullify(F_vgrid) 

      if (associated(F_keys1)) F_keys1 = RMN_OK
      if (associated(F_keys2)) F_keys2 = RMN_OK
      mydata = RMN_OK
      mydata(1:4) = (/F_datev, 0, 0, 0/)
      IF_MASTER: if (isiomaster_L) then
         datevfuzz = 0
         fuzzopr = FST_FIND_NEAR
         typvar_S = RMN_ANY_TYP
         lvltyp_S = ' '
         if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
         if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
         if (present(F_typvar_S)) typvar_S = F_typvar_S
         if (present(F_lvltyp_S) .or. present(F_vgrid)) then
            F_nkeys = fst_find(F_keys1, F_keys2, F_fileid, nomvar1_S, &
                 nomvar2_S, F_datev, F_ip1s, F_ip2, F_ip3, datevfuzz, fuzzopr, &
                 typvar_S, vgrid, lvltyp_S)
         else
            F_nkeys = fst_find(F_keys1, F_keys2, F_fileid, nomvar1_S, &
                 nomvar2_S, F_datev, F_ip1s, F_ip2, F_ip3, datevfuzz, fuzzopr, &
                 typvar_S)
         endif
         lvltyp = 0
         select case(lvltyp_S)
         case('SFC')
            lvltyp = 1
         case('M')
            lvltyp = 2
         case('T')
            lvltyp = 3
         end select
 !!$         if (F_nkeys > size(mydata)/2-4) error
         if (F_nkeys > 0 .and. associated(F_keys1) .and. associated(F_keys2) &
              .and. associated(F_ip1s)) then
            mydata(1:4) = (/F_datev, F_nkeys, size(F_ip1s), lvltyp/)
            mydata(5:5+F_nkeys-1) = F_keys1(1:F_nkeys)
            mydata(NMAX/3:NMAX/3+F_nkeys-1) = F_keys2(1:F_nkeys)
            if (size(F_ip1s) > 0) then
               mydata(2*NMAX/3:2*NMAX/3+size(F_ip1s)-1) = F_ip1s(:)
            endif
         else
            mydata(1:4) = (/F_datev, F_nkeys, 0, 0/)
         endif
      endif IF_MASTER

      mysize = size(mydata)
      call rpn_comm_bcast(mydata, mysize, RPN_COMM_INTEGER, comm_ipe_io_master, &
           trim(communicator_S), istat)

      IF_SLAVE: if (.not.isiomaster_L) then
         F_datev = mydata(1)
         F_nkeys = mydata(2)
         nip1s   = mydata(3)
         lvltyp_S = ' '
         select case(mydata(4))
         case(1)
            lvltyp_S = 'SFC'
         case(2)
            lvltyp_S = 'M'
         case(3)
            lvltyp_S = 'T'
         end select
 !!$         if (nkeys > size(mydata)/2-4) error
         if (F_nkeys > 0) then
            if (associated(F_keys1)) then
               if (size(F_keys1) < F_nkeys) then
                  call msg(MSG_WARNING, '(fstmpio) find - provided key list too small')
                  F_nkeys = RMN_ERR
               endif
            else
               allocate(F_keys1(F_nkeys), stat=istat)
            endif
            F_keys1(1:F_nkeys) = mydata(5:5+F_nkeys-1)
            if (associated(F_keys2)) then
               if (size(F_keys2) < F_nkeys) then
                  call msg(MSG_WARNING, '(fstmpio) find - provided key list too small')
                  F_nkeys = RMN_ERR
               endif
            else
               allocate(F_keys2(F_nkeys), stat=istat)
            endif
            F_keys2(1:F_nkeys) = mydata(NMAX/3:NMAX/3+F_nkeys-1)
         endif
         if (nip1s > 0) then
            if (.not.associated(F_ip1s)) then
               allocate(F_ip1s(nip1s), stat=istat)
            endif
            if (size(F_ip1s) < nip1s) then
               call msg(MSG_WARNING, '(fstmpio) find - provided ip1 list too small')
!!$               F_nkeys = RMN_ERR
            else
               F_ip1s(1:nip1s) = mydata(2*NMAX/3:2*NMAX/3+nip1s-1)
            endif
         endif
      endif IF_SLAVE

      if (present(F_lvltyp_S)) F_lvltyp_S = lvltyp_S
      if (present(F_vgrid) .and. any(lvltyp_S == (/'M', 'T'/)) .and. &
           associated(F_ip1s)) then
         itype = 0
         istat = vgrid_wb_bcast(vgrid, F_ip1s, itype, lvltyp_S, communicator_S, &
              comm_ipe_io_master)
         F_vgrid = vgrid
      endif

      ip1_S = '-1'
      if (associated(F_ip1s)) write(ip1_S, '(i10,a)') F_ip1s(1), ', ...'
      if (F_nkeys > 0) then
         write(msg_S, '(a, 2i10, a, i12, a)') '(fstmpio) Found: '//nomvar1_S(1:4)// &
              ' + '//nomvar2_S(1:4)//' [datev='//cmcdate_toprint(F_datev)// &
              '] [ip123='//trim(ip1_S), F_ip2, F_ip3, '] [key=', F_keys1(1), ']'
      else
         write(msg_S, '(a, 2i10, a, i12, a)') '(fstmpio) Not Found: '//nomvar1_S(1:4)// &
              ' + '//nomvar2_S(1:4)//' [datev='//cmcdate_toprint(F_datev)// &
              '] [ip123='//trim(ip1_S), F_ip2, F_ip3, '] [key=', -1, ']'
      endif
      call msg(MSG_DEBUG, msg_S)
      !---------------------------------------------------------------------
      return
   end function fstmpio_find_3d_vect


   !/@
   function fstmpio_rdhint_3d_r4(F_data1, F_status, F_keys1, F_hintlist_S, &
        F_fileids, F_rpncomm_gridid, F_outgridid, F_coregridid, F_realloc_L) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real, pointer :: F_data1(:,:,:)
      integer, intent(out) :: F_status(:)
      integer, intent(in) :: F_keys1(:)
      character(len=*), intent(in) :: F_hintlist_S(:)
      integer, intent(in) :: F_fileids(:)
      integer, intent(in) :: F_rpncomm_gridid
      integer, intent(in) :: F_outgridid
      integer, intent(in), optional :: F_coregridid
      logical,intent(in), optional :: F_realloc_L
      !#TODO: alternate itf with hgridis_S so we have min, max as well
      !@author
      !@return
      integer :: F_istat
      !@/
      real, target :: hdata0(1,1,1)
      real, pointer ::  hdata1(:,:,:)
      integer, pointer :: zlist(:), status(:)
      logical, pointer :: zlist_o(:)
      integer :: istat, coregridid, &
           nkeys, lnkeys, lnkeys2, nfids, irest, igkey, ilkey, kstart, kend, &
           gni, gnj, mini, maxi, minj, maxj, &
           starti(ptopo_grid_npex), counti(ptopo_grid_npex), &
           startj(ptopo_grid_npey), countj(ptopo_grid_npey), nhint, k, &
           comm_ipe_io_master
      character(len=64) :: communicator_S
      logical :: realloc_L, isiomaster_L, isiope_L
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fstmpio) rdhint_3d_r4 [BEGIN]')
      F_istat = RMN_ERR
      F_status = RMN_ERR

      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      realloc_L = .false.
      if (present(F_realloc_L)) realloc_L = F_realloc_L

      nkeys   = size(F_keys1)
      nhint   = size(F_hintlist_S)
      nfids   = size(F_fileids)
      lnkeys  = 0
      lnkeys2 = 0
      kstart  = 1

      if (size(F_status) < nkeys) then
         call msg(MSG_WARNING, '(fst) rdhint: status array too small')
         return
      endif

      F_istat = RMN_OK
      nullify(hdata1)
      IF_IOPE: if (isiope_L) then

         coregridid = F_outgridid
         if (present(F_coregridid)) coregridid = F_coregridid

         !#TODO: mv code below to s/r to be called from _vect as well
         lnkeys  = nkeys / ptopo_io_npe
         lnkeys2 = (nkeys + ptopo_io_npe - 1) / ptopo_io_npe
         irest   = nkeys - lnkeys * ptopo_io_npe
         igkey   = mod(ptopo_io_ipe, ptopo_io_npe)
         kstart  = 1 + lnkeys * igkey
         if (igkey < irest) then
            lnkeys = lnkeys + 1
            kstart = kstart + igkey
         else
            kstart = kstart + irest
         endif
         kend = min(kstart+lnkeys-1, nkeys)

         F_istat = fst_rdhint(hdata1, F_status(kstart:kend), &
              F_keys1(kstart:kend), &
              F_hintlist_S(min(kstart, nhint):min(kend, nhint)), &
              F_fileids(min(kstart, nfids):min(kend, nfids)), &
              F_outgridid, F_coregridid, realloc_L)

      else
         hdata1 => hdata0
         hdata1 = 0.
      endif IF_IOPE

      allocate(status(size(F_status)))
      status = F_status
      call rpn_comm_allreduce( &
           status, F_status, size(F_status), &
           RPN_COMM_INTEGER, RPN_COMM_MAX, RPN_COMM_GRID, istat)
      deallocate(status, stat=istat)
      F_istat = maxval(F_status)
      call collect_error(F_istat)
!!$      if (.not.RMN_IS_OK(F_istat)) return
      if (all(F_status < 0)) return

      !#TODO: mv code below to s/r to be called twice from _vect as well

      !#TODO: ptopo_io_grid_id, is it costly to call rpn_comm_create_2dgrid every time?
      ! (Not costly) Can we release a grid_id? (not sure how, MAX_GRIDS=64)
      istat = rpn_comm_get_2dgrid(F_rpncomm_gridid, &
           ptopo_grid_npex, ptopo_grid_npey, &
           gni, gnj, mini, maxi, minj, maxj, starti, counti, startj, countj)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in rpn_comm_get_2dgrid')
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif

      istat = fst_checkalloc2(F_data1, mini, maxi, minj, maxj, nkeys, &
           realloc_L)
      if (.not.RMN_IS_OK(istat)) then
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif

      allocate(zlist(max(1,lnkeys2)), zlist_o(nkeys))
      zlist = -1
      do ilkey = 1, lnkeys
         zlist(ilkey) =  kstart + ilkey - 1
      end do
      zlist_o = .false.
      call rpn_comm_bcast(lnkeys2, 1, RPN_COMM_INTEGER, &
           comm_ipe_io_master, trim(communicator_S), istat) !#TODO: reduce max?
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in rpn_comm_bcast')
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif
      F_data1 = 0. !# Init to 0 for case when maxij > lnij
      istat = rpn_comm_shuf_ezdist(ptopo_io_setno, F_rpncomm_gridid, &
           hdata1, lnkeys2, F_data1, nkeys, zlist, zlist_o)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in rpn_comm_shuf_ezdist')
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif
      do k=1, nkeys
         if (.not.zlist_o(k)) F_status(k) = RMN_ERR
      end do

      if (associated(hdata1) .and. .not.associated(hdata1, hdata0)) &
           deallocate(hdata1, stat=istat)
      if (associated(zlist)) deallocate(zlist, stat=istat)
      if (associated(zlist_o)) deallocate(zlist_o, stat=istat)
      call msg(MSG_DEBUG, '(fstmpio) rdhint_3d_r4 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpio_rdhint_3d_r4


   !/@
   function fstmpio_rdhint_3d_r4_vect(F_data1, F_data2, F_status, &
        F_keys1, F_keys2, F_hintlist_S, &
        F_fileids, F_rpncomm_gridid, F_outgridid, F_coregridid, F_realloc_L) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real, pointer :: F_data1(:,:,:), F_data2(:,:,:)
      integer, intent(out) :: F_status(:)
      integer, intent(in) :: F_keys1(:), F_keys2(:)
      character(len=*), intent(in) :: F_hintlist_S(:)
      integer, intent(in) :: F_fileids(:)
      integer, intent(in) :: F_rpncomm_gridid
      integer, intent(in) :: F_outgridid
      integer, intent(in), optional :: F_coregridid
      logical,intent(in), optional :: F_realloc_L
      !#TODO: alternate itf with hgridis_S so we have min, max as well
      !@author
      !@return
      integer :: F_istat
      !@/
      real, target :: hdata0(1,1,1)
      real, pointer :: hdata1(:,:,:), hdata2(:,:,:)
      integer, pointer :: zlist(:), status(:)
      logical, pointer :: zlist_o(:)
      integer :: istat, coregridid, &
           nkeys, lnkeys, lnkeys2, nfids, irest, igkey, ilkey, kstart, kend, &
           gni, gnj, mini, maxi, minj, maxj, &
           starti(ptopo_grid_npex), counti(ptopo_grid_npex), &
           startj(ptopo_grid_npey), countj(ptopo_grid_npey), nhint, k, &
           comm_ipe_io_master
      character(len=64) :: communicator_S
      logical :: realloc_L, isiomaster_L, isiope_L
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fstmpio) rdhint_3d_r4_vect [BEGIN]')
      F_istat = RMN_ERR
      F_status = RMN_ERR

      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      realloc_L = .false.
      if (present(F_realloc_L)) realloc_L = F_realloc_L

      nkeys   = size(F_keys1)
      nhint   = size(F_hintlist_S)
      nfids   = size(F_fileids)
      lnkeys  = 0
      lnkeys2 = 0
      kstart  = 1

      if (size(F_status) < nkeys) then
         call msg(MSG_WARNING, '(fst) rdhint: status array too small')
         return
      endif

      F_istat = RMN_OK
      nullify(hdata1, hdata2)
      IF_IOPE: if (isiope_L) then

         coregridid = F_outgridid
         if (present(F_coregridid)) coregridid = F_coregridid

         !#TODO: mv code below to s/r to be called from _vect as well
         lnkeys  = nkeys / ptopo_io_npe
         lnkeys2 = (nkeys + ptopo_io_npe - 1) / ptopo_io_npe
         irest   = nkeys - lnkeys * ptopo_io_npe
         igkey   = mod(ptopo_io_ipe, ptopo_io_npe)
         kstart  = 1 + lnkeys * igkey
         if (igkey < irest) then
            lnkeys = lnkeys + 1
            kstart = kstart + igkey
         else
            kstart = kstart + irest
         endif
         kend = min(kstart+lnkeys-1, nkeys)

         F_istat = fst_rdhint(hdata1, hdata2, F_status(kstart:kend), &
              F_keys1(kstart:kend), F_keys2(kstart:kend),&
              F_hintlist_S(min(kstart, nhint):min(kend, nhint)), &
              F_fileids(min(kstart, nfids):min(kend, nfids)), &
              F_outgridid, F_coregridid, realloc_L)

      else
         hdata1 => hdata0
         hdata2 => hdata0
         hdata1 = 0.
      endif IF_IOPE

      allocate(status(size(F_status)))
      status = F_status
      call rpn_comm_allreduce( &
           status, F_status, size(F_status), &
           RPN_COMM_INTEGER, RPN_COMM_MAX, RPN_COMM_GRID, istat)
      deallocate(status, stat=istat)
      F_istat = maxval(F_status)
      call collect_error(F_istat)
!!$      if (.not.RMN_IS_OK(F_istat)) return
      if (all(F_status < 0)) return

      !#TODO: mv code below to s/r to be called twice from _vect as well

      !#TODO: ptopo_io_grid_id, is it costly to call rpn_comm_create_2dgrid every time?
      ! (Not costly) Can we release a grid_id? (not sure how, MAX_GRIDS=64)
      istat = rpn_comm_get_2dgrid(F_rpncomm_gridid, &
           ptopo_grid_npex, ptopo_grid_npey, &
           gni, gnj, mini, maxi, minj, maxj, starti, counti, startj, countj)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in rpn_comm_get_2dgrid')
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif

      istat = fst_checkalloc2(F_data1, mini, maxi, minj, maxj, nkeys, &
           realloc_L)
      istat = min(fst_checkalloc2(F_data2, mini, maxi, minj, maxj, nkeys, &
           realloc_L), istat)
      if (.not.RMN_IS_OK(istat)) then
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif
 
      allocate(zlist(max(1,lnkeys2)), zlist_o(nkeys))
      zlist = -1
      do ilkey = 1, lnkeys
         zlist(ilkey) =  kstart + ilkey - 1
      end do

      call rpn_comm_bcast(lnkeys2, 1, RPN_COMM_INTEGER, &
           comm_ipe_io_master, trim(communicator_S), istat) !#TODO: reduce max?
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in rpn_comm_bcast')
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif

      zlist_o = .false.
      F_data1 = 0. !# Init to 0 for case when maxij > lnij
      istat = rpn_comm_shuf_ezdist(ptopo_io_setno, F_rpncomm_gridid, &
           hdata1, lnkeys2, F_data1, nkeys, zlist, zlist_o)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in rpn_comm_shuf_ezdist')
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif
      do k=1, nkeys
         if (.not.zlist_o(k)) F_status(k) = RMN_ERR
      end do

      zlist_o = .false.
      F_data2 = 0. !# Init to 0 for case when maxij > lnij
      istat = rpn_comm_shuf_ezdist(ptopo_io_setno, F_rpncomm_gridid, &
           hdata2, lnkeys2, F_data2, nkeys, zlist, zlist_o)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in rpn_comm_shuf_ezdist')
         F_istat = RMN_ERR
         F_status = RMN_ERR
      endif
      do k=1, nkeys
         if (.not.zlist_o(k)) F_status(k) = RMN_ERR
      end do

      if (associated(hdata1) .and. .not.associated(hdata1, hdata0)) &
           deallocate(hdata1, stat=istat)
      if (associated(hdata2) .and. .not.associated(hdata2, hdata0)) &
           deallocate(hdata2, stat=istat)
      if (associated(zlist)) deallocate(zlist, stat=istat)
      if (associated(zlist_o)) deallocate(zlist_o, stat=istat)
      call msg(MSG_DEBUG, '(fstmpio) rdhint_3d_r4_vect [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpio_rdhint_3d_r4_vect


   !/@
   function fstmpio_getmeta(F_key, F_nomvar_S, F_dateo, F_deet, F_npas, &
        F_ip1, F_ip2, F_ip3, F_etiket_S, F_typvar_S, &
        F_ni, F_nj, F_nk, F_nbits, F_datyp, &
        F_grtyp_S, F_ig1, F_ig2, F_ig3, F_ig4) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer, intent(in) :: F_key
      character(len=*), intent(out), optional :: F_nomvar_S, F_etiket_S, F_typvar_S, F_grtyp_S
      integer, intent(out), optional :: F_dateo, F_deet, F_npas, F_ip1, F_ip2, F_ip3, &
           F_ni, F_nj, F_nk, F_nbits, F_datyp, F_ig1, F_ig2, F_ig3, F_ig4
      !@author
      !@return
      integer :: F_istat
      !@/
      integer, parameter :: MYSTRLEN = 32
      integer, parameter :: NVAR_INT = 16
      integer, parameter :: NVAR_STR = 3
      logical :: isiomaster_L, isiope_L
      integer :: istat, ii, nn, mysize, comm_ipe_io_master, &
           buffer(NVAR_INT+NVAR_STR*MYSTRLEN/CHARPERBYTE)
      character(len=CHARPERBYTE) :: str4_S
      character(len=MYSTRLEN) :: typvar_S, nomvar_S, etiket_S, grtyp_S, communicator_S
      integer :: ni1, nj1, nk1, &
           dateo, deet, npas, nbits, datyp, ip1, ip2, ip3, &
           ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3
      ! ---------------------------------------------------------------------
      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      F_istat = RMN_OK
      buffer = 0
      typvar_S = ' '
      nomvar_S = ' '
      etiket_S = ' '
      grtyp_S  = ' '
      if (isiomaster_L) then
         if (F_key < 0) then
            F_istat = RMN_ERR
            buffer(:) = 0
            buffer(1) = F_istat
         else
            F_istat = fstprm(F_key, dateo, deet, npas, ni1, nj1, nk1, &
                 nbits, datyp, ip1, ip2, ip3, &
                 typvar_S, nomvar_S, etiket_S, &
                 grtyp_S, ig1, ig2, ig3, ig4, swa, lng, dltf, &
                 ubc, extra1, extra2, extra3)
            buffer(1) = F_istat
            buffer(2) = dateo
            buffer(3) = deet
            buffer(4) = npas
            buffer(5) = ip1
            buffer(6) = ip2
            buffer(7) = ip3
            buffer(8) = ni1
            buffer(9) = nj1
            buffer(10) = nk1
            buffer(11) = nbits
            buffer(12) = datyp
            buffer(13) = ig1
            buffer(14) = ig2
            buffer(15) = ig3
            buffer(16) = ig4
            ii = NVAR_INT
            mysize = MYSTRLEN/CHARPERBYTE
            do nn=1, mysize, CHARPERBYTE
               ii = ii + 1
               buffer(ii) = transfer(nomvar_S(nn:nn+CHARPERBYTE-1), istat)
            enddo
            do nn=1, mysize, CHARPERBYTE
               ii = ii + 1
               buffer(ii) = transfer(etiket_S(nn:nn+CHARPERBYTE-1), istat)
            enddo
            do nn=1, mysize, CHARPERBYTE
               ii = ii + 1
               buffer(ii) = transfer(typvar_S(nn:nn+CHARPERBYTE-1), istat)
            enddo
            do nn=1, mysize, CHARPERBYTE
               ii = ii + 1
               buffer(ii) = transfer(grtyp_S(nn:nn+CHARPERBYTE-1), istat)
            enddo
         endif
      endif
      mysize = size(buffer)
      call rpn_comm_bcast(buffer, mysize, RPN_COMM_INTEGER, comm_ipe_io_master, &
           trim(communicator_S), istat) !#TODO: review, costly
      F_istat = buffer(1)
      if (.not.RMN_IS_OK(F_istat)) return

      if (present(F_dateo)) F_dateo = buffer(2)
      if (present(F_deet)) F_deet = buffer(3)
      if (present(F_npas)) F_npas = buffer(4)
      if (present(F_ip1)) F_ip1 = buffer(5)
      if (present(F_ip2)) F_ip2 = buffer(6)
      if (present(F_ip3)) F_ip3 = buffer(7)
      if (present(F_ni)) F_ni = buffer(8)
      if (present(F_nj)) F_nj = buffer(9)
      if (present(F_nk)) F_nk = buffer(10)
      if (present(F_nbits)) F_nbits = buffer(11)
      if (present(F_datyp)) F_datyp = buffer(12)
      if (present(F_ig1)) F_ig1 = buffer(13)
      if (present(F_ig2)) F_ig2 = buffer(14)
      if (present(F_ig3)) F_ig3 = buffer(15)
      if (present(F_ig4)) F_ig4 = buffer(16)
      ii = NVAR_INT
      mysize = MYSTRLEN/CHARPERBYTE
      do nn=1, mysize, CHARPERBYTE
         ii = ii + 1
         nomvar_S(nn:nn+CHARPERBYTE-1) = transfer(buffer(ii), str4_S)
      enddo
      do nn=1, mysize, CHARPERBYTE
         ii = ii + 1
         etiket_S(nn:nn+CHARPERBYTE-1) = transfer(buffer(ii), str4_S)
      enddo
      do nn=1, mysize, CHARPERBYTE
         ii = ii + 1
         typvar_S(nn:nn+CHARPERBYTE-1) = transfer(buffer(ii), str4_S)
      enddo
      do nn=1, mysize, CHARPERBYTE
         ii = ii + 1
         grtyp_S(nn:nn+CHARPERBYTE-1) = transfer(buffer(ii), str4_S)
      enddo
      if (present(F_nomvar_S)) then
         F_nomvar_S = ' '
         F_nomvar_S = trim(nomvar_S)
      endif
      if (present(F_etiket_S)) then
         F_etiket_S = ' '
         F_etiket_S = trim(etiket_S)
      endif
      if (present(F_typvar_S)) then
         F_typvar_S = ' '
         F_typvar_S = trim(typvar_S)
      endif
      if (present(F_grtyp_S)) then
         F_grtyp_S = ' '
         F_grtyp_S = trim(grtyp_S)
      endif
      ! ---------------------------------------------------------------------
      return
   end function fstmpio_getmeta


   !/@
   function fstmpio_get_hgridid(F_fileid, F_key) result(F_gridid)
      implicit none
      !@objective
      !@arguments
      integer, intent(in) :: F_fileid, F_key
      !@author
      !@return
      integer :: F_gridid
      !@/
      logical :: isiomaster_L, isiope_L
      integer :: istat, comm_ipe_io_master
      character(len=64) :: communicator_S
      ! ---------------------------------------------------------------------
      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      F_gridid = RMN_ERR
      if (isiomaster_L) then
         F_gridid = fst_get_gridid(F_fileid, F_key)
      endif
      F_gridid = ezgrid_bcast(F_gridid, communicator_S)
      ! ---------------------------------------------------------------------
      return
   end function fstmpio_get_hgridid


   !/@
   function fstmpio_get_vgrid(F_fileid, F_key, F_vgrid, F_ip1s, F_lvltyp_S) &
        result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer, intent(in) :: F_fileid, F_key
      type(vgrid_descriptor), intent(out) :: F_vgrid
      integer, pointer :: F_ip1s(:)
      character(len=*), intent(out) :: F_lvltyp_S
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: isiomaster_L, isiope_L
      integer :: istat, itype, comm_ipe_io_master
      character(len=32) :: lvltyp_S, communicator_S
      ! ---------------------------------------------------------------------
      istat = ptopo_get_io_params(isiomaster_L, isiope_L, comm_ipe_io_master, &
           communicator_S)

      F_istat = RMN_OK
      itype = 0
      F_lvltyp_S = ' '
      lvltyp_S = ' '
      call vgrid_nullify(F_vgrid)
      if (isiomaster_L) then
         F_istat = fst_get_vgrid(F_fileid, F_key, F_vgrid, F_ip1s, lvltyp_S)
      endif
      call rpn_comm_bcast(F_istat, 1, RPN_COMM_INTEGER, comm_ipe_io_master, &
           trim(communicator_S), istat)
      if (.not.RMN_IS_OK(F_istat)) return
      F_istat = vgrid_wb_bcast(F_vgrid, F_ip1s, itype, lvltyp_S, communicator_S, &
           comm_ipe_io_master)
      F_lvltyp_S = lvltyp_S
      ! ---------------------------------------------------------------------
      return
   end function fstmpio_get_vgrid

end module fstmpio_rdhint_mod
