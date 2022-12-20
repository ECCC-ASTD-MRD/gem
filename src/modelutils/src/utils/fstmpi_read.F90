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
module fstmpi_read_mod
   use iso_c_binding
   use rpn_comm_itf_mod
   use vGrid_Descriptors
   use vgrid_ov, only: vgrid_nullify
   use ezgrid_mod
   use fst_mod
   use hinterp4yy_mod
   use ptopo_utils
   use vgrid_wb
   implicit none
   private
   !@objective
   !@author  Stephane Chamberland, 2011-06
   !@description
   ! Public functions
   public :: fstmpi_find, fstmpi_read, fstmpi_getmeta, fstmpi_get_gridid, &
        fstmpi_get_hgridid, fstmpi_get_vgrid, fstmpi_rdhint, &
        fstmpi_rdhint_3d_r4, fstmpi_rdhint_3d_r4_vect

   ! Public constants
   public :: FST_FIND_LT, FST_FIND_LE, FST_FIND_NEAR, FST_FIND_GE, FST_FIND_GT
!@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

   interface fstmpi_get_gridid !# Name kept for backward compatibility
      module procedure fstmpi_get_hgridid
   end interface

   interface fstmpi_read
      module procedure fstmpi_read_3d_r4
   end interface

   interface fstmpi_rdhint
      module procedure fstmpi_rdhint_3d_r4
      module procedure fstmpi_rdhint_3d_r4_vect
   end interface

   integer,parameter :: CHARPERBYTE = 4

contains

   !/@
   function fstmpi_find(F_fileid,F_nomvar,F_datev,F_ip1,F_ip2,F_ip3, &
        F_datevfuzz,F_fuzzopr,F_typvar_S) result(F_key)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      character(len=*),intent(in) :: F_nomvar
      integer,intent(inout) :: F_datev
      integer,intent(in) :: F_fileid,F_ip1,F_ip2,F_ip3
      integer,intent(in),optional :: F_datevfuzz,F_fuzzopr
      character(len=*),intent(in),optional :: F_typvar_S
      !@return
      integer :: F_key
      !@author
      !@revision
      !@/
      integer :: datevfuzz,fuzzopr,istat,mydata(2),mysize
      character(len=128) :: msg_S
      character(len=12) :: nomvar_S
      character(len=2) :: typvar_s
      !---------------------------------------------------------------------
      nomvar_S = F_nomvar
      write(msg_S,'(a,i12,a,3i10,a)') '(fstmpi) Find, looking for: '//nomvar_S(1:4)// &
           ' [datev=',F_datev,'] [ip123=',F_ip1,F_ip2,F_ip3,']'
      call msg(MSG_DEBUG,msg_S)
      call ptopo_init_var()
      F_key = RMN_OK
      mydata = (/F_key,F_datev/)
      if (ptopo_isblocmaster_L) then
         datevfuzz = 0
         fuzzopr = FST_FIND_NEAR
         typvar_S = RMN_ANY_TYP
         if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
         if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
         if (present(F_typvar_S)) typvar_S = F_typvar_S
         F_key = fst_find(F_fileid,F_nomvar,F_datev,F_ip1,F_ip2,F_ip3, &
              datevfuzz,fuzzopr,typvar_S)
         mydata = (/F_key,F_datev/)
      endif
      mysize = size(mydata)
      call rpn_comm_bcast(mydata, mysize, RPN_COMM_INTEGER, RPN_COMM_MASTER, RPN_COMM_BLOC_COMM,istat)
!!$      print *,'fstmpi_find from bcst',istat; call flush(6)
      if (.not.ptopo_isblocmaster_L) then
         F_key = min(mydata(1),RMN_OK)
         F_datev = mydata(2)
      endif
      if (RMN_IS_OK(F_key)) then
         write(msg_S,'(a,i12,a,3i10,a,i12,a)') '(fstmpi) Found: '//nomvar_S(1:4)// &
              ' [datev=',F_datev,'] [ip123=',F_ip1,F_ip2,F_ip3,'] [key=',F_key,']'
      else
         write(msg_S,'(a,i12,a,3i10,a,i12,a)') '(fstmpi) Not Found: '//nomvar_S(1:4)// &
              ' [datev=',F_datev,'] [ip123=',F_ip1,F_ip2,F_ip3,'] [key=',F_key,']'
      endif
      call msg(MSG_DEBUG,msg_S)
      !---------------------------------------------------------------------
      return
   end function fstmpi_find


   !/@
   function fstmpi_read_3d_r4(F_key, F_data, F_fileid, F_gridid, &
        F_realloc_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_key
      real,pointer :: F_data(:,:,:)
      integer,intent(in),optional :: F_fileid
      integer,intent(out),optional :: F_gridid
!!$      character(len=*),intent(out),optional :: F_nomvar_S, F_etiket_S, F_typvar_S
!!$      integer,intent(out),optional :: F_dateo,F_deet,F_npas,F_ip1, F_ip2, F_ip3
      logical,intent(in),optional :: F_realloc_L
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: realloc_L
      integer :: istat, nijk(4), mysize
      ! ---------------------------------------------------------------------
!!$      print *,'fstmpi_read',F_key; call flush(6)
      call ptopo_init_var()
      F_istat = RMN_ERR
      if (present(F_gridid)) F_gridid = RMN_ERR
      nijk = RMN_ERR
      realloc_L = .false.
      if (present(F_realloc_L)) realloc_L = F_realloc_L

      if (ptopo_isblocmaster_L) then
         if (present(F_fileid) .and. present(F_gridid)) then
            F_istat = fst_read(F_key,F_data,F_fileid,F_gridid, &
                 F_realloc_L=realloc_L)
         else
            F_istat = fst_read(F_key,F_data, F_realloc_L=realloc_L)
         endif
         nijk(4) = F_istat
         if (RMN_IS_OK(F_istat)) nijk(1:3) = shape(F_data)
         if (any(nijk(1:3) <= 0))  nijk(4) = RMN_ERR
      endif
      mysize = size(nijk)
      call rpn_comm_bcast(nijk, mysize, RPN_COMM_INTEGER, RPN_COMM_MASTER, RPN_COMM_BLOC_COMM,istat)
      F_istat = nijk(4)
      if (.not.RMN_IS_OK(F_istat)) return

      if (.not.ptopo_isblocmaster_L) then
         if (.not.associated(F_data)) then
            allocate(F_data(nijk(1),nijk(2),nijk(3)))
         endif
         if (size(F_data,1) /= nijk(1) .or. &
              size(F_data,2) /= nijk(2)  .or. &
              size(F_data,3) /= nijk(3)) then
            F_istat = RMN_ERR
         endif
      endif
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return

      mysize = size(F_data) !=nijk(1)*nijk(2)*nijk(3)
      call rpn_comm_bcast(F_data, mysize, RPN_COMM_REAL, RPN_COMM_MASTER, RPN_COMM_BLOC_COMM,istat)
      if (present(F_fileid) .and. present(F_gridid)) then
         F_gridid = ezgrid_bcast(F_gridid,RPN_COMM_BLOC_COMM)
         F_istat = min(F_gridid,F_istat)
      endif
!!$      print *,'fstmpi_read (end)',F_key,F_istat; call flush(6)

      ! ---------------------------------------------------------------------
      return
   end function fstmpi_read_3d_r4


  !/@
   function fstmpi_rdhint_3d_r4(F_data1, F_status, F_keys1, F_hintlist_S, &
        F_fileids, F_outgridid, F_coregridid, F_realloc_L) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real, pointer :: F_data1(:,:,:)
      integer, intent(out) :: F_status(:)
      integer, intent(in) :: F_keys1(:)
      character(len=*), intent(in) :: F_hintlist_S(:)
      integer, intent(in) :: F_fileids(:)
      integer, intent(in) :: F_outgridid
      integer, intent(in), optional :: F_coregridid
      logical,intent(in), optional :: F_realloc_L
      !#TODO: alternate itf with hgridis_S so we have min, max as well
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: realloc_L
      integer :: istat, coregridid, nkeys, nhint, nfids, nij(2), ikey, fid, ingridid
      character(len=12) :: nomvar_S
      real, pointer :: indata1(:,:,:)
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fstmpi) rdhint_3d_r4 [BEGIN]')
      F_istat = RMN_ERR
      F_status = RMN_ERR

      realloc_L = .false.
      coregridid = F_outgridid
      if (present(F_realloc_L)) realloc_L = F_realloc_L
      if (present(F_coregridid)) coregridid = F_coregridid

      nkeys   = size(F_keys1)
      nhint   = size(F_hintlist_S)
      nfids   = size(F_fileids)

      if (size(F_status) < nkeys) then
         call msg(MSG_WARNING, '(fst) rdhint: status array too small')
         return
      endif

      F_istat = ezgrid_params(F_outgridid, nij)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING, '(fst) rdhint: Problem getting grid params')
         return
      endif

      F_istat = fst_checkalloc(F_data1, nij(1), nij(2), nkeys, realloc_L)
      if (.not.RMN_IS_OK(F_istat)) return

      nullify(indata1)
      DO_NKEYS: do ikey = 1, nkeys
         if (.not.RMN_IS_OK(F_keys1(ikey))) cycle

         ingridid = -1
         fid = F_fileids(min(ikey, nfids))
         istat = fstmpi_read_3d_r4(F_keys1(ikey), indata1, fid, ingridid, &
              F_realloc_L=realloc_L)
         nomvar_S = ' ' !#TODO get nomvar from fst_read or getmeta
         if (.not.RMN_IS_OK(istat)) then
            F_status(ikey) = RMN_ERR
            call msg(MSG_WARNING, '(fstmpi) rdhint: Problem reading: '//trim(nomvar_S))
            cycle
         endif
         !#TODO: allow for optional stats of read field before interp

         F_status(ikey) = hinterp4yy2d(F_data1, indata1, ikey, ingridid, &
              F_outgridid, coregridid, F_hintlist_S(min(ikey,nhint)), &
              nomvar_S)
         if (.not.RMN_IS_OK(F_status(ikey))) then
            call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in hinterp4yy2d for: '//trim(nomvar_S))
            cycle
         endif

      enddo DO_NKEYS
      F_istat = maxval(F_status)

      call msg(MSG_DEBUG, '(fstmpi) rdhint_3d_r4 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_rdhint_3d_r4


   !/@
   function fstmpi_rdhint_3d_r4_vect(F_data1, F_data2, F_status, &
        F_keys1, F_keys2, F_hintlist_S, &
        F_fileids, F_outgridid, F_coregridid, F_realloc_L) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real, pointer :: F_data1(:,:,:), F_data2(:,:,:)
      integer, intent(out) :: F_status(:)
      integer, intent(in) :: F_keys1(:), F_keys2(:)
      character(len=*), intent(in) :: F_hintlist_S(:)
      integer, intent(in) :: F_fileids(:)
      integer, intent(in) :: F_outgridid
      integer, intent(in), optional :: F_coregridid
      logical,intent(in), optional :: F_realloc_L
      !#TODO: alternate itf with hgridis_S so we have min, max as well
      !@author
      !@return
      integer :: F_istat
      !@/
      logical :: realloc_L
      integer :: istat, coregridid, nkeys, nhint, nfids, nij(2), ikey, fid, ingridid
      character(len=12) :: nomvar1_S, nomvar2_S
      real, pointer :: indata1(:,:,:), indata2(:,:,:)
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fstmpi) rdhint_3d_r4_vect [BEGIN]')
      F_istat = RMN_ERR
      F_status = RMN_ERR

      realloc_L = .false.
      coregridid = F_outgridid
      if (present(F_realloc_L)) realloc_L = F_realloc_L
      if (present(F_coregridid)) coregridid = F_coregridid

      nkeys   = size(F_keys1)
      nhint   = size(F_hintlist_S)
      nfids   = size(F_fileids)

      if (nkeys /= size(F_keys2)) then
         call msg(MSG_WARNING, '(fst) rdhint: incompatible list size')
         F_istat = RMN_ERR
         return
      endif

      if (size(F_status) < nkeys) then
         call msg(MSG_WARNING, '(fst) rdhint: status array too small')
         return
      endif

      F_istat = ezgrid_params(F_outgridid, nij)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING, '(fst) rdhint: Problem getting grid params')
         return
      endif

      F_istat = fst_checkalloc(F_data1, nij(1), nij(2), nkeys, realloc_L)
      if (.not.RMN_IS_OK(F_istat)) return
      F_istat = fst_checkalloc(F_data2, nij(1), nij(2), nkeys, realloc_L)
      if (.not.RMN_IS_OK(F_istat)) return

      nullify(indata1, indata2)
      DO_NKEYS: do ikey = 1, nkeys
         if (.not.RMN_IS_OK(F_keys1(ikey))) cycle

         ingridid = -1
         fid = F_fileids(min(ikey, nfids))
         istat = fstmpi_read_3d_r4(F_keys1(ikey), indata1, fid, ingridid, &
              F_realloc_L=realloc_L)
         nomvar1_S = '  ' !#TODO get nomvar from fst_read or getmeta
         if (.not.RMN_IS_OK(istat)) then
            F_status(ikey) = RMN_ERR
            call msg(MSG_WARNING, '(fstmpi) rdhint: Problem reading: '//trim(nomvar1_S))
            cycle
         endif
         istat = fstmpi_read_3d_r4(F_keys2(ikey), indata2, fid, ingridid, &
              F_realloc_L=realloc_L)
         nomvar2_S = '  ' !#TODO get nomvar from fst_read or getmeta
         if (.not.RMN_IS_OK(istat)) then
            F_status(ikey) = RMN_ERR
            call msg(MSG_WARNING, '(fstmpi) rdhint: Problem reading: '//trim(nomvar2_S))
            cycle
         endif
         !#TODO: allow for optional stats of read field before interp

         F_status(ikey) = hinterp4yy2d(F_data1, F_data2, indata1, indata2, ikey, &
              ingridid, F_outgridid, coregridid, &
              F_hintlist_S(min(ikey, nhint)), nomvar1_S, nomvar2_S)
         if (.not.RMN_IS_OK(F_status(ikey))) then
            call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in hinterp4yy2d for: '//trim(nomvar1_S)//' + '//trim(nomvar2_S))
            cycle
         endif

      enddo DO_NKEYS
      F_istat = maxval(F_status)

      call msg(MSG_DEBUG, '(fstmpi) rdhint_3d_r4_vect [END]')
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_rdhint_3d_r4_vect


   !/@
   function fstmpi_getmeta(F_key,F_nomvar_S,F_dateo,F_deet,F_npas,&
        F_ip1, F_ip2, F_ip3, F_etiket_S,F_typvar_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_key
      character(len=*),intent(out),optional :: F_nomvar_S, F_etiket_S,F_typvar_S
      integer,intent(out),optional :: F_dateo,F_deet,F_npas,F_ip1, F_ip2, F_ip3
      !@author
      !@return
      integer :: F_istat
      !@/
      integer,parameter :: MYSTRLEN = 32
      integer :: istat,ii,nn,mysize,buffer(7+3*MYSTRLEN/CHARPERBYTE)
      character(len=CHARPERBYTE) :: str4_S
      character(len=MYSTRLEN) :: typvar_S,nomvar_S, etiket_S, grtyp_S
      integer :: ni1,nj1,nk1, &
           dateo,deet,npas,nbits,datyp,ip1,ip2,ip3,&
           ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3
      ! ---------------------------------------------------------------------
      call ptopo_init_var()
      F_istat = RMN_OK
      buffer = 0
      typvar_S = ' '
      nomvar_S = ' '
      etiket_S = ' '
      if (ptopo_isblocmaster_L) then
         if (F_key < 0) then
            F_istat = RMN_ERR
            buffer(1) = F_istat
            ii = 7
            do nn=1,3*mysize
               ii = ii + 1
               buffer(ii) = transfer(' ',istat)
            enddo
         else
            F_istat = fstprm(F_key, dateo,deet,npas, ni1,nj1,nk1, &
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
            ii = 7
            mysize = MYSTRLEN/CHARPERBYTE
            do nn=1,mysize,CHARPERBYTE
               ii = ii + 1
               buffer(ii) = transfer(nomvar_S(nn:nn+CHARPERBYTE-1),istat)
            enddo
            do nn=1,mysize,CHARPERBYTE
               ii = ii + 1
               buffer(ii) = transfer(etiket_S(nn:nn+CHARPERBYTE-1),istat)
            enddo
            do nn=1,mysize,CHARPERBYTE
               ii = ii + 1
               buffer(ii) = transfer(typvar_S(nn:nn+CHARPERBYTE-1),istat)
            enddo
         endif
     endif
      mysize = size(buffer)
      call rpn_comm_bcast(buffer, mysize, RPN_COMM_INTEGER, RPN_COMM_MASTER, RPN_COMM_BLOC_COMM,istat)
      F_istat = buffer(1)
      if (.not.RMN_IS_OK(F_istat)) return

      if (present(F_dateo)) F_dateo = buffer(2)
      if (present(F_deet)) F_deet = buffer(3)
      if (present(F_npas)) F_npas = buffer(4)
      if (present(F_ip1)) F_ip1 = buffer(5)
      if (present(F_ip2)) F_ip2 = buffer(6)
      if (present(F_ip3)) F_ip3 = buffer(7)
      ii = 7
      mysize = MYSTRLEN/CHARPERBYTE
      do nn=1,mysize,CHARPERBYTE
         ii = ii + 1
         nomvar_S(nn:nn+CHARPERBYTE-1) = transfer(buffer(ii),str4_S)
      enddo
      do nn=1,mysize,CHARPERBYTE
         ii = ii + 1
         etiket_S(nn:nn+CHARPERBYTE-1) = transfer(buffer(ii),str4_S)
      enddo
      do nn=1,mysize,CHARPERBYTE
         ii = ii + 1
         typvar_S(nn:nn+CHARPERBYTE-1) = transfer(buffer(ii),str4_S)
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
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_getmeta


   !/@
   function fstmpi_get_hgridid(F_fileid,F_key) result(F_gridid)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_fileid,F_key
      !@author
      !@return
      integer :: F_gridid
      !@/
      ! ---------------------------------------------------------------------
      call ptopo_init_var()
      F_gridid = RMN_ERR
      if (ptopo_isblocmaster_L) then
         F_gridid = fst_get_gridid(F_fileid,F_key)
      endif
      F_gridid = ezgrid_bcast(F_gridid,RPN_COMM_BLOC_COMM)
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_get_hgridid


   !/@
   function fstmpi_get_vgrid(F_fileid,F_key,F_vgrid,F_ip1list,F_lvltyp_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_fileid,F_key
      type(vgrid_descriptor),intent(out) :: F_vgrid
      integer,pointer :: F_ip1list(:)
      character(len=*),intent(out) :: F_lvltyp_S
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: itype
      ! ---------------------------------------------------------------------
      F_istat = RMN_OK
      call vgrid_nullify(F_vgrid)
      itype = 0
      if (ptopo_isblocmaster_L) then
         F_istat = fst_get_vgrid(F_fileid,F_key,F_vgrid,F_ip1list,F_lvltyp_S)
      endif
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return
      F_istat = vgrid_wb_bcast(F_vgrid,F_ip1list,itype,F_lvltyp_S,RPN_COMM_BLOC_COMM)
      ! ---------------------------------------------------------------------
      return
   end function fstmpi_get_vgrid

end module fstmpi_read_mod
