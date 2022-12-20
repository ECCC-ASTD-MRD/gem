!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

#include <rmn/msg.h>

!/@*
module bus_fold_mod
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use clib_itf_mod, only: clib_tolower
   use rmn_gmm
   implicit none
   private
   public :: bus_fold_set, bus_fold, bus_unfold
   !@objective Functions to trasfer data from bus/folded space to
   !           grided (2d/3d) space and back
   !@author Michel Desgagne  -  sping 2010
   !revision
   ! 2012-03 S.Chamberland: split from gemdyn
   !@description
   !TODO-later: describe bus structure
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   integer,parameter :: MAXBUS = 32
   logical,parameter :: UPDATE_L = .true.

   type :: bus_info
      character(len=32) :: name_S  !- Bus name
      real,pointer :: ptr(:,:)     !- ptr to bus
      integer :: bni,bnj           !- bus x,y dims
      integer :: bnikv             !- ubound(ptr,1)
      integer :: i0,in,j0,jn,k0,kn !- subset of 2d/3d data to fold (bus idx
                                   !- start at 1 always, thus k0=2 introduce a shift)
   end type bus_info

   type(bus_info),save :: buslist(MAXBUS)
   integer,save :: m_nbus = 0

   interface bus_fold
      module procedure bus_fold2d_1
      module procedure bus_fold3d_1
      module procedure bus_fold2d_2
      module procedure bus_fold3d_2
      module procedure bus_fold_gmm
   end interface

   interface bus_unfold
      module procedure bus_unfold2d_1
      module procedure bus_unfold3d_1
      module procedure bus_unfold2d_2
      module procedure bus_unfold3d_2
      module procedure bus_unfold_gmm
   end interface

contains

   function priv_bus_id(F_bus_S) result(F_id)
      implicit none
      character(len=*),intent(inout) :: F_bus_S
      integer :: F_id,id,istat
      !---------------------------------------------------------------
      F_id = RMN_ERR
      istat = clib_tolower(F_bus_S)
      do id = 1,m_nbus
         if (F_bus_S == buslist(id)%name_S) then
            F_id = id
            exit
         endif
      enddo
      !---------------------------------------------------------------
      return
   end function priv_bus_id


   !/@*
   function bus_fold_set(F_bus_S,F_bus,F_i0,F_in,F_j0,F_jn,F_k0,F_kn,F_bus_ni) result(F_istat)
      implicit none
      !@Objective Register a bus
      !@Arguments
      character(len=*),intent(in) :: F_bus_S !- Bus name
      integer,intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0,F_kn !- subset of 2d/3d data to fold
      integer,intent(in) :: F_bus_ni !- Bus x lenght
      real,pointer :: F_bus(:,:)
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      character(len=128) :: msg_S
      character(len=32) :: bus_S
      integer :: id
      !---------------------------------------------------------------
      F_istat = RMN_ERR
      bus_S = F_bus_S
      id = priv_bus_id(bus_S)
      if (RMN_IS_OK(id)) then
         call msg(MSG_ERROR,'(bus_fold) Bus name already reserved')
         return
      endif
      if (.not.associated(F_bus)) then
         call msg(MSG_ERROR,'(bus_fold) Bus pointer not associated')
         return
      endif
      if (m_nbus >= MAXBUS) then
         call msg(MSG_ERROR,'(bus_fold) Too many registered buses')
         return
      endif
      m_nbus = m_nbus + 1
      buslist(m_nbus)%name_S = bus_S
      buslist(m_nbus)%ptr => F_bus
      buslist(m_nbus)%bni = F_bus_ni
      buslist(m_nbus)%bnj = ubound(F_bus,2)
      buslist(m_nbus)%bnikv = ubound(F_bus,1)
      buslist(m_nbus)%i0 = F_i0
      buslist(m_nbus)%in = F_in
      buslist(m_nbus)%j0 = F_j0
      buslist(m_nbus)%jn = F_jn
      buslist(m_nbus)%k0 = F_k0
      buslist(m_nbus)%kn = F_kn
      write(msg_S,'(a,3i8,6(a,i6),a)') trim(F_bus_S)//"(", &
           buslist(m_nbus)%bni,buslist(m_nbus)%bnikv,buslist(m_nbus)%bnj, &
           ") <=> (",F_i0,":",F_in,",",F_j0,":",F_jn,",",F_k0,":",F_kn,")"
      call msg(MSG_INFO,'(bus_fold) Set: '//trim(msg_S))
      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function bus_fold_set


   !/@*
   function bus_fold2d_1(F_bus,F_data2d,F_busidx,F_i0,F_in,F_j0,F_jn,F_bus_ni) result(F_istat)
      implicit none
      !@Objective  Transfer data3d to bus/folded space
      !@Arguments
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      integer,intent(in) :: F_i0,F_in,F_j0,F_jn !- data3d subset to fold
      integer,intent(in) :: F_bus_ni !- Bus x lenght
      real,pointer :: F_data2d(:,:)
      real,pointer :: F_bus(:,:)
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      integer :: j,minxyz(3),maxxyz(3),k0,kn,bussize
      character(len=64) :: msg_S
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_fold2d_1 BEGIN')
      F_istat = RMN_ERR
      if (.not.associated(F_bus)) then
         call msg(MSG_WARNING,'(bus_fold) Bus pointer not associated')
         return
      endif
      if (F_busidx<1 .or. F_busidx-1 + F_bus_ni > ubound(F_bus,1)) then
         write(msg_S,'(i9," < 1 .or. ",i9," > ",i9)') F_busidx, F_busidx-1 + F_bus_ni,ubound(F_bus,1)
         call msg(MSG_WARNING,'(bus_fold) busidx out of bound:'//msg_s)
         return
      endif
      if (.not.associated(F_data2d)) then
         call msg(MSG_WARNING,'(bus_fold) Data pointer not associated')
         return
      endif
      F_istat = RMN_OK
      minxyz(1:2) = lbound(F_data2d) ; minxyz(3) = 1
      maxxyz(1:2) = ubound(F_data2d) ; maxxyz(3) = 1
      k0 = 1 ; kn = 1
      bussize = size(F_bus(F_busidx:,1))
      do j=1,ubound(F_bus,2)
         call bus_fold_row(F_bus(F_busidx:,j),F_data2d,j,F_i0,F_in,F_j0,F_jn, &
              k0,kn,F_bus_ni,bussize,minxyz(1),maxxyz(1),minxyz(2),maxxyz(2), &
              minxyz(3),maxxyz(3))
      enddo
      call msg(MSG_DEBUG,'(bus_fold) bus_fold2d_1 END')
      !---------------------------------------------------------------
      return
   end function bus_fold2d_1


   !/@*
   function bus_unfold2d_1(F_bus,F_data2d,F_busidx,F_i0,F_in,F_j0,F_jn,F_bus_ni) result(F_istat)
      implicit none
      !@Objective Transfer bus/folded to 3d space
      !@Arguments
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      integer,intent(in) :: F_i0,F_in,F_j0,F_jn !- data3d subset to fold
      integer,intent(in) :: F_bus_ni !- Bus x lenght
      real,pointer :: F_data2d(:,:)
      real,pointer :: F_bus(:,:)
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      integer :: j,minxyz(3),maxxyz(3),k0,kn,bussize
      character(len=64) :: msg_S
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold2d_1 BEGIN')
      F_istat = RMN_ERR
      if (.not.associated(F_bus)) then
         call msg(MSG_WARNING,'(bus_unfold) Bus pointer not associated')
         return
      endif
      if (F_busidx<1 .or. F_busidx-1 + F_bus_ni > ubound(F_bus,1)) then
         write(msg_S,'(i9," < 1 .or. ",i9," > ",i9)') F_busidx, F_busidx-1 + F_bus_ni,ubound(F_bus,1)
         call msg(MSG_WARNING,'(bus_unfold) busidx out of bound:'//msg_s)
         return
      endif
      if (.not.associated(F_data2d)) then
         call msg(MSG_WARNING,'(bus_unfold) Data pointer not associated')
         return
      endif
      F_istat = RMN_OK
      minxyz(1:2) = lbound(F_data2d) ; minxyz(3) = 1
      maxxyz(1:2) = ubound(F_data2d) ; maxxyz(3) = 1
      k0 = 1 ; kn = 1
      bussize = size(F_bus(F_busidx:,1))
      do j=1,ubound(F_bus,2)
         call bus_unfold_row(F_bus(F_busidx:,j),F_data2d,j,F_i0,F_in, &
              F_j0,F_jn,k0,kn,F_bus_ni,bussize,minxyz(1),maxxyz(1), &
              minxyz(2),maxxyz(2),minxyz(3),maxxyz(3))
      enddo
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold2d_1 END')
      !---------------------------------------------------------------
      return
   end function bus_unfold2d_1


   !/@*
   function bus_fold3d_1(F_bus,F_data3d,F_busidx,F_i0,F_in,F_j0,F_jn,F_k0,F_kn,F_bus_ni) result(F_istat)
      implicit none
      !@Objective  Transfer data3d to bus/folded space
      !@Arguments
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      integer,intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0,F_kn !- data3d subset to fold
      integer,intent(in) :: F_bus_ni !- Bus x lenght
      real,pointer :: F_data3d(:,:,:)
      real,pointer :: F_bus(:,:)
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      integer :: j,minxyz(3),maxxyz(3),bussize
      character(len=64) :: msg_S
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_fold3d_1 BEGIN')
      F_istat = RMN_ERR
      if (.not.associated(F_bus)) then
         call msg(MSG_WARNING,'(bus_fold) Bus pointer not associated')
         return
      endif
      if (F_busidx<1 .or. F_busidx-1 + F_bus_ni > ubound(F_bus,1)) then
         write(msg_S,'(i9," < 1 .or. ",i9," > ",i9)') F_busidx, F_busidx-1 + F_bus_ni,ubound(F_bus,1)
         call msg(MSG_WARNING,'(bus_fold) busidx out of bound:'//msg_s)
         return
      endif
      if (.not.associated(F_data3d)) then
         call msg(MSG_WARNING,'(bus_fold) Data pointer not associated')
         return
      endif
      F_istat = RMN_OK
      minxyz = lbound(F_data3d)
      maxxyz = ubound(F_data3d)
      bussize = size(F_bus(F_busidx:,1))
      do j=1,ubound(F_bus,2)
         call bus_fold_row(F_bus(F_busidx:,j),F_data3d,j,F_i0,F_in,F_j0,F_jn, &
              F_k0,F_kn,F_bus_ni,bussize,minxyz(1),maxxyz(1), &
              minxyz(2),maxxyz(2),minxyz(3),maxxyz(3))
      enddo
      call msg(MSG_DEBUG,'(bus_fold) bus_fold3d_1 END')
      !---------------------------------------------------------------
      return
   end function bus_fold3d_1


   !/@*
   function bus_unfold3d_1(F_bus,F_data3d,F_busidx,F_i0,F_in,F_j0,F_jn,F_k0,F_kn,F_bus_ni) result(F_istat)
      implicit none
      !@Objective Transfer bus/folded to 3d space
      !@Arguments
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      integer,intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0,F_kn !- data3d subset to fold
      integer,intent(in) :: F_bus_ni !- Bus x lenght
      real,pointer :: F_data3d(:,:,:)
      real,pointer :: F_bus(:,:)
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      integer :: j,minxyz(3),maxxyz(3),bussize
      character(len=64) :: msg_S
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold3d_1 BEGIN')
      F_istat = RMN_ERR
      if (.not.associated(F_bus)) then
         call msg(MSG_WARNING,'(bus_unfold) Bus pointer not associated')
         return
      endif
      if (F_busidx<1 .or. F_busidx-1 + F_bus_ni > ubound(F_bus,1)) then
         write(msg_S,'(i9," < 1 .or. ",i9," > ",i9)') F_busidx, F_busidx-1 + F_bus_ni,ubound(F_bus,1)
         call msg(MSG_WARNING,'(bus_unfold) busidx out of bound:'//msg_s)
         return
      endif
      if (.not.associated(F_data3d)) then
         call msg(MSG_WARNING,'(bus_unfold) Data pointer not associated')
         return
      endif
      F_istat = RMN_OK
      minxyz = lbound(F_data3d)
      maxxyz = ubound(F_data3d)
      bussize = size(F_bus(F_busidx:,1))
      do j=1,ubound(F_bus,2)
         call bus_unfold_row(F_bus(F_busidx:,j),F_data3d,j,F_i0,F_in,F_j0, &
              F_jn,F_k0,F_kn,F_bus_ni,bussize,minxyz(1),maxxyz(1),minxyz(2), &
              maxxyz(2),minxyz(3),maxxyz(3))
      enddo
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold3d_1 END')
      !---------------------------------------------------------------
      return
   end function bus_unfold3d_1


   !/@*
   function bus_fold2d_2(F_data2d,F_bus_S,F_busidx) result(F_istat)
      implicit none
      !@Objective  Transfer data3d to bus/folded space
      !@Arguments
      character(len=*),intent(in) :: F_bus_S
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      real,pointer :: F_data2d(:,:)
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      character(len=64) :: bus_S,msg_S
      integer :: j,minxyz(3),maxxyz(3),k0,kn,id,bussize
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_fold2d_2 BEGIN')
      F_istat = RMN_ERR
      bus_S = F_bus_S
      id = priv_bus_id(bus_S)
      if (.not.RMN_IS_OK(id)) then
         call msg(MSG_WARNING,'(bus_fold) Not such bus: '//trim(bus_S))
         return
      endif
      if (.not.associated(F_data2d)) then
         call msg(MSG_WARNING,'(bus_fold) Data pointer not associated')
         return
      endif
      if (F_busidx<1 .or. F_busidx-1 + buslist(id)%bni > buslist(id)%bnikv) then
         write(msg_S,'(i9," < 1 .or. ",i9," > ",i9)') F_busidx,F_busidx-1 + buslist(id)%bni,buslist(id)%bnikv
         call msg(MSG_WARNING,'(bus_fold) busidx out of bound: '//trim(bus_S)//' '//msg_s)
         return
      endif
      F_istat = RMN_OK
      minxyz(1:2) = lbound(F_data2d) ; minxyz(3) = 1
      maxxyz(1:2) = ubound(F_data2d) ; maxxyz(3) = 1
      k0 = 1 ; kn = 1

      bussize = size(buslist(id)%ptr(F_busidx:,1))
      write(msg_S,'(a,3i6,a,3i6)') trim(bus_S)//" ",minxyz," ; ",maxxyz
      call msg(MSG_DEBUG,'(bus_fold) '//trim(msg_S))
      do j=1,buslist(id)%bnj
         call bus_fold_row(buslist(id)%ptr(F_busidx:,j),F_data2d,j, &
              buslist(id)%i0,buslist(id)%in,buslist(id)%j0,buslist(id)%jn, &
              k0,kn,buslist(id)%bni,bussize, &
              minxyz(1),maxxyz(1),minxyz(2),maxxyz(2),minxyz(3),maxxyz(3))
      enddo
      call msg(MSG_DEBUG,'(bus_fold) bus_fold2d_2 END')
      !---------------------------------------------------------------
      return
   end function bus_fold2d_2


   !/@*
   function bus_unfold2d_2(F_data2d,F_bus_S,F_busidx) result(F_istat)
      implicit none
      !@Objective  Transfer bus/folded to 3d space
      !@Arguments
      character(len=*),intent(in) :: F_bus_S
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      real,pointer :: F_data2d(:,:)
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      character(len=64) :: bus_S,msg_S
      integer :: j,minxyz(3),maxxyz(3),k0,kn,id,bussize
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold2d_2 BEGIN')
      F_istat = RMN_ERR
      bus_S = F_bus_S
      id = priv_bus_id(bus_S)
      if (.not.RMN_IS_OK(id)) then
         call msg(MSG_WARNING,'(bus_unfold) Not such bus: '//trim(bus_S))
         return
      endif
      if (F_busidx<1 .or. F_busidx-1 + buslist(id)%bni > buslist(id)%bnikv) then
         write(msg_S,'(i9," < 1 .or. ",i9," > ",i9)') F_busidx,F_busidx-1 + buslist(id)%bni,buslist(id)%bnikv
         call msg(MSG_WARNING,'(bus_unfold) busidx out of bound: '//trim(bus_S)//' '//msg_s)
         return
      endif
      if (.not.associated(F_data2d)) then
         call msg(MSG_WARNING,'(bus_unfold) Data pointer not associated')
         return
      endif
      F_istat = RMN_OK
      minxyz(1:2) = lbound(F_data2d) ; minxyz(3) = 1
      maxxyz(1:2) = ubound(F_data2d) ; maxxyz(3) = 1
      k0 = 1 ; kn = 1
      bussize = size(buslist(id)%ptr(F_busidx:,1))
      write(msg_S,'(a,3i6,a,3i6)') trim(bus_S)//" ",minxyz," ; ",maxxyz
      call msg(MSG_DEBUG,'(bus_unfold) '//trim(msg_S))
      do j=1,buslist(id)%bnj
         call bus_unfold_row(buslist(id)%ptr(F_busidx:,j),F_data2d,j, &
              buslist(id)%i0,buslist(id)%in,buslist(id)%j0,buslist(id)%jn, &
              k0,kn,buslist(id)%bni,bussize, &
              minxyz(1),maxxyz(1),minxyz(2),maxxyz(2),minxyz(3),maxxyz(3))
      enddo
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold2d_2 END')
      !---------------------------------------------------------------
      return
   end function bus_unfold2d_2


   !/@*
   function bus_fold3d_2(F_data3d,F_bus_S,F_busidx,F_k0,F_kn) result(F_istat)
      implicit none
      !@Objective  Transfer data3d to bus/folded space
      !@Arguments
      character(len=*),intent(in) :: F_bus_S
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      real,pointer :: F_data3d(:,:,:)
      integer,intent(in),optional :: F_k0,F_kn !- subset of levels to fold
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      character(len=64) :: bus_S,msg_S
      integer :: minxyz(3),maxxyz(3),id,k0,kn
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_fold3d_2 BEGIN')
      F_istat = RMN_ERR
      bus_S = F_bus_S
      id = priv_bus_id(bus_S)
      if (.not.RMN_IS_OK(id)) then
         call msg(MSG_WARNING,'(bus_fold) Not such bus: '//trim(bus_S))
         return
      endif
      if (F_busidx<1 .or. F_busidx-1 + buslist(id)%bni > buslist(id)%bnikv) then
         write(msg_S,'(i9," < 1 .or. ",i9," > ",i9)') F_busidx,F_busidx-1 + buslist(id)%bni,buslist(id)%bnikv
         call msg(MSG_WARNING,'(bus_fold) busidx out of bound: '//trim(bus_S)//' '//msg_s)
         return
      endif
      if (.not.associated(F_data3d)) then
         call msg(MSG_WARNING,'(bus_fold) Data pointer not associated')
         return
      endif
      F_istat = RMN_OK
      minxyz = lbound(F_data3d)
      maxxyz = ubound(F_data3d)
      k0 = minxyz(3) ; if (present(F_k0)) k0 = F_k0
      kn = maxxyz(3) ; if (present(F_kn)) kn = F_kn
      k0 = min(maxval((/minxyz(3),buslist(id)%k0,k0/)),maxxyz(3))
      kn = max(k0,minval((/kn,buslist(id)%kn,maxxyz(3)/)))
      write(msg_S,'(a,3i6,a,3i6,a,2i4)') trim(bus_S)//" ",minxyz," ; ",maxxyz," ; ",k0,kn
      call msg(MSG_DEBUG,'(bus_fold) '//trim(msg_S))
!!$      bussize = size(buslist(id)%ptr(F_busidx:,1))
!!$      do j=1,ubound(buslist(id)%ptr,2)
!!$         call bus_fold_row(buslist(id)%ptr(F_busidx:,j),F_data3d,j, &
!!$              buslist(id)%i0,buslist(id)%in,buslist(id)%j0,buslist(id)%jn, &
!!$              k0,kn,buslist(id)%bni,bussize, &
!!$              minxyz(1),maxxyz(1),minxyz(2),maxxyz(2),minxyz(3),maxxyz(3))
!!$      enddo
      call bus_fold_allrows(buslist(id)%ptr,F_data3d,F_busidx, &
           buslist(id)%i0,buslist(id)%in,buslist(id)%j0,buslist(id)%jn, &
           k0,kn,buslist(id)%bni,buslist(id)%bnikv,buslist(id)%bnj, &
           minxyz(1),maxxyz(1),minxyz(2),maxxyz(2),minxyz(3),maxxyz(3))
      call msg(MSG_DEBUG,'(bus_fold) bus_fold3d_2 END')
      !---------------------------------------------------------------
      return
   end function bus_fold3d_2


   !/@*
   function bus_unfold3d_2(F_data3d,F_bus_S,F_busidx,F_k0,F_kn) result(F_istat)
      implicit none
      !@Objective  Transfer bus/folded space to 3d
      !@Arguments
      character(len=*),intent(in) :: F_bus_S
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      real,pointer :: F_data3d(:,:,:)
      integer,intent(in),optional :: F_k0,F_kn !- subset of levels to fold
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      character(len=64) :: bus_S,msg_S
      integer :: j,minxyz(3),maxxyz(3),id,bussize,k0,kn
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold3d_2 BEGIN')
      F_istat = RMN_ERR
      bus_S = F_bus_S
      id = priv_bus_id(bus_S)
      if (.not.RMN_IS_OK(id)) then
         call msg(MSG_WARNING,'(bus_unfold) Not such bus: '//trim(bus_S))
         return
      endif
      if (F_busidx<1 .or. F_busidx-1 + buslist(id)%bni > buslist(id)%bnikv) then
         write(msg_S,'(i9," < 1 .or. ",i9," > ",i9)') F_busidx,F_busidx-1 + buslist(id)%bni,buslist(id)%bnikv
         call msg(MSG_WARNING,'(bus_unfold) busidx out of bound: '//trim(bus_S)//' '//msg_s)
         return
      endif
      if (.not.associated(F_data3d)) then
         call msg(MSG_WARNING,'(bus_unfold) Data pointer not associated')
         return
      endif
      F_istat = RMN_OK
      minxyz = lbound(F_data3d)
      maxxyz = ubound(F_data3d)
      k0 = minxyz(3) ; if (present(F_k0)) k0 = F_k0
      kn = maxxyz(3) ; if (present(F_kn)) kn = F_kn
      k0 = min(maxval((/minxyz(3),buslist(id)%k0,k0/)),maxxyz(3))
      kn = max(k0,minval((/kn,buslist(id)%kn,maxxyz(3)/)))
      write(msg_S,'(a,3i6,a,3i6,a,2i4)') trim(bus_S)//" ",minxyz," ; ",maxxyz," ; ",k0,kn
      call msg(MSG_DEBUG,'(bus_unfold) '//trim(msg_S))
      bussize = size(buslist(id)%ptr(F_busidx:,1))
      do j=1,ubound(buslist(id)%ptr,2)
         call bus_unfold_row(buslist(id)%ptr(F_busidx:,j),F_data3d,j, &
              buslist(id)%i0,buslist(id)%in,buslist(id)%j0,buslist(id)%jn, &
              k0,kn,buslist(id)%bni,bussize, &
              minxyz(1),maxxyz(1),minxyz(2),maxxyz(2),minxyz(3),maxxyz(3))
      enddo
!!$      call bus_unfold_allrows(buslist(id)%ptr,F_data3d,F_busidx, &
!!$           buslist(id)%i0,buslist(id)%in,buslist(id)%j0,buslist(id)%jn, &
!!$           k0,kn,buslist(id)%bni,buslist(id)%bnikv,buslist(id)%bnj, &
!!$           minxyz(1),maxxyz(1),minxyz(2),maxxyz(2),minxyz(3),maxxyz(3))
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold3d_2 END')
      !---------------------------------------------------------------
      return
   end function bus_unfold3d_2


   !/@*
   function bus_fold_gmm(F_gmmname_S,F_bus_S,F_busidx,F_k0,F_kn) result(F_istat)
      implicit none
      !@Objective  Transfer 3d space to from bus/folded
      !@Arguments
      character(len=*),intent(in) :: F_gmmname_S,F_bus_S
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      integer,intent(in),optional :: F_k0,F_kn !- subset of levels to fold
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      integer :: istat,k0,kn
      real,pointer :: data2d(:,:),data3d(:,:,:)
      character(len=GMM_MAXNAMELENGTH) :: gmmname_S
      type(gmm_metadata) :: mymeta
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_fold_gmm BEGIN')
      F_istat = RMN_ERR
      call gmmx_name_ci(gmmname_S,F_gmmname_S,.not.UPDATE_L) !# Note: bug in GMM, NOT case independent
      istat = gmm_getmeta(gmmname_S,mymeta)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(bus_fold) Not such gmm field: '//trim(F_gmmname_S))
         return
      endif
      nullify(data2d,data3d)
      if (mymeta%l(3)%n == 0) then
         istat = gmm_get(gmmname_S,data2d)
         if (.not.(RMN_IS_OK(istat).and.associated(data2d))) then
            call msg(MSG_WARNING,'(bus_fold) Problem getting gmm field: '//trim(F_gmmname_S))
            return
         endif
      else
         istat = gmm_get(gmmname_S,data3d)
         if (.not.(RMN_IS_OK(istat).and.associated(data3d))) then
            call msg(MSG_WARNING,'(bus_fold) Problem getting gmm field: '//trim(F_gmmname_S))
            return
         endif
      endif
      if (associated(data2d)) then
         F_istat = bus_fold(data2d,F_bus_S,F_busidx)
      else
         k0 = lbound(data3d,3) !# max(1,lbound(data3d,3))
         kn = ubound(data3d,3)
         if (present(F_k0)) k0 = min(max(k0,F_k0),kn)
         if (present(F_kn)) kn = min(max(k0,F_kn),kn)
         F_istat = bus_fold(data3d,F_bus_S,F_busidx,k0,kn)
      endif
      call msg(MSG_DEBUG,'(bus_fold) bus_fold_gmm END')
      !---------------------------------------------------------------
      return
   end function bus_fold_gmm


   !/@*
   function bus_unfold_gmm(F_gmmname_S,F_bus_S,F_busidx,F_k0,F_kn) result(F_istat)
      implicit none
      !@Objective  Transfer from bus/folded to 3d space
      !@Arguments
      character(len=*),intent(in) :: F_gmmname_S,F_bus_S
      integer,intent(in) :: F_busidx !- Start index of bus var to fill
      integer,intent(in),optional :: F_k0,F_kn !- subset of levels to fold
      !@return
      integer :: F_istat
      !@author Michel Desgagne  -  sping 2010
      !revision
      ! 2012-03 S.Chamberland: split from gemdyn
      !*@/
      integer :: istat,k0,kn
      real,pointer :: data2d(:,:),data3d(:,:,:)
      character(len=GMM_MAXNAMELENGTH) :: gmmname_S
      type(gmm_metadata) :: mymeta
      !---------------------------------------------------------------
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold_gmm BEGIN')
      F_istat = RMN_ERR
      call gmmx_name_ci(gmmname_S,F_gmmname_S,.not.UPDATE_L) !# Note: bug in GMM, NOT case independent
      istat = gmm_getmeta(gmmname_S,mymeta)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING,'(bus_unfold) Not such gmm field: '//trim(F_gmmname_S))
         return
      endif
      nullify(data2d,data3d)
      if (mymeta%l(3)%n == 0) then
         istat = gmm_get(gmmname_S,data2d)
         if (.not.(RMN_IS_OK(istat).and.associated(data2d))) then
            call msg(MSG_WARNING,'(bus_unfold) Problem getting gmm field: '//trim(F_gmmname_S))
            return
         endif
      else
         istat = gmm_get(gmmname_S,data3d)
         if (.not.(RMN_IS_OK(istat).and.associated(data3d))) then
            call msg(MSG_WARNING,'(bus_unfold) Problem getting gmm field: '//trim(F_gmmname_S))
            return
         endif
      endif
      if (associated(data2d)) then
         F_istat = bus_unfold(data2d,F_bus_S,F_busidx)
      else
         k0 = lbound(data3d,3) !# max(1,lbound(data3d,3))
         kn = ubound(data3d,3)
         if (present(F_k0)) k0 = min(max(k0,F_k0),kn)
         if (present(F_kn)) kn = min(max(k0,F_kn),kn)
         F_istat = bus_unfold(data3d,F_bus_S,F_busidx,k0,kn)
      endif
      call msg(MSG_DEBUG,'(bus_fold) bus_unfold_gmm END')
      !---------------------------------------------------------------
      return
   end function bus_unfold_gmm

end module bus_fold_mod


!/@*
subroutine bus_fold_row(F_subbus,F_data3d,F_bus_j,F_i0,F_in,F_j0,F_jn,F_k0, &
     F_kn,F_bus_ni,F_bus_nik,F_minx,F_maxx,F_miny,F_maxy,F_minz,F_maxz)
   implicit none
   !@Objective  Transfer data3d to bus/folded space
   !@Arguments
   integer,intent(in) :: F_bus_j !- Bus row to fold
   integer,intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0,F_kn !- data3d subset to fold
   integer,intent(in) :: F_bus_ni !- Bus x lenght
   integer,intent(in) :: F_bus_nik,F_minx,F_maxx,F_miny,F_maxy,F_minz,F_maxz
   real,intent(in) :: F_data3d(F_minx:F_maxx,F_miny:F_maxy,F_minz:F_maxz)
   real,intent(inout) :: F_subbus(F_bus_nik) !- sub bus array starting at bus var idx
   !@author Michel Desgagne  -  sping 2010
   !revision
   ! 2012-03 S.Chamberland: split from gemdyn
   !*@/
!!!#include <arch_specific.hf>
   integer :: buszise,k0,kn,ij0,nip,njp,offi,offj,ijmax,k,idx0,i,idx,ijp,ip,jp
!!$   character(len=8) :: row_S
   !---------------------------------------------------------------
!!$   write(row_S,'(i8)') F_bus_j
!!$   call msg(MSG_DEBUG,'(bus_fold) bus_fold_row BEGIN '//trim(row_S))
   buszise = size(F_subbus)
   k0 = min(max(F_minz,F_k0),F_maxz)
   kn = min(max(k0,F_kn),F_maxz)
   ij0 = (F_bus_j-1) * F_bus_ni
   !TODO-later: check dim consistency
   nip = F_in - F_i0 + 1
   njp = F_jn - F_j0 + 1
   offi = F_i0 - 1
   offj = F_j0 - 1
   ijmax = nip*njp

!$omp parallel private(idx0,i,idx,ijp,ip,jp)shared(F_subbus,F_data3d)
!$omp do
   do k=k0,kn
      idx0 = (k-k0)*F_bus_ni
      do i=1,F_bus_ni
         idx = min(idx0+i,buszise) !TODO-later: error if maxed out
         ijp = min(ij0+i,ijmax)
         jp= ijp/nip + min(1,mod(ijp,nip))
         ip= ijp - (jp-1)*nip + offi
         jp= jp + offj
         ip = min(max(F_minx,ip),F_maxx) !TODO-later: error if maxed out
         jp = min(max(F_miny,jp),F_maxy) !TODO-later: error if maxed out
         F_subbus(idx) = F_data3d(ip,jp,k)
     end do
   end do
!$omp end do
!$omp end parallel
!!$   call msg(MSG_DEBUG,'(bus_fold) bus_fold_row END '//trim(row_S))
   !---------------------------------------------------------------
   return
end subroutine bus_fold_row


!/@*
subroutine bus_fold_allrows(F_bus,F_data3d,F_bus_idx,F_i0,F_in,F_j0,F_jn,F_k0, &
   F_kn,F_bus_ni,F_bus_nikv,F_bus_nj,F_minx,F_maxx,F_miny,F_maxy,F_minz,F_maxz)
   implicit none
   !@Objective  Transfer data3d to bus/folded space
   !@Arguments
   integer,intent(in) :: F_bus_idx !- Bus var idx to fold
   integer,intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0,F_kn !- data3d subset to fold
   integer,intent(in) :: F_bus_ni !- Bus x lenght
   integer,intent(in) :: F_bus_nikv,F_bus_nj,F_minx,F_maxx,F_miny,F_maxy,F_minz,F_maxz
   real,intent(in) :: F_data3d(F_minx:F_maxx,F_miny:F_maxy,F_minz:F_maxz)
   real,intent(inout) :: F_bus(F_bus_nikv,F_bus_nj) !- sub bus array starting at bus var idx
   !@author Michel Desgagne  -  sping 2010
   !revision
   ! 2012-03 S.Chamberland: split from gemdyn
   !*@/
!!!#include <arch_specific.hf>
   integer :: buszise,k0,kn,ij0,nip,njp,offi,offj,ijmax,k,idx0,i,idx,ijp,ip,jp,bus_off_idx,bus_j
   character(len=9) :: idx_S
   !---------------------------------------------------------------
   write(idx_S,'(i9)') F_bus_idx
   call msg(MSG_DEBUG,'(bus_fold) allrows [BEGIN] idx='//trim(idx_S))
   buszise = F_bus_nikv-F_bus_idx+1
   bus_off_idx = F_bus_idx - 1
   k0 = min(max(F_minz,F_k0),F_maxz)
   kn = min(max(k0,F_kn),F_maxz)

   !TODO-later: check dim consistency
   nip = F_in - F_i0 + 1
   njp = F_jn - F_j0 + 1
   offi = F_i0 - 1
   offj = F_j0 - 1
   ijmax = nip*njp

   DO_BUS_J: do bus_j=1,F_bus_nj
      ij0 = (bus_j-1) * F_bus_ni

!$omp parallel private(idx0,i,idx,ijp,ip,jp)shared(F_bus,F_data3d)
!$omp do
      do k=k0,kn
         idx0 = (k-k0)*F_bus_ni
         do i=1,F_bus_ni
            idx = min(idx0+i,buszise) !TODO-later: error if maxed out
            ijp = min(ij0+i,ijmax)
            jp= ijp/nip + min(1,mod(ijp,nip))
            ip= ijp - (jp-1)*nip + offi
            jp= jp + offj
            ip = min(max(F_minx,ip),F_maxx) !TODO-later: error if maxed out
            jp = min(max(F_miny,jp),F_maxy) !TODO-later: error if maxed out
            F_bus(idx+bus_off_idx,bus_j) = F_data3d(ip,jp,k)
         end do
      end do
!$omp end do
!$omp end parallel

   enddo DO_BUS_J
   call msg(MSG_DEBUG,'(bus_fold) allrows [END] idx='//trim(idx_S))
   !---------------------------------------------------------------
   return
end subroutine bus_fold_allrows


!/@*
subroutine bus_unfold_row(F_subbus,F_data3d,F_bus_j,F_i0,F_in,F_j0,F_jn,F_k0, &
   F_kn,F_bus_ni,F_bus_nik,F_minx,F_maxx,F_miny,F_maxy,F_minz,F_maxz)
   implicit none
   !@Objective  Transfer bus/folded to 3d space
   !@Arguments
   integer,intent(in) :: F_bus_j !- Bus row to fold
   integer,intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0,F_kn !- data3d subset to fold
   integer,intent(in) :: F_bus_ni !- Bus x lenght
   integer,intent(in) :: F_bus_nik,F_minx,F_maxx,F_miny,F_maxy,F_minz,F_maxz
   real,intent(inout) :: F_data3d(F_minx:F_maxx,F_miny:F_maxy,F_minz:F_maxz)
   real,intent(in) :: F_subbus(F_bus_nik)  !- sub bus array starting at bus var idx
   !@author Michel Desgagne  -  sping 2010
   !revision
   ! 2012-03 S.Chamberland: split from gemdyn
   !*@/
!!!#include <arch_specific.hf>
   integer :: buszise,k0,kn,ij0,nip,njp,offi,offj,ijmax,k,idx0,i,idx,ijp,ip,jp
!!$   character(len=8) :: row_S
   !---------------------------------------------------------------
!!$   write(row_S,'(i8)') F_bus_j
!!$   call msg(MSG_DEBUG,'(bus_fold) bus_unfold_row BEGIN '//trim(row_S))
   buszise = size(F_subbus)
   k0 = min(max(F_minz,F_k0),F_maxz)
   kn = min(max(k0,F_kn),F_maxz)
   ij0 = (F_bus_j-1) * F_bus_ni
   !TODO-later: check dim consistency
   nip = F_in - F_i0 + 1
   njp = F_jn - F_j0 + 1
   offi = F_i0 - 1
   offj = F_j0 - 1
   ijmax = nip*njp
!$omp parallel private(idx0,i,idx,ijp,ip,jp)shared(F_subbus,F_data3d)
!$omp do
   do k=k0,kn
      idx0 = (k-k0)*F_bus_ni
      do i=1,F_bus_ni
         idx = min(idx0+i,buszise) !TODO-later: error if maxed out
         ijp = min(ij0+i,ijmax)
         jp= ijp/nip + min(1,mod(ijp,nip))
         ip= ijp - (jp-1)*nip + offi
         jp= jp + offj
         ip = min(max(F_minx,ip),F_maxx) !TODO-later: error if maxed out
         jp = min(max(F_miny,jp),F_maxy) !TODO-later: error if maxed out
         F_data3d(ip,jp,k) = F_subbus(idx)
      end do
   end do
!$omp end do
!$omp end parallel
!!$   call msg(MSG_DEBUG,'(bus_fold) bus_unfold_row END '//trim(row_S))
     !---------------------------------------------------------------
   return
end subroutine bus_unfold_row

