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
module ptr_store
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland, 2012-02
   !@description
   ! Public functions
   public :: ptr_store_new,ptr_store_get,ptr_store_free
   ! Public constants
   integer,parameter,public :: PTR_STORE_NAME_LEN = 64
   !
!@/
#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>
#include <rmn/msg.h>

   interface ptr_store_new
      module procedure ptr_store_new_r4_3d
   end interface

   interface ptr_store_get
      module procedure ptr_store_get_r4_3d
      module procedure ptr_store_get_r4_3d_n
   end interface

   interface ptr_store_free
      module procedure ptr_store_free_r4_3d
   end interface

   type :: ptr_store_r4_3d
      real,pointer :: ptr(:,:,:)
      integer :: l_ijk(3),u_ijk(3)
      logical :: islock_L
      character(len=PTR_STORE_NAME_LEN) :: name_S
   end type ptr_store_r4_3d

   integer,parameter :: NMAX_PRT = 1024
   integer,save :: m_nbptr_r4_3d = 0
   type(ptr_store_r4_3d),save :: m_ptr_list_r4_3d(NMAX_PRT)

contains

   !/@*
   function ptr_store_new_r4_3d(F_name_S,F_data) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_name_S
      real,pointer :: F_data(:,:,:)
      integer :: F_istat
      !*@/
      integer :: ii,istat
      character(len=PTR_STORE_NAME_LEN) :: name_S
      real,pointer :: data(:,:,:)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (F_name_S == ' ' .or. .not.associated(F_data)) return
      name_S = F_name_S
      istat = clib_tolower(name_S)
      nullify(data)
      call ptr_store_get(F_name_S,data)
      if (associated(data)) then
         call msg(MSG_ERROR,'(ptr_store) New, Name already in use: '//trim(name_S))
         return
      endif
 
      do ii=1,m_nbptr_r4_3d
         if (associated(m_ptr_list_r4_3d(ii)%ptr,F_data)) then
            if (m_ptr_list_r4_3d(ii)%name_S /= ' ') then
               call msg(MSG_ERROR,'(ptr_store) New, Ptr already registered as '&
                    //trim(m_ptr_list_r4_3d(ii)%name_S)//', Cannot rename to: '&
                    //trim(name_S))
               return
            else
               !TODO-later: should we allow naming of locked ptr
               m_ptr_list_r4_3d(ii)%name_S = name_S
               m_ptr_list_r4_3d(ii)%islock_L = .true.
               call msg(MSG_INFO,'(ptr_store) New, Naming existing ptr: '//trim(name_S))
               F_istat = RMN_OK
               return
            endif
         endif
      enddo

      if (m_nbptr_r4_3d >= NMAX_PRT) then
         call msg(MSG_ERROR,'(ptr_store) New, Too many ptr, Cannot register: '//trim(name_S))
         return
      endif

      m_nbptr_r4_3d = m_nbptr_r4_3d + 1
      ii = m_nbptr_r4_3d
      m_ptr_list_r4_3d(ii)%ptr => F_data
      m_ptr_list_r4_3d(ii)%l_ijk = lbound(F_data)
      m_ptr_list_r4_3d(ii)%u_ijk = ubound(F_data)
      m_ptr_list_r4_3d(ii)%islock_L = .true.
      m_ptr_list_r4_3d(ii)%name_S = name_S
      call msg(MSG_INFO,'(ptr_store) New: '//trim(name_S))
      F_istat = RMN_OK
      ! ---------------------------------------------------------------------
      return
   end function ptr_store_new_r4_3d


   !/@*
   subroutine ptr_store_get_r4_3d_n(F_name_S,F_data)
      implicit none
      character(len=*),intent(in) :: F_name_S
      real,pointer :: F_data(:,:,:)
      !*@/
      integer :: ii,istat
      character(len=PTR_STORE_NAME_LEN) :: name_S
      ! ---------------------------------------------------------------------
      nullify(F_data)
      if (F_name_S == ' ') return
      name_S = F_name_S
      istat = clib_tolower(name_S)
      do ii=1,m_nbptr_r4_3d
         if (m_ptr_list_r4_3d(ii)%name_S == name_S) then
            F_data => m_ptr_list_r4_3d(ii)%ptr
            m_ptr_list_r4_3d(ii)%islock_L = .true.
            exit
         endif
      enddo
      ! ---------------------------------------------------------------------
      return
   end subroutine ptr_store_get_r4_3d_n


   !/@*
   subroutine ptr_store_get_r4_3d(F_data,l_ijk,u_ijk)
      implicit none
      real,pointer :: F_data(:,:,:)
      integer,intent(in) :: l_ijk(3),u_ijk(3)
      !*@/
      integer :: ii,istat,idx
      ! ---------------------------------------------------------------------
      nullify(F_data)
      idx = -1
      do ii=1,m_nbptr_r4_3d
         if (m_ptr_list_r4_3d(ii)%name_S /= ' ' .or. &
              m_ptr_list_r4_3d(ii)%islock_L) cycle

         if (all(m_ptr_list_r4_3d(ii)%l_ijk == l_ijk) .and. &
              all(m_ptr_list_r4_3d(ii)%u_ijk == u_ijk)) then
            if (associated(m_ptr_list_r4_3d(ii)%ptr)) then
               m_ptr_list_r4_3d(ii)%islock_L = .true.
               F_data => m_ptr_list_r4_3d(ii)%ptr
               return
            else
               idx = ii
               exit
            endif
         endif
      enddo

      if (idx < 1) then
         if (m_nbptr_r4_3d >= NMAX_PRT) return
         m_nbptr_r4_3d = m_nbptr_r4_3d + 1
         idx = m_nbptr_r4_3d
      endif

      allocate(m_ptr_list_r4_3d(idx)%ptr(l_ijk(1):u_ijk(1),l_ijk(2):u_ijk(2),l_ijk(3):u_ijk(3)),stat=istat)
      if (istat /= 0) then
         call msg(MSG_ERROR,'(ptr_store) Allocation problem')
         nullify(m_ptr_list_r4_3d(idx)%ptr)
         if (idx == m_nbptr_r4_3d) m_nbptr_r4_3d = m_nbptr_r4_3d - 1
         return
      endif
      F_data => m_ptr_list_r4_3d(idx)%ptr
      m_ptr_list_r4_3d(idx)%l_ijk = l_ijk
      m_ptr_list_r4_3d(idx)%u_ijk = u_ijk
      m_ptr_list_r4_3d(idx)%islock_L = .true.
      m_ptr_list_r4_3d(idx)%name_S = ' '
      ! ---------------------------------------------------------------------
      return
   end subroutine ptr_store_get_r4_3d


   !/@*
   subroutine ptr_store_free_r4_3d(F_data)
      implicit none
      real,pointer :: F_data(:,:,:)
      !*@/
      integer :: ii
      ! ---------------------------------------------------------------------
      if (.not.associated(F_data)) return
      do ii=1,m_nbptr_r4_3d
         if (m_ptr_list_r4_3d(ii)%islock_L) then
            if (associated(m_ptr_list_r4_3d(ii)%ptr,F_data)) then
               nullify(F_data)
               m_ptr_list_r4_3d(ii)%islock_L = .false.
               m_ptr_list_r4_3d(ii)%name_S = ' '
               return
            endif
         endif
      enddo
      ! ---------------------------------------------------------------------
      return
   end subroutine ptr_store_free_r4_3d


end module ptr_store
