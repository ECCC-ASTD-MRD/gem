! * libdescrip - Vertical grid descriptor library for FORTRAN programming
! * Copyright (C) 2016  Direction du developpement des previsions nationales
! *                     Centre meteorologique canadien
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.

module vgrid_utils
   
   implicit none
   private

   ! Public utilities
   public :: get_allocate                        !allocation of array values
   public :: get_error,put_error                 !get/put error messaging 
   public :: same_vec                            !check for equivalence of arrays
   public :: up                                  !convert string to upper-case

   ! Public vgrid_descriptor constants
#include "vgrid_descriptors.hf"
   
   ! Private constants
   integer, parameter :: LONG_STRING=1024        !number of characters in a long string
   
   ! Private variables
   
   interface get_allocate
      module procedure get_allocate_i1d
      module procedure get_allocate_r1d
      module procedure get_allocate_r3d
      module procedure get_allocate_r81d
      module procedure get_allocate_r83d
   end interface get_allocate
   
   interface same_vec
      module procedure same_vec_i
      module procedure same_vec_r
      module procedure same_vec_r8
      module procedure same_vec_r83d
   end interface same_vec

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate space for pointer returns

   integer function get_allocate_i1d(key_S,value,len,allow_reshape_L,msg_S) result(istat)
      ! Allocate space for the result value and report error
      implicit none
      character(len=*), intent(in) :: key_S
      integer, dimension(:), pointer :: value
      integer, intent(in) :: len
      logical :: allow_reshape_L
      character(len=*) :: msg_S
      !Local variables
      logical :: alloc_lev_L   
      external msg
      istat=-1
      alloc_lev_L=.false.
      if(.not.associated(value))then
         alloc_lev_L=.true.
      else
         if(size(value)/=len)then
            if(allow_reshape_L)then
               write(for_msg,*) 'reshaping 1D integer vector '//trim(msg_S)
               call msg(MSG_INFO,VGD_PRFX//for_msg)
               deallocate(value)
               alloc_lev_L=.true.
            else
               write(for_msg,*) '1D pointer already allocated with a different length, will not reallocate '//trim(msg_S)
               call msg(MSG_ERROR,VGD_PRFX//for_msg)
               return
            endif
         endif
      endif
      if(alloc_lev_L)then
         allocate(value(len),stat=istat)
         if (istat /= 0) then
            write(for_msg,*) 'unable to allocate space for '//trim(key_S)//' request '//trim(msg_S)
            call msg(MSG_CRITICAL,for_msg)
         endif
      else
         istat=0
      endif
   end function get_allocate_i1d
   
   integer function get_allocate_r1d(key_S,value,len,allow_reshape_L,msg_S) result(istat)
      ! Allocate space for the result value and report error
      implicit none
      character(len=*), intent(in) :: key_S
      real, dimension(:), pointer :: value
      integer, intent(in) :: len
      logical :: allow_reshape_L
      character(len=*) :: msg_S
      !Local variables
      logical :: alloc_lev_L
      external msg
      istat=-1
      alloc_lev_L=.false.
      if(.not.associated(value))then
         alloc_lev_L=.true.
      else
         if(size(value)/=len)then
            if(allow_reshape_L)then
               write(for_msg,*) 'reshaping 1D real vector '//trim(msg_S)
               call msg(MSG_INFO,VGD_PRFX//for_msg)
               deallocate(value)
               alloc_lev_L=.true.
            else
               write(for_msg,*) '1D pointer already allocated with a different length, will not reallocate '//trim(msg_S)
               call msg(MSG_ERROR,VGD_PRFX//for_msg)
               return
            endif
         endif
      endif
      if(alloc_lev_L)then
         allocate(value(len),stat=istat)
         if (istat /= 0) then
            write(for_msg,*) 'unable to allocate space for '//trim(key_S)//' request '//trim(msg_S)
            call msg(MSG_CRITICAL,for_msg)
         endif
      else
         istat=0
      endif
   end function get_allocate_r1d
   
   integer function get_allocate_r81d(key_S,value,len,allow_reshape_L,msg_S) result(istat)
      ! Allocate space for the result value and report error
      implicit none
      character(len=*), intent(in) :: key_S
      real(kind=8), dimension(:), pointer :: value
      integer, intent(in) :: len
      logical :: allow_reshape_L
      character(len=*) :: msg_S
      !Local variables
      logical :: alloc_lev_L
      external msg
      istat=-1
      alloc_lev_L=.false.
      if(.not.associated(value))then
         alloc_lev_L=.true.
      else
         if(size(value)/=len)then
            if(allow_reshape_L)then
               write(for_msg,*) 'reshaping 1D real(kind=8) vector '//trim(msg_S)
               call msg(MSG_INFO,VGD_PRFX//for_msg)
               deallocate(value)
               alloc_lev_L=.true.
            else
               write(for_msg,*) '1D pointer already allocated with a different length, will not reallocate '//trim(msg_S)
               call msg(MSG_ERROR,VGD_PRFX//for_msg)
               return
            endif
         endif
      endif
      if(alloc_lev_L)then
         allocate(value(len),stat=istat)
         if (istat /= 0) then
            write(for_msg,*) 'unable to allocate space for '//trim(key_S)//' request '//trim(msg_S)
            call msg(MSG_CRITICAL,for_msg)
         endif
      else
       istat=0
    endif
 end function get_allocate_r81d

  logical function same_vec_i(vec1,vec2) result(equal)
    ! Check for equality between a pair of pointer vectors
    implicit none
    integer, dimension(:), pointer :: vec1,vec2          !Pointer vectors to compare
    integer :: i
    equal = .false.
    if (associated(vec1)) then
       if (associated(vec2)) then
          if (size(vec1) == size(vec2)) then
             do i=1,size(vec1)
                if (vec1(i) /= vec2(i)) return
             enddo
          endif
       else
          return
       endif
    else
       if (associated(vec2)) return
    endif
    equal = .true.
    return
  end function same_vec_i

  logical function same_vec_r(vec1,vec2) result(equal)
    ! Check for equality between a pair of pointer vectors
    implicit none
    real, dimension(:), pointer :: vec1,vec2     !Pointer vectors to compare
    integer :: i
    equal = .false.
    if (associated(vec1)) then
       if (associated(vec2)) then
          if (size(vec1) == size(vec2)) then
             do i=1,size(vec1)
                if (vec1(i) /= vec2(i)) return
             enddo
          endif
       else
          return
       endif
    else
       if (associated(vec2)) return
    endif
    equal = .true.
    return
  end function same_vec_r
 
  logical function same_vec_r8(vec1,vec2) result(equal)
    ! Check for equality between a pair of pointer vectors
    implicit none
    real(kind=8), dimension(:), pointer :: vec1,vec2     !Pointer vectors to compare
    integer :: i
    equal = .false.
    if (associated(vec1)) then
       if (associated(vec2)) then
          if (size(vec1) == size(vec2)) then
             do i=1,size(vec1)
                if (vec1(i) /= vec2(i)) return
             enddo
          endif
       else
          return
       endif
    else
       if (associated(vec2)) return
    endif
    equal = .true.
    return
  end function same_vec_r8

  logical function same_vec_r83d(vec1,vec2) result(equal)
    ! Check for equality between a pair of pointer vectors
    implicit none
    real(kind=8), dimension(:,:,:), pointer :: vec1,vec2 !Pointer vectors to compare
    integer :: i,j,k
    equal = .false.
    if (associated(vec1)) then
       if (associated(vec2)) then
          if (size(vec1) == size(vec2)) then
             do k=1,size(vec1,dim=3)
                do j=1,size(vec1,dim=2)
                   do i=1,size(vec1,dim=1)
                      if (vec1(i,j,k) /= vec2(i,j,k)) return
                   enddo
                enddo
             enddo
          endif
       else
          return
       endif
    else
       if (associated(vec2)) return
    endif
    equal = .true.
    return
  end function same_vec_r83d

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Convert string to upper-case

  function up (string) result(upper_string)
    ! Convert a string to all upper-case
    implicit none
    character(len=*), intent(in) :: string      !Input string to upper-case
    character(len=LONG_STRING) :: upper_string  !Upper-cased result
    integer :: i
    external msg
    if (len_trim(string) > len(upper_string)) then
       write(for_msg,*) 'Long string truncated in up() ',trim(string)
       call msg(MSG_WARNING,for_msg)
    endif
    upper_string = string
    do i = 1,len_trim(string)
       if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
          upper_string(i:i) = achar(iachar(string(i:i)) - 32)
       endif
    enddo
    return
  end function up

 integer function get_allocate_r3d(key_S,value,len,allow_reshape_L,msg_S) result(istat)
    ! Allocate space for the result value and report error (len is result of 'shape()')
    implicit none
    character(len=*), intent(in) :: key_S
    real, dimension(:,:,:), pointer :: value
    integer, dimension(:), intent(in) :: len
    logical :: allow_reshape_L
    character(len=*) :: msg_S
    !Local variables
    logical :: alloc_lev_L
    external msg
    istat=-1
    if (size(len) < 3) then
       write(for_msg,*) 'wrong array shape specified for '//trim(key_S)
       call msg(MSG_CRITICAL,for_msg)
       return
    endif
    alloc_lev_L=.false.
    if(.not.associated(value))then
       alloc_lev_L=.true.
    else
       if(  size(value,1)/=len(1).or.&
            size(value,2)/=len(2).or.&
            size(value,3)/=len(3))then
          if(allow_reshape_L)then
             write(for_msg,*) 'reshaping 3D real table'//trim(msg_S)
             call msg(MSG_INFO,VGD_PRFX//for_msg)
             deallocate(value)
             alloc_lev_L=.true.
          else
             write(for_msg,*) '3D pointer already allocated with a different length, will not reallocate '//trim(msg_S)
             call msg(MSG_ERROR,VGD_PRFX//for_msg)
             return
          endif
       endif
    endif
    if(alloc_lev_L)then
       allocate(value(len(1),len(2),len(3)),stat=istat)
       if (istat /= 0) then
          write(for_msg,*) 'unable to allocate space for '//trim(key_S)//' request '//trim(msg_S)
          call msg(MSG_CRITICAL,for_msg)
       endif
    else
       istat=0
    endif
 end function get_allocate_r3d

 integer function get_allocate_r83d(key_S,value,len,allow_reshape_L,msg_S) result(istat)
    ! Allocate space for the result value and report error (len is result of 'shape()')
    implicit none
    character(len=*), intent(in) :: key_S
    real(kind=8), dimension(:,:,:), pointer :: value
    integer, dimension(:), intent(in) :: len
    logical :: allow_reshape_L
    character(len=*) :: msg_S
    !Local variables
    logical :: alloc_lev_L
    external msg
    istat=-1
    if (size(len) < 3) then
       write(for_msg,*) 'wrong array shape specified for '//trim(key_S)
       call msg(MSG_CRITICAL,for_msg)
       return
    endif
    alloc_lev_L=.false.
    if(.not.associated(value))then
       alloc_lev_L=.true.
    else
       if(  size(value,1)/=len(1).or.&
            size(value,2)/=len(2).or.&
            size(value,3)/=len(3))then
          if(allow_reshape_L)then
             write(for_msg,*) 'reshaping 3D real(kind=8) table'//trim(msg_S)
             call msg(MSG_INFO,VGD_PRFX//for_msg)
             deallocate(value)
             alloc_lev_L=.true.
          else
             write(for_msg,*) '3D pointer already allocated with a different length, will not reallocate '//trim(msg_S)
             call msg(MSG_ERROR,VGD_PRFX//for_msg)
             return
          endif
       endif
    endif
    if(alloc_lev_L)then
       allocate(value(len(1),len(2),len(3)),stat=istat)
       if (istat /= 0) then
          write(for_msg,*) 'unable to allocate space for '//trim(key_S)//' request '//trim(msg_S)
          call msg(MSG_CRITICAL,for_msg)
       endif
    else
       istat=0
    endif
 end function get_allocate_r83d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get/Put support functions
  
  real function get_error(key,quiet) result(value)
    ! Write error message and return a missing value     
    implicit none
    character(len=*), intent(in) :: key
    logical, optional, intent(in) :: quiet      !Do not print massages
    ! Local variables
    integer :: level_msg
    external msg
    level_msg=MSG_CRITICAL
    if (present(quiet)) then
       if(quiet)level_msg=MSG_QUIET    
    endif
    write(for_msg,*) 'Attempt to retrieve invalid key '//trim(key)//' returns VGD_MISSING'
    call msg(level_msg,for_msg)
    value = dble(VGD_MISSING)
    return
  end function get_error

  integer function put_error(key) result(error)
    character(len=*), intent(in) :: key
    external msg
    write(for_msg,*) 'WARNING: attempt to set useless value for '//trim(key)
    call msg(MSG_CRITICAL,for_msg)
    error = VGD_ERROR
    return
  end function put_error

end module vgrid_utils
