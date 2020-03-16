!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
!
! public object types used for user application interface
!
! all these objects must be considered OPAQUE, their contents are managed 
! by the appropriate RPN_COMM_... routines and MUST NOT be altered by the user
!
  integer, parameter :: RPNCOMM_types_version = 1

  type, bind(C)    :: rpncomm_context
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_context
  type(rpncomm_context), parameter :: NULL_rpncomm_context = rpncomm_context(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_window
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_window
  type(rpncomm_window), parameter :: NULL_rpncomm_window = rpncomm_window(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_communicator
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_communicator
  type(rpncomm_communicator), parameter :: NULL_rpncomm_communicator = rpncomm_communicator(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_info
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_info
  type(rpncomm_info), parameter :: NULL_rpncomm_info = rpncomm_info(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_file
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_file
  type(rpncomm_file), parameter :: NULL_rpncomm_file = rpncomm_file(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_request
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_request
  type(rpncomm_request), parameter :: NULL_rpncomm_request = rpncomm_request(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_datatype
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_datatype
  type(rpncomm_datatype), parameter :: NULL_rpncomm_datatype = rpncomm_datatype(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_operator
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_operator
  type(rpncomm_operator), parameter :: NULL_rpncomm_operator = rpncomm_operator(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_group
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_group
  type(rpncomm_group), parameter :: NULL_rpncomm_group = rpncomm_group(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_ptr
    type(c_ptr)    :: p     ! managed by rpn_comm routines
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_ptr
  type(rpncomm_ptr), parameter :: NULL_rpncomm_ptr = rpncomm_ptr(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_pattern
    type(c_ptr)    :: p     ! pointer to data description
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_pattern
  type(rpncomm_pattern), parameter :: NULL_rpncomm_pattern = rpncomm_pattern(C_NULL_PTR,0,0)

  type, bind(C)    :: rpncomm_array
    type(c_ptr)    :: p     ! pointer to data description
    integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
    integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_array
  type(rpncomm_array), parameter :: NULL_rpncomm_array = rpncomm_array(C_NULL_PTR,0,0)
!
! add interfaces to fortran <-> c conversion routines
!
include 'RPN_COMM_ftoc.inc'

