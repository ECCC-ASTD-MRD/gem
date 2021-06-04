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
! interface to the routines used to convert between Fortran and C 
! communicators, data types, operators, ....
!
interface RPN_COMM_c2f
  function type_c2f(c) result(ftn) bind(C,name='RPN_COMM_type_c2f')
    import :: rpncomm_datatype, C_INT
    implicit none
    type(rpncomm_datatype), intent(IN) :: c
    integer(C_INT) :: ftn
  end function type_c2f

  function group_c2f(c) result(ftn) bind(C,name='RPN_COMM_group_c2f')
    import :: rpncomm_group, C_INT
    implicit none
    type(rpncomm_group), intent(IN) :: c
    integer(C_INT) :: ftn
  end function group_c2f

  function win_c2f(c) result(ftn) bind(C,name='RPN_COMM_window_c2f')
    import :: rpncomm_window, C_INT
    implicit none
    type(rpncomm_window), intent(IN) :: c
    integer(C_INT) :: ftn
  end function win_c2f

  function op_c2f(c) result(ftn) bind(C,name='RPN_COMM_op_c2f')
    import :: rpncomm_operator, C_INT
    implicit none
    type(rpncomm_operator), intent(IN) :: c
    integer(C_INT) :: ftn
  end function op_c2f

  function file_c2f(c) result(ftn) bind(C,name='RPN_COMM_file_c2f')
    import :: rpncomm_file, C_INT
    implicit none
    type(rpncomm_file), intent(IN) :: c
    integer(C_INT) :: ftn
  end function file_c2f

  function comm_c2f(c) result(ftn) bind(C,name='RPN_COMM_comm_c2f')
    import :: rpncomm_communicator, C_INT
    implicit none
    type(rpncomm_communicator), intent(IN) :: c
    integer(C_INT) :: ftn
  end function comm_c2f

  function req_c2f(c) result(ftn) bind(C,name='RPN_COMM_request_c2f')
    import :: rpncomm_request, C_INT
    implicit none
    type(rpncomm_request), intent(IN) :: c
    integer(C_INT) :: ftn
  end function req_c2f

  function info_c2f(c) result(ftn) bind(C,name='RPN_COMM_info_c2f')
    import :: rpncomm_info, C_INT
    implicit none
    type(rpncomm_info), intent(IN) :: c
    integer(C_INT) :: ftn
  end function info_c2f
end interface

interface RPN_COMM_f2c
  function type_f2c(ftn,mpi_fval) result(c) bind(C,name='RPN_COMM_type_f2c')
    import :: C_PTR, rpncomm_datatype, C_INT
    implicit none
    type(rpncomm_datatype), intent(INOUT) :: ftn
    integer(C_INT), intent(IN), value :: mpi_fval
    type(C_PTR) :: c
  end function type_f2c

  function group_f2c(ftn,mpi_fval) result(c) bind(C,name='RPN_COMM_group_f2c')
    import :: C_PTR, rpncomm_group, C_INT
    implicit none
    type(rpncomm_group), intent(INOUT) :: ftn
    integer(C_INT), intent(IN), value :: mpi_fval
    type(C_PTR) :: c
  end function group_f2c

  function win_f2c(ftn,mpi_fval) result(c) bind(C,name='RPN_COMM_window_f2c')
    import :: C_PTR, rpncomm_window, C_INT
    implicit none
    type(rpncomm_window), intent(INOUT) :: ftn
    integer(C_INT), intent(IN), value :: mpi_fval
    type(C_PTR) :: c
  end function win_f2c

  function op_f2c(ftn,mpi_fval) result(c) bind(C,name='RPN_COMM_op_f2c')
    import :: C_PTR, rpncomm_operator, C_INT
    implicit none
    type(rpncomm_operator), intent(INOUT) :: ftn
    integer(C_INT), intent(IN), value :: mpi_fval
    type(C_PTR) :: c
  end function op_f2c

  function file_f2c(ftn,mpi_fval) result(c) bind(C,name='RPN_COMM_file_f2c')
    import :: C_PTR, rpncomm_file, C_INT
    implicit none
    type(rpncomm_file), intent(INOUT) :: ftn
    integer(C_INT), intent(IN), value :: mpi_fval
    type(C_PTR) :: c
  end function file_f2c

  function comm_f2c(ftn,mpi_fval) result(c) bind(C,name='RPN_COMM_comm_f2c')
    import :: C_PTR, rpncomm_communicator, C_INT
    implicit none
    type(rpncomm_communicator), intent(INOUT) :: ftn
    integer(C_INT), intent(IN), value :: mpi_fval
    type(C_PTR) :: c
  end function comm_f2c

  function req_f2c(ftn,mpi_fval) result(c) bind(C,name='RPN_COMM_request_f2c')
    import :: C_PTR, rpncomm_request, C_INT
    implicit none
    type(rpncomm_request), intent(INOUT) :: ftn
    integer(C_INT), intent(IN), value :: mpi_fval
    type(C_PTR) :: c
  end function req_f2c

  function info_f2c(ftn,mpi_fval) result(c) bind(C,name='RPN_COMM_info_f2c')
    import :: C_PTR, rpncomm_info, C_INT
    implicit none
    type(rpncomm_info), intent(INOUT) :: ftn
    integer(C_INT), intent(IN), value :: mpi_fval
    type(C_PTR) :: c
  end function info_f2c
end interface

