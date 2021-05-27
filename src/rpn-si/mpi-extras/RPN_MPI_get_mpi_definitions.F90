!/! RPN_MPI - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2020  Division de Recherche en Prevision Numerique
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
module RPN_MPI_mpi_layout
  use ISO_C_BINDING
  use rpn_mpi_mpif
  implicit none
#include <RPN_MPI_mpi_symbols.hf>
#include <RPN_MPI_mpi_layout.hf>
#include <RPN_MPI_nulls.hf>

  ! lr is statically initialized, lw will be initialized by RPN_MPI_init_mpi_layout
  type(mpi_layout_r), save, target  :: lr = mpi_layout_r_NULL
  type(mpi_layout_f), save, pointer :: lw => NULL()  ! wrapped version

  ! ml is statically initialized, mw will be initialized by RPN_MPI_init_mpi_layout
  type(mpi_layout_internal), save, target :: ml = &    ! RPN_MPI communicators, ranks, sizes
    mpi_layout_internal( &
      layout_version, &
      -1, -1, &                               ! host, numa node
      [-1, -1, -1], &                         ! colors
      grid_hierarchy( &                       ! communicators
        application(MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL), &
        application(MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL), &
        mpigrid(MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, &
                MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, &
                MPI_COMM_NULL, MPI_COMM_NULL), &
        mpigrid(MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, &
                MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL, &
                MPI_COMM_NULL, MPI_COMM_NULL), &
        subgrid(MPI_COMM_NULL, MPI_COMM_NULL, MPI_COMM_NULL) &
      ), &
      grid_hierarchy( &                       ! ranks
        application(-1, -1, -1), &
        application(-1, -1, -1), &
        mpigrid(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1), &
        mpigrid(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1), &
        subgrid(-1, -1, -1) &
      ), &
      grid_hierarchy( &                       ! sizes
        application(-1, -1, -1), &
        application(-1, -1, -1), &
        mpigrid(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1), &
        mpigrid(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1), &
        subgrid(-1, -1, -1) &
      ) &
    )
  type(mpi_layout), save, pointer         :: mw => NULL()    ! communicators, ranks, sizes, wrapped

  ! MPI constants, dr is statically initialized, dw will be initialized by RPN_MPI_init_mpi_layout
  type(RPN_MPI_mpi_definitions_raw), save, target :: dr    = RPN_MPI_mpi_definitions_raw( &
    mpi_symbols_version, &
    MPI_GROUP_NULL, MPI_REQUEST_NULL,MPI_ERRHANDLER_NULL, MPI_INFO_NULL,  MPI_WIN_NULL, &
    MPI_STATUS_SIZE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_SUCCESS, MPI_ERROR, &
    MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN, &
    MPI_COMM_NULL, MPI_COMM_WORLD, MPI_COMM_SELF, MPI_GROUP_EMPTY, MPI_COMM_TYPE_SHARED, &
    MPI_DATATYPE_NULL, &
    MPI_BYTE, MPI_PACKED, MPI_UB, MPI_LB, &
    MPI_CHARACTER, MPI_LOGICAL, MPI_INTEGER, MPI_INTEGER1, MPI_INTEGER2, MPI_INTEGER4, &
    MPI_INTEGER8, MPI_INTEGER16, MPI_REAL, MPI_REAL4, MPI_REAL8, MPI_REAL16, &
    MPI_DOUBLE_PRECISION, MPI_COMPLEX, MPI_COMPLEX8, MPI_COMPLEX16, MPI_COMPLEX32, &
    MPI_DOUBLE_COMPLEX, MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER, &
    MPI_OP_NULL, &
    MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_LAND, MPI_BAND, &
    MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MAXLOC, MPI_MINLOC, MPI_REPLACE, &
    MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE &
    )
  type(RPN_MPI_mpi_definitions), save, pointer :: dw => NULL()   ! MPI constants, "wrapped"
contains
!
! make "wrapped" definitions available by pointing them to the "raw" definitions
!
! - MUST BE CALLED ASAP by RPN_MPI_init
! - will get called by RPN_MPI_get_mpi_definitions.../RPN_MPI_get_mpi_layout...
!   if dw or mw is not "associated" at time of call
!
  subroutine RPN_MPI_init_mpi_layout
    implicit none
#include <RPN_MPI_system_interfaces.hf>
    integer :: cpu
    type(C_PTR) :: p

    p = C_LOC(dr)                                     ! definitions
    if(.not. associated(dw)) call C_F_POINTER(p, dw)  ! point to "raw" version

!     ml%version = layout_version
    ml%host = get_host_id()              ! get linux host id
    cpu = sched_get_my_cpu()             ! get logical core number
    ml%numa = numa_node(cpu)             ! get numa space number for this core

    if(lr%host == -1) then               ! hostid not initialized
!     lr%version = layout_version
      lr%host = get_host_id()            ! get linux host id
      cpu = sched_get_my_cpu()           ! get logical core number
      lr%numa = numa_node(cpu)           ! get numa space number for this core
    endif
    p = C_LOC(lr)                                     ! layout
    if(.not. associated(lw)) call C_F_POINTER(p, lw)  ! point to "raw" version

    p = C_LOC(ml)                                     ! layout
    if(.not. associated(mw)) call C_F_POINTER(p, mw)  ! point to "raw" version
    return
  end subroutine RPN_MPI_init_mpi_layout
end module RPN_MPI_mpi_layout

! get access to the most common MPI parameters, like MPI_COMM_WORLD, MPI_COMM_NULL, etc...
! useful to make calls to MPI library routines directly without having to include mpif.f or equivalent
! thereby creating a dependency for the code upon the MPI implementation (insulation layer)
!
! RPN_MPI_get_mpi_definitions_raw gets a copy of the "not wrapped, raw" parameters
! RPN_MPI_get_mpi_definitions     gets a copy of the "wrapped" parameters, useful for interface enforcement
!
 subroutine RPN_MPI_get_mpi_definitions_raw(what, ierr)          !InTf! ! get a copy of raw MPI definitions
  use ISO_C_BINDING
  use RPN_MPI_mpi_layout
  implicit none
!!  import :: RPN_MPI_mpi_definitions_raw                        !InTf! 
  type(RPN_MPI_mpi_definitions_raw), intent(INOUT) :: what      !InTf! 
  integer, intent(OUT) :: ierr                                  !InTf! 

  if(.not. associated(dw)) call RPN_MPI_init_mpi_layout

  ierr = MPI_ERROR
  if(what%version .ne. mpi_symbols_version) then
    write(0,*)'ERROR: (RPN_MPI_get_mpi_definitions_raw) version mismatch, expected :',mpi_symbols_version,' got :',what%version
    what%version = -1
    return
  endif
  ierr = MPI_SUCCESS
  what = dr
  return
 end subroutine RPN_MPI_get_mpi_definitions_raw          !InTf!

 subroutine RPN_MPI_get_mpi_definitions(what, ierr)          !InTf! ! get a copy of wrapped MPI definitions
  use ISO_C_BINDING
  use RPN_MPI_mpi_layout
  implicit none
!!  import :: RPN_MPI_mpi_definitions                        !InTf! 
  type(RPN_MPI_mpi_definitions), intent(INOUT) :: what      !InTf! 
  integer, intent(OUT) :: ierr                              !InTf! 

  ierr = MPI_ERROR
  if(what%version .ne. mpi_symbols_version) then
    write(0,*)'ERROR: (RPN_MPI_get_mpi_definitions) version mismatch, expected :',mpi_symbols_version,' got :',what%version
    what%version = -1
    return
  endif
  ierr = MPI_SUCCESS
  what = transfer(dr,dw)

  return
 end subroutine RPN_MPI_get_mpi_definitions                  !InTf! 

! RPN_MPI_get_mpi_layout_raw : same as RPN_MPI_get_mpi_layout 
! but advertized with type mpi_layout_internal rather than mpi_layout
!
! get a copy of RPN_MPI internal mpi layout (communicators, ranks, sizes)
!
! RPN_MPI_get_mpi_layout_raw gets a copy of the "not wrapped, raw" information (communicators)
! RPN_MPI_get_mpi_layout     gets a copy of the "wrapped" information, useful for interface enforcement
!
!! end interface                                             !InTf!
!! interface RPN_MPI_get_mpi_layout_raw                      !InTf! 
!
! OLD layout
!
 subroutine RPN_MPI_get_mpi_layout_raw1(what, ierr)          !InTf! 
  use RPN_MPI_mpi_layout
  implicit none
!!  import :: mpi_layout_internal, C_INT                     !InTf! 
  type(mpi_layout_internal), intent(INOUT) :: what           !InTf! 
  integer(C_INT), intent(OUT) :: ierr                        !InTf! 
  call RPN_MPI_get_mpi_layout1(what, ierr)
  return
 end subroutine RPN_MPI_get_mpi_layout_raw1                  !InTf! 
!
! NEW layout
!
 subroutine RPN_MPI_get_mpi_layout_raw2(what, ierr)          !InTf! 
  use RPN_MPI_mpi_layout
  implicit none
!!  import :: mpi_layout_r, C_INT                            !InTf! 
  type(mpi_layout_r), intent(INOUT) :: what                  !InTf! 
  integer(C_INT), intent(OUT) :: ierr                        !InTf! 
  call RPN_MPI_get_mpi_layout2(what, ierr)
  return
 end subroutine RPN_MPI_get_mpi_layout_raw2                  !InTf! 
!! end interface                                             !InTf!
!! interface                                                 !InTf!

! get a copy of RPN_MPI internal mpi layout (communicators, ranks, sizes)
!
!! end interface                                             !InTf!
!! interface RPN_MPI_get_mpi_layout                          !InTf! 
!
! OLD layout
!
 subroutine RPN_MPI_get_mpi_layout1(what, ierr)              !InTf! 
  use RPN_MPI_mpi_layout
  implicit none
!!  import :: mpi_layout, C_INT                              !InTf! 
! white lie in published interface, 
! "what" is treated as the "internal" type rather than the "wrapped" type
! but it is advertized to the user as "wrapped"
!!  type(mpi_layout), intent(INOUT) :: what                  !InTf! 
  type(mpi_layout_internal), intent(INOUT) :: what
  integer(C_INT), intent(OUT) :: ierr                        !InTf! 

  ierr = MPI_ERROR
  if(.not. associated(mw)) call RPN_MPI_init_mpi_layout

  if(what%version .ne. layout_version) then
    write(0,*)'ERROR: (RPN_MPI_get_mpi_layout) version mismatch, expected :',layout_version,' got :',what%version
    what%version = -1
    return
  endif
  ierr = MPI_SUCCESS

  what = ml
  what%comm%sgrd%row       = MPI_COMM_NULL   ! row is not defined for supergrids
  what%rank%sgrd%row       = -1
  what%size%sgrd%row       = -1
  what%comm%sgrd%column    = MPI_COMM_NULL   ! and neither is column
  what%rank%sgrd%column    = -1
  what%size%sgrd%column    = -1
  return
 end subroutine RPN_MPI_get_mpi_layout1                      !InTf! 
!
! NEW layout
!
 subroutine RPN_MPI_get_mpi_layout2(what, ierr)              !InTf! 
  use RPN_MPI_mpi_layout
  implicit none
!!  import :: mpi_layout_f, C_INT                            !InTf! 
! white lie in published interface, 
! "what" is treated as the "internal" type rather than the "wrapped" type
! but it is advertized to the user as "wrapped"
!!  type(mpi_layout_f), intent(INOUT) :: what                !InTf! 
  type(mpi_layout_r), intent(INOUT) :: what
  integer(C_INT), intent(OUT) :: ierr                        !InTf! 

  ierr = MPI_ERROR
  if(.not. associated(mw)) call RPN_MPI_init_mpi_layout

  if(what%version .ne. layout_version) then
    write(0,*)'ERROR: (RPN_MPI_get_mpi_layout) version mismatch, expected :',layout_version,' got :',what%version
    what%version = -1
    return
  endif
  ierr = MPI_SUCCESS

  what = lr
  what%sgrd%row%comm       = MPI_COMM_NULL   ! row is not defined for supergrids
  what%sgrd%row%rank       = -1
  what%sgrd%row%size       = -1
  what%sgrd%column%comm    = MPI_COMM_NULL   ! and neither is column
  what%sgrd%column%rank    = -1
  what%sgrd%column%size    = -1
  return
 end subroutine RPN_MPI_get_mpi_layout2                      !InTf! 
!! end interface                                             !InTf!
!! interface                                                 !InTf!


