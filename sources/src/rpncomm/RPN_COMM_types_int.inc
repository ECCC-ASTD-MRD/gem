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
!====================================================================================
! internal, private datatypes used by RPN_COMM routines
!
! a commented version of some public types is also present, with
! comments describing their internal usage
!====================================================================================
!   type, bind(C)    :: rpncomm_communicator ! (see RPN_COMM_comm.F90)
!     type(c_ptr)    :: p     ! signature ( to be changed to C MPI communicator)
!     integer(c_int) :: t1    ! t2 .xor. RPN_COMM_MAGIC
!     integer(c_int) :: t2    ! Fortran MPI communicator
!   end type rpncomm_communicator
!   type(rpncomm_communicator), parameter :: NULL_rpncomm_communicator = rpncomm_communicator(C_NULL_PTR,0,0)
!====================================================================================
!   type, bind(C)    :: rpncomm_datatype ! (see RPN_COMM_datyp.F90)
!     type(c_ptr)    :: p     ! signature ( to be changed to C MPI data type)
!     integer(c_int) :: t1    ! tag1 (index into data type table)
!     integer(c_int) :: t2    ! Fortran MPI data type
!   end type rpncomm_datatype
!   type(rpncomm_datatype), parameter :: NULL_rpncomm_datatype = rpncomm_datatype(C_NULL_PTR,0,0)
!====================================================================================
!   type, bind(C)    :: rpncomm_pattern  !  (see RPN_COMM_spread.F90)
!     type(c_ptr)    :: p     ! pointer to data description
!     integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
!     integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
!   end type rpncomm_pattern
!====================================================================================
!   type, bind(C)    :: rpncomm_array
!     type(c_ptr)    :: p     ! pointer to data description
!     integer(c_int) :: t1    ! tag1 (reserved for rpn_comm usage)
!     integer(c_int) :: t2    ! tag2 (reserved for rpn_comm usage)
!   end type rpncomm_array
!====================================================================================
!  type, bind(C)    :: rpncomm_window
!    type(c_ptr)    :: p     ! pointer to table entry  = C_LOC( win_tab(t2) )
!    integer(c_int) :: t1    ! t2 .xor. RPN_COMM_MAGIC
!    integer(c_int) :: t2    ! index into internal window table
!  end type rpncomm_window
!====================================================================================
!  type, bind(C)    :: rpncomm_operator
!    type(c_ptr)    :: p     ! signature ( to be changed to C MPI communicator)
!    integer(c_int) :: t1    ! t2 .xor. RPN_COMM_MAGIC
!    integer(c_int) :: t2    ! Fortran MPI operator
!  end type rpncomm_operator
!====================================================================================
  type :: rpncomm_windef
    integer, dimension(:), pointer :: remote   ! base translated as Fortran pointer to 1D array
    type(C_PTR) :: base    ! base address of array (C type pointer)
    integer :: indx        ! index of this entry in window table
    integer :: opr         ! MPI operator
    integer :: typ         ! MPI datatype
    integer :: ext         ! size of MPI datatype  (extent)
    integer :: win         ! MPI window
    integer :: com         ! MPI communicator
    integer :: grp         ! MPI group associated with communicator
    integer :: s_group     ! MPI send (put) group
    integer :: r_group     ! MPI receive (get) group
    integer :: siz         ! deliberately limiting window size to less than 2GB
    logical :: is_open     ! .true. if window is "exposed"
    logical :: is_user     ! .true. if array associated with window was user supplied rather than internal
    logical :: active_mode ! use active mode for this window's communications
  end type rpncomm_windef
  type(rpncomm_windef), parameter :: NULL_rpncomm_windef = rpncomm_windef(NULL(),C_NULL_PTR,0,0,0,0,0,0,0,0,0,0,.false.,.false.,.true.)
!====================================================================================
  type, bind(C)    :: rpncomm_bnds
    integer(c_int) :: lo, hi, lnp, gnp  ! low bound, high bound, local dimension, global dimension
  end type rpncomm_bnds
  type(rpncomm_bnds), parameter :: NULL_rpncomm_bnds = rpncomm_bnds(0,0,0,0)
!====================================================================================
  type, bind(C)    :: rpncomm_field     ! 1/2/3D array, represented as as 3D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: sz                ! base data element size in bytes
    type(rpncomm_bnds) :: x, y, z       ! dimension bounds (up to 3 dimensions)
    integer(c_int) :: hx, hy, hz        ! halos for exchange
    integer(c_int) :: ew3d, ns3d        ! user defined data types for halo exchange
  end type rpncomm_field
  type(rpncomm_field), parameter :: NULL_rpncomm_field = rpncomm_field(         &
                       C_NULL_PTR, C_NULL_PTR, 0,                               &
                       NULL_rpncomm_bnds, NULL_rpncomm_bnds, NULL_rpncomm_bnds, &
                       0, 0, 0, 0, 0 )
!====================================================================================
  type, bind(C)    :: rpncomm_vector    ! 1D vector
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x                 ! dimension
  end type rpncomm_vector
  type(rpncomm_vector), parameter :: NULL_rpncomm_vector = rpncomm_vector(C_NULL_PTR, C_NULL_PTR, 0)
!====================================================================================
  type, bind(C)    :: rpncomm_matrix     ! 2D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x, y              ! dimensions
  end type rpncomm_matrix
  type(rpncomm_matrix), parameter :: NULL_rpncomm_matrix = rpncomm_matrix(C_NULL_PTR, C_NULL_PTR, 0, 0)
!====================================================================================
  type, bind(C)    :: rpncomm_block     ! 3D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x, y, z           ! dimensions
  end type rpncomm_block
  type(rpncomm_block), parameter :: NULL_rpncomm_block = rpncomm_block(C_NULL_PTR, C_NULL_PTR, 0, 0, 0)
!====================================================================================
!
! Domain type
!
  type DOMM
     sequence
     character(len=12) nom
     character(len=1024) path
     integer npex, npey
  end type
!====================================================================================
!
!  symbol tables
!
  type :: SYMTAB
    integer :: number
    character (len=32) :: string
  end type
!====================================================================================
!
! Decomposition tables
!
  type :: DEC         ! decomposisition description
    integer :: id       ! "magic" pseudo unique identifier
    integer :: l1, l2   ! number of points per section
    integer :: ns       ! start of shorter tiles (1->np) (if ns .lt. np,  l2 must be l1-1)
    integer :: lg       ! total number of points
    integer :: np       ! number of tiles along axis
  end type DEC
  
  type :: DIST_1D     ! 1D type decomposition
    integer :: id       ! "magic" pseudo unique identifier
    integer :: a        ! DEC id
    integer    :: com   ! communicator (row, column, ...)
    integer    :: grp   ! group
  end type DIST_1D

  type :: DIST_2D     ! 2D type decomposition
    integer :: id       ! "magic" pseudo unique identifier
    integer :: idx,idy  ! E-W and N-S DIST_1D id's
    integer    :: com   ! communicator (grid, subgrid, ...)
    integer    :: grp   ! group
  end type DIST_2D
