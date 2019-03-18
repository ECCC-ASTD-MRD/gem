!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
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
! */
      include 'RPN_COMM_constants.inc'
!
!     integer,parameter :: RPN_COMM_MASTER = 0
!
!      character(len=*),parameter :: RPN_COMM_WORLD      = 'WORLD'      ! PLS do not use
!      character(len=*),parameter :: RPN_COMM_ALLGRIDS   = 'ALLGRIDS'   ! use rather than RPN_COMM_WORLD
!      character(len=*),parameter :: RPN_COMM_MULTIGRID  = 'MULTIGRID'
!      character(len=*),parameter :: RPN_COMM_GRID       = 'GRID'       ! use rather than RPN_COMM_DOMM
!      character(len=*),parameter :: RPN_COMM_DEFAULT    = 'DEFO'
!      character(len=*),parameter :: RPN_COMM_BLOCMASTER = 'BLOCMASTER'
!      character(len=*),parameter :: RPN_COMM_BLOC_COMM  = 'BLOC'
!      character(len=*),parameter :: RPN_COMM_EW         = 'EW'
!      character(len=*),parameter :: RPN_COMM_ROW        = 'EW'
!      character(len=*),parameter :: RPN_COMM_NS         = 'NS'
!      character(len=*),parameter :: RPN_COMM_COLUMN     = 'NS'
!
!      character(len=*),parameter :: RPN_COMM_INTEGER             =       &
!     &                             'MPI_INTEGER'
!      character(len=*),parameter :: RPN_COMM_INTEGER2            =       &
!     &                             'MPI_INTEGER2'
!      character(len=*),parameter :: RPN_COMM_COMPLEX             =       &
!     &                             'MPI_COMPLEX'
!      character(len=*),parameter :: RPN_COMM_DOUBLE_COMPLEX      =       &
!     &                             'MPI_DOUBLE_COMPLEX'
!      character(len=*),parameter :: RPN_COMM_DOUBLE_PRECISION    =       &
!     &                             'MPI_DOUBLE_PRECISION'
!      character(len=*),parameter :: RPN_COMM_REAL                =       &
!     &                             'MPI_REAL'
!      character(len=*),parameter :: RPN_COMM_REAL4               =       &
!     &                             'MPI_REAL4'
!      character(len=*),parameter :: RPN_COMM_REAL8               =       &
!     &                             'MPI_REAL8'
!      character(len=*),parameter :: RPN_COMM_LOGICAL             =       &
!     &                             'MPI_LOGICAL'
!      character(len=*),parameter :: RPN_COMM_CHARACTER           =       &
!     &                             'MPI_CHARACTER'
!
!      character(len=*),parameter :: RPN_COMM_SUM     = 'MPI_SUM'
!      character(len=*),parameter :: RPN_COMM_MIN     = 'MPI_MIN'
!      character(len=*),parameter :: RPN_COMM_MAX     = 'MPI_MAX'
!      character(len=*),parameter :: RPN_COMM_PROD    = 'MPI_PROD'
!      character(len=*),parameter :: RPN_COMM_LAND    = 'MPI_LAND'
!      character(len=*),parameter :: RPN_COMM_BAND    = 'MPI_BAND'
!      character(len=*),parameter :: RPN_COMM_LOR     = 'MPI_LOR'
!      character(len=*),parameter :: RPN_COMM_BOR     = 'MPI_BOR'
!      character(len=*),parameter :: RPN_COMM_LXOR    = 'MPI_LXOR'
!      character(len=*),parameter :: RPN_COMM_BXOR    = 'MPI_BXOR'
!      character(len=*),parameter :: RPN_COMM_OP_NULL = 'MPI_OP_NULL'
!      character(len=*),parameter :: RPN_COMM_MAXLOC  = 'MPI_MAXLOC'
!      character(len=*),parameter :: RPN_COMM_MINLOC  = 'MPI_MINLOC'
!
!      logical,parameter :: RPN_COMM_TOPO_X      = .true.
!      logical,parameter :: RPN_COMM_TOPO_Y      = .false.
!      logical,parameter :: RPN_COMM_TOPO_FILL   = .true.
!      logical,parameter :: RPN_COMM_TOPO_NOFILL = .false.
!
      interface
!
      subroutine rpn_comm_init(Userinit,Pelocal,Petotal,Pex,Pey)
! If not already started, initialize MPI mode and some common blocks 
! (internal use for RPN_COMM routines) and return the topology. 
! If pex and pey are greater than 0, those values will be used and 
! Userinit will not be called.
      external :: Userinit
      integer, intent(out)   :: Pelocal,Petotal
      integer, intent(inout) :: Pex,Pey
! Userinit: User specified subroutine which returns the process topology (subroutine userinit(pex,pey)
! Pelocal : rank of local process
! Petotal : total number of processes
! Pex     : number of PE along x axis: if 0, specified with Userinit
! Pey     : number of PE along y axis: if 0, specified with Userinit
      end subroutine rpn_comm_init
!
!
      integer function rpn_comm_init_multigrid                           &
     &   (Userinit,Pelocal,Petotal,Pex,Pey,ngrids)
! If not already started, initialize MPI mode and some common blocks 
! (internal use for RPN_COMM routines) and return the topology
! for a multi grid topology.
! If pex and pey are greater than 0, those values will be used and 
! Userinit will not be called.
      external :: Userinit
      integer, intent(out)   :: Pelocal,Petotal
      integer, intent(inout) :: Pex,Pey
      integer, intent(in)    :: ngrids
! Userinit: User specified subroutine which returns the process topology (subroutine userinit(pex,pey)
! Pelocal : rank of local process
! Petotal : total number of processes
! Pex     : number of PE along x axis: if 0, specified with Userinit
! Pey     : number of PE along y axis: if 0, specified with Userinit
! ngrids  : number of grids
      end function rpn_comm_init_multigrid
!
!
      integer function rpn_comm_mype(me,me_x,me_y)
! Return rank and position of the local process, in regards of his domain
      integer, intent(out) :: me,me_x,me_y
! me  : rank of process in its domain
! me_x: position along X axis in its domain
! me_y: position along Y axis in its domain
      end function rpn_comm_mype
!
!
      integer function rpn_comm_bloc(nblocx,nblocy)
! Creates the needed communicators to use blocks within a domain. 
! Those communicators are named "BLOC" for in-block communication and 
! "BLOCMASTER" for communication between blocks master PE. 
! Blocks are constructed by divinding the number of processes along 
! X by nblocx and the number of processes by nblocy. 
! Note that each block has to be identical, so it is impossible to 
! split 7 processes into 2 blocks (then RPN_COMM_bloc=-1 and 
! an error message appears).
      integer, intent(in) :: nblocx,nblocy
! nblocx,nblocy: Number of blocks in X and Y direction
      end function rpn_comm_bloc
!
!
      subroutine rpn_comm_carac(npex,npey,me,medomm,mex,mey,sizex,       &
     &                          sizey,ismaster, mymaster, mybloc,        &
     &                          myblocx,myblocy,blocme,domname)
! Returns local process info.
      integer, intent(out) :: npex,npey,me,medomm,mex,mey,sizex,sizey,   &
     &              ismaster, mymaster, mybloc, myblocx,myblocy,blocme
      character(len=*),intent(out) :: domname
! npex,npey: Number of PE along x and y
! me       : Rank in the context of communicator "ALL"
! medomm   : Rank in the context of communicator "DOMM"
! mex,mey  : x and y coordinates in domain
! sizex,sizey: Size of blocks along x and y axis
! ismaster : Equals 1 if PE is included in communicator "BLOCMASTER", 0 else
! mymaster : Rank (relative to the domain) of the master of the PE's block
! mybloc   : Rank of PE's block
! myblocx,myblocy: Coordinate of PE's block along x and y
! blocme   : Rank of PE relative to its block
! domname  : Domain name of the current PE
      end subroutine rpn_comm_carac
!
!
      subroutine rpn_comm_rank(domname,rank,ierr)
! Return rank of PE in specified domname/communicator
      character(len=*),intent(in) :: domname
      integer, intent(out) :: rank
      integer, intent(out) :: ierr
! domname : Domain/communicator name
! rank    : rank of process in domname
! ierr    : error code
      end subroutine rpn_comm_rank
!
!
      integer function rpn_comm_topo(nxg,minx,maxx,nxl,nxlmax,halox,     &
     &                               nx0,alongx,fill)
! Generate needed information about local tile along a specified axis. 
! The input is the total number of point to divide and the size of the halo. 
! The function will split the domain depending of the topology 
! induced by the PEs.
      integer, intent(in) :: nxg,halox
      logical, intent(in) :: alongx,fill
      integer, intent(out):: minx,maxx,nxl,nxlmax,nx0
! nxg    : number of points of global domain along one axis
! halox  : size of halo
! alongx : logical, true if the domain is divided along x
! fill   : logical, add one more point to maxx if nxlmax is even (for memory optimization purpose)
! minx maxx : dimensions of the local array needed to contain the local tile
! nxl    : number of local points along axis for this PE
! nxlmax : maximum of local points along axis for all PE
! nx0    : Position of the first element of the local PE over nxg
      end function rpn_comm_topo
!
!
      integer function rpn_comm_comm(domname)
! Makes the correspondance between RPN_COMM character communicators and MPI communicators... 
! i.e returns integer MPI communicator corresponding to RPN_COMM character communicator domname
! domname  : Domain name; RPN_COMM character communicator name
      character(len=*),intent(in) :: domname
      end function rpn_comm_comm
!
!
      integer function rpn_comm_datyp(datatype)
! Return the MPI integer datatype code
! corresponding to the RPN_COMM Datatype string
! datatype  :  RPN_COMM type name
      character(len=*),intent(in) :: datatype
      end function rpn_comm_datyp
!
!
      subroutine rpn_comm_barrier(comm, ierr)
! Stub for MPI_barrier routine using the same calling sequence, 
! except for datatypes, operators and communicators. 
! Those are replaced by an equivalent string, described in the table below:
! This RPN_COMM routine is a sync point for all PE's included in 
! the specified communicators. Its primary use is for debugging, 
! it may be useful in some other situations when a task has to be 
! terminated for everyone before starting a new one. 
      character(len=*), intent(in) :: comm
      integer, intent(out) :: ierr
! comm : communicator
! ierr : error code
      end subroutine rpn_comm_barrier
!
!
      subroutine rpn_comm_finalize( ierr)
!This RPN_COMM routine stops the MPI mode. It allows a MPI program to terminate gracefully.
      integer, intent(out) :: ierr
! ierr : error code
      end subroutine rpn_comm_finalize
!
!
      subroutine mpi_comm_split(old_comm, color, key, new_comm, ierr)
! Create a new communicator
      integer, intent(in)  :: old_comm, color, key
      integer, intent(out) :: new_comm, ierr
! old_comm : 
! color    :
! key      :
! new_comm :
! ierr     : error code
      end subroutine mpi_comm_split
!
!
      subroutine mpi_comm_rank(comm, rank, ierr)
! Get rank of PE in communicator comm
      integer, intent(in)  :: comm
      integer, intent(out) :: rank, ierr
! comm : communicator
! rank : rank of process in comm
! ierr : error code
      end subroutine mpi_comm_rank
!
        function RPN_COMM_spread_context(context,com,rootpe,pe,npts)    &
     &                                 result(status)
        use ISO_C_BINDING
        implicit none
          type(C_PTR), intent(OUT) :: context
          character (len=*), intent(IN) :: com             ! RPN_COMM communicator
          integer, intent(IN) :: npts                      ! number of data points
          integer, intent(IN) :: rootpe                    ! root PE for the spread operation
          integer, dimension(npts), intent(IN) :: pe       ! destination table, data point i will be sent to PE pe(i)
          integer :: status                                ! 0 if successful, non zero otherwise
        end function RPN_COMM_spread_context

        function RPN_COMM_spread(context,source,npts,ndata,dest)        &
     &                           result(status)
        use ISO_C_BINDING
        implicit none
          type(C_PTR), intent(IN) :: context
          integer, intent(IN) :: npts, ndata
          real, dimension(npts,ndata), intent(IN) :: source  ! source array, used only on root PE
          real, dimension(:,:), pointer, intent(INOUT) :: dest
          integer :: status
        end function RPN_COMM_spread

        function rpn_comm_shmget(comm,size) result(where)
        use ISO_C_BINDING
        implicit none
          integer(C_INT), intent(IN) :: comm               ! RPN_COMM communicator
          integer(C_INT), intent(IN) :: size               ! size in bytes of shared memory area
          type(C_PTR) :: where                             ! pointer to shared memory area
        end function rpn_comm_shmget
      end interface
