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
module numa
      use ISO_C_Binding
      use, intrinsic :: iso_fortran_env
      use clib_itf_mod
      use ptopo
      use lun
      implicit none
      public
      save
      include 'mpif.h'

      logical :: Numa_uniform_L
      integer :: Numa_sockcomm, nodecomm, Numa_peercomm, noderank, &
                 Numa_sockrank, Numa_peerrank
      integer :: Numa_cores_per_socket, Numa_active_cores_per_socket

contains

      subroutine numa_init
      implicit none

      character(len=64) fmt,ni
      integer ierr, isiz, ns, col, row
      integer, dimension(2,0:Ptopo_npey-1,0:Ptopo_npex-1) :: NuRNuP,GNU
!
!     ---------------------------------------------------------------
!      
      call RPN_COMM_split_by_socket ( Ptopo_intracomm, nodecomm     ,&
              Numa_sockcomm, Numa_peercomm, noderank, Numa_sockrank,&
              Numa_peerrank, isiz, ierr )
      call MPI_Comm_size(Numa_sockcomm, ns, ierr)
      call RPN_COMM_allreduce(ns,Numa_cores_per_socket,1,&
                      "MPI_INTEGER","MPI_MAX","grid",ierr)
      Numa_active_cores_per_socket= ns

!!! ns = Numa_cores_per_socket is not working at the moment
!!! so we will use the hard coded value 20 that we eventually
!!! will replaced by an env variable for external control
      
      ns = 20 ! on most of our current systems
      if (Ptopo_npey > ns) then
         Numa_uniform_L = mod(Ptopo_npey,ns)==0
      else
         Numa_uniform_L = mod(ns,Ptopo_npey)==0
      endif

      NuRNuP= 0.
      NuRNuP(1,Ptopo_myrow,Ptopo_mycol) = Numa_sockrank
      NuRNuP(2,Ptopo_myrow,Ptopo_mycol) = Numa_peerrank
      call RPN_COMM_allreduce(NuRNuP,GNU,2*Ptopo_numproc,&
                      "MPI_INTEGER","MPI_SUM","grid",ierr)

      if (Lun_debug_L) then
         write(lun_out,'(4(a25,i5/))') "Numa_sockrank:",Numa_sockrank,&
              "Numa_cores_per_socket:",Numa_cores_per_socket,&
              "Numa_peerrank:",Numa_peerrank,&
              "Numa_active_cores_per_socket:",Numa_active_cores_per_socket
         if (lun_out>0) then
         write(lun_out,'(a25,l)') "Numa_uniform_L:",Numa_uniform_L
         write(lun_out,*) Numa_sockrank,',',Numa_peerrank
         write (ni,'(i3)') Ptopo_npex
         fmt='('//trim(ni)//'("("i2,",",i3,")"))'
         do row=Ptopo_npey-1,0,-1
            write (Lun_out,trim(fmt)) (gnu(1,row,col),&
                     gnu(2,row,col),col=0,Ptopo_npex-1)
         end do
         endif
      endif
!
!     ---------------------------------------------------------------
!
      return
      end subroutine numa_init
      
      subroutine numa_space (F_pntr, F_msize, F_err)
      implicit none

      integer            , intent(OUT) :: F_err
      integer(kind=INT64), intent(IN ) :: F_msize
      integer, dimension(:), pointer, intent(INOUT) :: F_pntr

      type(C_PTR) :: baseptr
      integer(KIND=MPI_ADDRESS_KIND) :: wsiz
      integer :: rank,ierr,dispunit,win
!
!     ---------------------------------------------------------------
!
      F_err = -1
      call MPI_comm_rank ( Numa_sockcomm, rank, ierr )
      if (ierr .ne. MPI_SUCCESS) return
      if(rank == 0) then        ! everythng allocated by rank 0
         wsiz = F_msize
      else
         wsiz = 0
      endif
      dispunit = 4              ! words (integers/floats)
      wsiz = wsiz * dispunit    ! size in Bytes
      call MPI_win_allocate_shared (wsiz, dispunit, MPI_INFO_NULL,&
                                 Numa_sockcomm, baseptr, win, ierr)
      if (ierr .ne. MPI_SUCCESS) return
      call MPI_win_shared_query (win, MPI_PROC_NULL, wsiz, dispunit,&
                                 baseptr, F_err)
!call RPN_COMM_win_allocate_shared ( Numa_sockcomm, F_msize, win, &
!                                          baseptr, F_err )
      nullify(F_pntr)
      call c_f_pointer( baseptr, F_pntr, [1] )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine numa_space
      
end module numa
