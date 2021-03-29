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
                 Numa_sockrank, Numa_peerrank, win
      integer :: Numa_cores_per_socket, Numa_active_cores_per_socket
      type(C_PTR) :: baseptr

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

      if (Ptopo_npey > Numa_cores_per_socket) then
         Numa_uniform_L = mod(Ptopo_npey,Numa_cores_per_socket)==0
      else
         Numa_uniform_L = mod(Numa_cores_per_socket,Ptopo_npey)==0
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

      include 'mpif.h'
      include 'rpn_comm.inc'
      
      integer F_err
      integer(KIND=MPI_ADDRESS_KIND) :: F_msize
      integer, dimension(:), pointer :: F_pntr

      integer ierr
!
!     ---------------------------------------------------------------

      call RPN_COMM_win_allocate_shared ( Numa_sockcomm, F_msize, win, &
                                          baseptr, F_err )
      nullify(F_pntr)
      call c_f_pointer( baseptr, F_pntr, [1] )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine numa_space
      
end module numa
