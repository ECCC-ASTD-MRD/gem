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
module copy_and_open
      use iso_c_binding
      implicit none
      public
      save
      include 'mpif.h'
#include <rmnlib_basics.hf>

      interface
         function c_unlink(path) result(ok) BIND(C,name='unlink') ! remove file
         import :: C_CHAR, C_INT
         character(C_CHAR), dimension(*), intent(IN) :: path
         integer(C_INT) :: ok
         end function c_unlink
         function mkstemp (template) result(fd) BIND(c,name='mkstemp')
         import :: C_CHAR, C_INT
         character(C_CHAR), dimension (*), intent(IN) :: template
         integer(C_INT) :: fd
         end function mkstemp
      end interface

contains

      subroutine cp_n_open (F_unf, F_comm, F_file_S)
      implicit none
      character(len=*), intent(IN) :: F_file_S
      integer, intent(INOUT) :: F_unf
      integer, intent(IN   ) :: F_comm

      character(len=2048) tmp_name
      integer, parameter :: BUFSIZE=100000
      integer, dimension(BUFSIZE) :: buf
      integer :: rank, err
      logical, parameter :: cpornot = .false.
!
!-------------------------------------------------------------------
!
      if (cpornot) then
         call MPI_comm_rank ( F_comm, rank, err )
         if (rank == 0) then
            call array_from_file(buf,size(buf),F_file_S)
         endif
         call MPI_BCAST (buf, size(buf), MPI_INTEGER, 0, F_comm, err)
         tmp_name= '/tmp/gemtmpf_'//'XXXXXX'//C_NULL_CHAR
         err = mkstemp(tmp_name)
         call array_to_file (buf,size(buf),trim(tmp_name))
         F_unf= 0
         if (fnom (F_unf,trim(tmp_name), 'SEQ+OLD', 0) == 0) then
            err= c_unlink(trim(tmp_name)//achar(0))
         else
            F_unf= -1
         endif
      else
         F_unf= 0
         if (.not. (fnom (F_unf,trim(F_file_S), 'SEQ+OLD', 0) == 0)) F_unf= -1
      endif
         
!
!-------------------------------------------------------------------
!
      return
      end subroutine cp_n_open

end module copy_and_open
