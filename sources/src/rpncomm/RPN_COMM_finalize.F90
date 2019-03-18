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
!InTf!
      SUBROUTINE RPN_COMM_finalize(ierr) !InTf!

!	Luc Corbeil, 2000-11-21
!	mpi finalize
      use rpn_comm
      implicit none                      !InTf!
      integer, intent(OUT) ::  ierr      !InTf!
        
!*
!        include 'mpif.h'

!       if(allocated(pe_domains)) deallocate(pe_domains)
      if(allocated(pe_id))      deallocate(pe_id)
      if(allocated(pe_xtab))    deallocate(pe_xtab)
      if(allocated(pe_ytab))    deallocate(pe_ytab)
      if(allocated(ord_tab))    deallocate(ord_tab)


!       call rpn_comm_softbarrier(MPI_COMM_WORLD)
      call mpi_finalize(ierr)

      return
      end SUBROUTINE RPN_COMM_finalize         !InTf!
