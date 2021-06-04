!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
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
        module RPN_COMM_barrier_priv
          integer, save :: detail = 0   ! 0, 1, 2
          integer, save :: isize = -1   ! valid number of data points in times array
          integer, save :: irank = -1   ! rank of PE in communicator domain
          real *8, save, dimension(:), pointer :: times  ! array of timings (used if detail =2)
          logical, save :: valid = .false. ! if .true. , times array is allocated
          real *8, save :: t0 = -1.0    ! max time around barrier (used if detail = 1)
        end module RPN_COMM_barrier_priv
!InTf!
        SUBROUTINE RPN_COMM_barrier(com,ierr)                    !InTf!

!	Luc Corbeil, 2000-11-21
!	mpi barrier
        use RPN_COMM_barrier_priv
        implicit none                                            !InTf!
        integer, intent(OUT) :: ierr                             !InTf!
        character(len=*), intent(IN) ::  com                     !InTf!
!*
        integer :: comm
        include 'mpif.h'

        integer, external :: RPN_COMM_comm
        logical RPN_COMM_grank
        real *8 t1,t2

	comm=rpn_comm_comm(com)

        if(.not.RPN_COMM_grank(com)) return
        if(detail > 0) then
          t1=mpi_wtime()
        endif
        call mpi_barrier(comm,ierr)
        if(detail > 0) then
          t0 = -1.0    ! indicate that no valid data is available
          isize = -1
          irank = -1
          t2=mpi_wtime()-t1
          if(detail==1) then  ! just get max barrier wait
            call mpi_reduce(t2,t0,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm,ierr)
          endif
          if(detail==2) then  ! get all barrier waits
            call mpi_comm_size(comm,isize,ierr)
            call mpi_comm_rank(comm,irank,ierr)
            if(irank==0) then
              if(valid)then
                if(size(times)<isize) then
                  deallocate(times)
                  allocate(times(isize))
                endif
              else
                allocate(times(isize))
              endif
              valid = .true.
            endif
            call mpi_gather(t2,1,MPI_DOUBLE_PRECISION,times,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
          endif
        endif
	return
	end SUBROUTINE RPN_COMM_barrier                                !InTf!
!InTf!
        integer function RPN_COMM_barrier_data(level,values,nvalues)                                !InTf!
        use RPN_COMM_barrier_priv
        implicit none
        integer , intent(IN) :: level  ! >= 0, set detail level,  <0 get data from last barrier call !InTf!
        integer , intent(IN) :: nvalues                                                              !InTf!
        real(kind=kind(1.d0)), dimension(nvalues), intent(OUT) :: values                                           !InTf!

        RPN_COMM_barrier_data = -1  ! precondition to error

        if(level>=0) then  ! set detail level, ignore other arguments
          detail = level
          RPN_COMM_barrier_data = 0
          return
        endif

        if(detail==0)return     ! no data available
        if(irank /= 0) return   ! only root of communicator has valid data
        irank = -1              ! invalidate data

        if(nvalues==1.and.t0>0)then  ! return max time value
          values(1)=t0
          RPN_COMM_barrier_data = 1
          t0 = -1.0
        endif

        if(nvalues>1.and.isize>0)then  ! return the whole array of timings
          values(1:min(isize,nvalues))=times(1:min(isize,nvalues))
          RPN_COMM_barrier_data = min(isize,nvalues)
          isize = -1
        endif

        return ! function values is the number of valid data points in values
        end function RPN_COMM_barrier_data                                            !InTf!
