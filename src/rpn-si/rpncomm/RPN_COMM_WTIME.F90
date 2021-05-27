!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2013 Division de Recherche en Prevision
! Numerique
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
!InTf!
        real(kind=kind(1.d0)) function RPN_COMM_wtime()            !InTf!
        implicit none                               !InTf!
!
        include 'mpif.h'
!
        integer :: myrank, ierr
        logical,save :: display = .TRUE.
!
        if (display) then
           call MPI_COMM_RANK(MPI_COMM_WORLD, myrank,ierr)
           if (myrank == 0) then
              write (6,100) MPI_Wtick()
           endif
           display = .FALSE.
        endif

        RPN_COMM_wtime = MPI_Wtime()

 100    format ('RPN_COMM_wtime timer resolution:',ES9.2,' s')
        return
        end function RPN_COMM_wtime                  !InTf!

