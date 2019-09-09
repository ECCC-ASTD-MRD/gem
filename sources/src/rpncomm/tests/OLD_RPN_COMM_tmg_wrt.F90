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
!InTf!
      SUBROUTINE RPN_COMM_tmg_wrt(Iun)           !InTf! 
      use rpn_comm
      implicit none                              !InTf!
      integer,intent(IN) :: Iun                  !InTf!
!arguments
!  I	Iun	file to write timings into
!*
!        include 'rpn_comm.h'
!        include 'mpif.h'
!
      external r8ivecc, r8ifpec, r8irtc
      real *8 r8ivecc, r8ifpec
      complex *16 r8irtc
      integer i,j,ip,nwds,ipe,ierr,tag,me,totpe
      integer status(MPI_STATUS_SIZE)
!
!
!	temporary patch to allow programs using MSG to
!	collect statistics, information is avalilable
!	in pe_me and pe_tot
!
      call MPI_COMM_RANK(pe_defcomm,me,ierr)
      call MPI_COMM_SIZE(pe_defcomm,totpe,ierr)
!
      nwds = MAXTMGELEM*MAXTMG
!
!	PE 0 will write its own timings,
!	collect timings from other processors
!	and write the collected timings.
!	All messages are sent at FULL size, but tmg_tbl(1,1)
!	contains the useful number of timings collected by a PE
!
      if(me .eq. 0) then  ! I am the ROOT (PE = 0)
        write(Iun)totpe,MAXTMGELEM,MAXTMG
        if(tmg_indx.gt.1) write(Iun)0,MAXTMGELEM,tmg_indx-1, &
                   ((tmg_tbl(i,j),i=1,MAXTMGELEM),j=2,tmg_indx)
        do ip = 1,totpe-1  ! loop over other PEs
          ipe = ip
          tag = ipe
          call MPI_RECV(tmg_tbl,nwds,MPI_DOUBLE_PRECISION,ipe, &
               tag,pe_defcomm,status,ierr)
          tmg_indx = tmg_tbl(1,1)
!	    DO NOT WRITE empty timing records
          if(tmg_indx.gt.1) write(Iun)ipe,MAXTMGELEM,tmg_indx-1, &
                      ((tmg_tbl(i,j),i=1,MAXTMGELEM),j=2,tmg_indx)
        enddo
      else  ! I am an ordinary PE, send timings to ROOT (PE = 0)
        ipe = 0
        tag = me
        tmg_tbl(1,1) = tmg_indx
        call MPI_SSEND(tmg_tbl,nwds,MPI_DOUBLE_PRECISION,ipe, &
                tag,pe_defcomm,ierr)
      endif
      return
!
!**ENTRY RPN_COMM_tmg
      entry RPN_COMM_tmg
!notes
!	main timing point, timings are collected and counter
!	is incremented
!*

1	tmg_indx=max(1,min(MAXTMG,tmg_indx+1))
      tmg_new(TMG_CLOCK) = MPI_wtime()
      tmg_new(TMG_CPU)   = dimag(r8irtc())
      tmg_new(TMG_VEC)   = r8ivecc()
      tmg_new(TMG_FLOP ) = r8ifpec()
      do i=1,MAXTMGELEM
        tmg_tbl(i,tmg_indx) = tmg_new(i) - tmg_old(i)
        tmg_old(i) = tmg_new(i)
        tmg_new(i) = 0
      enddo

      return
!**ENTRY RPN_COMM_tmg_in
      entry RPN_COMM_tmg_in
!
!notes
!	special entry called when a communication routine is entered
!	it it used to collect overhead statistics
!*

      tmg_old(TMG_CPU_C) = dimag(r8irtc())
      tmg_old(TMG_CLOCK_C) = real(r8irtc())

      return
!**ENTRY RPN_COMM_tmg_in
      entry RPN_COMM_tmg_out
!
!notes
!	special entry called when a communication routine is exited
!	it it used to collect overhead statistics
!*

      tmg_new(TMG_CPU_C) = tmg_new(TMG_CPU_C) + dimag(r8irtc())-tmg_old(TMG_CPU_C)
      tmg_new(TMG_CLOCK_C) = tmg_new(TMG_CLOCK_C) + real(r8irtc())-tmg_old(TMG_CLOCK_C)
      tmg_old(TMG_CPU_C) = 0
      tmg_old(TMG_CLOCK_C) = 0

      return

      end SUBROUTINE RPN_COMM_tmg_wrt         !InTf!
