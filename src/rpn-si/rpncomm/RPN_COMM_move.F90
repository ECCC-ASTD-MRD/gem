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

        SUBROUTINE RPN_COMM_move( sendbuf, sendcount, mtyp, dest,&
     &     recvbuf, recvcount,periodicity, ierr)

        use rpn_comm
        implicit none
!        include 'rpn_comm.h'
!        include 'mpif.h'

        integer sendbuf, sendcount, recvbuf, recvcount
        integer sendtag, recvtag, ierr
        character(len=*), intent(IN) ::  mtyp
        character(len=*), intent(IN) ::  dest
        character(len=1) dest2
        integer sendt, recvt, icomm, icomm2, idest, irecv
        logical lsend, lrecv,lsendrecv, borders, borderr
        logical periodx,periody,periodicity
        integer status(MPI_STATUS_SIZE)

        integer  RPN_COMM_datyp, RPN_COMM_comm
        external RPN_COMM_datyp, RPN_COMM_comm

        periodx=periodicity
        periody=periodicity
        ierr = -1
        sendt=rpn_comm_datyp(mtyp)
        recvt=rpn_comm_datyp(mtyp)
        icomm=pe_defcomm
        borders = .false.
        borderr = .false.
        lsend = .false.
        lrecv= .false.
        lsendrecv = .false.

        call rpn_comm_low2up(dest,dest2)

        if(dest2.eq.'N') then
           idest = pe_id(pe_mex,pe_mey+1)
           irecv = pe_id(pe_mex,pe_mey-1)
           borders = (pe_mey.eq.(pe_ny-1)).and..not.periody
           borderr = (pe_mey.eq.0).and..not.periody
        else if (dest2.eq.'S') then
           idest = pe_id(pe_mex,pe_mey-1)
           irecv = pe_id(pe_mex,pe_mey+1)
           borderr = (pe_mey.eq.(pe_ny-1)).and..not.periody
           borders = (pe_mey.eq.0).and..not.periody
        else if (dest2.eq.'E') then
           idest = pe_id(pe_mex+1,pe_mey)
           irecv = pe_id(pe_mex-1,pe_mey)
           borders = (pe_mex.eq.(pe_nx-1)).and..not.periodx
           borderr = (pe_mex.eq.0).and..not.periodx
        else if (dest2.eq.'W') then
           idest = pe_id(pe_mex-1,pe_mey)
           irecv = pe_id(pe_mex+1,pe_mey)
           borderr = (pe_mex.eq.(pe_nx-1)).and..not.periodx
           borders = (pe_mex.eq.0).and..not.periodx
        else
           write(rpn_u,*) 'RPN_COMM_move ERROR: invalid destination'
           write(rpn_u,*) "Valid destinations: 'N','S','E','W' "
           return
        endif

        sendtag = pe_me
        recvtag = irecv
        if((sendcount.gt.0).and.(recvcount.gt.0).and..not.borders&
     &       .and..not.borderr) then
           call MPI_sendrecv( sendbuf, sendcount, sendt,&
     &          idest, sendtag, recvbuf, recvcount, recvt, irecv, &
     &          recvtag, icomm, status, ierr )
           if(ierr.ne.0) then
              write(rpn_u,*) 'RPN_COMM_move ERROR: sendrecv failed'
              return
           endif
        else 
           if(sendcount.gt.0.and..not.borders) then
!              write(rpn_u,*) 'J''envoie!',pe_me,sendbuf,sendcount
              call MPI_send(sendbuf,sendcount,sendt,idest,&
     &             sendtag, icomm, ierr)
              if(ierr.ne.0) then
                 write(rpn_u,*) 'RPN_COMM_move ERROR: send failed'
                 return
              endif
              
           else if(recvcount.gt.0.and..not.borderr) then
!               write(rpn_u,*) 'Je recois!',pe_me
             call RPN_COMM_recv(recvbuf,recvcount,recvt,irecv,&
     &             recvtag, icomm, status, ierr)
              if(ierr.ne.0) then
                 write(rpn_u,*) 'RPN_COMM_move ERROR: send failed'
                 return
              endif
           endif
        endif

        return
        end SUBROUTINE RPN_COMM_move

           
