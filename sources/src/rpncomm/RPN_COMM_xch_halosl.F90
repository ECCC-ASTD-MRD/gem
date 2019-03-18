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

      SUBROUTINE RPN_COMM_xch_halosl(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,gni,npol_row,nimax)
!     nj, haloy, periody, npol_row : NOT USED
      use rpn_comm
      implicit none
!      include 'mpif.h'
!     
!     exchange a halo with neighbours, semi-lagrangian fashion for E-W part of exchange
!     (really an allgather, G needs to be large enough for the entire E-W dimension)
!
      integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
      integer gni,npol_row
      logical periodx,periody
      integer, intent(INOUT) :: g(minx:maxx,miny:maxy,nk)
      integer nimax
      integer tempo(nimax,miny:maxy,nk,0:1)  ! first dimension must be able to accomodate largest tile

      integer i,j,k, proc, nijk
      integer gmin, gmax
      integer count(pe_nx + pe_ny)  ! lazy shortcut for  max(pe_nx,pe_ny)
      integer depl(pe_nx + pe_ny)   ! lazy shortcut for  max(pe_nx,pe_ny)
      

      integer status(MPI_STATUS_SIZE)
      integer sendtag, recvtag, sendpe, recvpe, sendi, recvi
      integer procx, ierr, RPN_COMM_limit, offset, temp


      ierr = RPN_COMM_limit(pe_mex,pe_nx,1,gni,gmin,gmax,count,depl)
      if(pe_nx.eq.1) goto 100

!     We put our own piece of the array into its place in the global arrray
!     and copy it into tempo buffer
      if(pe_mex.gt.0) then  ! this IF not really needed
         offset = depl(pe_mex+1) - 1
         do k=1,nk
            do j=miny,maxy
               do i=1,ni
!!                temp = g(i,j,k)
!!                tempo(i,j,k,0) = temp
!!                g(offset+i,j,k) = temp
                  g(i+pe_mex*nimax,j,k) = g(i,j,k)  !! to be fixed when distribution rules change
               enddo
            enddo
         enddo
      endif

! put own piece into tempo buffer (following nest of loops can be suppressed when previous IF suppressed)
      do k=1,nk
         do j=miny,maxy
            do i=1,ni
               tempo(i,j,k,0) = g(i+pe_mex*nimax,j,k)  !! to be fixed when distribution rules change
            enddo
         enddo
      enddo
! pass the bucket around, all the data will make the rounds in a ring fashion
      sendi = 0
      recvi = 1
      sendpe = pe_id(pe_mex-1,pe_mey)  ! send to west PE (with periodicity)
      recvpe = pe_id(pe_mex+1,pe_mey)  ! get from east PE (with periodicity)
      procx = mod(pe_mex + 1, pe_nx)   ! data in for this round will be from east PE
      nijk=nimax*(maxy-miny+1)*nk      ! largest amount of data to send or receive
      sendtag= pe_me
      recvtag= recvpe

      do proc = 1,pe_nx-1 
!	 call tmg_start(93,'COMM SEMI-LAG')
         call mpi_sendrecv(tempo(1,miny,1,sendi),nijk,MPI_integer,   &
                           sendpe,sendtag,                           &
                           tempo(1,miny,1,recvi),nijk,MPI_integer,   &
                           recvpe,recvtag, PE_DEFCOMM,status,ierr)
!	 call tmg_stop(93)
         offset = depl(procx+1) - 1
         do k=1,nk  ! data coming in was from PE numbered procx
            do j=miny,maxy
               do i=1,count(procx+1)
                  g(i+procx*nimax,j,k)=tempo(i,j,k,recvi)  !! to be fixed when distribution rules change
!!                g(offset+i,j,k) = tempo(i,j,k,recvi)
               enddo
            enddo
         enddo
         procx = mod(procx + 1, pe_nx) ! data in will be from PE one further east (with wrap around)
         sendi = recvi     ! will send what i got, swap send and receive buffers
         recvi = mod(sendi+1,2)
!!         recvi = 1 - sendi
      enddo
! E-W periodicity if needed
 100  continue
      if(periodx) then
         do k=1,nk
            do j = miny,maxy
               do i=1,halox
                  g(-halox+i,j,k) = g(gni-halox+i,j,k) 
                  g(gni+i,j,k) = g(i,j,k)
               enddo
            enddo
         enddo
      endif
      return
      end
