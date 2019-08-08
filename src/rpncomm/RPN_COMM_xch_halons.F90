!COMP_ARCH=xlf13 -suppress='-O[ ]*2' -add='-O 0'
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

      SUBROUTINE RPN_COMM_xch_halons(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody)
      use rpn_comm
      implicit none
!
!	exchange a halo with neighbours
!

      integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
      logical periodx,periody
!	integer *8 mem_time, exch_time, ewtime
      integer g(minx:maxx,miny:maxy,nk)
!
!	include 'mpif.h'
!
      integer halo_to_east(halox,nj,nk),halo_to_west(halox,nj,nk)
      integer halo_from_east(halox,nj,nk),halo_from_west(halox,nj,nk)
      integer halo_to_north(1-halox:ni+halox,haloy,nk)
      integer halo_to_south(1-halox:ni+halox,haloy,nk)
      integer halo_from_north(1-halox:ni+halox,haloy,nk)
      integer halo_from_south(1-halox:ni+halox,haloy,nk)
!
      integer i, j, k, m
      integer nwds
      integer sendtag, gettag, ierr
      integer status(MPI_STATUS_SIZE)
      logical east,west,north,south
      integer eastpe,westpe,northpe,southpe

      integer tag_2s, tag_2n   ! tags for north_to_south and south_to_north moves
        integer, dimension(4) :: requests            ! table of requests
        integer, dimension(MPI_STATUS_SIZE,4) :: statuses  ! table of statuses
        integer messages ! number of pending asynchronous messages

! print *,'NS halo entered',pe_me
      east=(bnd_east) .and. (.not.periodx)
      eastpe=pe_id(pe_mex+1,pe_mey)
      west=(bnd_west) .and. (.not.periodx)
      westpe=pe_id(pe_mex-1,pe_mey)
      north=(bnd_north) .and. (.not.periody)
      northpe=pe_id(pe_mex,pe_mey+1)
      south=(bnd_south) .and. (.not.periody)
      southpe=pe_id(pe_mex,pe_mey-1)

        if ( pe_ny.eq.1 ) then 
           if(periody)then
              do k=1,nk
            do m=1,haloy
            do i=1-halox,ni+halox
                 g(i,nj+m,k)=g(i,m,k)
                 g(i,m-haloy,k)=g(i,nj+m-haloy,k)
              enddo
              enddo
            enddo
           endif
           goto 9999   ! return
        endif

      do k=1,nk
      do m=1,haloy
      do i=1-halox,ni+halox
        halo_to_south(i,m,k  )=g(i,m,k  )
 	  halo_to_north(i,m,k  )=g(i,nj+m-haloy,k  )
      enddo
      enddo
      enddo
! print *,'NS halo peeled',pe_me
      nwds=nk*haloy*(2*halox+ni)
      sendtag=pe_medomm
      gettag=northpe
      if(north) then
!	  send to south_neighbor
        if(.not.south)then
          call MPI_SEND(halo_to_south,nwds,MPI_INTEGER,southpe,sendtag,PE_DEFCOMM,ierr)
        endif
      else if(south) then
!	  receive from north_neighbor
        call MPI_RECV(halo_from_north,nwds,MPI_INTEGER,northpe,gettag,PE_DEFCOMM,status,ierr)
      else
!	  send to south_neighbor and receive from north_neighbor
        call MPI_SENDRECV(                                      &
              halo_to_south,nwds,MPI_INTEGER,southpe,sendtag,   &
              halo_from_north,nwds,MPI_INTEGER,northpe,gettag,  &
              PE_DEFCOMM,status,ierr)
      endif
      sendtag=pe_medomm
      gettag=southpe
      if(south) then
!	send to north neighbor   
        if(.not.north)then
          call MPI_SEND(halo_to_north,nwds,MPI_INTEGER,northpe,sendtag,PE_DEFCOMM,ierr)
        endif
      else if(north) then
!	  receive from south_neighbor
        call MPI_RECV(halo_from_south,nwds,MPI_INTEGER,southpe,gettag,PE_DEFCOMM,status,ierr)
      else
!	  send to north_neighbor and receive from south_neighbor
        call MPI_SENDRECV(                                       &
               halo_to_north,nwds,MPI_INTEGER,northpe,sendtag,   &
      	       halo_from_south,nwds,MPI_INTEGER,southpe,gettag,  &
      	       PE_DEFCOMM,status,ierr)
      endif        
! print *,'NS halo exchanged',pe_me
!
 	if(.not.north)then
 	do k=1,nk
 	do m=1,haloy
 	do i=1-halox,ni+halox
 	  g(i,nj+m,k  )=halo_from_north(i,m,k  )
 	enddo
 	enddo
 	enddo
 	endif
! print *,'NS north halo inserted',pe_me
      

      if(.not.south)then
      do k=1,nk
      do m=1,haloy
      do i=1-halox,ni+halox
        g(i,m-haloy,k  )=halo_from_south(i,m,k  )
      enddo
      enddo
      enddo
      endif
! print *,'NS south halo inserted',pe_me
9999    continue
      return
        end
