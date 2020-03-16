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
      module RPN_COMM_bloc_mgt ! block distribution table
      type :: block_entry
        integer :: hash
        integer :: pe_medomm, pe_indomm, pe_mex, pe_mey
        integer :: nblocx, nblocy, pe_nx, pe_ny, pe_myrow, pe_mycol
        integer :: pe_bloc, pe_gr_bloc, pe_blocmaster, pe_gr_blocmaster
        integer :: BLOC_SIZEX, BLOC_SIZEY, BLOC_MASTER
        logical :: BLOC_EXIST
        integer :: BLOC_mybloc, BLOC_myblocx, BLOC_myblocy
        integer :: BLOC_me, BLOC_corner
        integer :: BLOC_comm_world, bloc_comm_row, bloc_comm_col
      end type
      integer, parameter :: MAX_ENTRIES=64
      integer, save :: valid_entries=0
      type(block_entry),dimension(MAX_ENTRIES),target,save :: btab
      end module RPN_COMM_bloc_mgt
!InTf!
      integer function RPN_COMM_bloc(nblocx,nblocy) ! switch to nblocx by nblocy block distribution !InTf!
      implicit none                                                               !InTf!
      integer, intent(IN) :: nblocx, nblocy                                       !InTf!
      integer, external :: RPN_COMM_bloc_create, RPN_COMM_bloc_find
      integer :: index
!     old code
!      RPN_COMM_bloc = RPN_COMM_bloc_create(nblocx,nblocy)
!      return
!     new code
      index = RPN_COMM_bloc_find(nblocx,nblocy,.true.) ! does the required block distribution exist ?
      if(index > 0)then
        RPN_COMM_bloc = 0   ! success
      else
        RPN_COMM_bloc = RPN_COMM_bloc_create(nblocx,nblocy) ! no preexisting answer exists, try to create one
      endif
      return
      end function RPN_COMM_bloc                                                   !InTf!
!InTf!
      integer function RPN_COMM_bloc_find(nblocx,nblocy,set) ! find and optionally use an already defined block distribution !InTf!
      use rpn_comm
      use RPN_COMM_bloc_mgt
      use rpncomm_com
      implicit none                                                                !InTf!
      integer, intent(IN) :: nblocx, nblocy                                        !InTf!
      logical, intent(IN) :: set   ! if block distribution found, apply it         !InTf!
      type(block_entry), pointer :: t
      integer :: i, hash

      RPN_COMM_bloc_find = -1  !  precondition for failure
      if(valid_entries < 1) return  ! no entries in table
      hash=nblocx+nblocy+pe_nx+pe_ny+pe_mex+pe_mey+pe_medomm+pe_indomm+pe_myrow+pe_mycol
      do i=1,valid_entries  ! all ingredients used to compute a block distribution must match
        t=>btab(i)
        if(t%hash      /= hash)      cycle ! quick inexpensive shortcut test
        if(t%nblocx    /= nblocx)    cycle
        if(t%nblocy    /= nblocy)    cycle
        if(t%pe_nx     /= pe_nx)     cycle
        if(t%pe_ny     /= pe_ny)     cycle
        if(t%pe_mex    /= pe_mex)    cycle
        if(t%pe_mey    /= pe_mey)    cycle
        if(t%pe_medomm /= pe_medomm) cycle
        if(t%pe_indomm /= pe_indomm) cycle
        if(t%pe_myrow  /= pe_myrow)  cycle
        if(t%pe_mycol  /= pe_mycol)  cycle
        RPN_COMM_bloc_find = i
        if(set) then  ! set value of variables in module rpn_comm from block distribution table entry 
          if(.not. associated(com_tab)) call init_com_tab
          BLOC_master     = t%BLOC_master
          BLOC_exist      = t%BLOC_exist
          BLOC_SIZEX      = t%BLOC_SIZEX
          BLOC_SIZEY      = t%BLOC_SIZEY
          BLOC_myblocx    = t%BLOC_myblocx
          BLOC_myblocy    = t%BLOC_myblocy
          BLOC_mybloc     = t%BLOC_mybloc
          BLOC_me         = t%BLOC_me
          BLOC_corner     = t%BLOC_corner
          BLOC_comm_world = t%BLOC_comm_world
          BLOC_comm_row   = t%BLOC_comm_row
          BLOC_comm_col   = t%BLOC_comm_col
          pe_bloc         = t%pe_bloc
          com_tab(13)%number = pe_bloc
!    write(rpn_u,*) 'INFO: found pe_bloc=',pe_bloc
          pe_gr_bloc      = t%pe_gr_bloc
          pe_blocmaster   = t%pe_blocmaster
          com_tab(12)%number = pe_blocmaster
          pe_gr_blocmaster = t%pe_gr_blocmaster
!          write(rpn_u,*) 'INFO: using valid block distribution for',nblocx,' by',nblocy
!        else
!          write(rpn_u,*) 'INFO: found valid block distribution for',nblocx,' by',nblocy
        endif
        exit       ! all ingredients match, exit loop
      enddo
      return
      end function RPN_COMM_bloc_find                                               !InTf!
!InTf!
      integer function RPN_COMM_bloc_create(nblocx,nblocy)                          !InTf!
      use rpn_comm
      use RPN_COMM_bloc_mgt
      use rpncomm_com
      implicit none                                                                 !InTf!
      integer, intent(IN) :: nblocx, nblocy                                         !InTf!
!arguments
!     I nblocx, nblocy: number of blocks on the subgrid in x-y direction
!     O RPN_COMM_bloc_create : error status (-1 if error, else 0)
!
      integer nblocs, ierr, indices(nblocx*nblocy)
      integer longx, longy, i,j,n
      integer mybloc
      type(block_entry), pointer :: t
!     
      RPN_COMM_bloc_create = -1
      nblocs=nblocx*nblocy

      if(.not. associated(com_tab)) call init_com_tab
!
!      if(nblocs.eq.1.and.(BLOC_SIZEX*BLOC_SIZEY.eq.1)) then
!         RPN_COMM_bloc_create = 0
!         return
!      endif
!
!     test if dimensions are suitable and table is not full
!
      if(mod(pe_nx,nblocx).ne.0) then
         write(rpn_u,*) 'ERROR: (RPN_COMM_bloc_create) mod(pe_nx,blocx).ne.0'
         return
      endif
      if(mod(pe_ny,nblocy).ne.0) then
         write(rpn_u,*) 'ERROR: (RPN_COMM_bloc_create) mod(pe_ny,blocy).ne.0'
         return
      endif
      if(valid_entries>=MAX_ENTRIES) then ! OOPS no space left in table
        write(rpn_u,*) 'ERROR: (RPN_COMM_bloc_create) block distribution table full'
        return
      endif

      valid_entries=valid_entries+1
      t=>btab(valid_entries)  ! store ingredients necessary to compute table entry
      t%hash = nblocx+nblocy+pe_nx+pe_ny+pe_mex+pe_mey+pe_medomm+pe_indomm+pe_myrow+pe_mycol
      t%nblocx = nblocx
      t%nblocy = nblocy
      t%pe_nx = pe_nx
      t%pe_ny = pe_ny
      t%pe_mex = pe_mex
      t%pe_mey = pe_mey
      t%pe_medomm = pe_medomm
      t%pe_indomm = pe_indomm
      t%pe_myrow = pe_myrow
      t%pe_mycol = pe_mycol

      longx=pe_nx/nblocx
      longy=pe_ny/nblocy

      BLOC_master = 0
      BLOC_EXIST = .true.
      BLOC_SIZEX   = nblocx
      BLOC_SIZEY   = nblocy
      BLOC_myblocx = pe_mex / longx
      BLOC_myblocy = pe_mey / longy
      BLOC_mybloc  = BLOC_myblocx + nblocx*BLOC_myblocy
      BLOC_me      = pe_mex-BLOC_myblocx*longx + longx*(pe_mey-BLOC_myblocy*longy)
      BLOC_corner  = pe_medomm-(pe_mex-BLOC_myblocx*longx) -pe_nx*(pe_mey-BLOC_myblocy*longy)
      BLOC_comm_world = pe_indomm
      BLOC_comm_row   = pe_myrow
      BLOC_comm_col   = pe_mycol
      
      if(BLOC_corner==pe_medomm) BLOC_master=1
      
      pe_bloc = mpi_comm_null
      call MPI_COMM_SPLIT(BLOC_comm_world,BLOC_mybloc,BLOC_me,pe_bloc,ierr)  ! new communicator for blockpeers
      
      n=1
      do j=1,nblocy
         do i=1,nblocx
!            indices(n) = pe_pe0+(i-1)*longx+pe_nx*longy*(j-1)
            indices(n) = (i-1)*longx+pe_nx*longy*(j-1)
            n=n+1
         enddo
      enddo
      call MPI_Group_incl(pe_gr_indomm, nblocs, indices, pe_gr_blocmaster, ierr)  ! group created for blockmasters
!!      call MPI_Group_rank(pe_gr_bloc,mybloc,ierr) ! mybloc is not used
      
      call MPI_Comm_group(pe_bloc,pe_gr_bloc,ierr)  ! new group for blockpeers
      call MPI_Comm_create(BLOC_comm_world, pe_gr_blocmaster, pe_blocmaster, ierr)  ! communicator created for blockmasters

      t%BLOC_master      = BLOC_master  ! store description of block distribution
      t%BLOC_EXIST       = BLOC_EXIST
      t%BLOC_SIZEX       = BLOC_SIZEX
      t%BLOC_SIZEY       = BLOC_SIZEY
      t%BLOC_myblocx     = BLOC_myblocx
      t%BLOC_myblocy     = BLOC_myblocy
      t%BLOC_mybloc      = BLOC_mybloc
      t%BLOC_me          = BLOC_me
      t%BLOC_corner      = BLOC_corner
      t%BLOC_comm_world  = BLOC_comm_world
      t%BLOC_comm_row    = BLOC_comm_row
      t%BLOC_comm_col    = BLOC_comm_col
      t%pe_bloc          = pe_bloc
      com_tab(13)%number = pe_bloc
      t%pe_gr_bloc       = pe_gr_bloc
      t%pe_blocmaster    = pe_blocmaster
      com_tab(12)%number = pe_blocmaster
      t%pe_gr_blocmaster = pe_gr_blocmaster
!  write(rpn_u,*) 'INFO: created pe_bloc=',pe_bloc
      if(pe_me==0) write(rpn_u,*) 'INFO: (RPN_COMM_bloc_create) created block distribution for',nblocx,' by',nblocy
      RPN_COMM_bloc_create = 0  ! success

        
      end function RPN_COMM_bloc_create                                            !InTf!

      
