!/* RPNCOMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 2012  Centre ESCER U.Q.A.M
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
      module rpn_comm_localdist
      implicit none
      SAVE
!
!     dimension of clients is (4,1:nclients)
!         clients(:,N) = exchange number, PE number, number of values (for Nth client), data offset
!
!     with this PE as root of transfer, there will be nclients transfers
!     clients(3,N) values will be sent to / received from  PE clients(2,N),
!     data will be read from / written into 
!     starting at buffer position clients(4,N)
!
!     dimension of sources is (3,1:nsources)
!         sources(:,N) = exchange number, PE number, number of values (for Nth source), data offset
!
!     with this PE as client for transfer, there will be nsources transfers
!     sources(3,N) values will be sent to / received from  PE sources(2,N),
!     data will be read from / written into 
!     starting at buffer position sources(4,N)
!
!     a "root" PE has multiple clients to send data to or receive data from
!     an "exchange" is a "root" sending to / receiving from one or more "client(s)"
!           root PE <--->  client PE 1 ,  client PE 2 , ..... , client PE N
!
!     function RPN_COMM_add_dist_entry receives a multi row array where each row
!     describes such an exchange (number of values, root PE, client PE 1, ...., client PE n)
!     a value of -1 for the client PE means a non existent PE (not all exchanges have the same
!     number of "clients" and all rows must have the same dimension)
!     this function adds an entry into the local PE's distribution table and returns the
!     table index of this entry to the caller
!
!     this table index, a data buffer, and a metadata buffer will then be passed to 
!     function RPN_COMM_do_dist that will perform the transfer and fill the metadata 
!     buffer with pertinent information
!
      type :: dist_table_entry                        ! exchange table entry 
        integer, pointer, dimension(:,:) :: clients   ! list of exchange clients for this PE
        integer, pointer, dimension(:,:) :: sources   ! list of exchange "roots" for this PE
        integer, pointer, dimension(:)   :: offsets   ! offset table into data for all rows
        integer :: nclients                           ! number of exchanges where I am the root
        integer :: nsources                           ! number  of exchanges where I am a client
        integer :: ndata                              ! data buffer dimension
        integer :: comm                               ! mpi communicator
      end type dist_table_entry
      type(dist_table_entry), pointer, dimension(:) :: dist_table
      integer :: dist_table_size=-1
      integer :: dist_table_entries=-1

      end module rpn_comm_localdist
!
! ========================================================================================
!
      integer function RPN_COMM_alloc_dist_table(nentries)
!     create the distribution table with nentries entries 
!     (one needed per group of exchange patterns)
!     function returns -1 if failure
!     function returns nentries upon success
      use rpn_comm_localdist
      implicit none
      integer, intent(IN) :: nentries

      RPN_COMM_alloc_dist_table = -1       ! precondition to FAILED
      if(dist_table_size == -1 .and. nentries > 0)then
        allocate(dist_table(1:nentries))
        dist_table_size           = nentries   ! number of possible entries in table
        RPN_COMM_alloc_dist_table = nentries   ! SUCCESS
        dist_table_entries        = 0          ! table is EMPTY
      endif
      return
      end
!
! ========================================================================================
!
      function RPN_COMM_dist_offsets(pattern,offsets,noffsets)
      use rpn_comm_localdist
      implicit none
      integer RPN_COMM_dist_offsets
      integer, intent(IN) :: pattern, noffsets
      integer, dimension(noffsets), intent(OUT) :: offsets

      RPN_COMM_dist_offsets = -1

      if(pattern > dist_table_entries) return  ! out of table
      if(pattern <= 0) return                  ! out of table
      if(noffsets < size(dist_table(pattern)%offsets)) return   ! array offsets is too small

      offsets(1:size(dist_table(pattern)%offsets)) = dist_table(pattern)%offsets
      RPN_COMM_dist_offsets = size(dist_table(pattern)%offsets)
      return
      end
!
! ========================================================================================
!
      integer function RPN_COMM_dist_pattern(pattern)  ! is pattern valid ?
      use rpn_comm_localdist
      implicit none
      integer, intent(IN) :: pattern

      RPN_COMM_dist_pattern = -1  ! precondition to failure
      if(pattern > dist_table_entries) return  ! out of table
      if(pattern <= 0) return                  ! out of table

      RPN_COMM_dist_pattern = pattern  ! pattern is valid
      return
      end
!
! ========================================================================================
!
      integer function RPN_COMM_get_dist_meta(pattern,metadata,nmeta)
!
!     get metadata for an exchange pattern
!
!     metadata : describes what will be sent/received
!         metadata(1,N) is offset of data in array values
!         metadata(2,N) is number of data elements to be collected
!         metadata(3,N) is the exchange number (<0 means this PE is a client)
!         metadata(4,N) is the PE ordinal of the other end of the exchange
!     nmeta has been previously obtained from RPN_COMM_add_dist_entry
!
      use rpn_comm_localdist
      implicit none
      integer, intent(IN) :: pattern, nmeta
      integer, intent(OUT), dimension(4,nmeta) :: metadata

      integer, pointer, dimension(:,:) :: clients, sources
      integer nclients, nsources, i, row, length, offset, request
      integer the_pe

      RPN_COMM_get_dist_meta = -1
      if(pattern > dist_table_entries) return  ! out of table
      if(pattern <= 0) return                  ! out of table

      clients => dist_table(pattern)%clients  ! array describing my clients
      sources => dist_table(pattern)%sources  ! array describing my root data sources
      nclients = dist_table(pattern)%nclients ! number of clients
      nsources = dist_table(pattern)%nsources ! number of root data sources
      nclients = dist_table(pattern)%nclients

      if(nmeta /= nclients + nsources) return  ! inconsistency in metadata size
      request = 0
      do i = 1 , nclients
         request = request + 1
         row     = clients(1,i)
         the_pe  = clients(2,i)
         length  = clients(3,i)
         offset  = clients(4,i)
         metadata(1,request) = offset
         metadata(2,request) = length
         metadata(3,request) = row
         metadata(4,request) = the_pe
      enddo
      do i = 1 , nsources
         request = request + 1
         row     = sources(1,i)
         the_pe  = sources(2,i)
         length  = sources(3,i)
         offset  = sources(4,i)
         metadata(1,request) = offset
         metadata(2,request) = length
         metadata(3,request) = -row
         metadata(4,request) = the_pe
      enddo
      RPN_COMM_get_dist_meta = 0
      return
      end
!
! ========================================================================================
!
      subroutine RPN_COMM_get_dist_dims(pattern,ndata,nmeta)
!
!     get needed sizes for data and metadata buffers when performing
!     exchange with table index = pattern performed by function RPN_COMM_do_dist
!
      use rpn_comm_localdist
      implicit none
      integer, intent(IN) :: pattern
      integer, intent(OUT) :: ndata, nmeta

      integer nclients, nsources, i

      ndata = -1
      nmeta = -1
      if(pattern > dist_table_entries) return  ! out of table
      if(pattern <= 0) return                  ! out of table

      ndata    = 0
      nclients = dist_table(pattern)%nclients
      nsources = dist_table(pattern)%nsources
      do i = 1, nclients        ! where I am root
         ndata = ndata + dist_table(pattern)%clients(3,i)  ! add nb of values in this exchange
      enddo
      do i = 1, nsources        ! where I am a client
         ndata = ndata + dist_table(pattern)%sources(3,i)  ! add nb of values in this exchange
      enddo
      nmeta = nclients + nsources
      end
!
! ========================================================================================
!
      integer function RPN_COMM_add_dist_entry(table,max_clients,nrows,communicator,ndata,nmeta)
!
!     build a table entry describing a group of data exchanges described in array table
!     and add it to the exchange table
!
!     table(-2,j) = empty, will receive offset for this row in data buffer (OUTPUT from this function)
!     table(-1,j) = number of values for this exchange (must be > 0)
!     table( 0,j) = ordinal of root in communicator for this exchange (must be >= 0)
!     table( i,j) = ordinal of a PE in communicator for this exchange (-1 means entry is void)
!     table(1:max_clients,j) = ordinals of all client PEs for this exchange
!     information will be sent from root PE to client PEs or gathered on root PE from client PEs
!     max_client: maximum number of clients in exchanges
!     nrows : number of exchanges described in table
!     communicator : RPN COMM communicator for this exchange ('GRID','SUPERGRID','ALLGRIDS',....)
!     ndata : size of the data buffer needed for the call to RPN_COMM_do_dist  (OUTPUT)
!     nmeta : second dimension of metadata table needed for the call to RPN_COMM_get_dist_meta  (OUTPUT)
!
!     the function will return a pattern index to be used for the actual exchange 
!     or -1 in case of failure
!
!     this function builds the exchange tables for the current PE from the
!     global table
!
      use rpn_comm_localdist

      implicit none
      integer, intent(IN) :: max_clients,nrows
      integer, dimension(-2:max_clients,nrows), intent(INOUT) :: table
      character *(*), intent(IN) :: communicator
      integer, intent(OUT) :: ndata, nmeta

      integer RPN_COMM_comm, RPN_COMM_alloc_dist_table
      external RPN_COMM_comm, RPN_COMM_alloc_dist_table
      integer the_comm, this_pe
      integer my_clients, my_sources, i, j, ierr
      integer, pointer, dimension(:,:) :: clients, sources
      integer no_client, no_source, offset

      table(-2,:) = 999999999
      RPN_COMM_add_dist_entry = -1    ! precondition to FAILED
      if(dist_table_size == -1) then  ! no table allocated, allocate 16 entries by default
         ierr = RPN_COMM_alloc_dist_table(16)
      endif
      if(dist_table_entries >= dist_table_size) return    ! OUCH, table is full

      the_comm = RPN_COMM_comm(communicator)  ! translate communicator from RPN_COMM to MPI
      call RPN_COMM_rank( communicator, this_pe ,ierr )     ! my rank in this communicator
      my_clients = 0
      my_sources = 0
      ndata = 0
      dist_table_entries = dist_table_entries + 1         ! new table entry
      do j =  1 , nrows        ! find number of clients and sources fort this set of exchanges
         if(table(0,j) == this_pe) then  ! I am the root for this exchange
            ndata = ndata + table(-1,j)  ! count needed data buffer space
            do i = 1, max_clients                         ! count valid clients in row
               if(table(i,j) >= 0 .and. table(i,j) /= this_pe) then   ! I have a client other than myself
                  my_clients = my_clients + 1
               endif
            enddo
         else                             ! I am not the root, am i a client ?
            do i = 1, max_clients
               if(table(i,j) == this_pe) then ! I am a client
                  ndata = ndata + table(-1,j) ! count needed data buffer space
                  my_sources = my_sources + 1
                  exit       ! I can only be client once in an exchange
               endif
            enddo
         endif
      enddo

      nmeta = my_clients + my_sources
      allocate(dist_table(dist_table_entries)%clients(4,1:my_clients))    ! allocate client table
      clients => dist_table(dist_table_entries)%clients
      allocate(dist_table(dist_table_entries)%sources(4,1:my_sources))    ! allocate sources table
      sources => dist_table(dist_table_entries)%sources
      allocate(dist_table(dist_table_entries)%offsets(nrows)) ! table of offsets
!      allocate(dist_table(dist_table_entries)%data(ndata))    ! allocate data buffer
      dist_table(dist_table_entries)%clients = -999
      dist_table(dist_table_entries)%nclients = my_clients    ! number of clients
      dist_table(dist_table_entries)%sources = -999
      dist_table(dist_table_entries)%nsources = my_sources    ! number of root sources
      dist_table(dist_table_entries)%offsets = 0              ! 0 means offset is not valid
!      dist_table(dist_table_entries)%data = -1
      dist_table(dist_table_entries)%ndata = ndata
      dist_table(dist_table_entries)%comm = the_comm          ! grid communicator
      no_client = 0
      no_source = 0
!
!     rows in clients will contain the exchange number, client PE ordinal, and data length
!     rows in sources will contain the exchange number,  root  PE ordinal, and data length
!
      offset = 1
      do j =  1 , nrows
        if(table(0,j) == this_pe) then  ! I am the root for this exchange
          do i = 1, max_clients                                  ! count valid clients in row
            if(table(i,j) >= 0 .and. table(i,j) /= this_pe) then         ! I have a client other than myself
               no_client = no_client + 1
               clients(1,no_client) = j           ! exchange number
               clients(2,no_client) = table(i,j)  ! PE ordinal
               clients(3,no_client) = table(-1,j) ! data length
               clients(4,no_client) = offset      ! offset
            endif
          enddo
          dist_table(dist_table_entries)%offsets(j) = offset  ! positive offset means a root
          offset = offset + table(-1,j)      ! bump offset by data length
        else                                        ! I am not the root, am i a client ?
           do i = 1, max_clients
             if(table(i,j) == this_pe) then ! I am a client
               no_source = no_source + 1
               sources(1,no_source) = j   ! row number
               sources(2,no_source) = table(0,j)  ! PE ordinal
               sources(3,no_source) = table(-1,j) ! data length
               sources(4,no_source) = -offset      ! offset
               dist_table(dist_table_entries)%offsets(j) = -offset  ! negative offset means a client
               offset = offset + table(-1,j)      ! bump offset by data length
               exit        ! I can only be client once in an exchange
             endif
           enddo
         endif
      enddo
      table(-2,:) = dist_table(dist_table_entries)%offsets
      RPN_COMM_add_dist_entry = dist_table_entries   !  SUCCESS
      return
      end
!
! ========================================================================================
!
      function RPN_COMM_do_dist(pattern,from_root,values,nvalues)
!
!     pattern is the pattern number returned by RPN_COMM_add_dist_entry for this data exchange
!     values buffer must have been filled according to what was obtained when calling 
!     RPN_COMM_add_dist_entry (data offsets are found in table(-2,:) )
!     (root data if From_root, client data if >not from_root)
!
!     values : 1D array containing data to send / receive
!     nvalues : dimension of array values, max number of values to be distributed/collected
!     from_root : distribution direction flag
!          .true.  : data flows from root PE to client PEs
!          .false. : data flows from client PEs to root PE
!
      use rpn_comm_localdist
      implicit none
      integer :: RPN_COMM_do_dist
      integer, intent(IN) :: pattern, nvalues
      logical, intent(IN) :: from_root
      integer, dimension(nvalues), intent(INOUT) :: values

      include 'mpif.h'

      integer row, j, request, offset, length
      integer, pointer, dimension(:,:) :: clients, sources
      integer, allocatable, dimension(:,:) :: statuses
      integer, allocatable, dimension(:) :: requests
      integer nclients, nsources, the_comm, ierr, the_pe

      RPN_COMM_do_dist = -1   ! precondition to FAILED
      if(pattern > dist_table_entries) return  ! out of table ?
      if(pattern <= 0) return

      clients => dist_table(pattern)%clients  ! array describing my clients
      sources => dist_table(pattern)%sources  ! array describing my root data sources
      nclients = dist_table(pattern)%nclients ! number of clients
      nsources = dist_table(pattern)%nsources ! number of root data sources
      the_comm = dist_table(pattern)%comm     ! MPI communicator
      allocate( statuses(MPI_STATUS_SIZE,nclients+nsources) ) ! asynchronous messages status
      allocate( requests(nclients+nsources) )

      request = 0
      do j = 1 , nclients  ! as root send to / receive from clients
         request = request + 1
         row     = clients(1,j)
         the_pe  = clients(2,j)
         length  = clients(3,j)
         offset  = clients(4,j)
         if(from_root) then ! send  to client, use row as communication tag
            call MPI_isend(values(offset),length,MPI_INTEGER,the_pe,row,the_comm,requests(request),ierr)
         else               ! receive from client, use row as communication tag
            call MPI_irecv(values(offset),length,MPI_INTEGER,the_pe,row,the_comm,requests(request),ierr)
         endif
      enddo
      do j = 1 , nsources  ! as client receive from / send to roots
         request = request + 1
         row     = sources(1,j)
         the_pe  = sources(2,j)
         length  = sources(3,j)
         offset  = abs(sources(4,j))
         if(from_root) then ! receive from root, use row as communication tag
            call MPI_irecv(values(offset),length,MPI_INTEGER,the_pe,row,the_comm,requests(request),ierr)
         else               ! send  to root, use row as communication tag
            call MPI_isend(values(offset),length,MPI_INTEGER,the_pe,row,the_comm,requests(request),ierr)
         endif
      enddo
!
      call MPI_waitall(request,requests,statuses,ierr)    ! wait for all transfers to complete
!
      RPN_COMM_do_dist = 0    ! SUCCESS
!
      deallocate(statuses)
      deallocate(requests)
      return
      end
!
! ========================================================================================
!
      integer function RPN_COMM_dist_test(npetot)
      use rpn_comm_localdist
      implicit none
      integer, intent(IN) :: npetot

      integer, parameter :: MAX_CLIENTS=3

      integer RPN_COMM_dist_offsets, RPN_COMM_do_dist
      external RPN_COMM_dist_offsets, RPN_COMM_do_dist
      integer :: RPN_COMM_add_dist_entry, RPN_COMM_alloc_dist_table
      external :: RPN_COMM_add_dist_entry, RPN_COMM_alloc_dist_table
      integer RPN_COMM_get_dist_meta
      external RPN_COMM_get_dist_meta
      integer, dimension(-2:MAX_CLIENTS,npetot) :: my_table
      integer :: maxpe, i, mype, ierr, pattern, j, low, up
      integer :: ndata, nmeta, nrows
      integer, dimension(20) :: localdata
      integer, pointer, dimension(:,:) :: metadata
      integer, pointer, dimension(:)   :: values
      integer, pointer, dimension(:)   :: offsets
      logical debug

      debug = .false.
      maxpe = npetot -1
      nrows = npetot
      my_table = -1   ! void all entries
      RPN_COMM_dist_test = -1
      call RPN_COMM_rank( 'GRID', mype ,ierr )
!
!     step 0, allocate a communication table (optional if size <16)
!
      ierr = RPN_COMM_alloc_dist_table(10)
!
!     step 1, fill pattern table and create communication pattern
!             then allocate data buffer
!
      my_table = -1
      do i = 1, npetot
        my_table(-1,i) = 1 + mod(i,3)
        my_table(-1,i) = 2
        my_table( 0,i) = maxpe / 2   ! middle PE will be the root
        if(mod(i,maxpe+1) /= my_table( 0,i)) my_table( 1,i)=mod(i,maxpe+1)
        if(mod(i+1,maxpe+1) /= my_table( 0,i)) my_table( 2,i)=mod(i+1,maxpe+1)
      enddo
      if(mype == 0) then
        print *,'ROOT PE will be PE no', maxpe / 2 
      endif
      if(mype == 0 .and. debug) then
        do i = 1, npetot
          print 100,i,'/=/',my_table(-1:MAX_CLIENTS,i)
        enddo
      endif
      pattern = RPN_COMM_add_dist_entry(my_table,MAX_CLIENTS,nrows,'GRID',ndata,nmeta)
100   format(I3,A,20I4)
      if(mype == maxpe / 2) then
        print 100,mype,': CLIENTS',dist_table(pattern)%clients(2,:)  ! clients
        print 100,mype,': EXCHNG ',dist_table(pattern)%clients(1,:)
        print 100,mype,': LENGTH ',dist_table(pattern)%clients(3,:)
        print 100,mype,': OFFSET ',dist_table(pattern)%clients(4,:)
      else
        print 100,mype,': ROOTPE ',my_table(0,:)  ! clients
        print 100,mype,': OFFSETS',my_table(-2,:)  ! offsets table
      endif
       if(debug) print *,mype,' : pattern,ndata,nmeta= ',pattern,ndata,nmeta
102   format(A,I3,A,20I4)
      if(dist_table(1)%nclients > 0 .and. debug) then
          print 102,'++',mype,'//',dist_table(1)%clients(1,:)
          print 102,'++',mype,'//',dist_table(1)%clients(2,:)
          print 102,'++',mype,'//',dist_table(1)%clients(3,:)
      endif
      if(dist_table(1)%nsources > 0 .and. debug) then
          print 102,'//',mype,'//',dist_table(1)%sources(1,:)
          print 102,'//',mype,'//',dist_table(1)%sources(2,:)
          print 102,'//',mype,'//',dist_table(1)%sources(3,:)
      endif
      allocate(metadata(4,nmeta))
      allocate(values(ndata))
      values = -1
!
!     step 3, put data into the proper places
!
      do i = 1, npetot
         if(my_table(-2,i) /= 0) then     ! this row is active
           if(my_table(0,i) == mype) then  ! I am the root for this exchange
             low = abs(my_table(-2,i))
             up = low + my_table(-1,i) - 1
             values(low:up) = mype + 100*i
             print 101,'<=',mype,' =>',low,up,values(low:up)
101          format(A,I2,A,20I4)
           endif
         endif
      enddo
!
!     step 4, perform the exchange
!
      ierr = RPN_COMM_do_dist(pattern,.true.,values,ndata)
      do i = 1, npetot
         if(my_table(-2,i) /= 0) then     ! this row is active
           if(my_table(0,i) /= mype) then  ! I am a client for this exchange
             low = abs(my_table(-2,i))
             up = low + my_table(-1,i) - 1
             print 101,'==',mype,' ==',low,up,values(low:up),values(low:up) - ((maxpe / 2)+100*i)
           endif
         endif
      enddo
      RPN_COMM_dist_test = 0
      return
      end
