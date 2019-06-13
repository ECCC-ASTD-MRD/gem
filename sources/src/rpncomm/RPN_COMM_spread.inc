!
! spread metadata (private information for RPN_COMM_spread and associates)
!                 (application code gets a bind C pointer to this information)
!
  type, bind(C)    :: spread_context
    type(c_ptr)    :: order      ! reordering table (useful only on root PE)
    type(c_ptr)    :: nlocal     ! number of samples to send to each PE (useful only on root PE)
    type(c_ptr)    :: offset     ! offset table for scatterv (useful only on root PE)
    integer(C_INT) :: ntotal     ! total number of samples (useful only on root PE)
    integer(C_INT) :: npoints    ! number of samples expected on this PE
    integer(C_INT) :: n_pes      ! number pf PEs (useful only on root PE)
    integer(C_INT) :: comm       ! MPI communicator for this "spread" operation
    integer(C_INT) :: myrank     ! rank ot this PE in comm
    integer(C_INT) :: rootpe     ! rank of root pe for "spread" operration
    integer(C_INT) :: shiftcnt   ! shift count for PE / point number extraction
  end type spread_context
