!
! helper routine used by create_ioset
! for now only method 0 is supported for IO PE dispersion
!
  subroutine RPN_COMM_make_io_pe_list(x,y,npes,pe_nx,pe_ny,method)  !InTf!
    implicit none
    integer, intent(IN) :: npes                 !InTf!   ! number of PEs in set
    integer, dimension(npes), intent(OUT) :: x  !InTf!   ! x coordinates of PEs in set
    integer, dimension(npes), intent(OUT) :: y  !InTf!   ! y coordinates of PEs in set
    integer, intent(IN) :: pe_nx, pe_ny         !InTf!   ! number of PEs along x and y in PE grid
    integer, intent(IN) :: method               !InTf!   ! fill method
    integer :: i
    integer :: deltax, deltay, pe_nxy
    integer, save :: scramblex = 0
    integer, save :: scrambley = 0
    integer, save :: pe_nx_old = 0
    integer, save :: pe_ny_old = 0
    logical :: modulo_x_y
    integer, dimension(16), save :: primes = [ &
         5,      7,      11,     13,   &
        17,     19,      23,     29,   &
        31,     37,      41,     43,   &
        47,     53,      59,     61  ]
!
    if(scramblex == 0 .or. pe_nx_old .ne. pe_nx .or. pe_ny_old .ne. pe_ny) then    ! (re)initialize
      scramblex = 1
      scrambley = 1
      pe_nx_old = pe_nx
      pe_ny_old = pe_ny
      if(pe_nx > pe_ny) then   ! scramblex = lowest number that is prime with respect to pe_nx
        do i = 1 , size(primes)
          scramblex = primes(i)
          if(mod(pe_nx,primes(i)) .ne. 0) exit
        enddo
      else                     ! scrambley = lowest number that is prime with respect to pe_ny
        do i = 1 , size(primes)
          scrambley = primes(i)
          if(mod(pe_ny,primes(i)) .ne. 0) exit
        enddo
      endif
#if defined(SELF_TEST)
      print *,"DEBUG: scramblex, scrambley =",scramblex, scrambley
#endif
    endif
    x = -1
    y = -1
    if(method .ne. 0) then
      print *,"ERROR: method MUST be zero for the time being"
      return      ! method 0 is the only one supported for the time being
    endif
    deltax = 1
    deltay = 1
    pe_nxy = min(pe_nx,pe_ny)
    if(method == 0) then
      deltax = scramblex
      deltay = scrambley
    endif
!     if( npes > pe_nxy * pe_nxy ) then
!       print *,"ERROR: too many PEs requested in set (",npes," ),max permitted:",pe_nxy * pe_nxy
!       return
!     endif
!     modulo_x_y = mod(pe_ny,pe_nx) == 0 .or. mod(pe_nx,pe_ny) == 0 ! one is a multiple of the other
    modulo_x_y = .true.   ! forced true for the time being
    do i = 0 , npes-1
!       if(npes > pe_nxy) then
!         if(pe_nx > pe_ny) then
!           x(i+1) = mod( i * deltax , pe_nx )
!           y(i+1) = mod( mod( i , pe_ny) + i / pe_nx , pe_ny)
!         else
!           x(i+1) = mod( mod( i , pe_nx ) + i / pe_ny , pe_nx)
!           y(i+1) = mod( i * deltay , pe_ny)
!         endif
!       else
        if(modulo_x_y) then
          if(pe_nx < pe_ny) then
            x(i+1) = mod( i , pe_nx)
            y(i+1) = mod(x(i+1)+(i/pe_nx)*scrambley,pe_ny)
          else
            y(i+1) = mod( i , pe_ny)
            x(i+1) = mod(y(i+1)+(i/pe_ny)*scramblex,pe_nx)
          endif
!           y(i+1) = mod(y(i+1)+i/pe_ny,pe_ny)
        else
          x(i+1) = mod( i * deltax , pe_nx )
          y(i+1) = mod( i * deltay , pe_ny )
        endif
!       endif
    enddo
  end subroutine RPN_COMM_make_io_pe_list  !InTf!
!
!
! check that this PE set is valid (no duplicates) and print IO PE map
! also check that no column has 2 IO PEs in a group and neither has any row
! used to validate what create_ioset has produced
!
  function RPN_COMM_check_ioset(newset, x ,y, npes, pe_nx, pe_ny, pe_me, diag) result(status)  !InTf!
    implicit none
    integer, intent(IN) :: newset, npes        !InTf!   ! set number, number of PEs in set
    integer, intent(IN) :: pe_nx, pe_ny        !InTf!   ! number of PEs along x and y
    integer, intent(IN) :: pe_me               !InTf!   ! ordinal in grid of this PE
    logical, intent(IN) :: diag                !InTf!   ! if .true., PE 0 prints the IO PE map for this set
    integer, intent(IN), dimension(npes) :: x  !InTf!   ! x coordinates of IO PEs
    integer, intent(IN), dimension(npes) :: y  !InTf!   ! y coordinates of IO PEs
    integer :: status                          !InTf!   ! 0 if set is OK, -1 otherwise
!
    integer*1, dimension(0:pe_nx-1,0:pe_ny-1) :: safe
    integer*1, dimension(0:pe_ny-1) :: row      ! there are pe_ny rows
    integer*1, dimension(0:pe_nx-1) :: col      ! there are pe_nx columns
    integer :: i, j, groupsize, low, high
!
    status = -1
    safe = 0
    groupsize = min(pe_nx, pe_ny)
    do i = 1, npes, groupsize  ! loop over groups
      row = 0
      col = 0
      do j = i , min(i+groupsize-1,npes)      ! for PE at coordinates ( x(j), y(j) )
        row(y(j)) = row(y(j)) + 1             ! add one to count for this row
        col(x(j)) = col(x(j)) + 1             ! add one to count for this column
      enddo
      if(any(row > 1) .or. any(col > 1) ) then
        print *,"ERROR: problem creating IO set, there are 2 or more PEs on row or column"
        print *,"ERROR: x = ",x(i:min(i+groupsize-1,npes))
        print *,"ERROR: row = ",row
        print *,"ERROR: y = ",y(i:min(i+groupsize-1,npes))
        print *,"ERROR: col = ",col
        print *,"ERROR: group, group limits = ",(i-1)/groupsize+1,(j,j=i,min(i+groupsize-1,npes))
        status = -1
        return
      endif
    enddo
    do i = 1 , npes
      if(safe(x(i),y(i)) .ne. 0) then  ! OOPS, 2 PEs in this set are the same
        print *,"ERROR: problem creating IO set, there are duplicates",x(i),y(i),safe(x(i),y(i))
        print *,"ERROR: x = ",x(i:min(i+groupsize-1,npes))
        print *,"ERROR: y = ",y(i:min(i+groupsize-1,npes))
        status = -1
        return
      else
        safe(x(i),y(i)) = 1 + mod(  ( (i-1) / min(pe_nx, pe_ny) ) , 9)  ! group number, 9 if group number > 9
      endif
    enddo
    if(pe_me == 0 .and. diag)then
      print 101,"===== IO PE map for set",newset," (",(npes+min(pe_nx, pe_ny)-1)/min(pe_nx, pe_ny)," groups) ====="
      print 102,"INFO: x =",x(1:npes)
      print 102,"INFO: y =",y(1:npes)
      do j = pe_ny-1 , 0 , -1
        print 100,j,safe(:,j)
      enddo
    endif
100   format(1X,I4,1x,128I1)
101   format(A,I3,A,I3,A)
102   format(A10,20I4(/10X,20I4))
    status = 0
    return
  end function RPN_COMM_check_ioset  !InTf!
#if defined(SELF_TEST)
program test_valid_set
  implicit none
  integer, parameter :: MAXPES=1024
  integer, dimension(MAXPES) :: x    ! x coordinates of PEs in set
  integer, dimension(MAXPES) :: y    ! y coordinates of PEs in set
  integer :: npes                   ! number of PEs in set
  integer :: penx, peny             ! number of PEs along x and y in PE grid
  integer :: method                 ! fill method
  logical :: diag                   ! if .true. print IO PE map and diagnostics
! Notes
!   some error messages are printed in case of error
!   diagnostics include a map of IO PEs
!******

  integer, external :: RPN_COMM_check_ioset

  diag = .true.
  method = 0
  npes = 1
  do while(npes >0)
    print *,'npes,penx,peny ?'
    read(5,*)npes,penx,peny
    if(npes <=0) stop
    call RPN_COMM_make_io_pe_list(x,y,npes,penx,peny,method)
    if(x(1) == -1) stop   ! miserable failure at creation
    if( RPN_COMM_check_ioset(0, x ,y, npes, penx, peny, 0, diag) .ne. 0 ) stop
  enddo
  stop
end program test_valid_set 
#endif