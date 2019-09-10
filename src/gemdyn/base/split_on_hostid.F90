subroutine split_on_hostid(comm,me,howmany,newcomm)
  use iso_c_binding
  implicit none
!  include 'mpif.h'
  integer, intent(IN) :: comm
  integer, intent(OUT) :: me, newcomm, howmany
  integer :: my_color, me_in_comm, ierr
  interface
    integer(C_INT) function f_gethostid()BIND(C,name='gethostid')
      import :: C_INT
    end function f_gethostid
  end interface
  my_color = abs(f_gethostid())
  call MPI_COMM_RANK(comm,me_in_comm,ierr)
  call MPI_COMM_SPLIT(comm,my_color,me_in_comm,newcomm,ierr)
  call MPI_COMM_RANK(newcomm,me,ierr)
  call MPI_COMM_SIZE(newcomm,howmany,ierr)
  return
end subroutine split_on_hostid
