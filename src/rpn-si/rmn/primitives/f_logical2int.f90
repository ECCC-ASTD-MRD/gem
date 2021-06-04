!**s/r f_logical2int
subroutine f_logical2int(dest,src,n)
  implicit none
!arguments
  integer :: n
  integer, dimension(n) :: dest
  logical, dimension(n) :: src
!object
!  handler function call in C to store/retreive fortran logical
!** end of rpn-doc sections
  dest = 0
  where(src) dest = 1
  return
end subroutine f_logical2int


!**s/r f_int2logical
subroutine f_int2logical(dest,src,n)
  implicit none
!arguments
  integer :: n
  integer, dimension(n) :: src
  logical, dimension(n) :: dest
!object
!  handler function call in C to store/retreive fortran logical
!** end of rpn-doc sections
  dest = src.ne.0
  return
end subroutine f_int2logical


!**s/r f_logical_move
subroutine f_logical_move(dest,src,n)
  implicit none
!arguments
  integer :: n
  logical, dimension(n) :: dest,src
!object
!  handler function call in C to store/retreive fortran logical
!** end of rpn-doc sections
  dest=src
  return
end subroutine