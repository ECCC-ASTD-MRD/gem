!!end interface     !InTf!
!! interface rpn_comm_ptr     !InTf!
#define IN_RPN_COMM_ptr
  subroutine RPN_COMM_ptr_i4_1d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  integer(KIND=4), dimension(1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1))
  return
  end subroutine RPN_COMM_ptr_i4_1d    !InTf!

  subroutine RPN_COMM_ptr_i4_2d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  integer(KIND=4), dimension(1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1))
  return
  end subroutine RPN_COMM_ptr_i4_2d    !InTf!

  subroutine RPN_COMM_ptr_i4_3d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  integer(KIND=4), dimension(1,1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1,1))
  return
  end subroutine RPN_COMM_ptr_i4_3d    !InTf!

  subroutine RPN_COMM_ptr_i4_4d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  integer(KIND=4), dimension(1,1,1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1,1,1))
  return
  end subroutine RPN_COMM_ptr_i4_4d    !InTf!

  subroutine RPN_COMM_ptr_i8_1d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  integer(KIND=8), dimension(1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1))
  return
  end subroutine RPN_COMM_ptr_i8_1d    !InTf!

  subroutine RPN_COMM_ptr_i8_2d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  integer(KIND=8), dimension(1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1))
  return
  end subroutine RPN_COMM_ptr_i8_2d    !InTf!

  subroutine RPN_COMM_ptr_i8_3d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  integer(KIND=8), dimension(1,1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1,1))
  return
  end subroutine RPN_COMM_ptr_i8_3d    !InTf!

  subroutine RPN_COMM_ptr_i8_4d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  integer(KIND=8), dimension(1,1,1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1,1,1))
  return
  end subroutine RPN_COMM_ptr_i8_4d    !InTf!

  subroutine RPN_COMM_ptr_r4_1d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  real(KIND=4), dimension(1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1))
  return
  end subroutine RPN_COMM_ptr_r4_1d    !InTf!

  subroutine RPN_COMM_ptr_r4_2d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  real(KIND=4), dimension(1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1))
  return
  end subroutine RPN_COMM_ptr_r4_2d    !InTf!

  subroutine RPN_COMM_ptr_r4_3d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  real(KIND=4), dimension(1,1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1,1))
  return
  end subroutine RPN_COMM_ptr_r4_3d    !InTf!

  subroutine RPN_COMM_ptr_r4_4d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  real(KIND=4), dimension(1,1,1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1,1,1))
  return
  end subroutine RPN_COMM_ptr_r4_4d    !InTf!

  subroutine RPN_COMM_ptr_r8_1d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  real(KIND=8), dimension(1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1))
  return
  end subroutine RPN_COMM_ptr_r8_1d    !InTf!

  subroutine RPN_COMM_ptr_r8_2d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  real(KIND=8), dimension(1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1))
  return
  end subroutine RPN_COMM_ptr_r8_2d    !InTf!

  subroutine RPN_COMM_ptr_r8_3d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  real(KIND=8), dimension(1,1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1,1))
  return
  end subroutine RPN_COMM_ptr_r8_3d    !InTf!

  subroutine RPN_COMM_ptr_r8_4d(what,ptr)    !InTf!
  use ISO_C_BINDING    !InTf!
!!  import :: rpncomm_ptr      !InTf!
  implicit none    !InTf!
  include 'RPN_COMM_types.inc'
  real(KIND=8), dimension(1,1,1,1),intent(IN), target :: what    !InTf!
  type(rpncomm_ptr), intent(OUT) :: ptr    !InTf!
  ptr%p = c_loc(what(1,1,1,1))
  return
  end subroutine RPN_COMM_ptr_r8_4d    !InTf!

!!end interface rpn_comm_ptr      !InTf!
!!interface     !InTf!
