subroutine rpn_comm_test_010
use ISO_C_BINDING
implicit none
#include "RPN_COMM_int.inc"

integer :: a1(10), a2(10,10), a3(10,10,10), a4(10,10,10,10)
integer, pointer, dimension(:) :: pp
integer, allocatable, dimension(:) :: b1
integer, allocatable, dimension(:,:) :: b2
integer, allocatable, dimension(:,:,:) :: b3
integer, allocatable, dimension(:,:,:,:) :: b4

type(rpncomm_ptr) :: p1, p2, p3, p4

allocate(b1(10))
allocate(b2(10,10))
allocate(b3(10,10,10))
allocate(b4(10,10,10,10))

call rpn_comm_ptr(a1,p1)
call c_f_pointer(p1%p,pp,[10])
call ploc(pp,a1)
call rpn_comm_ptr(a2,p2)
call c_f_pointer(p2%p,pp,[10])
call ploc(pp,a2)
call rpn_comm_ptr(a3,p3)
call c_f_pointer(p3%p,pp,[10])
call ploc(pp,a3)
call rpn_comm_ptr(a4,p4)
call c_f_pointer(p4%p,pp,[10])
call ploc(pp,a4)

call rpn_comm_ptr(b1,p1)
call c_f_pointer(p1%p,pp,[10])
call ploc(pp,b1)
call rpn_comm_ptr(b2,p2)
call c_f_pointer(p2%p,pp,[10])
call ploc(pp,b2)
call rpn_comm_ptr(b3,p3)
call c_f_pointer(p3%p,pp,[10])
call ploc(pp,b3)
call rpn_comm_ptr(b4,p4)
call c_f_pointer(p4%p,pp,[10])
call ploc(pp,b4)


stop
end
subroutine ploc(p,ref)
integer, dimension(*) :: p,ref
print *,loc(p),loc(ref),loc(p)==loc(ref)
return
end

