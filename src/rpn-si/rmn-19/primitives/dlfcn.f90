! File dlfcn.f90   
FUNCTION C_F_STRING_CONVERT(CPTR) RESULT(FPTR)
   USE ISO_C_BINDING
   IMPLICIT NONE
  ! Convert a null-terminated C string into a Fortran character array pointer
  TYPE(C_PTR), INTENT(IN) :: CPTR ! The C address
  CHARACTER(KIND=C_CHAR), DIMENSION(:), POINTER :: FPTR
  CHARACTER(C_CHAR), DIMENSION(1), SAVE, TARGET :: dummy_string="?"

  INTERFACE ! strlen is a standard C function from <string.h>
      ! int strlen(char *string)
      FUNCTION strlen(string) RESULT(len) BIND(C,NAME="strlen")
        USE ISO_C_BINDING
        INTEGER(C_INT) :: len
        TYPE(C_PTR), VALUE :: string ! A C pointer
      END FUNCTION
  END INTERFACE

  IF(C_ASSOCIATED(CPTR)) THEN
      CALL C_F_POINTER(FPTR=FPTR, CPTR=CPTR, SHAPE=[strlen(CPTR)])
  ELSE
      ! To avoid segfaults, associate FPTR with a dummy target:
      FPTR=>dummy_string
  END IF

END FUNCTION C_F_STRING_CONVERT

SUBROUTINE DLFCN_Test_Shared_Lib
   USE ISO_C_BINDING
   IMPLICIT NONE
   include 'dlfcn.inc'

   ! Local variables:
   CHARACTER(KIND=C_CHAR,LEN=1024) :: dll_name, sub_name
   TYPE(C_PTR) :: handle=C_NULL_PTR
   TYPE(C_FUNPTR) :: funptr=C_NULL_FUNPTR
   INTEGER(C_INT) :: status
   integer :: iname
   real*8 :: one

   ! The dynamic subroutine has a simple interface:
   ABSTRACT INTERFACE
      SUBROUTINE MySub(x) BIND(C)
         USE ISO_C_BINDING
         REAL(C_DOUBLE), VALUE :: x
      END SUBROUTINE
   END INTERFACE
   PROCEDURE(MySub), POINTER :: dll_sub ! Dynamically-linked procedure

   call register_dl_routines
!   WRITE(*,*) "Enter the name of the DL and the name of the DL subroutine:"
!   READ(*,"(A)") dll_name ! Enter "shared.so"
!   READ(*,"(A)") sub_name ! Enter "MySub"
   dll_name="./libshared.so"

   ! Open the DL:
   handle=DL_Open(TRIM(dll_name)//C_NULL_CHAR, IOR(RTLD_NOW, RTLD_GLOBAL))
      ! The use of IOR is not really proper...wait till Fortran 2008
   IF(.NOT.C_ASSOCIATED(handle)) THEN
      WRITE(*,*) "Error in dlopen: ", C_F_STRING_CONVERT(DL_Error())
      STOP
   END IF
   one = 1.0_c_double
   do iname=1,3

   if(iname==1) sub_name='MySub1'
   if(iname==2) sub_name='MySub2'
   if(iname==3) sub_name='MySub3'

   ! Find the subroutine in the DL:
   funptr=DL_Sym(handle,TRIM(sub_name)//C_NULL_CHAR)
   IF(.NOT.C_ASSOCIATED(funptr)) THEN
      WRITE(*,*) "Error in dlsym: ", C_F_STRING_CONVERT(DL_Error())
      STOP
   END IF
   ! Now convert the C function pointer to a Fortran procedure pointer
   CALL C_F_PROCPOINTER(CPTR=funptr, FPTR=dll_sub)
   ! Finally, invoke the dynamically-linked subroutine:
   CALL dll_sub(one)
   one = one + .001

   enddo

   ! Now close the DL:
   status=DL_Close(handle)
   IF(status/=0) THEN
      WRITE(*,*) "Error in dlclose: ", C_F_STRING_CONVERT(DL_Error())
      STOP
   END IF

END SUBROUTINE

