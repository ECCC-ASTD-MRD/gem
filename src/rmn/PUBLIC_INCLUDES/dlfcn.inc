INTEGER, PARAMETER :: RTLD_LAZY=1, RTLD_NOW=2, RTLD_GLOBAL=256, RTLD_LOCAL=0
!
INTERFACE ! All we need is interfaces for the prototypes in <dlfcn.h>
  FUNCTION strlen(string) RESULT(len) BIND(C,NAME="strlen")
      ! int strlen(char *string)
      USE ISO_C_BINDING
      INTEGER(C_INT) :: len
      TYPE(C_PTR), VALUE :: string ! A C pointer
  END FUNCTION
  FUNCTION DL_Open(file,mode) RESULT(handle) BIND(C,NAME="DlOpen")
      ! void *dlopen(const char *file, int mode);
      USE ISO_C_BINDING
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: file
        ! C strings should be declared as character arrays
      INTEGER(C_INT), VALUE :: mode
      TYPE(C_PTR) :: handle
  END FUNCTION
  FUNCTION DL_Sym(handle,name) RESULT(funptr) BIND(C,NAME="DlSym")
      ! void *dlsym(void *handle, const char *name);
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: name
      TYPE(C_FUNPTR) :: funptr ! A function pointer
  END FUNCTION
  FUNCTION DL_Close(handle) RESULT(status) BIND(C,NAME="DlClose")
      ! int dlclose(void *handle);
      USE ISO_C_BINDING
      TYPE(C_PTR), VALUE :: handle
      INTEGER(C_INT) :: status
  END FUNCTION
  FUNCTION DL_Error() RESULT(error) BIND(C,NAME="DlError")
      ! char *dlerror(void);
      USE ISO_C_BINDING
      TYPE(C_PTR) :: error
  END FUNCTION
  FUNCTION C_F_STRING_CONVERT(CPTR) RESULT(FPTR)
    USE ISO_C_BINDING
    IMPLICIT NONE
    ! Convert a null-terminated C string into a Fortran character array pointer
    TYPE(C_PTR), INTENT(IN) :: CPTR ! The C address
    CHARACTER(KIND=C_CHAR), DIMENSION(:), POINTER :: FPTR
  END FUNCTION C_F_STRING_CONVERT
END INTERFACE
