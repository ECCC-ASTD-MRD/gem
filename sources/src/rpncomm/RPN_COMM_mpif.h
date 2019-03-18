!
!     FAKE mpif.h that contains enough definitions to compile RPN_COMM
!     and its stubs
!
!     NULL "handles"
!
      integer MPI_GROUP_NULL, MPI_COMM_NULL, MPI_DATATYPE_NULL
      integer MPI_REQUEST_NULL, MPI_OP_NULL, MPI_ERRHANDLER_NULL
      integer MPI_INFO_NULL, MPI_WIN_NULL

      parameter (MPI_GROUP_NULL=0)
      parameter (MPI_COMM_NULL=2)
      parameter (MPI_DATATYPE_NULL=0)
      parameter (MPI_REQUEST_NULL=0)
      parameter (MPI_OP_NULL=0)
      parameter (MPI_ERRHANDLER_NULL=0)
      parameter (MPI_INFO_NULL=0)
      parameter (MPI_WIN_NULL=0)
!
!     varia
!
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=10)
      integer MPI_ANY_SOURCE, MPI_ANY_TAG
      parameter (MPI_ANY_SOURCE=-1)
      parameter (MPI_ANY_TAG=-1)
      integer MPI_SUCCESS
      parameter( MPI_SUCCESS = 0)

      integer MPI_COMM_WORLD, MPI_COMM_SELF
      integer MPI_GROUP_EMPTY
      integer MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN

      parameter (MPI_COMM_WORLD=0)
      parameter (MPI_COMM_SELF=1)
      parameter (MPI_GROUP_EMPTY=1)
      parameter (MPI_ERRORS_ARE_FATAL=1)
      parameter (MPI_ERRORS_RETURN=2)

!
!     data types
!
      integer MPI_BYTE, MPI_PACKED, MPI_UB, MPI_LB
      integer MPI_CHARACTER, MPI_LOGICAL
      integer MPI_INTEGER, MPI_INTEGER1, MPI_INTEGER2, MPI_INTEGER4
      integer MPI_INTEGER8, MPI_INTEGER16
      integer MPI_REAL, MPI_REAL2, MPI_REAL4, MPI_REAL8, MPI_REAL16
      integer MPI_DOUBLE_PRECISION
      integer MPI_COMPLEX, MPI_COMPLEX8, MPI_COMPLEX16, MPI_COMPLEX32
      integer MPI_DOUBLE_COMPLEX
      integer MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
      integer MPI_2COMPLEX, MPI_2DOUBLE_COMPLEX
      integer MPI_LOGICAL1, MPI_LOGICAL2, MPI_LOGICAL4, MPI_LOGICAL8

      parameter (MPI_BYTE=1)
      parameter (MPI_PACKED=2)
      parameter (MPI_UB=3)
      parameter (MPI_LB=4)
      parameter (MPI_CHARACTER=5)
      parameter (MPI_LOGICAL=6)
      parameter (MPI_INTEGER=7)
      parameter (MPI_INTEGER1=8)
      parameter (MPI_INTEGER2=9)
      parameter (MPI_INTEGER4=10)
      parameter (MPI_INTEGER8=11)
      parameter (MPI_INTEGER16=12)
      parameter (MPI_REAL=13)
      parameter (MPI_REAL4=14)
      parameter (MPI_REAL8=15)
      parameter (MPI_REAL16=16)
      parameter (MPI_DOUBLE_PRECISION=17)
      parameter (MPI_COMPLEX=18)
      parameter (MPI_COMPLEX8=19)
      parameter (MPI_COMPLEX16=20)
      parameter (MPI_COMPLEX32=21)
      parameter (MPI_DOUBLE_COMPLEX=22)
      parameter (MPI_2REAL=23)
      parameter (MPI_2DOUBLE_PRECISION=24)
      parameter (MPI_2INTEGER=25)
      parameter (MPI_2COMPLEX=26)
      parameter (MPI_2DOUBLE_COMPLEX=27)
      parameter (MPI_REAL2=28)
      parameter (MPI_LOGICAL1=29)
      parameter (MPI_LOGICAL2=30)
      parameter (MPI_LOGICAL4=31)
      parameter (MPI_LOGICAL8=32)
!
!     operators
!
      integer MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_LAND
      integer MPI_BAND, MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR
      integer MPI_MAXLOC, MPI_MINLOC, MPI_REPLACE

      parameter (MPI_MAX=1)
      parameter (MPI_MIN=2)
      parameter (MPI_SUM=3)
      parameter (MPI_PROD=4)
      parameter (MPI_LAND=5)
      parameter (MPI_BAND=6)
      parameter (MPI_LOR=7)
      parameter (MPI_BOR=8)
      parameter (MPI_LXOR=9)
      parameter (MPI_BXOR=10)
      parameter (MPI_MAXLOC=11)
      parameter (MPI_MINLOC=12)
      parameter (MPI_REPLACE=13)
