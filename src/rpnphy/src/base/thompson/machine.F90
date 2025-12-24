      module machine
         use, intrinsic :: iso_fortran_env, only: REAL32, REAL64, INT32
!! \section arg_table_machine
!! \htmlinclude machine.html
!!

      implicit none

#ifndef SINGLE_PREC
      integer, parameter :: kind_io4  = REAL32, kind_io8  = REAL64 , kind_ior = REAL64 &
     &,                     kind_evod = REAL64, kind_dbl_prec = REAL64            &
#ifdef __PGI
     &,                     kind_qdt_prec = REAL64                           &
#else
     &,                     kind_qdt_prec = 16                          &
#endif
     &,                     kind_rad  = REAL64                               &
     &,                     kind_phys = REAL32     ,kind_taum=REAL64              &
     &,                     kind_grid = REAL64                               &
     &,                     kind_REAL = REAL64                               &! used in cmp_comm
     &,                     kind_LOGICAL = 4                            &
     &,                     kind_INTEGER = INT32                             ! -,,-

#else
      integer, parameter :: kind_io4  = REAL32, kind_io8  = REAL64 , kind_ior = REAL64 &
     &,                     kind_evod = REAL32, kind_dbl_prec = REAL64            &
#ifdef __PGI
     &,                     kind_qdt_prec = REAL64                           &
#else
     &,                     kind_qdt_prec = 16                          &
#endif
     &,                     kind_rad  = REAL32                               &
     &,                     kind_phys = REAL32     ,kind_taum=REAL32              &
     &,                     kind_grid = REAL32                               &
     &,                     kind_REAL = REAL32                               &! used in cmp_comm
     &,                     kind_LOGICAL = 4                            &
     &,                     kind_INTEGER = INT32                             ! -,,-

#endif

#ifdef OVERLOAD_R4
      integer, parameter :: kind_dyn  = REAL32
#else
      integer, parameter :: kind_dyn  = REAL64
#endif

!
      real(kind=kind_evod), parameter :: mprec = 1.e-12           ! machine precision to restrict dep
      real(kind=kind_evod), parameter :: grib_undef = 9.99e20     ! grib undefine value
!
      end module machine
