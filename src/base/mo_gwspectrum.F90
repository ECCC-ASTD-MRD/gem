
module mo_gwspectrum

   use, intrinsic :: iso_fortran_env, only: REAL32

   real(REAL32), parameter :: kstar = 5.E-5
   integer, parameter :: naz = 8
   real(REAL32), parameter :: rnaz = 1./float(naz)
   real(REAL32), parameter :: m_min = 5.23598E-4
   real(REAL32), parameter :: slope = 1.
   real(REAL32), parameter :: f1 = 1.5
   real(REAL32), parameter :: f2 = 0.3
   real(REAL32), parameter :: f3 = 1.0
   real(REAL32), parameter :: f5 = 2.0
   real(REAL32), parameter :: f6 = 1.0
   real(REAL32), parameter :: alt_cutoff = 105.E3
   integer, parameter :: icutoff = 0
   real(REAL32), parameter :: smco = 2.0
   real(REAL32), parameter :: rsum_wts = 1. / (smco + 2.)
   integer, parameter :: nsmax = 5
   logical, save :: hasspp = .false.

end module mo_gwspectrum
