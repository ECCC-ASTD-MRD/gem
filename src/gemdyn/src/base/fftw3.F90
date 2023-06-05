! FFTW3 module, as per directions in section 7.7 of the FFTW3 manual.

! Installing FFTW as a module helps segregate these routines from the
! remainder of the GEM code

module FFTW3
   use, intrinsic :: iso_c_binding
#include "fftw3.f03"
end module
