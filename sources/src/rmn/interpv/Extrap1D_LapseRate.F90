!!!s/r Extrap1D_LapseRate Interfaces
!
!AUTHOR
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - new code
! v1_3    Blezius J.W. OCT 2003 - add the extended interface
!
!OBJECT
!        To provide a means to easily compile the contained routine twice, once
!        with single-precision arguments, and once with double-precision
!        arguments.  The unextended interface is likewise compiled twice.
!
!NOTES
!
!!

!version that accepts real(single) arguments, with the extended interface
subroutine Extrap1D_LapseRate_X  &
#define real48 single
#include "Extrap1D_LapseRate_Body.inc"
#undef real48

!version that accepts real(double) arguments, with the extended interface
subroutine Extrap1D_LapseRate_X8  &
#define real48 double
#include "Extrap1D_LapseRate_Body.inc"
#undef real48




!version that accepts real(single) arguments
subroutine Extrap1D_LapseRate  &
#define real48 single
#define DestnFunc Extrap1D_LapseRate_X
#include "Interp1D_PreX_Shell.inc"
#undef DestnFunc
#undef real48

!version that accepts real(double) arguments
subroutine Extrap1D_LapseRate8  &
#define real48 double
#define DestnFunc Extrap1D_LapseRate_X8
#include "Interp1D_PreX_Shell.inc"
#undef DestnFunc
#undef real48
