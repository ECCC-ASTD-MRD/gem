!!!s/r Extrap1D_Abort Interfaces
!
!AUTHOR
!     J.W. Blezius MAY 2014 First revision for this file
!
!REVISION
!         Blezius J.W. MAY 2014 - add abortion on extrapolation (this file)
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
subroutine Extrap1D_Abort_X  &
#define real48 single
#include "Extrap1D_Abort_Body.inc"
#undef real48

!version that accepts real(double) arguments, with the extended interface
subroutine Extrap1D_Abort_X8  &
#define real48 double
#include "Extrap1D_Abort_Body.inc"
#undef real48




!version that accepts real(single) arguments
subroutine Extrap1D_Abort  &
#define real48 single
#define DestnFunc Extrap1D_Abort_X
#include "Interp1D_PreX_Shell.inc"
#undef DestnFunc
#undef real48

!version that accepts real(double) arguments
subroutine Extrap1D_Abort8  &
#define real48 double
#define DestnFunc Extrap1D_Abort_X8
#include "Interp1D_PreX_Shell.inc"
#undef DestnFunc
#undef real48
