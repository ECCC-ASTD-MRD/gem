!!!s/r Interp1D_FindPos Interfaces
!
!AUTHOR
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - new code
! v1_3    Blezius J.W. OCT 2003 - add this header
!
!OBJECT
!        To provide a means to easily compile the contained routine twice, once
!        with single-precision arguments, and once with double-precision
!        arguments.
!
!NOTES
!        Unlike the other routines, the interface to this routine was not
!        extended in revision 1.3; there is not 'extended' interface.
!
!!

!version that accepts real(single) arguments
subroutine Interp1D_FindPos  &
#define real48 single
#include "Interp1D_FindPos_Body.inc"
#undef real48

!version that accepts real(double) arguments
subroutine Interp1D_FindPos8  &
#define real48 double
#include "Interp1D_FindPos_Body.inc"
#undef real48
