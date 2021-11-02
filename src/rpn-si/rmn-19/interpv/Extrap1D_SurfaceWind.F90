!!!s/r Extrap1D_SurfaceWind Interfaces
!
!AUTHOR
!     J.W. Blezius OCT 2003 new kind of extrapolation
!
!REVISION
! v1_3    Blezius J.W. OCT 2003 - new code
!
!OBJECT
!        To provide a means to easily compile the contained routine twice, once
!        with single-precision arguments, and once with double-precision
!        arguments.
!
!NOTES
!        Because this routine did not exist prior to v1.3, the 'unextended'
!        interface is not supported.  More to the point, this routine cannot even
!        function with the 'unextended' interface; this is the reason that the
!        interface was extended.
!
!!

!version that accepts real(single) arguments, with the extended interface
subroutine Extrap1D_SurfaceWind_X  &
#define real48 single
#include "Extrap1D_SurfaceWind_Body.inc"
#undef real48

!version that accepts real(double) arguments, with the extended interface
subroutine Extrap1D_SurfaceWind_X8  &
#define real48 double
#include "Extrap1D_SurfaceWind_Body.inc"
#undef real48
