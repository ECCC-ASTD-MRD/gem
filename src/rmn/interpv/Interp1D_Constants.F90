!!!mod Interp1D_Constants - Constants for 1-D interpolation package
module Interp1D_Constants
!
!AUTHOR
!     J.W. Blezius MAY 2002 first library to replace duplicate interpolation
!                           routines
!
!REVISION
! v1_0    Blezius J.W.          - initial version
!
!OBJECT
!        This file gathers together into one place the constants that control
!        the 1-D interpolation package.
!
!!
! variable precision -- *** COMPILER-DEPENDENT VALUES ***
  public
  integer, parameter :: single=4, &     ! 4-bit reals, used for interpolation
                        double=8        ! 8-bit reals, used for interpolation

end module Interp1D_Constants 
