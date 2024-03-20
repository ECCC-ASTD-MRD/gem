!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

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
