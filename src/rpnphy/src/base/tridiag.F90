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

module tridiag
  implicit none
  private

  ! API functions
  public :: td_solve

contains

  subroutine td_solve(F_fld, F_a, F_b, F_c, F_d, F_ni, F_nk)
    implicit none

    ! Argument description
    integer, intent(in) :: F_ni                         !horizontal dimension
    integer, intent(in) :: F_nk                         !vertical dimension
    real, dimension(F_ni,F_nk), intent(in) :: F_a       !matrix subdiagonal
    real, dimension(F_ni,F_nk), intent(in) :: F_b       !matrix diagonal
    real, dimension(F_ni,F_nk), intent(in) :: F_c       !matrix superdiagonal
    real, dimension(F_ni,F_nk), intent(in) :: F_d       !rhs forcing
    real, dimension(F_ni,F_nk), intent(out) :: F_fld    !updated atmospheric field

    !@Author R. McTaggart-Cowan and A. Zadra (Spring 2024)
    !@Revisions
    !@Object
    !    Solve the matrix problem Mx = F_d, where x is F_fld and M is the
    !    tridiagonal matrix consisting of vectors F_a, F_b and F_c.
    !    This solver returns the updated F_fld value using a standard
    !    Thomas algorithm.
    !*@/

    ! Internal declarations
    integer :: i,k
    real, dimension(F_nk) :: cprime
    
    ! Solve to update field value
    do i=1,F_ni
       ! Compute new coefficients
       cprime(1) = F_c(i,1) / F_b(i,1)
       F_fld(i,1) = F_d(i,1) / F_b(i,1)
       do k=2,F_nk
          cprime(k) = F_c(i,k) / (F_b(i,k) - F_a(i,k)*cprime(k-1))
          F_fld(i,k) = (F_d(i,k) - F_a(i,k) * F_fld(i,k-1)) / &
               (F_b(i,k) - F_a(i,k)*cprime(k-1))
       enddo
       ! Update field via back-substitution
       do k=F_nk-1,1,-1
          F_fld(i,k) = F_fld(i,k) - cprime(k) * F_fld(i,k+1)
       enddo
    enddo

  end subroutine td_solve

end module tridiag
