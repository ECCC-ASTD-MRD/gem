!-------------------------------------- LICENCE BEGIN --------------------------
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
   subroutine linoz_xcol( x_vmr, p_low, p_up, x_over, ni,nkm1)

!-----------------------------------------------------------------------
! Calculate overhead columns *at* levels:  cst_o3 * SUM( vmr * delta_p ) 
! where delta_p is between levels and vmr are at mid-levels
!-----------------------------------------------------------------------

      use linoz_param, only: avno, mwt_air, du_o3
      use tdpack_const, only: GRAV

      implicit none
#include <arch_specific.hf>
!-----------------------------------------------------------------------
!  input x_col : overhead column (D.U.) for level above
! output x_over : overhead column (D.U.) for current level
!-----------------------------------------------------------------------

    integer,    intent(in)    :: ni,nkm1

    real, intent(in)    :: x_vmr(ni, nkm1)         ! mole /mole vmr
    real, intent(in)    :: p_low(ni, nkm1), p_up(ni, nkm1)   ! pressure (Pa) at layer interfaces
    real, intent(inout) :: x_over(ni, nkm1)

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------

! avno      = 0.6022000000000e+24 ! avogadro's num        atoms mol-1
! mwt_air   = 0.2897000000000e+02 ! mol. wgt. of air
! du_o3     = 0.2687000000000e+20 !1 DU = 2.687 E-19 molecules of O3 per square metre
!

    integer :: i, k
    real :: cst_o3, x_col
    real :: du1_o3    = 0.2687000000000e+1 !1 DU = 2.687 E-19 molecules of O3 per square metre
!   real :: avno1     = 0.6022000000000e+5 !avogadro's num   0.6022e+24   atoms mol-1 
    real :: avno2     = 0.6022000000000e+1 !avogadro's num   0.6022e+24   atoms mol-1

!      cst_o3 = (1.0E-4/grav) * ( avno / (1.0E-3 * mwt_air) )* 1.0E+03 / du_o3  ! molec/cm2 -> DU
!      cst_o3 = (1.0E-4/grav) * ( avno1 / mwt_air) / du1_o3  ! molec/cm2 -> DU

    cst_o3 = (1.0E+6/grav) * ( avno2 / mwt_air) / du1_o3  ! molec/cm2 -> DU

    do i = 1, ni

        x_col = 0.

    do k = 1, nkm1

      if( x_col == 0. ) then
         x_col = cst_o3 * x_vmr(i,k) * ( p_up(i,k) - 0. )
       else
         x_col = x_col + cst_o3 * ( p_low(i,k) - p_up(i,k) ) * x_vmr(i,k)
      end if
      x_over(i,k) = x_col

    enddo
    enddo
         
      !----------------------------------------------------------------
      return
   end subroutine linoz_xcol
