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
   subroutine ghg_xcol( x_vmr, p_low, p_up, x_col)

!-----------------------------------------------------------------------
! Calculate overhead columns *at* levels:  cst_x * SUM( vmr * delta_p ) 
! where delta_p is between levels and vmr are at mid-levels
!-----------------------------------------------------------------------

      use linoz_mod, only: avno, mwt_air, local_dbg
      use tdpack_const, only: GRAV
      implicit none
#include <arch_specific.hf>

!-----------------------------------------------------------------------
!  input x_vmr : concentration in mole/mole vmr
! output x_col : overhead column (molec /cm2) for current level
!-----------------------------------------------------------------------

      real, intent(in)    :: x_vmr                      ! mole/mole vmr 
      real, intent(in)    :: p_low, p_up                ! pressure (Pa) at layer interfaces
      real, intent(inout) :: x_col                      ! 1E+15 (Peta) molecules /cm2

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------

! avno      = 0.6022000000000e+24 ! avogadro's num        atoms mol-1
! mwt_air   = 0.2897000000000e+02 ! mol. wgt. of air  g mol-1
!
      real :: cst_x, x_nmr, dpgm, x_pcol

! Column X
!
! Col X (molec/cm2) = Y [mole/mole] /mwt_air * avno * 1.E+3 (g/kg) * 1E-4 (m2/cm2)
!
      ! 1E-4 m2/cm2; 1E+3 g/kg => 1E-1
      cst_x = grav /mwt_air * avno * 1.0E+1    ! molec/cm2 

      x_nmr = x_vmr * avno                     ! molecules (x) mole^-1
      dpgm = 1.0E-1*(p_low-p_up)/grav/mwt_air    ! mole (air) cm-2
      x_pcol = 1.0E-15* x_nmr * dpgm           ! 1E+15 (Peta) molecules cm-2
      x_col = x_col + x_pcol                   ! 1E+15 (Peta) molecules cm-2

      !if( x_col == 0. ) then
      !   x_col = cst_x * x_vmr * ( p_up - 0. )
      ! else
      !  x_col = x_col + cst_x * ( p_low - p_up ) * x_vmr
      ! end if

!      if (local_dbg) then
!         write(lun_out, *) "in ghg_xcol: dpgm, vmr, pcol, col=", dpgm, x_vmr, x_pcol, x_col
!      end if

      !----------------------------------------------------------------
      return
   end subroutine ghg_xcol
