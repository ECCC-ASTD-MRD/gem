module linoz_xcol
   implicit none
   private
   public :: linoz_xcol1, ghg_xcol1

contains

   subroutine linoz_xcol1(x_vmr, p_low, p_up, x_over, ni,nkm1)
      use linoz_param, only: avno, mwt_air, du_o3
      use tdpack_const, only: GRAV
      implicit none

      !-----------------------------------------------------------------------
      ! Calculate overhead columns *at* levels:  cst_o3 * SUM( vmr * delta_p ) 
      ! where delta_p is between levels and vmr are at mid-levels
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      !  input x_col : overhead column (D.U.) for level above
      ! output x_over : overhead column (D.U.) for current level
      !-----------------------------------------------------------------------

      integer, intent(in) :: ni, nkm1
      real, intent(in)    :: x_vmr(ni,nkm1)                 ! mole /mole vmr
      real, intent(in)    :: p_low(ni,nkm1), p_up(ni,nkm1)  ! pressure (Pa) at layer interfaces
      real, intent(inout) :: x_over(ni,nkm1)

      !-----------------------------------------------------------------------
      !  Local variables
      !-----------------------------------------------------------------------

      ! avno      = 0.6022000000000e+24 ! avogadro's num        atoms mol-1
      ! mwt_air   = 0.2897000000000e+02 ! mol. wgt. of air
      ! du_o3     = 0.2687000000000e+20 !1 DU = 2.687 E-19 molecules of O3 per square metre

      integer :: i, k
      real :: x_col(ni)

      real, parameter :: du1_o3 = 0.2687000000000e+1  !1 DU = 2.687 E-19 molecules of O3 per square metre
      real, parameter :: avno2  = 0.6022000000000e+1  !avogadro's num   0.6022e+24   atoms mol-1
      real, parameter :: cst_o3 = (1.0E+6/grav) * (avno2 / mwt_air) / du1_o3  ! molec/cm2 -> DU    

      !----------------------------------------------------------------
      x_col = 0.
      do k = 1, nkm1
         do i = 1, ni
            if( x_col(i) == 0. ) then
               x_col(i) = cst_o3 * x_vmr(i,k) * p_up(i,k)
            else
               x_col(i) = x_col(i) + cst_o3*(p_low(i,k)-p_up(i,k))*x_vmr(i,k)
            end if
            x_over(i,k) = x_col(i)
         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine linoz_xcol1

   
   subroutine ghg_xcol1(x_vmr, p_low, p_up, x_over, ni, nkm1, hu_vmr)
      use linoz_param, only: avno, mwt_air, hu_linoz_tropo
      use tdpack_const, only: GRAV
      implicit none

      !-----------------------------------------------------------------------
      ! Calculate overhead columns *at* levels:  cst_x * SUM( vmr * delta_p ) 
      ! where delta_p is between levels and vmr are at mid-levels
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      !  input x_vmr : concentration in mole/mole vmr
      ! output x_over : overhead column (molec /cm2) for current level
      !-----------------------------------------------------------------------

      integer, intent(in) :: ni, nkm1

      real, intent(in)    :: x_vmr(ni,nkm1)                 ! mole/mole vmr 
      real, intent(in)    :: p_low(ni,nkm1), p_up(ni,nkm1)  ! pressure (Pa) at layer interfaces
      real, intent(inout) :: x_over(ni,nkm1)                ! 1E+15 (Peta) molecules /cm2
      real, intent(in), optional :: hu_vmr(ni,nkm1)         ! ppmv vmr

      !-----------------------------------------------------------------------
      !  Local variables
      !-----------------------------------------------------------------------

      integer :: i, k
      real :: x_nmr, dpgm, x_pcol, x_col(ni)
      logical :: mask(ni,nkm1)
      !----------------------------------------------------------------
      if (present(hu_vmr)) then
         mask = (hu_vmr >= hu_linoz_tropo)
      else
         mask = .true.
      endif
       
      x_col = 0.
      do k = 1, nkm1
         do i = 1, ni
            if (mask(i,k)) then
               x_nmr = x_vmr(i,k) * avno          ! molecules (x) mole^-1
               dpgm = 1.0E-1*(p_low(i,k)-p_up(i,k))/grav/mwt_air  ! mole (air) cm-2
               x_pcol = 1.0E-15 * x_nmr * dpgm    ! 1E+15 (Peta) molecules cm-2
               x_col(i) = x_col(i) + x_pcol       ! 1E+15 (Peta) molecules cm-2

               x_over(i,k) = x_col(i)
            endif
         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine ghg_xcol1

end module linoz_xcol
