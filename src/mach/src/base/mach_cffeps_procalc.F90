!============================================================================!
!
! Projet / Project : GEM-MACH
! Fichier / File   : mach_cffeps_procalc.ftn90
! Creation         : Kerry Anderson
!                  : (C To Fortran)
!                    A. Akingunola, J. Chen, and P. Makar - Fall 2018
! Description      : FBP 91 fire shape parameter calculations
! Notices          : Copyright (c) 1993, Forestry Canada, All Rights Reserved
!
! Extras:
!      From one non-zero fire shape input parameter (area, perimeter, forward
!      spread distance, time1, time2), PROCalc calculates all other fire shape
!      parameters.  Accelleration effects may be included.

!      Note that all desired outputs must be filled with 0 when called. If all
!      input parameters are 0, the outputs are zero.

!      If the headfire rate of spread (ROS) is zero, then the second elapsed time
!      is set to zero and all remaining parameters are untouched (presumably
!      all but one have already been set to zero by the calling routine).
!============================================================================
!!if_on
 subroutine mach_cffeps_procalc(ifuel, fbp, t1, t2, a, p, d)
   use mach_cffeps_mod,  only: fbp_type, cffeps_fbp_accel
!!if_off

   implicit none
!!if_on
   integer(kind=4), intent   (in) :: ifuel
   type(fbp_type),  intent   (in) :: fbp
   real(kind=4),    intent(inout) :: t1, t2, a, p, d
!!if_off
   real(kind=4) :: alpha, lt, lt1, lt2, la, lb, lm, lw, lp
   real(kind=4), parameter :: pi = 3.1415926

   if (ifuel == 1 .or. fbp%cfb == 0.0) then  ! case ("C1")
      alpha = 0.115
   else
      alpha = (0.115 - 18.8 * fbp%cfb**2.5 * exp(-8.0 * fbp%cfb)) !/* 72 */
   end if

   if (fbp%ros /= 0.0) then    ! /* rate of spread */
      la = (fbp%ros + fbp%bros) * 0.5
      lb = fbp%fros
      lm = (la - lb) / (la + lb)
      lw = ((fbp%ros + fbp%bros) * 0.5 + fbp%fros) * (1.0 + lm * lm * 0.25)
      lp = pi * 0.5 * (fbp%ros + fbp%bros) * fbp%fros
   end if

   if (fbp%ros == 0.0) then     ! /* no rate of spread */
      t2 = 999.9
   else if (t2 > 0.0) then ! /* hrs */
      if (cffeps_fbp_accel) then
         lt2 = t2 * 60.0
         lt1 = t1 * 60.0
      else
         lt2 = t2 * 60.0 + (exp(-alpha * t2 * 60.0) - 1.0) / alpha  !/* 71a */
         lt1 = t1 * 60.0 + (exp(-alpha * t1 * 60.0) - 1.0) / alpha  !/* 71a */
      end if

!/* NB: values calculated are the changes in distance/area/perimeter
!       over the time interval */

      d = fbp%ros * (lt2 - lt1) * 1.0e-3              ! /* 71 */
      a = lp * (lt2 * lt2 - lt1 * lt1) * 1.0e-4       ! /* 84 */
      p = pi * lw * (lt2 - lt1) * 1.0e-3              ! /* 85 */
   else if (a > 0.0) then ! /* ha  */
      lt = sqrt(a / lp) * 100.0                       ! /* 84 */
      d = (fbp%ros * lt) * 1.0e-3                     ! /* 71 */
      p = (pi * lt * lw) * 1.0e-3                     ! /* 85 */
      t2 = converge(alpha, lt)
      t1 = 0.0
   else if (p > 0.0) then  ! /* km  */
      lt = p / (pi * lw) * 1.0e3                      ! /* 85 */
      d = (fbp%ros * lt) * 1.0e-3                     ! /* 71 */
      a = (lp * lt * lt) * 1.0e-4                     ! /* 84 */
      t2 = converge(alpha, lt)
      t1 = 0.0
   else if(d > 0.0) then  ! /* km  */
      lt = (d / fbp%ros) * 1.0e3                      ! /* 71 */
      a = (lp * lt * lt) * 1.0e-4                     ! /* 84 */
      p = (pi * lt * lw) * 1.0e-3                     ! /* 85 */
      t2 = converge(alpha, lt)
      t1 = 0.0
   end if

   return

   contains
      real(kind=4) function converge(alpha, lt)
         implicit none

         real(kind=4), intent(in) :: alpha, lt

         real(kind=4)    :: llt
         integer(kind=4) :: i

         llt = lt
         if (cffeps_fbp_accel) then
            do i = 1, 10
               llt = lt - exp(-alpha * llt) / alpha + 1.0 / alpha
            end do
         end if
         converge = llt / 60.0

         return
      end function converge

 end subroutine mach_cffeps_procalc
