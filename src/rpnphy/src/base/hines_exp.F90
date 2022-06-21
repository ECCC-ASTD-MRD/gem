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
!-------------------------------------- LICENCE END --------------------------

!/@*
subroutine hines_exp(darr, data_zmax, alt, alt_exp,     &
     &             il1, il2, lev1, lev2, nlons, nlevs)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   !@Author  aug. 13/95 - c. mclandress
   !@Object
   !  This routine exponentially damps a longitude by altitude array
   !  of darr above a specified altitude.
   !@Arguments
   !            - Output -
   ! darr       modified data array.
   !            - Input arguments:
   ! darr       original data array.
   ! alt        altitudes.
   ! alt_exp    altitude above which exponential decay applied.
   ! il1        first longitudinal index to use (il1 >= 1).
   ! il2        last longitudinal index to use (il1 <= il2 <= nlons).
   ! lev1       first altitude level to use (lev1 >=1).
   ! lev2       last altitude level to use (lev1 < lev2 <= nlevs).
   ! nlons      number of longitudes.
   ! nlevs      number of vertical
   ! data_zmax  data values just above altitude alt_exp (work array)

   integer  il1, il2, lev1, lev2, nlons, nlevs
   real(REAL64) :: alt_exp
   real(REAL64) :: darr(nlons,nlevs), data_zmax(nlons), alt(nlons,nlevs)
   !*@/

   integer  levbot, levtop, lincr, i, l
   real(REAL64) :: rhscale, exp_fac
   !-----------------------------------------------------------------------
    rhscale = 1./5.e3
   

    !  index of lowest altitude level (bottom of drag calculation).
 
    levbot = lev2
    levtop = lev1
    lincr  = 1
    if (levbot > levtop)  then
       levbot = lev1
       levtop = lev2
       lincr  = -1
    else
       call physeterror('hines_exp', 'level index not increasing downward')
       return
    end if

    !  data values at first level above alt_exp.

    do i = il1,il2
       do l = levtop,levbot,lincr
          if (alt(i,l) .ge. alt_exp)  then
             data_zmax(i) = darr(i,l)
          end if
       end do
    end do

    !  exponentially damp field above alt_exp to model top at l=1.

    do l = 1,lev2
       do i = il1,il2
          if (alt(i,l) .ge. alt_exp)  then
             exp_fac = exp((alt_exp-alt(i,l))*rhscale)
             darr(i,l) = data_zmax(i) * exp_fac
          end if
       end do
    end do

    return
    !-----------------------------------------------------------------------
 end subroutine hines_exp
