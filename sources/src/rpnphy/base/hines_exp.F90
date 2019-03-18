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
SUBROUTINE hines_exp (darr, data_zmax, alt, alt_exp,     &
       &                il1, il2, lev1, lev2, nlons, nlevs)

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

    implicit none
#include <arch_specific.hf>

    INTEGER  il1, il2, lev1, lev2, nlons, nlevs
    REAL*8  alt_exp
    REAL*8  darr(nlons,nlevs), data_zmax(nlons), alt(nlons,nlevs)
    REAL*8  exp_fac(nlons,nlevs)
!*@/

    INTEGER  levbot, levtop, lincr, i, l
    REAL*8  hscale
    !-----------------------------------------------------------------------
    hscale = 5.e3

    !  index of lowest altitude level (bottom of drag calculation).
    !
    levbot = lev2
    levtop = lev1
    lincr  = 1
    IF (levbot > levtop)  THEN
       levbot = lev1
       levtop = lev2
       lincr  = -1
    ELSE
       call physeterror('hines_exp', 'level index not increasing downward')
       return
    END IF
    !
    !  data values at first level above alt_exp.
    !
    DO i = il1,il2
       DO l = levtop,levbot,lincr
          IF (alt(i,l) .GE. alt_exp)  THEN
             data_zmax(i) = darr(i,l)
          END IF
       END DO
    END DO
    !
    !  exponentially damp field above alt_exp to model top at l=1.
    !
    DO l = 1,lev2
       DO i = il1,il2
         exp_fac(i,l)=(alt_exp-alt(i,l))/hscale
       END DO
    END DO
    call vexp(exp_fac,exp_fac,lev2*(il2-il1+1))

    DO l = 1,lev2
       DO i = il1,il2
          IF (alt(i,l) .GE. alt_exp)  THEN
             darr(i,l) = data_zmax(i) * exp_fac(i,l)
          END IF
       END DO
    END DO
    !
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_exp
