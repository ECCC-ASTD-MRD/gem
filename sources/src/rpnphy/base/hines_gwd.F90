!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!*** S/P VERT_SMOOTH
  SUBROUTINE vert_smooth (darr, work, coeff, nsmooth,         &
       &                  il1, il2, lev1, lev2, nlons, nlevs)
!  subroutine arguements.
    !
    implicit none
#include <arch_specific.hf>

    INTEGER  nsmooth, il1, il2, lev1, lev2, nlons, nlevs
    REAL*8  coeff
    REAL*8  darr(nlons,nlevs), work(nlons,nlevs)
!Author
!
!  aug. 3/95 - c. mclandress
!
!Reviison
!
!Object
!
!  Smooth a longitude by altitude array in the vertical over a
!  specified number of levels using a three point smoother. 
!
!Arguments
!
!  note: input array darr is modified on output!
!
!          - Output/Input -
! darr     smoothed array (on output),unsmoothed array of data (on input)
!
!          - Input
! darr     unsmoothed array of data (on input).
! work     work array of same dimension as darr.
! coeff    smoothing coefficient for a 1:coeff:1 stencil.
!          (e.g., coeff = 2 will result in a smoother which
!          weights the level l gridpoint by two and the two 
!          adjecent levels (l+1 and l-1) by one).
! nsmooth  number of times to smooth in vertical.
!          (e.g., nsmooth=1 means smoothed only once, 
!          nsmooth=2 means smoothing repeated twice, etc.)
! il1      first longitudinal index to use (il1 >= 1).
! il2      last longitudinal index to use (il1 <= il2 <= nlons).
! lev1     first altitude level to use (lev1 >=1). 
! lev2     last altitude level to use (lev1 < lev2 <= nlevs).
! nlons    number of longitudes.
! nlevs    number of vertical levels.
!

    !
    !  internal variables.
    !
    INTEGER  i, l, ns, lev1p, lev2m
    REAL*8  sum_wts
    !-----------------------------------------------------------------------     
    !
    !  calculate sum of weights.
    !
    sum_wts = coeff + 2.
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    !
    !  smooth nsmooth times
    !
    DO ns = 1,nsmooth
       !
       !  copy darr into work array.
       !
       DO l = lev1,lev2
          DO i = il1,il2
             work(i,l) = darr(i,l)
          END DO
       END DO
       !
       !  smooth array work in vertical direction and put into darr.
       !
       DO l = lev1p,lev2m
          DO i = il1,il2
             darr(i,l) = (work(i,l+1)+coeff*work(i,l)+work(i,l-1) ) / sum_wts 
          END DO
       END DO
    END DO

    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE vert_smooth
