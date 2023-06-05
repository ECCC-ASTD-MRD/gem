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


subroutine secajus(tconv, t    , s    , ps,   niter, &
     conv , cdt1 , ni   , nk)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>

   integer ni,nk, niter
   real cdt1, conv

   real ps(ni)
   real s(ni,nk), t(ni,nk), tconv(ni,nk)

   !@Author Alain Patoine

   !@Revision
   !
   ! 001    B. Bilodeau (Jan 1997) - Adaptation from 3D to 2D;
   !        add tendencies calculations and dynamic allocation.
   ! 002    C. Girard (Mar 1997) - Conserve theta rather than T.
   ! 003    B. Bilodeau (Jan 2001) - Automatic arrays

   !@Object perform a dry convective adjustment

   !@Notes
   !     The algorith is the same as the one used in an old version of the rfe
   !     model. Examples of the original code were provided by both C. Beaudoin
   !     and B. Bilodeau

   !@Arguments
   !          - Output -
   ! tconv    temperature tendency due to dry convective adjustement
   ! niter    number of iterations
   !          - Input -
   ! t        temperature field
   ! s        sigma levels
   ! conv     convergence criteria
   ! cdt1     timestep
   ! ni       field dimension in x-direction
   ! nk       field dimension in z-direction

   integer, parameter :: nitmax=25

   logical adj

   integer i, k, nkm, nkp

   real, dimension(ni)    :: wrk11,wrk12,wrk14
   real, dimension(ni,nk) :: wrk21,wrk22,wrk23,p

   !***********************************************************************
 
   !     calcul de la pression
   do k=1,nk
      do i=1,ni
         p(i,k) = s(i,k) * ps(i)
      enddo
   enddo

   do i=1,ni
      wrk23(i,1)   = (p(i,2)   - p(i,1)  )   * 0.5
   enddo

   do k=2,nk-1
      do i=1,ni
         wrk23(i,k)   = (p(i,k+1) - p(i,k-1))   * 0.5
      enddo
   enddo

   do i=1,ni
      wrk23(i,nk) = (p(i,nk) - p(i,nk-1)) * 0.5
   enddo

   do k=1,nk-1
      do i=1,ni
         wrk21(i,k) = (p(i,k)/p(i,k+1))**cappa
         wrk23(i,k) = wrk23(i,k+1) / wrk23(i,k)
         wrk22(i,k) = 1. / (1.+wrk23(i,k))
      enddo
   enddo

   !---------------------------------------------------------------------------

   DOITER: do niter=1,nitmax
      adj = .false.

      do i=1,ni
         wrk11(i) = t(i,1)
      enddo

      do k=1,nk-1
         nkm = k-1
         nkp = k+1

         if (k .gt. 1) then
            do i=1,ni
               tconv(i,nkm) = (wrk11(i) - t(i,nkm))/cdt1
               wrk11(i)     = wrk12(i)
            enddo
         endif

         do i=1,ni
            wrk12(i) = t(i,nkp)
            wrk14(i) = wrk11(i)-wrk12(i)*wrk21(i,k)

            wrk14(i) = max(0.,-wrk14(i))*wrk22(i,k)

            wrk11(i) = wrk11(i)+wrk14(i)*wrk23(i,k)
            wrk12(i) = wrk12(i)-wrk14(i)/wrk21(i,k)

            if (abs(wrk14(i)) .gt. conv) adj=.true.
         enddo

      enddo

      do i=1,ni
         tconv(i,nk-1) = (wrk11(i) - t(i,nk-1))/cdt1
         tconv(i,nk  ) = (wrk12(i) - t(i,nk  ))/cdt1
      enddo

   enddo DOITER

   return
end subroutine secajus
