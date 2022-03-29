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

subroutine lightning2(foudre_rt, zp0_plus, zsigm, ztplus, zwplus, q_grpl, iiwc, ni, nk)
   use tdpack_const, only: GRAV, TCDK
   implicit none

!!!#include <arch_specific.hf>

   !@Arguments
   integer, intent(in)                  :: ni,nk      !horizontal/vertical dimensions
   real, dimension(ni),    intent(out)  :: foudre_rt  !lightning flash rate density,            # m-2 s-1
   real, dimension(ni),    intent(in)   :: zp0_plus   !surface pressure                         Pa
   real, dimension(ni,nk), intent(in)   :: zsigm      !sigma
   real, dimension(ni,nk), intent(in)   :: ztplus     !air temperature                          K
   real, dimension(ni,nk), intent(in)   :: zwplus     !vertical motion
   real, dimension(ni,nk-1), intent(in) :: q_grpl     !mixing ratio of graupel+hail             kg kg-1
   real, dimension(ni,nk-1), intent(in) :: iiwc       !mixing ratio of sum of all ice species   kg kg-1

   !@Author Anna Glazer   (Sept 2014)
   !@Revisions
   ! Anna Glazer (Avril 2015) - Pointers for version 4.8-a1
   ! Anna Glazer (Mai 2015)   - Adapted to 4.8.a5
   !@Object
   !  Compute the lightning threat expressed as
   !  flash-rate density (number of flashes per sec and m2)
   !  F3 = r1*F1 + F2*F2   ! f3 = f(foudre)
   !
   !  ref: McCaul et al., Wea. Forecasting 2009, vol. 40, pp. 709-729


   !Local variables and parameters:
   integer :: i, k, k_T15
   real    :: dsg, dpsg, i_conv, i_grav, iiwp, F1, F2, minDiff, tmp

   real, parameter :: CONV  = 3.e+8       !to convert f3 unit from 5minlcn_p3i1)
   real, parameter :: FMIN  = 1.e-12      !s-1 m-2, flash-rate density threshold for output
   real, parameter :: SEUIL = 1.e-6       !kg kg-1, mixing ratio threshold for graupel flux calculation
   real, parameter :: T15   = TCDK - 15.  ! = 258.15 (K)
   real, parameter :: R1    = 0.95        !tunable parameter in McCaul el al. (2009) lighting diagnostic
   real, parameter :: R2    = 0.05        !tunable parameter
   real, parameter :: K1    = 0.042       !tunable parameter
   real, parameter :: K2    = 0.20        !tunable parameter

   i_conv = 1./CONV
   i_grav = 1./GRAV

   I_LOOP: do i = 1,ni

      !find lowest level k where T(k) is closest to -15 C:
      !  (note: starts at nk-1, not nk, since we need a common level with T (and w) and q_grpl)
      k_T15   = 1
      minDiff = 999.
      do k = nk-1,1,-1
         tmp = abs(ztplus(i,k)-T15)
         if (tmp < minDiff) then
            minDiff = tmp
            k_T15 = k
            if (minDiff < 1.) exit   !assume level is found (1 deg differenc is close enough)
         endif
      enddo

      !compute F1 (function of graupel flux at -15 C):
      if (zwplus(i,k_T15) > 0. .and. q_grpl(i,k_T15) > SEUIL) then
         F1 = 1.e+3*K1*zwplus(i,k_T15)*q_grpl(i,k_T15)
      else
         F1 = 0.
      endif

      !compute F2 (function of ice-water path):
      iiwp = 0.
      do k = 2,nk-1
         dsg  = 0.5 * ( zsigm(i,k+1) - zsigm(i,k-1) )
         dpsg = zp0_plus(i)*dsg*i_grav
         iiwp = iiwp + max(iiwc(i,k), 0.) * dpsg
      enddo
      F2 = K2*iiwp

      !combine F1 + F2 to compute lighting threat (flash density rate; F3 in McCaul et al. 2008):
      foudre_rt(i) = i_conv*(R1*F1 + R2*F2)
      if (foudre_rt(i) < FMIN) foudre_rt(i) = 0.

   enddo I_LOOP

   return
end subroutine lightning2

