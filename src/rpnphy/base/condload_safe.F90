!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html

!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

!> Perform the precipitation fallout
!>
!> @param[inout] qliq Liquid water in the cloud layer as input
!>    liquid water after the precipitation fallout as output
!> @param[inout] qice Cloud ice in the cloud layer as input
!>    cloud ice after the precipitation fallout as output
!> @param[inout] wtw W*W of the updraft as input new vertical
!>    velocity with the condensation load effect as output
!>
!> @param[in] dz Row number
!> @param[in] buoyancy Buoyancy term for the updraft
!> @param[in] enterm Entrainment term for the updraft
!> @param[in] rate Rate of precipitation fallout production in the
!>          updraft (= 0.01)
!> @param[in] qnewlq Fresh liquid water condensate in the updraft
!> @param[in] qnewic Fresh cloud ice in the updraft
!>
!> @param[out] qlqout Precipitation fallout (liquid water)
!> @param[out] qicout Precipitation fallout (cloud ice)
!>
!>
!>  9/18/88...This precipitation fallout scheme is based on the scheme
!>  used by Ogura and Cho (1973). Liquid water fallout from a parcel
!>  is calculated using the equation dq=-rate*q*dt, but to simulate a
!>  quasi-continuous process, and to eliminate a dependency on vertical
!>  resolution, this is expressed as q=q*exp(-rate*dz).
!>
!> @date 1998
!> @author Kain and Fritsch
!> @author Stephane Belair
!> @author A-M Leduc
subroutine condload_safe(qliq, qice, wtw, dz, buoyancy, enterm, rate, &
   qnewlq, qnewic, qlqout, qicout, grav)

#ifdef __INTEL_COMPILER
   use ifport
#endif
   implicit none

   real, intent(inout) :: qliq, qice, wtw
   real, intent(in) :: dz, buoyancy, enterm, rate, qnewlq, qnewic
   real, intent(out) :: qlqout, qicout

   real :: g1, qnew, qest, qtot, wavg
   ! These varaibles have values assigned, but never used
   ! I guess they exist to test the FP operations status
   real :: conv, ratio3, ratio4
   integer(2) :: control, control_nodivzero, status

   qtot = qliq + qice
   qnew = qnewlq + qnewic

   ! Estimate the vertical velocity so that
   ! an average vertical velocity can be
   ! calculated to estimate the time required
   ! ascent between model levels.

   qest = 0.5 * (qtot + qnew)

   ! Solve for the vertical motion.

   g1 = wtw + buoyancy - enterm - 2.0 * grav * dz * qest / 1.5
   if (g1 .lt. 0.0) g1 = 0.0
   wavg = (sqrt(wtw) + sqrt(g1)) / 2.0

   ! Conversion
#ifdef __INTEL_COMPILER
   call getcontrolfpqq(control)
   control_nodivzero = ior(control, FPCW$ZERODIVIDE)
   call setcontrolfpqq(control_nodivzero)
   ! Run all possible div-by-zero calculations.
   conv = rate * dz / wavg
   ratio3 = qnewlq / (qnew + 1.E-10)
   qtot = qtot + 0.6 * qnew
   ratio4 = (0.6 * qnewlq + qliq) / (qtot + 1.E-10)
   call getstatusfpqq(status)
   call clearstatusfpqq()
   call setcontrolfpqq(control)
   if (iand(status, FPSW$ZERODIVIDE) > 0) then
      call condload_infconv(qliq, qice, wtw, dz, buoyancy, enterm, &
         rate,qnewlq, qnewic, qlqout, qicout, grav)
   else
      call condload(qliq, qice, wtw, dz, buoyancy, enterm, rate, &
         qnewlq, qnewic, qlqout, qicout, grav)
   endif
#else
   call condload(qliq, qice, wtw, dz, buoyancy, enterm, rate, qnewlq, &
      qnewic, qlqout, qicout, grav)
#endif

   return
end subroutine condload_safe
