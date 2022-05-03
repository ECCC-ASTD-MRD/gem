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

subroutine update4(TS, T2, WG, W2, WF, WL, WR, WS, &
     ALPHAS, RHOS, RHOSNO, VMOD, CD, RHOA, &
     HFLUX, EFLUX, TST, T2T, WGT, W2T, &
     WFT, WLT, WRT, WST, &
     ALPHAST, RHOST, &
     BM, FQ, ALFAT, ALFAQ, &
     CTU, TH, HU, &
     N)
   use tdpack_const, only: CPD, RAUW
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>

   integer N
   real TS(N), T2(N), WG(N), W2(N), WR(N), WS(N)
   real WF(N), WL(N)
   real ALPHAS(N), RHOS(N), RHOSNO(N)
   real VMOD(N), CD(N), RHOA(N)
   real HFLUX(N), EFLUX(N)
   real TST(N), T2T(N), WGT(N), W2T(N), WRT(N), WST(N)
   real WFT(N), WLT(N)
   real ALPHAST(N), RHOST(N)
   real BM(N), FQ(N)
   real ALFAT(N), ALFAQ(N)
   real CTU(N), TH(N), HU(N)

   !@Author S. Belair (January 1997)
   !@Revisions
   ! 001      J. Mailhot (March 1998) - Addition of sea ice surface
   ! 002      J. Mailhot (Oct 1998) - Remove AQ and add ALFAQ
   !                                  Change name from UPDATE1 to UPDATE2
   ! 003      S. Belair (Jan 1999) - New variables WF and WL
   ! 004      S. Belair (Feb 1999) - Change formulation for the vertical
   !                                 diffusion boundary terms, and cleanup
   ! 005      S. Belair (March 1999)
   !                                Remove temporal filtering for ISBA's
   !                                prognostic variables
   ! 006      J.-F. Mahfouf (Spring 2003) -
   !              Add implicit boundary condition option for vert. diff.
   !@Object Update the prognostic variables
   !
   !@Arguments
   !          - Input/Output -
   ! TS
   ! T2
   ! WG
   ! W2       prognostic variables at time -
   ! WF
   ! WL
   ! WR
   ! WS
   ! ALPHAS
   ! RHOS
   ! RHOSNO   density of snow (kg/m3) for output only
   !          - Input -
   ! VMOD     module of the low-level wind
   ! CD       surface transfer coefficient for momentum
   ! RHOA     low-level air density
   ! HFLUX    surface flux of sensible heat
   ! EFLUX    water vapor flux
   !           - Output -
   ! TST
   ! T2T
   ! WGT
   ! W2T       prognostic variables at time +
   ! WFT
   ! WLT
   ! WRT
   ! WST
   ! ALPHAST
   ! RHOST
   ! BM        homogeneous boundary condition term in the
   !           diffusion equation for U and V
   ! ALFAT     inhomogeneous boundary term in the diffusion equation for Theta
   ! ALFAQ     inhomogeneous boundary term in the diffusion equation for Q
   ! FQ        surface momentum flux

   integer :: I

   do I=1,N

      TS(I)     = TST(I)
      T2(I)     = T2T(I)
      WG(I)     = WGT(I)
      W2(I)     = W2T(I)
      WF(I)     = WFT(I)
      WL(I)     = WLT(I)
      WR(I)     = WRT(I)
      WS(I)     = WST(I)
      ALPHAS(I) = ALPHAST(I)
      RHOS(I)   = RHOST(I)
      RHOSNO(I) = RHOST(I)*RAUW

      WG(I)     = max( WG(I)   , 0.001  )
      W2(I)     = max( W2(I)   , 0.001  )

   end do

   ! Feedback on the vertical diffusion

   do I=1,N
      BM(I)      =  VMOD(I)*CD(I)
      FQ(I)      =  RHOA(I)*CD(I)*VMOD(I)*VMOD(I)
      ALFAT(I)   =  -HFLUX(I) / (CPD*RHOA(I))
      ALFAQ(I)   =  -EFLUX(I)
      if (IMPFLX) then
         ALFAT(I)   =   ALFAT(I) - CTU(I)* TH(I)
         ALFAQ(I)   =   ALFAQ(I) - CTU(I)* HU(I)
      endif
   end do

   return
end subroutine update4
