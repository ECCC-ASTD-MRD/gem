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
!------------------

subroutine WINDGUST(WGE, WGMAX, WGMIN, THVE, EN, U, V, UD, VD, &
     ZE, H, N, NK)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>

      integer, intent(in) :: N,NK

      real WGE(N), WGMAX(N), WGMIN(N)
      real THVE(N,NK), EN(N,NK), U(N,NK), V(N,NK), ZE(N,NK)
      real UD(N), VD(N), H(N)
!
!
!
!Author
!          J. Mailhot (September 2008)
!
!Revision
!
!Object
!           Calculates an estimate of surface wind gusts
!           (based on method of Brasseur 2001)
!
!Arguments
!                        -Output-
!
! WGE       wind gust estimate
! WGMAX     wind gust maximum (upper bound on gust estimate)
! WGMIN     wind gust minimum (lower bound on gust estimate)
!
!                         -Input-
!
! THVE      virtual potential temperature (on 'E' levels)
! EN        turbulent kinetic energy (on 'E' levels)
! U         east-west component of wind (on full levels)
! V         north-south component of wind (on full levels)
! UD        east-west component of wind on diagnostic level (10m)
! VD        north-south component of wind on diagnostic levels (10m)
! ZE        height of the sigma levels (on 'E' levels)
! H         height of the boundary layer
! N         horizontal dimension
! NK        vertical dimension
!
!Notes
!           This is based on Brasseur (2001, MWR 129, 5-25).
!
      integer J, K, KI
!
      real RATIO
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      integer, dimension(N,NK   ) :: KIK
      integer, dimension(N,NK   ) :: KIK1
!
      real, dimension(N      ) :: WINSPD
      real, dimension(N,NK   ) :: ENLOCAL
      real, dimension(N,NK   ) :: ENSUM
      real, dimension(N,NK   ) :: BUOYSUM
!
!***********************************************************************

      RATIO = 2.5/11.0
!
      do K=1,NK
      do J=1,N
        ENLOCAL(J,K) = min( EN(J,K), 4. )
        if (ZE(J,K).gt.H(J))  ENLOCAL(J,K) = 0.
        ENSUM(J,K)=0.0
        BUOYSUM(J,K) = 0.0
      end do
      end do
!
!
!                      --------- computes the layer-averaged turbulent kinetic energy
!
      do K=NK-1,1,-1
      do J=1,N
        ENSUM(J,K) = ENSUM(J,K+1) + &
             0.5*( ZE(J,K) - ZE(J,K+1) ) * &
             ( ENLOCAL(J,K) + ENLOCAL(J,K+1) )
      end do
      end do
!
      do K=1,NK-1
      do J=1,N
        ENSUM(J,K) = ENSUM(J,K) /( ZE(J,K) - ZE(J,NK) )
      end do
      end do
!
!
!                      --------- restricts to the boundary layer depth
!
      do K=1,NK
      do J=1,N
        if (ZE(J,K) .gt. H(J))   ENSUM(J,K)=0.0
      end do
      end do
!
!
!
!                      --------- computes the buoyant energy
!
      do KI=1,NK-1
      do K=KI,NK-1
      do J=1,N
        BUOYSUM(J,KI) = BUOYSUM(J,KI) + &
             ( GRAV/(THVE(J,K) + THVE(J,K+1)) )* &
             ( ZE(J,K) - ZE(J,K+1) ) * &
             ( 2.*THVE(J,KI) - THVE(J,K) - THVE(J,K+1) )
      end do
      end do
      end do
!
!
!                      --------- computes the downward parcel displacement
!
      do K=1,NK-1
      do J=1,N
        KIK(J,K) = 1
        KIK1(J,K) = 1
        if (ZE(J,K) .le. H(J)) then
          if (ENSUM(J,K).ge.BUOYSUM(J,K)) &
                                           KIK(J,K) = NK
          if (RATIO*ENLOCAL(J,K).ge.BUOYSUM(J,K)) &
                                           KIK1(J,K) = NK
        endif
      end do
      end do
!
!
!
!                      --------- computes the wind gusts
!
      do J=1,N
        WINSPD(J) = sqrt( UD(J)*UD(J) + VD(J)*VD(J) )
        WGE(J) = WINSPD(J)
        WGMAX(J) = WINSPD(J)
        WGMIN(J) = WINSPD(J)
      end do
!
!
      do K=1,NK-1
      do J=1,N
        WINSPD(J) = 0.5*sqrt( (U(J,K)+U(J,K+1)) &
                             *(U(J,K)+U(J,K+1)) &
                            + (V(J,K)+V(J,K+1)) &
                             *(V(J,K)+V(J,K+1)) )
        if (KIK(J,K).eq.NK) WGE(J) = max( WGE(J) , WINSPD(J) )
        if (KIK1(J,K).eq.NK) WGMIN(J) = max( WGMIN(J) , WINSPD(J) )
        if (ZE(J,K).le.H(J)) WGMAX(J) = max( WGMAX(J) , WINSPD(J) )
      end do
      end do
!
      do J=1,N
        WINSPD(J) = sqrt( UD(J)*UD(J) + VD(J)*VD(J) )
        WGMIN(J) = max( WGMIN(J) , WINSPD(J) )
        WGE  (J) = max( WGE(J)   , WGMIN(J) )
        WGMAX(J) = max( WGMAX(J) , WGE(J) )
      end do
!
!
      return
      end
