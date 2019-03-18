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
!**S/P WINDGUST
!
      SUBROUTINE WINDGUST( WGE, WGMAX, WGMIN, THVE, EN, U, V, UD, VD, &
                           ZE, H, N, NK)
!
!
      implicit none
#include <arch_specific.hf>
!
      INTEGER N,NK
!
      REAL WGE(N), WGMAX(N), WGMIN(N)
      REAL THVE(N,NK), EN(N,NK), U(N,NK), V(N,NK), ZE(N,NK)
      REAL UD(N), VD(N), H(N)
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
      INTEGER J, K, KI
!
      REAL RATIO
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      INTEGER, dimension(N,NK   ) :: KIK
      INTEGER, dimension(N,NK   ) :: KIK1
!
      REAL, dimension(N      ) :: WINSPD
      REAL, dimension(N,NK   ) :: ENLOCAL
      REAL, dimension(N,NK   ) :: ENSUM
      REAL, dimension(N,NK   ) :: BUOYSUM
!
!***********************************************************************
!
#include "tdpack_const.hf"
!
!
!
      RATIO = 2.5/11.0
!
      DO K=1,NK
      DO J=1,N
        ENLOCAL(J,K) = MIN( EN(J,K), 4. )
        IF (ZE(J,K).GT.H(J))  ENLOCAL(J,K) = 0.
        ENSUM(J,K)=0.0
        BUOYSUM(J,K) = 0.0
      END DO
      END DO
!
!
!                      --------- computes the layer-averaged turbulent kinetic energy
!
      DO K=NK-1,1,-1
      DO J=1,N
        ENSUM(J,K) = ENSUM(J,K+1) + &
             0.5*( ZE(J,K) - ZE(J,K+1) ) * &
             ( ENLOCAL(J,K) + ENLOCAL(J,K+1) )
      END DO
      END DO
!
      DO K=1,NK-1
      DO J=1,N
        ENSUM(J,K) = ENSUM(J,K) /( ZE(J,K) - ZE(J,NK) )
      END DO
      END DO
!
!
!                      --------- restricts to the boundary layer depth
!
      DO K=1,NK
      DO J=1,N
        IF (ZE(J,K) .GT. H(J))   ENSUM(J,K)=0.0
      END DO
      END DO
!
!
!
!                      --------- computes the buoyant energy
!
      DO KI=1,NK-1
      DO K=KI,NK-1
      DO J=1,N
        BUOYSUM(J,KI) = BUOYSUM(J,KI) + &
             ( GRAV/(THVE(J,K) + THVE(J,K+1)) )* &
             ( ZE(J,K) - ZE(J,K+1) ) * &
             ( 2.*THVE(J,KI) - THVE(J,K) - THVE(J,K+1) )
      END DO
      END DO
      END DO
!
!
!                      --------- computes the downward parcel displacement
!
      DO K=1,NK-1
      DO J=1,N
        KIK(J,K) = 1
        KIK1(J,K) = 1
        IF (ZE(J,K) .LE. H(J)) THEN
          IF (ENSUM(J,K).GE.BUOYSUM(J,K)) &
                                           KIK(J,K) = NK
          IF (RATIO*ENLOCAL(J,K).GE.BUOYSUM(J,K)) &
                                           KIK1(J,K) = NK
        ENDIF
      END DO
      END DO
!
!
!
!                      --------- computes the wind gusts
!
      DO J=1,N
        WINSPD(J) = SQRT( UD(J)*UD(J) + VD(J)*VD(J) )
        WGE(J) = WINSPD(J)
        WGMAX(J) = WINSPD(J)
        WGMIN(J) = WINSPD(J)
      END DO
!
!
      DO K=1,NK-1
      DO J=1,N
        WINSPD(J) = 0.5*SQRT( (U(J,K)+U(J,K+1)) &
                             *(U(J,K)+U(J,K+1)) &
                            + (V(J,K)+V(J,K+1)) &
                             *(V(J,K)+V(J,K+1)) )
        IF (KIK(J,K).EQ.NK) WGE(J) = MAX( WGE(J) , WINSPD(J) )
        IF (KIK1(J,K).EQ.NK) WGMIN(J) = MAX( WGMIN(J) , WINSPD(J) )
        IF (ZE(J,K).LE.H(J)) WGMAX(J) = MAX( WGMAX(J) , WINSPD(J) )
      END DO
      END DO
!
      DO J=1,N
        WINSPD(J) = SQRT( UD(J)*UD(J) + VD(J)*VD(J) )
        WGMIN(J) = MAX( WGMIN(J) , WINSPD(J) )
        WGE  (J) = MAX( WGE(J)   , WGMIN(J) )
        WGMAX(J) = MAX( WGMAX(J) , WGE(J) )
      END DO
!
!
      RETURN
      END
