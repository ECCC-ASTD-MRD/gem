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
!**S/P TWIND
!
      SUBROUTINE TWIND( WGE, WGMAX, WGMIN, SDTSWS, SDTSWD, TVE, &
                        EN, U, V, UD, VD, SE, ZE, H, USTAR, WSTAR, N, NK)
!
!
      implicit none
#include <arch_specific.hf>
!
      INTEGER N,NK
!
      REAL WGE(N), WGMAX(N), WGMIN(N), SDTSWS(N), SDTSWD(N)
      REAL TVE(N,NK), EN(N,NK), U(N,NK), V(N,NK), ZE(N,NK), SE(N,NK)
      REAL UD(N), VD(N)
      REAL H(N), USTAR(N), WSTAR(N)
!
!
!
!Author
!          J. Mailhot (September 2008)
!
!Revision
!
!Object
!           Interface routine for turbulent wind computations:
!           standard deviation of turbulent surface wind (speed and direction)
!           and wind gust estimates (based on method of Brasseur 2001)
!
!Arguments
!                        -Output-
!
! WGE       wind gust estimate
! WGMAX     wind gust maximum (upper bound on gust estimate)
! WGMIN     wind gust minimum (lower bound on gust estimate)
! SDTSWS    standard deviation of turbulent surface wind speed
! SDTSWD        "        "           "         "      "  direction (in degrees)
!
!                         -Input-
!
! TVE       virtual temperature (on 'E' levels)
! EN        turbulent kinetic energy (on 'E' levels)
! U         east-west component of wind (on full levels)
! V         north-south component of wind (on full levels)
! UD        east-west component of wind on diagnostic level (10m)
! VD        north-south component of wind on diagnostic level (10m)
! SE        sigma levels (on 'E' levels)
! ZE        height of the sigma levels (on 'E' levels)
! H         height of the boundary layer
! USTAR     friction velocity
! WSTAR     convective velocity scale
! N         horizontal dimension
! NK        vertical dimension
!
!Notes
!           Standard deviation calculations are based on values
!            from Wyngaard and Cote (1974, BLM 7, 289-308).
!           Wind gust estimates are based on Brasseur (2001, MWR 129, 5-25).
!
      INTEGER J, K
      REAL RAD2DEG
!
      REAL CU,CV,CW
      SAVE CU,CV,CW
      DATA CU, CV, CW  / 4.0, 1.75, 0.2 /
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      REAL, dimension(N      ) :: WINSPD
      REAL, dimension(N      ) :: SDX
      REAL, dimension(N,NK   ) :: THVE
      REAL, dimension(N,NK   ) :: WORK
      REAL, dimension(N,2    ) :: W1
!
!***********************************************************************
!
#include "tdpack_const.hf"
!
!
      RAD2DEG = 180.0/PI
!                                Standard deviation
!                                 (speed and direction in degrees)
      DO J=1,N
        W1(J,1) = USTAR(J)**2
        W1(J,2) = CW*WSTAR(J)**2
        SDTSWS(J) = CU*W1(J,1) + W1(J,2)
        SDX(J) = CV*W1(J,1) + W1(J,2)
        WINSPD(J) = UD(J)**2 + VD(J)**2
      END DO
!
      CALL VSSQRT (SDTSWS, SDTSWS, N)
      CALL VSSQRT (SDX   , SDX   , N)
      CALL VSSQRT (WINSPD, WINSPD, N)
      CALL VSATAN2 (SDTSWD, SDX, WINSPD, N)
!
      DO J=1,N
        SDTSWD(J) = RAD2DEG*ABS( SDTSWD(J) )
      END DO
!
!
!                                Virtual potential temperature (THVE)
      CALL VSPOWN1 (WORK,SE,-CAPPA,N*NK)
!
      DO K=1,NK
      DO J=1,N
        THVE(J,K) = TVE(J,K)*WORK(J,K)
      END DO
      END DO
!
!                                Wind gust estimates
!
      CALL WINDGUST( WGE, WGMAX, WGMIN, THVE, EN, U, V, UD, VD, &
                     ZE, H, N, NK)
!
!
      RETURN
      END
