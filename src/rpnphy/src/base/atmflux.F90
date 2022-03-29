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

      SUBROUTINE ATMFLUX3(R,VAR,TVAR,COEFK,GAMMA,FLUX_NG, &
                         ALPHA,BETA,PS,T,Q,TAU,SG,SELOC, &
                         VCOEF,C,D,TYPVAR,N,NK,TRNCH)
      use tdpack
      use series_mod, only: series_xst
      implicit none
!!!#include <arch_specific.hf>

      INTEGER N,NK,TRNCH
      REAL R(N,NK+1),VAR(N,NK),TVAR(N,NK)
      REAL GAMMA(N,NK),FLUX_NG(N,NK)
      REAL ALPHA(N), BETA(N), PS(N), TAU
      REAL, DIMENSION(N,NK), TARGET :: SG,SELOC,T,Q,COEFK
      REAL C(N,NK), D(N,NK)
      REAL VCOEF(N,NK,2)
      INTEGER TYPVAR

!@Author S. Belair (February 1996)
!@Revision
! 001      S. Belair (Oct 1996) - Include the countergradient term
! 002      L. Spacek (Dec 2007) - add "vertical staggering" option
!                                 change the name to atmflux1
! 003      J. Mailhot/A. Lock (Aug 2012) Revisions for non-local scaling
!                      Change calling sequence and rename ATMFLUX2
!@Object Calculate the atmospheric fluxes for heat, vapour, and momentum.
!@Arguments
!                        -Output-
! R         Resulting atmospheric flux
!           For U and V:  rho (w'v') = rho Km dv/ds
!           For Theta:    rho cp (w'theta')
!           For qv:       rho L (w'qv')
!                         -Input-
! VAR       Variable at t (U,V,Theta, or qv)
! TVAR      Time tendency of the variable
! COEFK     Vertical diffusion coefficient in sigma form
! GAMMA     Counter-gradient term (non-zero only for theta and hu)
! FLUX_NG   Non-gradient flux
! ALPHA     inhomogeneous bottom boundary condition
! BETA      homogeneous bottom boundary condition
! PS        Surface pressure
! T         Temperature
! Q         Specific humidity
! TAU       Timestep
! SG        Sigma levels
! SELOC     Staggered sigma levels
! C         Work field
! D         Work field
! TYPVAR    Type of variable to treat
!           '0'  --->  U
!           '1'  --->  V
!           '2'  --->  Q
!           '3'  --->  Theta
! N,NK      Horizontal and vertical dimensions
!*********************************************************************

      real, dimension(n,nk), target :: qmom,tmom,kmom
      real, dimension(n,nk+1) :: sigstag
      real, dimension(:,:), pointer :: tstag,qstag,kcoef
      INTEGER J,K, nksfc
      REAL A
      logical :: momentum_levels

!                                  Calculate staggered Q and T
!
      call vint_thermo2mom(tmom, t, vcoef, n, nk)
      call vint_thermo2mom(qmom, q, vcoef, n, nk)
      call vint_thermo2mom(kmom, coefk, vcoef, n, nk)
      if (typvar < 2) then
         momentum_levels = .true.
         nksfc = nk
         sigstag(:,1:nk) = seloc
         sigstag(:,nk+1) = sigstag(:,nk)
         tstag => t
         qstag => q
         kcoef => coefk
      else
         momentum_levels = .false.
         nksfc = nk+1
         sigstag(:,1:nk) = sg
         sigstag(:,nksfc) = 1.
         tstag => tmom
         qstag => qmom
         kcoef => kmom
      endif
!
!
!                                  Calculate VAR(t+1) (put into C)
      DO K=1,NK
        DO J=1,N
          C(J,K) = VAR(J,K) + TVAR(J,K)*TAU
        END DO
      END DO
!
!                                  Vertical derivative of VAR(t+1)
!                                  (Values on staggered levels).
!                                  Put into D.
!
      if (momentum_levels) then
         do k=1,nk-1
            d(:,k) = (c(:,k+1)-c(:,k)) / (sg(:,k+1)-sg(:,k))
         enddo
      else
         do k=2,nk
            d(:,k) = (c(:,k)-c(:,k-1)) / (seloc(:,k)-seloc(:,k-1))
         enddo
         d(:,1) = d(:,2)
      endif

!
!                                  Product K * dVAR/ds (into R).
!                                  Values on staggered levels.
!                                  CAREFUL:  We must divide by
!                                  a factor A = g sigma / R T
!                                  (on staggered levels again)
!
!
      DO K=1,nksfc-1
        DO J=1,N
          A = ( RGASD*tstag(J,K) ) / ( GRAV*sigstag(J,K) )
          R(J,K) = ( kcoef(J,K) * D(J,K) + FLUX_NG(J,K) ) * A
!
!
!
!                                  Add the countergradient part of the
!                                  flux:
!                                         = K * gamma / A**2
!
          IF (TYPVAR.EQ.2) THEN
             R(J,K) = R(J,K) + kcoef(J,K) * GAMMA(J,K) * A
             GAMMA(J,K) = GAMMA(J,K) / A
          END IF
          IF (TYPVAR.EQ.3) &
             R(J,K) = R(J,K) + kcoef(J,K) * &
                      GAMMA(J,K) * A * A
        END DO
      END DO
!
!                                  The same product for the
!                                  lowest level NK (actually,
!                                  NK is one less than the number
!                                  of model levels).
!                                  Km*dVAR/ds = alfa + beta*u(t+1)
!                                  Note:  Need to multiply by A also.
!                                  Assume FLUX_NG(J,NK)=0
!
      DO J=1,N
        A = ( RGASD*tstag(J,NK) ) / ( GRAV*sigstag(J,NK) )
        R(J,nksfc)   = (ALPHA(J) + BETA(J)*C(J,NK) ) * A
        IF (TYPVAR.EQ.2) &
              GAMMA(J,NK) = GAMMA(J,NK) / A
      END DO
!
!                                  Multiply by the air density
!                                  rho = p / (R Tv) = s ps / RTv
!                                  CAREFUL:  the calculations must
!                                            be on the SELOC levels,
!                                            but the variables T and
!                                            Q are defined on SG levels.
!
!                                  For all the levels except NK, NK+1
      DO K=1,nksfc-1
        DO J=1,N
          R(J,K) = sigstag(J,K)*PS(J) / &
                   ( RGASD * FOTVT(tstag(J,K),qstag(J,K)) ) &
                 * R(J,K)
        END DO
      END DO
!
!                                  For NK and NK+1
!
      DO J=1,N
          R(J,nksfc) = sigstag(J,NK)*PS(J) / &
                   ( RGASD * FOTVT(tstag(J,NK),qstag(J,NK)) ) &
                 *   R(J,nksfc)
      END DO
      if (momentum_levels) R(:,nk+1) = R(:,nk)



!
!
!                                  At this point, R contains
!                                  rho (w'var').  For the temperature
!                                  and the specific humidity, however,
!                                  we need to multiply by cp and L,
!                                  respectively.
!
!                                  Furthermore, a factor (T/theta)
!                                  is kept for the fluxes of sensible
!                                  heat (because we want w'T', and
!                                  not w'theta')
!
      DO K=1,nksfc
        DO J=1,N
          IF (TYPVAR.EQ.2) R(J,K) = CHLC * R(J,K)
          IF (TYPVAR.EQ.3) &
            R(J,K) = CPD  * R(J,K) * &
                   ( sigstag(J,K)*PS(J) / 100000. )**0.286
        END DO
      END DO
!
      if (momentum_levels) R(:,nk+1) = R(:,nk)

      ! Time series for the fluxes
      if (typvar == 0) call series_xst(r, 'f5', trnch)
      if (typvar == 1) call series_xst(r, 'f6', trnch)
      if (typvar == 2) call series_xst(r, 'f3', trnch)
      if (typvar == 3) call series_xst(r, 'f4', trnch)

      RETURN
      END
