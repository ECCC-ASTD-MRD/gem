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

      SUBROUTINE CONRES1(RI,GAMA,GAMAQ,FN,T,TVE,Q,QE,PS,HOL,SIGMA,SE,S, &
                            QC,N,M,NK)
      use tdpack
      use phy_options
      implicit none
!!!#include <arch_specific.hf>

      INTEGER N,M,NK,j,k,stat
      INTEGER, DIMENSION(N) :: klev
      INTEGER, EXTERNAL :: neark

      REAL RI(N,NK),GAMA(N,NK),GAMAQ(N,NK),FN(N,NK)
      REAL T(M,NK),TVE(N,NK),Q(M,NK),QE(N,NK)
      REAL PS(N),HOL(N)
      REAL SIGMA(n,NK),SE(n,NK),S(N,NK)
      REAL QC(N)
      REAL QFIX(N)

      REAL X,TVK,TVN,TVBAR,GAMAA,GAMAV,GAMAVS,FAC,DZ
      REAL BETA,GAMAS,AA,RGASD_OV_GRAV,GRAV_OV_CPD
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      REAL, dimension(N,NK) :: WORK
      REAL, dimension(N,NK) :: TE
      REAL, dimension(N,NK) :: QSAT
      REAL, dimension(N,NK) :: PRES
!
!Author
!        C .Girard RPN March 1993
!
!Revisions
! 001    R. Benoit (August 93) - Local sigma (Sigma, SE (2D))
! 002    C .Girard (March 96) - Clean-up and New shallow convection
! 003    M. Lepine (March 2003) -  CVMG... Replacements
! 004    A. Plante (May 2003) - IBM conversion
!           - replace CVMG* by if/else statements
!           - calls to vslog routine (from massvp4 library)
!           - calls to optimized routin mfoqst3
!           - constantes precomputations
! 005    B. Bilodeau (May 2005) - New comdeck fintern
! 006    A. Zadra (June 2014) - Replace use of Q(NK) by 
!            Q near a fixed pressure level above surface, 
!            to reduce sensitivity to position of lowest level
!
!Object
!        Parameterization of certain effects of shallow convection:
!        -Estimate a convective cloud fraction which will interact
!         with radiation schemes
!        -Modify the stability through Ri for all diffused variables
!         including wind
!        -Modify the equilibrium gradients of temperature and moisture
!
!
!Arguments
!
!          -Input/Output-
! RI       nombre de Richardson
! GAMA     Gradient at equilibrium for temperature
! GAMAQ    Gradient at equilibrium for moisture
! FN       Cloud fraction related to shallow convection
!
!          -Input-
! T        Temperature (M,NK)
! TVE      Virtual Temperature at level 'E' (N,NK)
! Q        Specific humidity(M,NK)
! QE       Specific humidity at level 'E' (N,NK)
! PS       surface pressure (N)
! HOL      Indicator of stability in the boundary layer (N)
!                 (Unstable when negative)
! SIGMA    Sigma level (n,NK)
! SE       Sigma level for 'E' (n,NK)
! N        Horizontal dimension
! M        1st dimension of T and Q in the calling program
! NK       vertical dimension
!
!          -Work-
! QC       specific humidity of cloud (N)

      RGASD_OV_GRAV=RGASD/GRAV
      GRAV_OV_CPD=GRAV/CPD

!     Find Q at nearest level 10mb above the surface
      stat = neark(se,ps,1000.,n,nk,klev)
 
      DO j=1,N
         QFIX(j) = Q(j,klev(j))
!
!           POSSIBILITE DE PRESENCE D'UN NUAGE CONVECTIF:
!             QC=QFIX    -AU-DESSUS D'UNE COUCHE LIMITE INSTABLE
!
         if (-HOL(j) .ge. 0.) then
            QC(j) = QFIX(j)
         else
            QC(j) = 0.
         endif
!
      END DO
!
!     Precomputations for optimisation on IBM
      DO k=1,NK
         DO j=1,N
            WORK(j,k)=log(SIGMA(j,NK)/SIGMA(j,k))
            TE(j,k)=FOTTV(TVE(j,k),QE(j,k))
            PRES(j,k)=SE(j,k)*PS(j)
         END DO
      END DO
!
      CALL mfoqst3(QSAT,TE,PRES,N,NK,N)
!
      DO k=NK-1,2,-1
         DO j=1,N
!
!           CALCUL DE TVK, TVN ET TVBAR
!
            TVK = FOTVT(T(j,k),Q(j,k))
            TVN = FOTVT(T(j,NK),Q(j,NK))
            TVBAR = 0.5*(TVK+TVN)
!
!           CALCUL APPROX. DE GAMAV
!
            DZ = RGASD_OV_GRAV*TVBAR*WORK(J,K)
            GAMAV = (TVK-TVN)/DZ + GRAV_OV_CPD
!
!           CALCUL APPROX. DE GAMAVS
!
            BETA=1.35E7*QSAT(j,k)/(TE(j,k)*TVE(j,k))
            GAMAS=GRAV_OV_CPD*(1.-6.46E-4*TE(j,k))*BETA/(1.+BETA)
!
            GAMAVS=(TVE(j,k)/TE(j,k))*GAMAS
!
!           CALCUL DU RAPPORT: -GAMAV/GAMAVS
!
            X = - GAMAV/GAMAVS
!
!           DOIT ETRE ENTRE -1. ET 0.
!
            X = MAX(-1.,MIN(0.,X))
!
!           CALCUL DE LA FRACTION NUAGEUSE: FN
!            (variante de la formule de Bjerknes)
!
            FN(j,k) = (1+X)/(fnnmod+X)
!
!           QC S'EVANOUIT LORSQUE FN=0
!
            if (FN(j,k) .eq. 0.)then
               QC(j) = 0.
            endif
!
!           FN S'EVANOUIT A SON TOUR SI QC.LT.QSAT
!
            if (QC(j).LT.QSAT(j,k))then
               FN(j,k) = 0.
            endif
!
!           CALCUL DE GAMA = AA x GAMAS
!
            AA=SQRT(FN(j,k))
            GAMAA=AA*GAMAS
!
!           MODIFICATION DE RI
!
            FAC = GRAV / ( TE(j,k) * (max(S(j,k),1.e-10)) )
            RI(j,k) =  RI(j,k) - FAC * GAMAA
!
            GAMA(j,k) = GAMA(j,k) + GAMAA
!
            GAMAQ(j,k) = GAMAQ(j,k)
!
         END DO
      END DO
!
      RETURN
      END
