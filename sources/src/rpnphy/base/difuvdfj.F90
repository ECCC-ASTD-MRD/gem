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

!/@*
SUBROUTINE DIFUVDFj1(TU, U, KU, GU, JNG, R, ALFA, BETA, S, SK, &
                           TAU, TYPE, F, NU, NR, N, NK)
   implicit none
#include <arch_specific.hf>
      INTEGER NU, NR, N, NK
      REAL TU(NU, NK), U(NU, NK), KU(NR, NK), GU(NR, NK), R(NR,NK)
      REAL JNG(NR, NK)
      REAL ALFA(N), BETA(N), S(n,NK), SK(n,NK), TAU, F
      INTEGER TYPE
!
!Author
!          R. Benoit (Mar 89)
!
!Revisions
! 001      R. Benoit (Aug 93) -Local sigma: s and sk become 2D
! 002      B. Bilodeau (Dec 94) - "IF" tests on integer
!          instead of character.
! 003      J. Mailhot (Sept 00) - Add type 4='EB'
! 004      A. PLante (June 2003) - IBM conversion
!             - calls to vrec routine (from massvp4 library)
! 005      J. Mailhot/L. Spacek (Dec 07) - Add type 5='ET' and cleanup
! 006      J. Mailhot (Aug 12) - Add argument and modifications (non-gradient flux,
!                                geometric averaging) and change name to difuvdfj1
!
!Object
!          to solve a vertical diffusion equation by finite
!          differences
!
!Arguments
!
!          - Output -
! TU       U tendency (D/DT U) due to the vertical diffusion and to
!          term R
!
!          - Input -
! U        variable to diffuse (U,V,T,Q,E)
! KU       diffusion coefficient
! GU       optional countergradient term
! JNG      optional non-gradient flux term
! R        optional inhomogeneous term
! ALFA     inhomogeneous term for the surface flux (for type 1='U', 2='UT' or 5='ET'or 6='ST')
!          surface boundary condition (for type 4='EB')
! BETA     homogeneous term for the surface flux
! S        sigma coordinates of full levels
! SK       sigma coordinates of diffusion coefficient
!          (or staggered variables) levels
! TAU      length of timestep
! TYPE     type of variable to diffuse (1='U',2='UT',3='E',4='EB' or 5='ET' or 6='ST')
! F        weighting factor for time 'N+1'
! NU       1st dimension of TU and U
! NR       1st dimension of KU, GU, JNG and R
! N        number of columns to process
! NK       vertical dimension
!
!Notes
!          D/DT U = D(U) + R
!          D(U) = D/DS ( J(U) + JNG )
!          J(U) = KU*(D/DS U + GU)
!          Limiting Conditions where S=ST: J=0(for 'U'/'ET'), D=0(for 'UT'
!          and ST=1)
!**        U=0(for 'E' and 'EB')
!          J=0(for 'E' and 'EB')
!          Limiting Conditions where S=SB: J=ALFA+BETA*U(S(NK))(for
!          'U'/'UT'/'ET'), J=0(for 'E'), U=ALFA(for 'EB')
!          ST = S(1)-1/2 (S(2)-S(1)) (except for 'TU')
!          SB = SK(NK) = 1.
!
!*@/
!
      INTEGER I, K, NKX
      REAL ST, SB, HM, HP, KUM, KUP, SCK1
      REAL, DIMENSION(N) :: HD
      REAL, DIMENSION(N,NK) :: VHM,VHP,A,B,C,D
      REAL(KIND=8), DIMENSION(N,NK) :: RHD,RHMD,RHPD
      LOGICAL :: SFCFLUX
      character(len=16) :: msg_S
      EXTERNAL DIFUVD1, DIFUVD2
!
      st(i)=s(i,1)-0.5*(s(i,2)-s(i,1))
      sb(i)=1.
!
      IF (TYPE.LE.2) THEN
         NKX=NK
         SCK1=1
         IF (TYPE.EQ.2) THEN
               SCK1=0
         ENDIF
      ELSE IF (TYPE.EQ.3 .OR. TYPE.EQ.4) THEN
         NKX=NK-1
      ELSE IF (TYPE.EQ.5 .OR. TYPE.EQ.6) THEN
         NKX=NK
      ELSE
         write(msg_S, *), type
         call physeterror('difuvdfj', 'Type inconnu: '//trim(msg_S))
         return
      ENDIF
!
! (1) CONSTRUIRE L'OPERATEUR TRIDIAGONAL DE DIFFUSION N=(A,B,C)
!                ET LE TERME CONTRE-GRADIENT (DANS D)
!
      IF (TYPE.LE.2) THEN
!
!     K=1
!
         HM=0
         DO 10 I=1,N
            HP=S(i,2)-S(i,1)
            HD(I)=SK(i,1)-ST(i)
            A(I,1)=0
            C(I,1)=SCK1*KU(I,1)/(HP*HD(I))
            B(I,1)=-A(I,1)-C(I,1)
10          D(I,1)=SCK1*(KU(I,1)*GU(I,1)+JNG(I,1))/HD(I)
!
!     K=2...NK-1
!
         DO K=2,NK-1,1
            DO I=1,N
!              THE FOLLOWING LHS ARE IN REAL
               VHM(I,K)=S(I,K)-S(I,K-1)
               VHP(I,K)=S(I,K+1)-S(I,K)
               HD(I)=SK(I,K)-SK(I,K-1)
!	       THE FOLLOWING LHS ARE IN REAL*8
               RHD(I,K)=HD(I)
               RHMD(I,K)=VHM(I,K)*HD(I)
               RHPD(I,K)=VHP(I,K)*HD(I)
            ENDDO
         ENDDO
         CALL VREC(RHD (1,2), RHD(1,2),N*(NK-2))
         CALL VREC(RHMD(1,2),RHMD(1,2),N*(NK-2))
         CALL VREC(RHPD(1,2),RHPD(1,2),N*(NK-2))
         DO K=2,NK-1,1
            DO I=1,N
               A(I,K)=KU(I,K-1)*RHMD(I,K)
               C(I,K)=KU(I,K)*RHPD(I,K)
               B(I,K)=-A(I,K)-C(I,K)
               D(I,K)=( KU(I,K)*GU(I,K)-KU(I,K-1)*GU(I,K-1) &
                             +JNG(I,K)-JNG(I,K-1) )*RHD(I,K)
            ENDDO
         ENDDO
!
!     K=NK
!
         HP=0
         DO 12 I=1,N
            HM=S(i,NK)-S(i,NK-1)
            HD(I)=SB(i)-SK(i,NK-1)
            A(I,NK)=KU(I,NK-1)/(HM*HD(I))
            C(I,NK)=0
            B(I,NK)=-A(I,NK)-C(I,NK)
12          D(I,NK)=(0-KU(I,NK-1)*GU(I,NK-1)-JNG(I,NK-1))/HD(I)
!
      ELSE IF (TYPE.EQ.3 .OR. TYPE.EQ.4 .OR. TYPE.EQ.5 .OR. TYPE.EQ.6) THEN
!
!     TYPE='E' or 'EB' or 'ET'
!
!     K=1
!
         DO 13 I=1,N
            HM=SK(i,1)-ST(i)
            HP=SK(i,2)-SK(i,1)
            HD(I)=S(i,2)-S(i,1)
!        Limiting Conditions at S=ST: U=0(for 'E' or 'EB')`
!**         KUM=0.5*KU(I,1)
!        Limiting Conditions at S=S(1): J=0(for 'E' or 'EB' or 'ET')
            KUM=0
            KUP=0.5*(KU(I,1)+KU(I,2))
            A(I,1)=KUM/(HM*HD(I))
            C(I,1)=KUP/(HP*HD(I))
            B(I,1)=-A(I,1)-C(I,1)
13          D(I,1)=(KUP*(GU(I,1)+GU(I,2))-KUM*GU(I,1) &
                        +(JNG(I,1)+JNG(I,2)) )/(2.*HD(I))
!
!     K=2...NKX-1
!
         DO K=2,NKX-1,1
            DO I=1,N
!              THE FOLLOWING LHS ARE IN REAL
               VHM(I,K)=SK(I,K)-SK(I,K-1)
               VHP(I,K)=SK(I,K+1)-SK(I,K)
               HD(I)=S(I,K+1)-S(I,K)
            ENDDO
            IF (K==NK-1 .AND. TYPE==6) THEN !VIRTUAL LEVEL FOR TYPE='ST'
               DO I=1,N
                  HD(I) = 0.5*(SK(I,K+1)+SK(I,K))-S(I,K)
               ENDDO
            ENDIF
            DO I=1,N
!	       THE FOLLOWING LHS ARE IN REAL*8
               RHD(I,K)=HD(I)
               RHMD(I,K)=VHM(I,K)*HD(I)
               RHPD(I,K)=VHP(I,K)*HD(I)
            ENDDO
         ENDDO
         CALL VREC( RHD(1,2), RHD(1,2),N*(NKX-2))
         CALL VREC(RHMD(1,2),RHMD(1,2),N*(NKX-2))
         CALL VREC(RHPD(1,2),RHPD(1,2),N*(NKX-2))
         DO K=2,NKX-1,1
            DO I=1,N
               KUM=0.5*(KU(I,K-1)+KU(I,K))
               KUP=0.5*(KU(I,K+1)+KU(I,K))
               A(I,K)=KUM*RHMD(I,K)
               C(I,K)=KUP*RHPD(I,K)
               B(I,K)=-A(I,K)-C(I,K)
               D(I,K)=.5*(KUP*(GU(I,K)+GU(I,K+1)) &
                      -KUM*(GU(I,K-1)+GU(I,K)) &
                       +(JNG(I,K)+JNG(I,K+1)) &
                       -(JNG(I,K-1)+JNG(I,K)))*RHD(I,K)
            ENDDO
         ENDDO
!
!     K=NKX
!
        IF (TYPE.EQ.3 .OR. TYPE.EQ.5 .OR. TYPE.EQ.6) THEN
!
!       TYPE='E' or 'ET' or 'ST'
!
           IF (TYPE.EQ.6) THEN !virtual level for TYPE='ST'
              DO I=1,N
                 HD(I)=SB(i)-0.5*(SK(i,NKX)+SK(i,NKX-1))
              ENDDO
           ELSE
              DO I=1,N
                 HD(I)=SB(i)-S(i,NKX)
              ENDDO
           ENDIF
           DO I=1,N
              HM=SK(i,NKX)-SK(i,NKX-1)
              KUM=0.5*(KU(I,NKX)+KU(I,NKX-1))
              KUP=0
              A(I,NKX)=KUM/(HM*HD(I))
              C(I,NKX)=0
              B(I,NKX)=-A(I,NKX)-C(I,NKX)
              D(I,NKX)=(0-KUM*(GU(I,NKX)+GU(I,NKX-1)) &
                           -(JNG(I,NKX)+JNG(I,NKX-1)) )/(2.*HD(I))
           ENDDO
!
        ELSE IF (TYPE.EQ.4) THEN
!
!       TYPE='EB'
!
           DO I=1,N
              HM=SK(i,NK-1)-SK(i,NK-2)
              HP=SB(i)-SK(i,NK-1)
              HD(I)=S(i,NK)-S(i,NK-1)
              KUM=0.5*(KU(I,NK-1)+KU(I,NK-2))
              KUP=0.5*(KU(I,NK)+KU(I,NK-1))
              A(I,NKX)=KUM/(HM*HD(I))
              B(I,NKX)=-A(I,NKX) -KUP/(HP*HD(I))
              C(I,NKX)=0
              D(I,NKX)=(KUP*(GU(I,NK)+GU(I,NK-1)) &
                       -KUM*(GU(I,NK-1)+GU(I,NK-2)) &
                       +(JNG(I,NK)+JNG(I,NK-1)) &
                       -(JNG(I,NK-1)*JNG(I,NK-2)))/(2.*HD(I)) &
                       +KUP*ALFA(I)/(HD(I)*HP)
           ENDDO
!
        ENDIF
!
      ENDIF
!
!
! (2) CALCULER LE COTE DROIT D=TAU*(N(U)+R+D/DS(KU*GU+JNG))
!
      CALL DIFUVD1 (D, 1., A, B, C, U, D, N, NU, NKX)
      DO 20 K=1,NKX
         DO 20 I=1,N
20       D(I,K)=TAU*(D(I,K)+R(I,K))
!
! (3) CALCULER OPERATEUR DU COTE GAUCHE
!
      DO 30 K=1,NKX
         DO 30 I=1,N
            A(I,K)= -F*TAU*A(I,K)
            B(I,K)=1-F*TAU*B(I,K)
30          C(I,K)= -F*TAU*C(I,K)
!
! (4) AJOUTER TERME DE FLUX DE SURFACE
!
      SFCFLUX = .TRUE.
      SELECT CASE (TYPE)
      CASE (:2) !TYPE='U'/'UT'
         DO I=1,N
            HD(I) = SB(i) - SK(i,NK-1)
         ENDDO
      CASE (5) !TYPE='ET'
         DO I=1,N
            HD(I) = SB(i) - S(i,NKX)
         ENDDO
      CASE (6) !TYPE='ST'
         DO I=1,N
            HD(I) = SB(i) - 0.5*(SK(i,NKX)+SK(i,NKX-1))
         ENDDO
      CASE DEFAULT
         SFCFLUX = .FALSE.
      END SELECT
      IF (SFCFLUX) THEN
         DO I=1,N
            B(I,NKX)=B(I,NKX)-TAU*BETA(I)/HD(I)
            D(I,NKX)=D(I,NKX)+(ALFA(I)+BETA(I)*U(I,NKX))*TAU/HD(I)
         ENDDO
      ENDIF
!
! (5) RESOUDRE SYSTEME TRIDIAGONAL [A,B,C] X = D. METTRE X DANS TU.
!
      CALL DIFUVD2 (TU, A, B, C, D, D, NU, N, NKX)
!
! (6) OBTENIR TENDANCE
!
      DO 60 K=1,NKX
         DO 60 I=1,N
60       TU(I,K)=TU(I,K)/TAU
!     K=NKX+1..NK
      DO 70 K=NKX+1,NK
         DO 70 I=1,N
70       TU(I,K)=0
!
      RETURN
      END
