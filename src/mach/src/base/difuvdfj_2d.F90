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
!!if_on
 subroutine DIFUVDFj_2D(TU, U, KU, GU, JNG, R, ALFA, BETA, S, SK, &
                        BND, TAU, type, F, NU, NR, N, NK)
!!if_off
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
!!if_on
   integer, intent(in)  :: NU, NR, N, NK
   real,    intent(out) :: TU(NU, NK)
   real,    intent(in)  :: U(NU, NK), KU(NR, NK), GU(NR, NK), R(NR,NK)
   real,    intent(in)  :: JNG(NR, NK), BND(NR, NK)
   real,    intent(in)  :: ALFA(N), BETA(N), S(n,NK), SK(n,NK), TAU, F
   integer, intent(in)  :: type
!!if_off

   !@Author R. Benoit (Mar 89)
   !@Revisions
   ! 001      R. Benoit (Aug 93) -Local sigma: s and sk become 2D
   ! 002      B. Bilodeau (Dec 94) - "IF" tests on integer
   !          instead of character.
   ! 003      J. Mailhot (Sept 00) - Add type 4='EB'
   ! 004      A. PLante (June 2003) - IBM conversion
   !             - calls to vrec routine (from massvp4 library)
   ! 005      J. Mailhot/L. Spacek (Dec 07) - Add type 5='ET' and cleanup
   ! 006      J. Mailhot (Aug 12) - Add argument and modifications (non-gradient flux,
   !                                geometric averaging) and change name to difuvdfj1
   !@Object
   !          to solve a vertical diffusion equation by finite
   !          differences
   !@Arguments
   !          - Output -
   ! TU       U tendency (D/DT U) due to the vertical diffusion and to
   !          term R
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
   !@Notes
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
   !*@/

   integer I, K, NKX
   real ST, SB, HM, HP, KUM, KUP, SCK1
   real, dimension(N) :: HD
   real, dimension(N,NK) :: VHM,VHP,A,B,C,D,HD2
   real(KIND=8), dimension(N,NK) :: RHD,RHMD,RHPD
   logical :: SFCFLUX
   character(len=16) :: msg_S
   external DIFUVD1, DIFUVD2
   external physeterror, VREC

   st(i)=s(i,1)-0.5*(s(i,2)-s(i,1))
   sb(i)=1.

   if (type.le.2) then
      NKX=NK
      SCK1=1
      if (type.eq.2) then
         SCK1=0
      endif
   else if (type.eq.3 .or. type.eq.4) then
      NKX=NK-1
   else if (type.eq.5 .or. type.eq.6) then
      NKX=NK
   else
      write(msg_S, *) type
      call physeterror('difuvdfj', 'Type inconnu: '//trim(msg_S))
      return
   endif

   ! (1) CONSTRUIRE L'OPERATEUR TRIDIAGONAL DE DIFFUSION N=(A,B,C)
   !                ET LE TERME CONTRE-GRADIENT (DANS D)

   if (type.le.2) then

      !     K=1

      HM=0
      do I=1,N
         HP=S(i,2)-S(i,1)
         HD(I)=SK(i,1)-ST(i)
         A(I,1)=0
         C(I,1)=SCK1*KU(I,1)/(HP*HD(I))
         B(I,1)=-A(I,1)-C(I,1)
         D(I,1)=SCK1*(KU(I,1)*GU(I,1)+JNG(I,1))/HD(I)
      enddo

      !     K=2...NK-1

      do K=2,NK-1,1
         do I=1,N
            !              THE FOLLOWING LHS ARE IN REAL
            VHM(I,K)=S(I,K)-S(I,K-1)
            VHP(I,K)=S(I,K+1)-S(I,K)
            HD(I)=SK(I,K)-SK(I,K-1)
            !           THE FOLLOWING LHS ARE IN real(REAL64)
            RHD(I,K)=1.d0/HD(I)
            RHMD(I,K)=1.d0/(VHM(I,K)*HD(I))
            RHPD(I,K)=1.d0/(VHP(I,K)*HD(I))

            A(I,K)=KU(I,K-1)*RHMD(I,K)
            C(I,K)=KU(I,K)*RHPD(I,K)
            B(I,K)=-A(I,K)-C(I,K)
            D(I,K)=( KU(I,K)*GU(I,K)-KU(I,K-1)*GU(I,K-1) &
                 +JNG(I,K)-JNG(I,K-1) )*RHD(I,K)
         enddo
      enddo

      !     K=NK

      HP=0
      do I=1,N
         HM=S(i,NK)-S(i,NK-1)
         HD(I)=SB(i)-SK(i,NK-1)
         A(I,NK)=KU(I,NK-1)/(HM*HD(I))
         C(I,NK)=0
         B(I,NK)=-A(I,NK)-C(I,NK)
         D(I,NK)=(0-KU(I,NK-1)*GU(I,NK-1)-JNG(I,NK-1))/HD(I)
      enddo

   else if (type.eq.3 .or. type.eq.4 .or. type.eq.5 .or. type.eq.6) then

      !     TYPE='E' or 'EB' or 'ET'

      !     K=1

      do I=1,N
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
         D(I,1)=(KUP*(GU(I,1)+GU(I,2))-KUM*GU(I,1) &
              +(JNG(I,1)+JNG(I,2)) )/(2.*HD(I))
      enddo

      !     K=2...NKX-1

      do K=2,NKX-1,1
         do I=1,N
            !              THE FOLLOWING LHS ARE IN REAL
            VHM(I,K)=SK(I,K)-SK(I,K-1)
            VHP(I,K)=SK(I,K+1)-SK(I,K)
            HD(I)=S(I,K+1)-S(I,K)
         enddo
         if (K==NK-1 .and. type==6) then !VIRTUAL LEVEL FOR TYPE='ST'
            do I=1,N
               HD(I) = 0.5*(SK(I,K+1)+SK(I,K))-S(I,K)
            enddo
         endif
         do I=1,N
            RHD(I,K)=1.d0/HD(I)
            RHMD(I,K)=1.d0/(VHM(I,K)*HD(I))
            RHPD(I,K)=1.d0/(VHP(I,K)*HD(I))

            KUM=0.5*(KU(I,K-1)+KU(I,K))
            KUP=0.5*(KU(I,K+1)+KU(I,K))
            A(I,K)=KUM*RHMD(I,K)
            C(I,K)=KUP*RHPD(I,K)
            B(I,K)=-A(I,K)-C(I,K)
            D(I,K)=.5*(KUP*(GU(I,K)+GU(I,K+1)) &
                 -KUM*(GU(I,K-1)+GU(I,K)) &
                 +(JNG(I,K)+JNG(I,K+1)) &
                 -(JNG(I,K-1)+JNG(I,K)))*RHD(I,K)
         enddo
      enddo

      !     K=NKX

      if (type.eq.3 .or. type.eq.5 .or. type.eq.6) then

         !       TYPE='E' or 'ET' or 'ST'

         if (type.eq.6) then !virtual level for TYPE='ST'
            do I=1,N
               HD(I)=SB(i)-0.5*(SK(i,NKX)+SK(i,NKX-1))
            enddo
         else
            do I=1,N
               HD(I)=SB(i)-S(i,NKX)
            enddo
         endif
         do I=1,N
            HM=SK(i,NKX)-SK(i,NKX-1)
            KUM=0.5*(KU(I,NKX)+KU(I,NKX-1))
            KUP=0
            A(I,NKX)=KUM/(HM*HD(I))
            C(I,NKX)=0
            B(I,NKX)=-A(I,NKX)-C(I,NKX)
            D(I,NKX)=(0-KUM*(GU(I,NKX)+GU(I,NKX-1)) &
                 -(JNG(I,NKX)+JNG(I,NKX-1)) )/(2.*HD(I))
         enddo

      else if (type.eq.4) then

         !       TYPE='EB'

         do I=1,N
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
         enddo

      endif

   endif


   ! (2) CALCULER LE COTE DROIT D=TAU*(N(U)+R+D/DS(KU*GU+JNG))

   call DIFUVD1(D, 1., A, B, C, U, D, N, NU, NKX)
   do K=1,NKX
      do I=1,N
         D(I,K)=TAU*(D(I,K)+R(I,K))
      enddo
   enddo

   ! (3) CALCULER OPERATEUR DU COTE GAUCHE

   do K=1,NKX
      do I=1,N
         A(I,K)= -F*TAU*A(I,K)
         B(I,K)=1-F*TAU*B(I,K)
         C(I,K)= -F*TAU*C(I,K)
      enddo
   enddo

   ! (4) AJOUTER TERME DE FLUX DE SURFACE

   SFCFLUX = .true.
   select case (type)
   case (:2) !TYPE='U'/'UT'
      do I=1,N
         HD(I) = SB(i) - SK(i,NK-1)
      enddo
   case (5) !TYPE='ET'
      do I=1,N
         HD(I) = SB(i) - S(i,NKX)
         HD2(I,1) = S(I,1)
         HD2(I,NKX) = HD(I)
      enddo
      do K=2,NKX-1
         do I=1,N
            HD2(I,K) = S(I,K) - S(I,K-1)
         enddo
      enddo
   case (6) !TYPE='ST'
      do I=1,N
         HD(I) = SB(i) - 0.5*(SK(i,NKX)+SK(i,NKX-1))
         HD2(I,1) = 0.5 * (SK(I,1) + SK(I,2))
         HD2(I,NKX) = HD(I)
      enddo
      do K=2,NKX-1
         do I=1,N
            HD2(I,K) = 0.5 * (SK(I,K+1) - SK(I,K-1))
         enddo
      enddo
   case DEFAULT
      SFCFLUX = .false.
   end select
   if (SFCFLUX) then
      do I=1,N
         B(I,NKX)=B(I,NKX)-TAU*BETA(I)/HD(I)
         D(I,NKX)=D(I,NKX)+(ALFA(I)+BETA(I)*U(I,NKX))*TAU/HD(I)
      enddo
      do K=1,NKX-1
         do I=1,N
            D(I,K) = D(I,K) + BND(I,K)*TAU/HD2(I,K)
         enddo
      enddo
   endif

   ! (5) RESOUDRE SYSTEME TRIDIAGONAL [A,B,C] X = D. METTRE X DANS TU.

   call DIFUVD2(TU, A, B, C, D, D, NU, N, NKX)

   ! (6) OBTENIR TENDANCE

   do K=1,NKX
      do I=1,N
         TU(I,K)=TU(I,K)/TAU
      enddo
   enddo
   !     K=NKX+1..NK
   do K=NKX+1,NK
      do I=1,N
         TU(I,K)=0.0
      enddo
   enddo

   return
 end subroutine DIFUVDFj_2D
