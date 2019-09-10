!-------------------------------------- LICENCE BEGIN ------------------------
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

subroutine WFLUX(ROS,G,WC1,TO1,RE1,TR1,RMUE,RE2,TR2,LMX,ILG)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer ILG,I,LMX
   real TOP(LMX)
   real(REAL64) :: RI0(LMX),RI0D(LMX),RI1(LMX),RI1D(LMX),RK(LMX),RM2(LMX),WCP(LMX),WM(LMX)
   real(REAL64) :: X1(LMX),X2(LMX),XM2P(LMX),XP2P(LMX)
   real(REAL64) :: A11(LMX),A12(LMX),A13(LMX),A21(LMX),A22(LMX),A23(LMX),ALPHA(LMX)
   real(REAL64) :: AM2B(LMX),AP2B(LMX),BETA(LMX)
   real(REAL64) :: C1(LMX),C2(LMX),D(LMX),DT(LMX),EXKM(LMX),EXKP(LMX),EXMUO(LMX),FF(LMX)
   real(REAL64) :: GP(LMX),RP(LMX)
   real(REAL64) :: TMP1EXP(LMX),TMP2EXP(LMX)
   real ROS(LMX),G(LMX),WC1(LMX),TO1(LMX),RE1(LMX),TR1(LMX), &
        RMUE(LMX),RE2(LMX),TR2(LMX)

   !@Author L.Garand (1989)
   !@Revision
   ! 001      G.Pellerin(Mar90)Standard documentation
   ! 002      D. Talbot (June 2003) - IBM conversion
   !               - calls to vexp routine (from massvp4 library)
   !               - divisions replaced by reciprocals
   !@Object
   !          to estimate layer transmission and reflexion using
   !          the Delta Eddington Approximation(vectorized version)
   !@Arguments
   !          - Input -
   ! ROS      reflectivity of underlaying layer
   ! G        asymmetry factor
   ! WC1      omega parameter
   ! TO1      optical thickness / omega
   !          - Output -
   ! RE1      reflected radiation without reflexion from underlying
   !          layer
   ! TR1      transmitted radiation without reflexion from underlying
   !          layer
   !          - Input -
   ! RMUE     cosine of equivalent sun zenith angle
   !          - Output -
   ! RE2      reflected radiation with reflexion from underlying layer
   ! TR2      transmitted radiation with reflexion from underlying
   !          layer
   !          - Input -
   ! LMX      maximum number of profiles that can be requested
   ! ILG      actual number of profiles requested

   real REC_RMUE(lmx)
   real(REAL64) :: REC_X2(LMX),REC_D(LMX)
   real REC_3

   REC_3=1./3.

   do I=1,ILG

      FF(I)=G(I)*G(I)
      GP(I)=G(I)/(1.+G(I))
      TOP(I)=(1.-WC1(I)*FF(I))*TO1(I)
      WCP(I)=(1-FF(I))*WC1(I)/(1.-WC1(I)*FF(I))
      DT(I)=2.*REC_3
      X1(I)=1.-WCP(I)*GP(I)
      WM(1)=1.-WCP(I)
      REC_RMUE(I)=1.0/rmue(i)
      RM2(I)=RMUE(I)*RMUE(I)
      !     apres plusieurs essais avec vsqrt ca ne valide pas
      RK(I)=sqrt(3.*WM(1)*X1(I))
      X2(I)=4.*(1.-RK(I)*RK(I)*RM2(I))
      !     ENDDO

      !     ce call est tel que ca ne valide plus
      !     CALL VREC (REC_X2,X2,ILG)

      !     DO I=1,ILG
      REC_X2(I)=1.0d0/X2(I)
      RP(I)=sqrt(3.*WM(1)/X1(I))
      ALPHA(I)=3.*WCP(I)*RM2(I)*(1.+GP(I)*WM(1))*REC_X2(I)
      BETA(I)=3.*WCP(I)*RMUE(I)*(1.+3.*GP(I)*RM2(I)*WM(1))*REC_X2(I)
      TMP1EXP(I)=-TOP(I)*REC_RMUE(I)
      TMP2EXP(I)=RK(I)*TOP(I)
   enddo

   call VEXP (EXMUO,TMP1EXP,ILG)
   call VEXP (EXKP,TMP2EXP,ILG)
   call VREC (EXKM,EXKP,ILG)

   do I=1,ILG
      XP2P(I)=1.+DT(I)*RP(I)
      XM2P(I)=1.-DT(I)*RP(I)
      AP2B(I)=ALPHA(I)+DT(I)*BETA(I)
      AM2B(I)=ALPHA(I)-DT(I)*BETA(I)

      !  WITHOUT REFLEXION FROM THE UNDERLYING LAYER

      A11(I)=XP2P(I)
      A12(I)=XM2P(I)
      A13(I)=AP2B(I)
      A22(I)=XP2P(I)*EXKP(I)
      A21(I)=XM2P(I)*EXKM(I)
      A23(I)=AM2B(I)*EXMUO(I)
      D(I)=A11(I)*A22(I)-A21(I)*A12(I)
   enddo

   call VREC (REC_D,D,ILG)

   do I=1,ILG
      C1(I)=(A22(I)*A13(I)-A12(I)*A23(I))*REC_D(I)
      C2(I)=(A11(I)*A23(I)-A21(I)*A13(I))*REC_D(I)
      RI0(I)=C1(I)+C2(I)-ALPHA(I)
      RI1(I)=RP(I)*(C1(I)-C2(I))-BETA(I)
      RE1(I)=(RI0(I)-DT(I)*RI1(I))*REC_RMUE(I)
      RI0D(I)=C1(I)*EXKM(I)+C2(I)*EXKP(I)-ALPHA(I)*EXMUO(I)
      RI1D(I)=RP(I)*(C1(I)*EXKM(I)-C2(I)*EXKP(I))-BETA(I)*EXMUO(I)
      TR1(I)=EXMUO(I)+(RI0D(I)+DT(I)*RI1D(I))*REC_RMUE(I)

      !  WITH REFLEXION FROM THE UNDERLYING LAYER

      A21(I)=A21(I)-ROS(I)*XP2P(I)*EXKM(I)
      A22(I)=A22(I)-ROS(I)*XM2P(I)*EXKP(I)
      A23(I)=A23(I)-ROS(I)*EXMUO(I)*(AP2B(I)-RMUE(I))
      D(I)=A11(I)*A22(I)-A21(I)*A12(I)
   enddo

   call VREC (REC_D,D,ILG)

   do I=1,ILG
      C1(I)=(A22(I)*A13(I)-A12(I)*A23(I))*REC_D(I)
      C2(I)=(A11(I)*A23(I)-A21(I)*A13(I))*REC_D(I)
      RI0(I)=C1(I)+C2(I)-ALPHA(I)
      RI1(I)=RP(I)*(C1(I)-C2(I))-BETA(I)
      RE2(I)=(RI0(I)-DT(I)*RI1(I))*REC_RMUE(I)
      RI0D(I)=C1(I)*EXKM(I)+C2(I)*EXKP(I)-ALPHA(I)*EXMUO(I)
      RI1D(I)=RP(I)*(C1(I)*EXKM(I)-C2(I)*EXKP(I))-BETA(I)*EXMUO(I)
      TR2(I)=EXMUO(I)+(RI0D(I)+DT(I)*RI1D(I))*REC_RMUE(I)

   enddo
   return
end subroutine WFLUX
