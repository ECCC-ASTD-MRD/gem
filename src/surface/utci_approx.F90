
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
module MODI_UTCI_APPROX

   !#TODO: no need for an interface here, put the funtion UTCI_APPROX in the module to avoid non consistency between module and function interface

   interface
      function UTCI_APPROX(PTA,PEHPA,PTMRT,PVA) result(PUTCI_APPROX)
         real, dimension(:), intent(IN) :: PTA
         real, dimension(:), intent(IN) :: PEHPA
         real, dimension(:), intent(IN) :: PTMRT
         real, dimension(:), intent(IN) :: PVA
         real, dimension(size(PTA)) :: PUTCI_APPROX
      end function UTCI_APPROX
   end interface
end module MODI_UTCI_APPROX


!~ *************************************************************
function UTCI_APPROX(PTA,PEHPA,PTMRT,PVA) result(PUTCI_APPROX)
   !~ *************************************************************
   !    PURPOSE       : Computes the Universal Thermal and Climate Index Equivalent temperature
   !    AUTHOR        : P. Broede
   !    REFERENCE     : www.utci.org
   !    MODIFICATIONS :  
   !    METHOD        :  UTCI method ~ computed by a 6th order approximating polynomial from the 4 Input paramters 
   !                     ~  UTCI_approx, Version a 0.002, October 2009
   !                     ~  Copyright (C) 2009  Peter Broede
   !~ DOUBLE PRECISION Function value is the UTCI in degree Celsius
   !~ Input parameters (all of type DOUBLE PRECISION)
   !~ - Ta       : air temperature, degree Celsius
   !~ - ehPa    : water vapour presure, hPa=hecto Pascal
   !~ - Tmrt   : mean radiant temperature, degree Celsius
   !~ - va10m  : wind speed 10 m above ground level in m)s
   !-------------------------------------------------------------------------------
   !
   implicit none
!!!#include <arch_specific.hf>
   !~ type of input of the argument list
   real, dimension(:), intent(IN) :: PTA
   real, dimension(:), intent(IN) :: PEHPA
   real, dimension(:), intent(IN) :: PTMRT
   real, dimension(:), intent(IN) :: PVA
   real, dimension(size(PTA)) :: PUTCI_APPROX

   !local variables
   real, dimension(size(PTA)) :: ZPA
   real, dimension(size(PTA)) :: ZD_TMRT

   real, dimension(7,7,7,7) :: ZZ

   real :: Z0, Z1, Z2, Z3, Z4, Z5, Z6, ZF, ZS
   ZS(Z0,Z1,Z2,Z3,Z4,Z5,Z6,ZF) = Z0 + Z1*ZF + Z2*ZF**2 + Z3*ZF**3 + Z4*ZF**4 + Z5*ZF**5 + Z6*ZF**6
   real, dimension(size(PTA),7) :: ZC_TA, ZC_VA, ZC_TMRT, ZC_PA
   integer :: J1, J2, J3, J4, JJ
 
   !REAL(KIND=JPRB) :: ZHOOK_HANDLE

   !IF (LHOOK) CALL DR_HOOK('UTCI_APPROX',0,ZHOOK_HANDLE)
 
   ZZ(:,:,:,:)=0.
   ! va
   !                  Ta**0,         Ta**1,           Ta**2,         Ta**3,          Ta**4,           Ta**5,          Ta**6
   ZZ(1,1,1,1:7) = (/ 6.07562052D-01,-2.27712343D-02, 8.06470249D-04,-1.54271372D-04,-3.24651735D-06,&
        7.32602852D-08, 1.35959073D-09/) ! va**0
   ZZ(1,1,2,1:6) = (/-2.25836520D+00, 8.80326035D-02, 2.16844454D-03,-1.53347087D-05,-5.72983704D-07,&
        -2.55090145D-09/) ! va**1
   ZZ(1,1,3,1:5) = (/-7.51269505D-01,-4.08350271D-03,-5.21670675D-05, 1.94544667D-06, 1.14099531D-08/) ! va**2
   ZZ(1,1,4,1:4) = (/ 1.58137256D-01,-6.57263143D-05, 2.22697524D-07,-4.16117031D-08/)                 ! va**3
   ZZ(1,1,5,1:3) = (/-1.27762753D-02, 9.66891875D-06, 2.52785852D-09/)                                 ! va**4
   ZZ(1,1,6,1:2) = (/ 4.56306672D-04,-1.74202546D-07/)                                                 ! va**5
   ZZ(1,1,7,1)   = -5.91491269D-06                                                                     ! va**6
   ! D_Tmrt / va
   ZZ(1,2,1,1:6) = (/ 3.98374029D-01, 1.83945314D-04,-1.73754510D-04,-7.60781159D-07, 3.77830287D-08, 5.43079673D-10/)
   ZZ(1,2,2,1:5) = (/-2.00518269D-02, 8.92859837D-04, 3.45433048D-06,-3.77925774D-07,-1.69699377D-09/)
   ZZ(1,2,3,1:4) = (/ 1.69992415D-04,-4.99204314D-05, 2.47417178D-07, 1.07596466D-08/)
   ZZ(1,2,4,1:3) = (/ 8.49242932D-05, 1.35191328D-06,-6.21531254D-09/)
   ZZ(1,2,5,1:2) = (/-4.99410301D-06,-1.89489258D-08/)
   ZZ(1,2,6,1)   = 8.15300114D-08
   ! D_Tmrt**2 / va
   ZZ(1,3,1,1:5) = (/ 7.55043090D-04,-5.65095215D-05,-4.52166564D-07, 2.46688878D-08, 2.42674348D-10/)
   ZZ(1,3,2,1:4) = (/ 1.54547250D-04, 5.24110970D-06,-8.75874982D-08,-1.50743064D-09/)
   ZZ(1,3,3,1:3) = (/-1.56236307D-05,-1.33895614D-07, 2.49709824D-09/)
   ZZ(1,3,4,1:2) = (/ 6.51711721D-07, 1.94960053D-09/)
   ZZ(1,3,5,1)   = -1.00361113D-08
   !D_Tmrt**3 / va
   ZZ(1,4,1,1:4) = (/-1.21206673D-05,-2.18203660D-07, 7.51269482D-09, 9.79063848D-11/)
   ZZ(1,4,2,1:3) = (/ 1.25006734D-06,-1.81584736D-09,-3.52197671D-10/)
   ZZ(1,4,3,1:2) = (/-3.36514630D-08, 1.35908359D-10/)
   ZZ(1,4,4,1)   = 4.17032620D-10
   !D_Tmrt**4 / va
   ZZ(1,5,1,1:3) = (/-1.30369025D-09, 4.13908461D-10, 9.22652254D-12/)
   ZZ(1,5,2,1:2) = (/-5.08220384D-09,-2.24730961D-11/)
   ZZ(1,5,3,1)   = 1.17139133D-10
   !D_Tmrt**5 / va
   ZZ(1,6,1,1:2) = (/6.62154879D-10, 4.03863260D-13/)
   ZZ(1,6,2,1)   = 1.95087203D-12
   !D_Tmrt**6
   ZZ(1,7,1,1)   = -4.73602469D-12
   ! Pa / va
   ZZ(2,1,1,1:6) = (/ 5.12733497D+00,-3.12788561D-01,-1.96701861D-02, 9.99690870D-04, 9.51738512D-06,-4.66426341D-07/)
   ZZ(2,1,2,1:5) = (/ 5.48050612D-01,-3.30552823D-03,-1.64119440D-03,-5.16670694D-06, 9.52692432D-07/)
   ZZ(2,1,3,1:4) = (/-4.29223622D-02, 5.00845667D-03, 1.00601257D-06,-1.81748644D-06/)
   ZZ(2,1,4,1:3) = (/-1.25813502D-03,-1.79330391D-04, 2.34994441D-06/)
   ZZ(2,1,5,1:2) = (/ 1.29735808D-04, 1.29064870D-06/)
   ZZ(2,1,6,1)   = -2.28558686D-06
   ! Pa / D_Tmrt / va
   ZZ(2,2,1,1:5) = (/-3.69476348D-02, 1.62325322D-03,-3.14279680D-05, 2.59835559D-06,-4.77136523D-08/)
   ZZ(2,2,2,1:4) = (/ 8.64203390D-03,-6.87405181D-04,-9.13863872D-06, 5.15916806D-07/)
   ZZ(2,2,3,1:3) = (/-3.59217476D-05, 3.28696511D-05,-7.10542454D-07/)
   ZZ(2,2,4,1:2) = (/-1.24382300D-05,-7.38584400D-09/)
   ZZ(2,2,5,1)   = 2.20609296D-07
   ! Pa / D_Tmrt**2 / va
   ZZ(2,3,1,1:4) = (/-7.32469180D-04,-1.87381964D-05, 4.80925239D-06,-8.75492040D-08/)
   ZZ(2,3,2,1:3) = (/ 2.77862930D-05,-5.06004592D-06, 1.14325367D-07/)
   ZZ(2,3,3,1:2) = (/ 2.53016723D-06,-1.72857035D-08/)
   ZZ(2,3,4,1)   = -3.95079398D-08
   ! Pa / D_Tmrt**3 / va
   ZZ(2,4,1,1:3) = (/-3.59413173D-07, 7.04388046D-07,-1.89309167D-08/)
   ZZ(2,4,2,1:2) = (/-4.79768731D-07, 7.96079978D-09/)
   ZZ(2,4,3,1)   = 1.62897058D-09
   ! Pa / D_Tmrt**4 / va
   ZZ(2,5,1,1:2) = (/ 3.94367674D-08,-1.18566247D-09/)
   ZZ(2,5,2,1)   = 3.34678041D-10
   ! Pa / D_Tmrt**5
   ZZ(2,6,1,1)   = -1.15606447D-10 
   ! Pa**2 / va
   ZZ(3,1,1,1:5) = (/-2.80626406D+00, 5.48712484D-01,-3.99428410D-03,-9.54009191D-04, 1.93090978D-05/)
   ZZ(3,1,2,1:4) = (/-3.08806365D-01, 1.16952364D-02, 4.95271903D-04,-1.90710882D-05/)
   ZZ(3,1,3,1:3) = (/ 2.10787756D-03,-6.98445738D-04, 2.30109073D-05/)
   ZZ(3,1,4,1:2) = (/ 4.17856590D-04,-1.27043871D-05/)
   ZZ(3,1,5,1)   = -3.04620472D-06
   ! Pa**2 / D_Tmrt / va
   ZZ(3,2,1,1:4) = (/ 5.14507424D-02,-4.32510997D-03, 8.99281156D-05,-7.14663943D-07/)
   ZZ(3,2,2,1:3) = (/-2.66016305D-04, 2.63789586D-04,-7.01199003D-06/)
   ZZ(3,2,3,1:2) = (/-1.06823306D-04, 3.61341136D-06/)
   ZZ(3,2,4,1)   = 2.29748967D-07
   ! Pa**2 / D_Tmrt**2 / va
   ZZ(3,3,1,1:3) = (/3.04788893D-04,-6.42070836D-05, 1.16257971D-06/)
   ZZ(3,3,2,1:2) = (/7.68023384D-06,-5.47446896D-07/)
   ZZ(3,3,3,1)   = -3.59937910D-08
   ! Pa**2 / D_Tmrt**3 / va
   ZZ(3,4,1,1:2) = (/-4.36497725D-06, 1.68737969D-07/)
   ZZ(3,4,2,1)   = 2.67489271D-08
   ! Pa**2 / D_Tmrt**4
   ZZ(3,5,1,1) =  3.23926897D-09
   ! Pa**3 / va
   ZZ(4,1,1,1:4) = (/-3.53874123D-02, -2.21201190D-01, 1.55126038D-02, -2.63917279D-04/)
   ZZ(4,1,2,1:3) = (/4.53433455D-02, -4.32943862D-03, 1.45389826D-04/)
   ZZ(4,1,3,1:2) = (/2.17508610D-04, -6.66724702D-05/)
   ZZ(4,1,4,1)   = 3.33217140D-05
   ! Pa**3 / D_Tmrt / va
   ZZ(4,2,1,1:3) = (/-2.26921615D-03, 3.80261982D-04, -5.45314314D-09/)
   ZZ(4,2,2,1:2) = (/-7.96355448D-04, 2.53458034D-05/)
   ZZ(4,2,3,1)   = -6.31223658D-06
   ! Pa**3 / D_Tmrt**2 / va
   ZZ(4,3,1,1:2) = (/3.02122035D-04, -4.77403547D-06/)
   ZZ(4,3,2,1)   = 1.73825715D-06
   ! Pa**3 / D_Tmrt**3
   ZZ(4,4,1,1) = -4.09087898D-07
   !  Pa**4 / va
   ZZ(5,1,1,1:3) = (/6.14155345D-01, -6.16755931D-02, 1.33374846D-03/)
   ZZ(5,1,2,1:2) = (/3.55375387D-03, -5.13027851D-04/)
   ZZ(5,1,3,1)   = 1.02449757D-04
   !  Pa**4 / D_Tmrt / va
   ZZ(5,2,1,1:2) = (/-1.48526421D-03, -4.11469183D-05/)
   ZZ(5,2,2,1)   = -6.80434415D-06
   !  Pa**4 / D_Tmrt**2 / va
   ZZ(5,3,1,1) = -9.77675906D-06
   !  Pa**5 / va
   ZZ(6,1,1,1:2) = (/8.82773108D-02, -3.01859306D-03/)
   ZZ(6,1,2,1)   = 1.04452989D-03
   !  Pa**5 / D_Tmrt
   ZZ(6,2,1,1) = 2.47090539D-04
   ! Pa**6
   ZZ(7,1,1,1) = 1.48348065D-03
   !
   ZD_TMRT = PTMRT - PTA
   ZPA = PEHPA / 10.0; !~ use vapour pressure in kPa
   !
   ZC_TA(:,:) = 0.
   ZC_VA(:,:) = 0.
   ZC_TMRT(:,:) = 0.
   ZC_PA(:,:) = 0.
   !
   do J4 = 1,7
      do J3 = 1,7
         do J2 = 1,7
            do J1 = 1,7
               ZC_TA(:,J1) = ZZ(J4,J3,J2,J1)
            enddo
            do JJ=1,size(PTA)
               !ZC_VA(JJ,J2) = ZC_TA(JJ,1)+ZC_TA(JJ,2)*PTA(JJ)+ZC_TA(JJ,3)*PTA(JJ)**2+ZC_TA(JJ,4)*PTA(JJ)**3+ZC_TA(JJ,5)*PTA(JJ)**4 &
               !               +ZC_TA(JJ,6)*PTA(JJ)**5+ZC_TA(JJ,7)*PTA(JJ)**6
               ZC_VA(JJ,J2) = ZS(ZC_TA(JJ,1),ZC_TA(JJ,2),ZC_TA(JJ,3),ZC_TA(JJ,4),ZC_TA(JJ,5),ZC_TA(JJ,6),ZC_TA(JJ,7),PTA(JJ))
            enddo
         enddo
         do JJ=1,size(PTA)
            !ZC_TMRT(JJ,J3) = ZC_VA(JJ,1)+ZC_VA(JJ,2)*PVA(JJ)+ZC_VA(JJ,3)*PVA(JJ)**2+ZC_VA(JJ,4)*PVA(JJ)**3+ZC_VA(JJ,5)*PVA(JJ)**4 &
            !               +ZC_VA(JJ,6)*PVA(JJ)**5+ZC_VA(JJ,7)*PVA(JJ)**6
            ZC_TMRT(JJ,J3) = ZS(ZC_VA(JJ,1),ZC_VA(JJ,2),ZC_VA(JJ,3),ZC_VA(JJ,4),ZC_VA(JJ,5),ZC_VA(JJ,6),ZC_VA(JJ,7),PVA(JJ))
         enddo
      enddo
      do JJ=1,size(PTA)
         !ZC_PA(JJ,J4) = ZC_TMRT(JJ,1)+ZC_TMRT(JJ,2)*ZD_TMRT(JJ)+ZC_TMRT(JJ,3)*ZD_TMRT(JJ)**2+ZC_TMRT(JJ,4)*ZD_TMRT(JJ)**3 &
         !               +ZC_TMRT(JJ,5)*ZD_TMRT(JJ)**4 +ZC_TMRT(JJ,6)*ZD_TMRT(JJ)**5+ZC_TMRT(JJ,7)*ZD_TMRT(JJ)**6
         ZC_PA(JJ,J4) = ZS(ZC_TMRT(JJ,1),ZC_TMRT(JJ,2),ZC_TMRT(JJ,3),ZC_TMRT(JJ,4),ZC_TMRT(JJ,5),ZC_TMRT(JJ,6),ZC_TMRT(JJ,7),ZD_TMRT(JJ))
      enddo
   enddo
   !
   do JJ=1,size(PTA)
      !PUTCI_APPROX(JJ) = PTA(JJ) + ZC_PA(JJ,1)+ZC_PA(JJ,2)*ZPA(JJ)+ZC_PA(JJ,3)*ZPA(JJ)**2+ZC_PA(JJ,4)*ZPA(JJ)**3 &
      !                  +ZC_PA(JJ,5)*ZPA(JJ)**4 + ZC_PA(JJ,6)*ZPA(JJ)**5+ZC_PA(JJ,7)*ZPA(JJ)**6
      PUTCI_APPROX(JJ) = PTA(JJ) + ZS(ZC_PA(JJ,1),ZC_PA(JJ,2),ZC_PA(JJ,3),ZC_PA(JJ,4),ZC_PA(JJ,5),ZC_PA(JJ,6),ZC_PA(JJ,7),ZPA(JJ))
   enddo
   !
   !IF (LHOOK) CALL DR_HOOK('UTCI_APPROX',1,ZHOOK_HANDLE)
   !
end function UTCI_APPROX

