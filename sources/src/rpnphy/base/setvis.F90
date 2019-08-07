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

subroutine setvis4(DSIG,DSH,DSC,DZ,RMUO,QOZSTP, &
     SIG, T, PSOL, QOZ, XLAT, XLON, HZ, DAYOFYEAR, &
     LMX, LEV, M)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>

   real DAYOFYEAR
   integer LMX,LEV,KMX
   integer I,J,K,M
   real DSIG(lmx,LEV),DSH(lmx,LEV),DSC(lmx,LEV), &
        PSOL(LMX),QOZSTP(LMX,LEV),RMUO(LMX),DZ(LMX,LEV)
   real SIG(lmx,LEV),T(M,LEV),XLAT(LMX),XLON(LMX),QOZ(LMX,LEV)
   real A1,A2,HZ

   !@Author
   !          L.Garand (1989)

   !@Revision
   ! 001      G.Pellerin(Mar90)Standard documentation
   ! 002      N. Brunet  (May91)
   !             New version of thermodynamic functions
   !             and file of constants
   ! 003      R. Benoit (Aug 93) Local Sigma
   ! 004      L. Garand (Apr 95) Reduction mode
   ! 005      L. Garand (April 95) Aerosols and liquid water
   !          no more calculated here
   ! 006      J. Toviessi (July 2009) changed the call from suncos to
   !                                  suncos1


   !@Object
   !          to calculate the inputs for solar radiation

   !@Arguments
   !          - Output -
   ! DSH      sigma thickness with exponent 1.9
   ! DSC      sigma thickness with exponent 1.75
   ! DZ       thickness of layers in metres
   ! RMUO     cosines of the solar angle
   ! QOZSTP   ozone in CMSTP ( cm at standard pressure )
   !          for each layer
   !          - Input -
   ! DSIG     sigma thickness of layers
   ! SIG      sigma levels
   ! T        temperature at that level
   ! PSOL     surface pressure
   ! QOZ      ozone in kg/m2 for subroutine RADFACE
   ! XLAT     latitude in radians
   ! XLON     longitude in radians
   ! HZ       Greenwich hour
   ! DAYOFYEAR  day of year (1..365)
   ! LMX      number of points to process
   ! LEV      number of layers
   ! M        1st dimension of T

   !@Notes
   !          1st dimension of T can be different from N

   real DUMMY1(lmx),DUMMY2(lmx),DUMMY3(lmx),DUMMY4(lmx)

   KMX=LEV-1
   !    CALCUL DES EPAISSEURS

   do i=1,lmx
      A2=(SIG(i,1)+SIG(i,2))/2.
      A1=2.*SIG(i,1)-A2
      !        A1=AMAX1(A1,1.E-8)
      !        la ligne suivante peut etre substituee a la ligne precedente
      A1=AMAX1(A1,sig(i,1)/2.)
      !        DSIG(i,1)=A2-A1
      DSH(i,1)=A2**1.9-A1**1.9
      DSC(i,1)=A2**1.75-A1**1.75
   enddo
   do K=2,KMX
      do i=1,lmx
         A1=(SIG(i,K)+SIG(i,K-1))/2.
         A2=(SIG(i,K)+SIG(i,K+1))/2.
         !           DSIG(i,K)=A2-A1
         DSH(i,K)= A2**1.9 - A1**1.9
         DSC(i,K)= A2**1.75-A1**1.75
      enddo
   enddo

   do i=1,lmx
      A2=(SIG(i,KMX)+SIG(i,KMX+1))/2.
      !        DSIG(i,LEV)=1.- A2
      DSH(i,LEV)=1.- A2 **1.9
      DSC(i,LEV)=1.- A2 **1.75
   enddo

   do J=1,LEV
      do I=1,LMX
         DZ(I,J) = DSIG(i,J)*RGASD*T(I,J)/GRAV/SIG(i,J)
      enddo
   enddo

   !     COSINUS ANGLE SOLAIRE
   call SUNCOS2(RMUO,DUMMY1,DUMMY2,DUMMY3,DUMMY4, &
        LMX,XLAT,XLON,HZ,DAYOFYEAR,.false.)


   !     OZONE EN CMSTP
   do J=1,LEV
      do I=1,LMX
         !   INVERSE DE CONVERSION DANS RADFACE
         QOZSTP(I,J) = QOZ(I,J)*PSOL(I)*DSIG(i,J)/GRAV/2.144E-2
      enddo
   enddo

   return
end subroutine setvis4
