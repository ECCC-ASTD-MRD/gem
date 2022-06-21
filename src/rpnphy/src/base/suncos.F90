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

subroutine SUNCOS2(SCOS,SSIN,STAN,BSIN,BCOS,LMX,XLAT,XLON,HZ,DAYOFYEAR, slope_L)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>
   integer I,LMX
   real XLAT(LMX),XLON(LMX),SCOS(LMX),HZ,DAYOFYEAR
   real DH,SDEC,CDEC,RDEC,AJOUR,SSIN(LMX),STAN(LMX),BSIN(LMX),BCOS(LMX)
   real A,EOT
   real DAYRAD
   logical slope_L

!Author
!          L.Garand (1989)
!
!Revision
! 001      G.Pellerin(Mar90)Standard documentation
! 002      N. Brunet  (May91)
!                New version of thermodynamic functions
!                and file of constants
! 003      L. Garand (Fev95) Add equation of time
! 004      J.P. Toviessi (June 2003) - IBM conversion
!               - calls to vscos, vssin routine (from massvp4 library)
!               - unnecessary calculations removed
! 005      J. P. Toviessi (July 2009) added modifications for radslope
!
! 006      P. Vaillancourt, I. Paunova (Oct 2009) Correct calculation of solar declination
!
!Object
!          to calculate the cosines of the solar angle for LMX
!          points
!
!Arguments
!
!          - Output -
! SCOS     cosines of the solar angle
! SSIN     sine of the solar angle
! STAN     tangents of the solar angle
! BSIN     sines of beta
! BCOS     cosines of beta
!
!          - Input -
! LMX      number of points
! XLAT     latitude in radians
! XLON     longitude in radians
! HZ       Greenwich hour (0 to 24)
! DAYOFYEAR  day of year (0 to 366) (real number)
!
!*

      real, dimension(LMX) :: tmcos, tmsin

      AJOUR=1.
       if(DAYOFYEAR.ne.0.) AJOUR=DAYOFYEAR

! Declinaision solaire de Robertson et Russelo 1968
      dayrad=AJOUR*2.*PI/365
      rdec=.3964 + 3.631*sin(dayrad) - 22.97*cos(dayrad) + .03838*sin(2.*dayrad) -0.3885*cos(2.*dayrad) + &
                 .07659*sin(3.*dayrad) -0.1587*cos(3.*dayrad)-0.01021*cos(4.*dayrad)

      rdec=rdec*PI/180.

! Declinaison solaire: approximation qui suppose que l'orbite est circulaire
!      RDEC=0.412*COS((AJOUR+10.)*2.*PI/365.25 -PI)

      SDEC=sin(RDEC)
      CDEC=cos(RDEC)
! correction for "equation of time"
      A = DAYOFYEAR/365.*2.*PI
! in minutes
      EOT = .002733 -7.343*sin(a)+ .5519*cos(a) -9.47*sin(2.*a) &
        -3.02*cos(2.*a) -.3289*sin(3.*a) -.07581*cos(3.*a) &
       -.1935*sin(4.*a) -.1245*cos(4.*a)
! express in a fraction of hour
      EOT=EOT/60.
! express in radians
      EOT=EOT*15.*PI/180.

      do I=1,LMX
      DH=HZ*PI/12. +XLON(I) -PI + EOT
      tmcos(I)=cos(XLAT(I))
      tmsin(I)=sin(XLAT(I))
      SCOS(I)=AMAX1(tmsin(I)*SDEC + tmcos(I)*CDEC*cos(DH), 0.00001)
      SCOS(I)=AMIN1(SCOS(I), 1.0)
      enddo

!----------------------------------------------------------

      if(slope_L) then

         do I=1,LMX

            DH=HZ*PI/12. +XLON(I) -PI + EOT

!           To calculate Sin Z

            SSIN(I)=amin1(sqrt(1.-(SCOS(I)*SCOS(I))), 1.)
            SSIN(I)=amax1(SSIN(I), 0.00001)

!           to calculate Tan Z

            STAN(I)=SSIN(I)/SCOS(I)

!           to calculate Sin B

            BSIN(I)=amin1((CDEC*sin(DH))/SSIN(I), 1.)
            BSIN(I)=amax1(BSIN(I), -1.)

!           to calulate Cos B

            BCOS(I)=(SCOS(I)*tmsin(I)-SDEC)/(SSIN(I)*amax1(tmcos(I),0.00001))
            BCOS(I)=amin1(BCOS(I), 1.)
            BCOS(I)=amax1(BCOS(I), -1.)

         enddo

      endif

      return
      end

