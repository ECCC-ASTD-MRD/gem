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
!**S/P SUNCOS1
!
!!if_on
 SUBROUTINE MACH_SUNCOS(SCOS,LMX,XLAT,XLON,HZ,DATE)
!!if_off
      use chm_consphychm_mod, only: pi
      implicit none
!! Jack C. (application for mach_adom2_messy.ftn90 in GEMMACH
!!if_on
          INTEGER(kind=4), intent(in):: LMX
          REAL(kind=4), intent(out)  :: SCOS(LMX)
          REAL(kind=4), intent(in)   :: XLAT(LMX),XLON(LMX),HZ,DATE
!!if_off
          EXTERNAL VSCOS,VSSIN
! local variables
          INTEGER I
          REAL DH,SDEC,CDEC,RDEC,AJOUR
          REAL A,EOT
          REAL DAYRAD

!!! taken from GEM phys GEM/x/4.8.rc8 for gemmach photolysis module
!!! "suncos1.f90" --> "mach_suncos.ftn90"
!!! this is necessary for COS(SZA)< 0 or SZA > 90deg.
!!!  where in the GEM application, minval is 0.00001 (see below)
!
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
!
!          - Input -
! LMX      number of points
! XLAT     latitude in radians
! XLON     longitude in radians
! HZ       Greenwich hour (0 to 24)
! DATE     julian day (0 to 366) (real number)
!
!*
!
      REAL, dimension(LMX) :: tmcos
      REAL, dimension(LMX) :: tmsin
!
       IF(DATE.NE.0.) AJOUR=DATE

! Declinaision solaire de Robertson et Russelo 1968
      dayrad=AJOUR*2.*PI/365
      rdec=.3964 + 3.631*sin(dayrad) - 22.97*cos(dayrad) + .03838*sin(2.*dayrad) -0.3885*cos(2.*dayrad) + &
                 .07659*sin(3.*dayrad) -0.1587*cos(3.*dayrad)-0.01021*cos(4.*dayrad)

      rdec=rdec*PI/180.

! Declinaison solaire: approximation qui suppose que l'orbite est circulaire
!      RDEC=0.412*COS((AJOUR+10.)*2.*PI/365.25 -PI)

      SDEC=SIN(RDEC)
      CDEC=COS(RDEC)
! correction for "equation of time"
      A = DATE/365.*2.*PI
! in minutes
      EOT = .002733 -7.343*sin(a)+ .5519*cos(a) -9.47*sin(2.*a) &
        -3.02*cos(2.*a) -.3289*sin(3.*a) -.07581*cos(3.*a) &
       -.1935*sin(4.*a) -.1245*cos(4.*a)
! express in a fraction of hour
      EOT=EOT/60.
! express in radians
      EOT=EOT*15.*PI/180.
!
      call VSCOS(tmcos(1),XLAT(1),LMX)
      call VSSIN(tmsin(1),XLAT(1),LMX)
!
      DO I=1,LMX
      DH=HZ*PI/12. +XLON(I) -PI + EOT

! modified below for photolysis calculation in GEMMACH - need SZA>90deg.
! in GEM:
!!    SCOS(I)=AMAX1(tmsin(I)*SDEC + tmcos(I)*CDEC*COS(DH), 0.00001)
! for GEMMACH:
      SCOS(I)=AMAX1(tmsin(I)*SDEC + tmcos(I)*CDEC*COS(DH), -1.0)
      SCOS(I)=AMIN1(SCOS(I), 1.0)
      ENDDO

      RETURN
 end SUBROUTINE MACH_SUNCOS
