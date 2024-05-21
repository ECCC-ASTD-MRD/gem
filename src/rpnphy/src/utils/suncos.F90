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

module suncos
   private
   public :: suncos1, suncos3

contains

   !/@*
   subroutine suncos1(scos, lmx, xlat, xlon, hz, dayofyear, &
        F_eot, F_sdec, F_cdec, F_tmcos, F_tmsin)
      use tdpack_const, only: PI
      implicit none
      !@object Calculate the cosines of the solar angle (scos = [-1., 1.])
      !@Author L.Garand (1989)
      !@arguments
      !          - Output -
      ! SCOS     cosines of the solar angle
      ! F_eot, F_sdec, F_cdec, F_tmcos, F_tmsin: optional internal vars for use in suncos2 with slope=.T.
      !          - Input -
      ! LMX      number of points
      ! XLAT     latitude in radians
      ! XLON     longitude in radians
      ! HZ       Greenwich hour (0 to 24)
      ! DAYOFYEAR  day of year (0 to 366) (real number)
      integer, intent(in) :: lmx
      real, intent(in) :: hz, dayofyear
      real, intent(in), dimension(lmx) :: xlat, xlon
      real, intent(out), dimension(lmx) :: scos
      real, intent(out), optional :: F_eot, F_sdec, F_cdec
      real, intent(out), dimension(lmx), optional :: F_tmcos, F_tmsin
      !*@/
      integer :: i
      real :: ajour, dayrad, rdec, sdec, cdec, a, eot, dh
      real, dimension(lmx) :: tmcos, tmsin
      !----------------------------------------------------------------
      ajour = 1.
      if (dayofyear /= 0.) ajour = dayofyear

      ! declinaision solaire de robertson et russelo 1968
      dayrad = ajour*2.*PI/365
      rdec = .3964 + 3.631*sin(dayrad) - 22.97*cos(dayrad) + .03838*sin(2.*dayrad) -0.3885*cos(2.*dayrad) + &
           .07659*sin(3.*dayrad) -0.1587*cos(3.*dayrad)-0.01021*cos(4.*dayrad)

      rdec = rdec*PI/180.

      ! declinaison solaire: approximation qui suppose que l'orbite est circulaire
      !      rdec=0.412*cos((ajour+10.)*2.*PI/365.25 -PI)

      sdec = sin(rdec)
      cdec = cos(rdec)
      ! correction for "equation of time"
      a = dayofyear/365.*2.*PI
      ! in minutes
      eot = .002733 - 7.343*sin(a) + .5519*cos(a) - 9.47*sin(2.*a) &
           - 3.02*cos(2.*a) - .3289*sin(3.*a) - .07581*cos(3.*a) &
           - .1935*sin(4.*a) - .1245*cos(4.*a)
      ! express in a fraction of hour
      eot = eot/60.
      ! express in radians
      eot = eot*15.*PI/180.

      do i=1,lmx
         dh = hz*PI/12. + xlon(i) - PI + eot
         tmcos(i) = cos(xlat(i))
         tmsin(i) = sin(xlat(i))
         scos(i) = amax1(-1.0, tmsin(i)*sdec + tmcos(i)*cdec*cos(dh))
         scos(i) = amin1(scos(i), 1.0)
      enddo
      if (present(F_eot)) F_eot = eot
      if (present(F_sdec)) F_sdec = sdec
      if (present(F_cdec)) F_cdec = cdec
      if (present(F_tmcos)) F_tmcos = tmcos
      if (present(F_tmsin)) F_tmsin = tmsin
      !----------------------------------------------------------------
      return
   end subroutine suncos1

   
   !/@*
   subroutine suncos3(scos, lmx, xlat, xlon, hz, dayofyear)
      use tdpack_const
      implicit none
!!!#include <arch_specific.hf>
      !@object Calculate the cosines of the solar angle (scos = ]0., 1.])
      !@arguments
      !          - Output -
      ! SCOS     cosines of the solar angle
      !          - Input -
      ! LMX      number of points
      ! XLAT     latitude in radians
      ! XLON     longitude in radians
      ! HZ       Greenwich hour (0 to 24)
      ! DAYOFYEAR  day of year (0 to 366) (real number)
      integer, intent(in) :: lmx
      real, intent(in) :: hz, dayofyear
      real, intent(in), dimension(lmx) :: xlat, xlon
      real, intent(out), dimension(lmx) :: scos
      !@Author L.Garand (1989)
      !*@/
      !----------------------------------------------------------------
      call suncos1(scos, lmx, xlat, xlon, hz, dayofyear)
      scos = amax1(0.00001, scos)
      !----------------------------------------------------------------
      return
   end subroutine suncos3

end module suncos


!/@*
subroutine suncos2(scos, ssin, stan, bsin, bcos, lmx, xlat, xlon, hz, dayofyear, slope_l)
   use suncos, only: suncos1
   use tdpack_const, only: PI
   implicit none
!!!#include <arch_specific.hf>
   !@object Calculate the cosines of the solar angle (scos = ]0., 1.])
   !@arguments
   !          - Output -
   ! SCOS     cosines of the solar angle
   ! SSIN     sine of the solar angle     (if slope_l=.T.)
   ! STAN     tangents of the solar angle (if slope_l=.T.)
   ! BSIN     sines of beta               (if slope_l=.T.)
   ! BCOS     cosines of beta             (if slope_l=.T.)
   !          - Input -
   ! LMX      number of points
   ! XLAT     latitude in radians
   ! XLON     longitude in radians
   ! HZ       Greenwich hour (0 to 24)
   ! DAYOFYEAR  day of year (0 to 366) (real number)
   logical, intent(in) :: slope_L
   integer, intent(in) :: lmx
   real, intent(in) :: hz, dayofyear
   real, intent(in), dimension(lmx) :: xlat, xlon
   real, intent(out), dimension(lmx) :: scos, ssin, stan, bsin, bcos
   !@Author L.Garand (1989)
   !@ Revision
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
   !*@/
   integer :: i
   real :: dh, eot, sdec, cdec
   real, dimension(lmx) :: tmcos, tmsin
   !----------------------------------------------------------------
   if (slope_L) then
      call suncos1(scos, lmx, xlat, xlon, hz, dayofyear, eot, sdec, cdec, tmcos, tmsin)
   else
      call suncos1(scos, lmx, xlat, xlon, hz, dayofyear)
   endif
   scos = amax1(0.00001, scos)

   if(slope_l) then

      do i=1,lmx

         dh = hz*PI/12. + xlon(i) - PI + eot

         !           to calculate sin z

         ssin(i) = amin1(sqrt(1.-(scos(i)*scos(i))), 1.)
         ssin(i) = amax1(ssin(i), 0.00001)

         !           to calculate tan z

         stan(i) = ssin(i)/scos(i)

         !           to calculate sin b

         bsin(i) = amin1((cdec*sin(dh))/ssin(i), 1.)
         bsin(i) = amax1(bsin(i), -1.)

         !           to calulate cos b

         bcos(i) = (scos(i) * tmsin(i) - sdec) / (ssin(i) * amax1(tmcos(i), 0.00001))
         bcos(i) = amin1(bcos(i), 1.)
         bcos(i) = amax1(bcos(i), -1.)

      enddo
   else
      ssin = 0.
      stan = 0.
      bsin = 0.
      bcos = 0.         
   endif
   !----------------------------------------------------------------
   return
end subroutine suncos2
