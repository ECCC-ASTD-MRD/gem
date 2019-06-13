!-------------------------------------- LICENCE BEGIN --------------------------
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
!-------------------------------------- LICENCE END ----------------------------

module surf_precip
   use tdpack_const, only: TCDK
   use phy_options, only: stcond
   implicit none
   private
   public :: surf_precip1, surf_precip3

   real, parameter :: EPSILON_4 = 1.E-10

contains
 
   !/@*
   subroutine surf_precip1(t, tlc, tls, tsc, tss, rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>

      !@Arguments
      !          - Input -
      ! T        low-level air temperature
      ! tlc      total liquid "convective" precipitation rate
      ! tls      total liquid "stratiform" precipitation rate
      ! tcs      total solid  "convective" precipitation rate
      ! tss      total solid  "stratiform" precipitation rate
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: t, tlc, tls, tsc, tss
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Author S. Belair (November 1998)
      !@Revisions
      ! 001      B. Bilodeau (January 2001)
      !          Automatic arrays
      ! 002      B. Bilodeau and S. Belair (July 2001)
      !          Change conditions for option 0
      ! 003      S. Belair (February 2004)
      !          Bug correction related to partition of total
      !          precipitation rate into liquid and solid.
      ! 004      L. Spacek (Aug 2004) - cloud clean-up
      !          elimination of ISTCOND=2,6,7,8 ICONVEC=4
      ! 005      D. Talbot (March 2005)
      !          Modification to the previous correction
      !          so that it does not apply to consun
      !
      !@Object Determine the phase of precipitation reaching the surface
      !*@/

      integer :: I, option
      real, dimension(ni) :: totrate, expli
      !----------------------------------------------------------------

      ! Sum of the precipitation rate at
      ! the surface

      do I=1,ni
         rainrate(I) = 0.
         snowrate(I) = 0.
         totrate (I) = tlc(I)+tls(I)+tsc(I)+tss(I)
      end do

      ! OPTIONS:  Depending on the explicit scheme
      !           used during the integration,
      !           we have three options for specifying
      !           the phase of precipitation reaching
      !           the ground

      if (stcond == 'NIL')  then
         option = 0
      else if (stcond == 'NEWSUND' .or. stcond == 'CONSUN') then
         option = 1
      else
         option = 2
      endif

      IF_OPTION: if (option == 0) then

         ! FIRST OPTION:
         ! The phase of the precipitation reaching
         ! the surface is simply determined from
         ! the low-level air temperature:
         !  T > 0.  then rain
         !  T < 0.  then snow

         do I=1,ni
            if (T(I) >= TCDK) then
               rainrate(I) = totrate(I)
               snowrate(I) = 0.
            else if (T(I) < TCDK) then
               rainrate(I) = 0.
               snowrate(I) = totrate(I)
            end if
         end do

      else if (option == 1) then  !IF_OPTION

         ! SECOND OPTION:
         ! The phase of precipitation at the surface
         ! is determined by results from the "explicit"
         ! condensation scheme (NEWSUND and CONSUN).
         ! It was noted in GEM15 that OPTION 2 was not a viable option
         ! for consun since the scheme produces a non zero tls
         ! in cold temperature in winter causing an erronous
         ! warm surface temperature bias especially in the arctic.

         do I=1,ni

            ! If "stratiform" (stable condensation) is greater
            ! than 0, then the explicit scheme determines if
            ! the precipitation is rain or snow

            if (tss(I) >= EPSILON_4) then
               rainrate(I) = 0.0
               snowrate(I) = totrate(I)
            else if (tls(I) >= EPSILON_4) then
               rainrate(I) = totrate(I)
               snowrate(I) = 0.0

               ! If "stratiform" precipitation is null, then
               ! we use low-level air temperature to specify
               ! the phase of precip.

            else if (T(I) >= TCDK) then
               rainrate(I) = totrate(I)
               snowrate(I) = 0.0
            else if (T(I) < TCDK) then
               rainrate(I) = 0.
               snowrate(I) = totrate(I)
            end if


         end do

      else  !IF_OPTION

         ! THIRD OPTION:
         ! The phase of precipitation at the surface
         ! is determined by results from the "explicit"
         ! condensation scheme (MIXED-PHASE, EXMO, KONG-YAU).
         ! It was noted in the LAM2.5 that KONG-YAU
         ! could produce a non zero tss over Quebec in summer
         ! and therefore OPTION 1 was no longer appropriate

         do I=1,ni
            expli(I) = tls(I) + tss(I)

            ! If "stratiform" (stable condensation) is greater
            ! than 0, then the explicit scheme determines if
            ! the precipitation is rain or snow

            if (expli(I) >= EPSILON_4) then
               rainrate(I) = tls(I)/expli(I) * totrate(I)
               snowrate(I) = tss(I)/expli(I) * totrate(I)

               ! If "stratiform" precipitation is null, then
               ! we use low-level air temperature to specify
               ! the phase of precip.

            else if (T(I) >= TCDK) then
               rainrate(I) = totrate(I)
               snowrate(I) = 0.0
            else if (T(I) < TCDK) then
               rainrate(I) = 0.
               snowrate(I) = totrate(I)
            end if

         end do

      end if IF_OPTION

      !----------------------------------------------------------------
      return
   end subroutine surf_precip1
 
 
   !/@*
   subroutine surf_precip3(t, tlc, tls, tsc, tss, &
        rainrate, snowrate, fneige, fip,ni)
      implicit none
!!!#include <arch_specific.hf>

      !@Arguments
      !          - Input -
      ! T        low-level air temperature
      ! tlc      total liquid "convective" precipitation rate
      ! tls      total liquid "stratiform" precipitation rate
      ! tsc      total solid  "convective" precipitation rate
      ! tss      total solid  "stratiform" precipitation rate
      ! FNEIGE   snow fraction from Bourge1 subroutine
      ! FIP      Fraction of liquid precipitation in the form of ice pellets.
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: T, tlc, tls, tsc, tss, FNEIGE, FIP
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Author S. Belair (November 1998)
      !@Revisions
      ! 001      B. Bilodeau (January 2001)
      !          Automatic arrays
      ! 002      B. Bilodeau and S. Belair (July 2001)
      !          Change conditions for option 0
      ! 003      S. Belair (February 2004)
      !          Bug correction related to partition of total
      !          precipitation rate into liquid and solid.
      ! 004      L. Spacek (Aug 2004) - cloud clean-up
      !          elimination of ISTCOND=2,6,7,8 ICONVEC=4
      ! 005      D. Talbot (March 2005)
      !          Modification to the previous correction
      !          so that it does not apply to consun
      ! 006      A-M.  Leduc (April 2006). Add tlcs and tscs to the totrate.
      !                      Change the option 1 to use FNEIGE. CHANGE name to
      !                      SURF_PRECIP2.
      ! 007     A-M. Leduc. A. Plante. (Nov 2007) add FIP. change s/r name to surf_precip3
      !                change calculation of snowrate to include all solid precipitation
      !                with f_solid.
      !@Object Determine the phase of precipitation reaching the surface
      !*@/

      integer :: I, option
      real :: F_SOLID
      real, dimension(ni) :: totrate, expli
      !----------------------------------------------------------------

      ! Sum of the precipitation rate at
      ! the surface

      do I=1,ni
         rainrate(I) = 0.
         snowrate(I) = 0.
         totrate (I) = tlc(I)+tls(I)+tsc(I)+tss(I)
      end do

      ! OPTIONS:  Depending on the explicit scheme
      !           used during the integration,
      !           we have three options for specifying
      !           the phase of precipitation reaching
      !           the ground

      if (stcond == 'NIL')  then
         option=0
      else if (stcond == 'NEWSUND' .or. stcond == 'CONSUN') then
         option=1
      else
         option=2
      endif

      IF_OPTION: if (option == 0) then

         ! FIRST OPTION:
         ! The phase of the precipitation reaching
         ! the surface is simply determined from
         ! the low-level air temperature:

         !  T > 0.  then rain
         !  T < 0.  then snow

         do I=1,ni
            if (T(I) >= TCDK) then
               rainrate(I) = totrate(I)
               snowrate(I) = 0.
            else if (T(I) < TCDK) then
               rainrate(I) = 0.
               snowrate(I) = totrate(I)
            end if
         end do

      else if (option == 1) then  !IF_OPTION

         ! SECOND OPTION:
         ! The phase of precipitation at the surface
         ! is determined using the bourge snow fraction
         ! fneige and ice pellet fraction fip.

         do I=1,ni

            F_SOLID     = FNEIGE(I) + ( 1.-FNEIGE(I) )* FIP(I)
            snowrate(I) = totrate(I) * F_SOLID
            rainrate(I) = totrate(I) * (1.-F_SOLID)

         end do


      else  !IF_OPTION

         ! THIRD OPTION:
         ! The phase of precipitation at the surface
         ! is determined by results from the "explicit"
         ! condensation scheme (MIXED-PHASE, EXMO, KONG-YAU).
         ! It was noted in the LAM2.5 that KONG-YAU
         ! could produce a non zero tss over Quebec in summer
         ! and therefore OPTION 1 was no longer appropriate

         do I=1,ni
            expli(I) = tls(I) + tss(I)

            ! If "stratiform" (stable condensation) is greater
            ! than 0, then the explicit scheme determines if
            ! the precipitation is rain or snow

            if (expli(I) >= EPSILON_4) then
               rainrate(I) = tls(I)/expli(I) * totrate(I)
               snowrate(I) = tss(I)/expli(I) * totrate(I)

               ! If "stratiform" precipitation is null, then
               ! we use low-level air temperature to specify
               ! the phase of precip.

            else if (T(I) >= TCDK) then
               rainrate(I) = totrate(I)
               snowrate(I) = 0.0
            else if (T(I) < TCDK) then
               rainrate(I) = 0.
               snowrate(I) = totrate(I)
            end if

         end do

      end if IF_OPTION
      !----------------------------------------------------------------
      return
   end subroutine surf_precip3

end module surf_precip
