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
   use phy_options, only: stcond, pcptype, mpdiag_for_sfc
   use modi_wetbulbt, only: wetbulbt
   use modi_hydromet_temp, only: hydromet_temp
   implicit none
   private
   public :: surf_precip1, surf_precip3

   real, parameter :: EPSILON_4 = 1.E-10
   integer, parameter :: OPT_NIL = 0
   integer, parameter :: OPT_CONSUN = 1
   integer, parameter :: OPT_OTHER = 2
   integer, parameter :: OPT_SPS_W19 = 3
   integer, parameter :: OPT_SPS_FRC = 4
   integer, parameter :: OPT_SPS_H13 = 5

contains
 
   !/@*
   subroutine surf_precip1(tt,  hu, pp, tlc, tls, tsc, tss, &
        rainfrac, snowfrac, frfrac, pefrac, &
        rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! tt       low-level air temperature
      ! hu       low-level air specific humidity
      ! pp       low-level air pressure
      ! tlc      total liquid "convective" precipitation rate
      ! tls      total liquid "stratiform" precipitation rate
      ! tcs      total solid  "convective" precipitation rate
      ! tss      total solid  "stratiform" precipitation rate
      ! rainfrac fraction of precip. falling as rain (only used for SPS_FRC)
      ! snowfrac fraction of precip. falling as snow (only used for SPS_FRC)
      ! frfrac fraction of precip. falling as freezing rain (only used for SPS_FRC)
      ! pefrac fraction of precip. falling as ice pellets (only used for SPS_FRC)
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: tt, tlc, tls, tsc, tss, hu, pp
      real, dimension(:), pointer :: rainfrac, snowfrac, frfrac, pefrac  !# intent(in)
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
      ! 006      V. Vionnet  (July 2021)
      !          Add new options using near-surface 
      !          met variables for SPS
      !
      !*@/
 
      integer :: i, option
      real, dimension(ni) :: totrate
      !----------------------------------------------------------------

      ! Sum of the precipitation rate at the surface
      do i=1,ni
         rainrate(i) = 0.
         snowrate(i) = 0.
         totrate (i) = tlc(i)+tls(i)+tsc(i)+tss(i)
      end do

      ! OPTiONS:  Depending on the explicit scheme used during the integration,
      !           we have five options for specifying the phase of 
      !           precipitation reaching the ground

      if (stcond == 'NIL')  then
         option = OPT_NIL
         if (pcptype == 'SPS_W19') option = OPT_SPS_W19
         if (pcptype == 'SPS_FRC') option = OPT_SPS_FRC
         if (pcptype == 'SPS_H13') option = OPT_SPS_H13
      else if (stcond == 'CONSUN') then
         option = OPT_CONSUN
      else
         option = OPT_OTHER
      endif

      select case (option)
      case (OPT_NIL) 
         call surf_precip_nil(tt, totrate, rainrate, snowrate, ni)
      case (OPT_CONSUN)
         call surf_precip_consun(tt, tls, tss, totrate, rainrate, snowrate, ni)
      case (OPT_OTHER)
         call surf_precip_other(tt, tls, tss, totrate, rainrate, snowrate, ni)
      case (OPT_SPS_W19)
         call surf_precip_sps_w19(tt, hu, pp, totrate, rainrate, snowrate, ni)
      case (OPT_SPS_H13)
         call surf_precip_sps_h13(tt, hu, pp, totrate, rainrate, snowrate, ni)
      case (OPT_SPS_FRC)
         call surf_precip_sps_frc(tt, hu, pp, totrate, &
              rainfrac, snowfrac, frfrac, pefrac, &
              rainrate, snowrate, ni)
      end select
      !----------------------------------------------------------------
      return
   end subroutine surf_precip1
 
 
   !/@*
   subroutine surf_precip3(tt, tlc, tls, tsc, tss, &
        rainrate, snowrate, fneige, fip,ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! tt       low-level air temperature
      ! tlc      total liquid "convective" precipitation rate
      ! tls      total liquid "stratiform" precipitation rate
      ! tsc      total solid  "convective" precipitation rate
      ! tss      total solid  "stratiform" precipitation rate
      ! fneige   snow fraction from Bourge1 subroutine
      ! fip      Fraction of liquid precipitation in the form of ice pellets.
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: tt, tlc, tls, tsc, tss, fneige, fip
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
      !          elimination of iSTCOND=2,6,7,8 iCONVEC=4
      ! 005      D. Talbot (March 2005)
      !          Modification to the previous correction
      !          so that it does not apply to consun
      ! 006      A-M.  Leduc (April 2006). Add tlcs and tscs to the totrate.
      !                      Change the option 1 to use FNEiGE. CHANGE name to
      !                      SURF_PRECiP2.
      ! 007     A-M. Leduc. A. Plante. (Nov 2007) add FiP. change s/r name to surf_precip3
      !                change calculation of snowrate to include all solid precipitation
      !                with f_solid.
      !*@/

      integer :: i
      real, dimension(ni) :: totrate
      !----------------------------------------------------------------

      ! Sum of the precipitation rate at the surface
      do i=1,ni
         rainrate(i) = 0.
         snowrate(i) = 0.
         totrate (i) = tlc(i)+tls(i)+tsc(i)+tss(i)
      end do

      ! OPTiONS:  Depending on the explicit scheme used during the integration,
      !           we have three options for specifying the phase of 
      !           precipitation reaching the ground
      select case (stcond(1:3))
      case ('NIL') 
         call surf_precip_nil(tt, totrate, rainrate, snowrate, ni)
      case ('CON')
         call surf_precip3_consun(totrate, fneige, fip, rainrate, snowrate, ni)
      case ('MP_')
         if (mpdiag_for_sfc) then
            call surf_precip3_consun(totrate, fneige, fip, rainrate, snowrate, ni)
         else
            call surf_precip_other(tt, tls, tss, totrate, rainrate, snowrate, ni)
         endif
      case DEFAULT
         call surf_precip_other(tt, tls, tss, totrate, rainrate, snowrate, ni)
      end select
      !----------------------------------------------------------------
      return
   end subroutine surf_precip3

   !---- Private functions --------------------------------------------


   !/@*
   subroutine surf_precip_nil(tt, totrate, rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! tt       low-level air temperature
      ! totrate  total precipitation rate
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: tt, totrate
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Notes
      !  The phase of the precipitation reaching
      !  the surface is simply determined from
      !  the low-level air temperature:
      !    tt > 0.  then rain
      !    tt < 0.  then snow
      !*@/
      integer :: i
      !----------------------------------------------------------------
      do i=1,ni
         if (tt(i) >= TCDK) then
            rainrate(i) = totrate(i)
            snowrate(i) = 0.
         else if (tt(i) < TCDK) then
            rainrate(i) = 0.
            snowrate(i) = totrate(i)
         end if
      end do
      !----------------------------------------------------------------
      return
   end subroutine surf_precip_nil


   !/@*
   subroutine surf_precip_consun(tt, tls, tss, totrate, rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! tt       low-level air temperature
      ! tls      total liquid "stratiform" precipitation rate
      ! tss      total solid  "stratiform" precipitation rate
      ! totrate  total precipitation rate
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: tt, tls, tss, totrate
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Notes
      !   The phase of precipitation at the surface
      !   is determined by results from the "explicit"
      !   condensation scheme (CONSUN).
      !   it was noted in GEM15 that OPTiON 2 was not a viable option
      !   for consun since the scheme produces a non zero tls
      !   in cold temperature in winter causing an erronous
      !   warm surface temperature bias especially in the arctic.
      !*@/
      integer :: i
      !----------------------------------------------------------------
      do i=1,ni

         ! if "stratiform" (stable condensation) is greater
         ! than 0, then the explicit scheme determines if
         ! the precipitation is rain or snow

         if (tss(i) >= EPSiLON_4) then
            rainrate(i) = 0.0
            snowrate(i) = totrate(i)
         else if (tls(i) >= EPSiLON_4) then
            rainrate(i) = totrate(i)
            snowrate(i) = 0.0

            ! if "stratiform" precipitation is null, then
            ! we use low-level air temperature to specify
            ! the phase of precip.

         else if (tt(i) >= TCDK) then
            rainrate(i) = totrate(i)
            snowrate(i) = 0.0
         else if (tt(i) < TCDK) then
            rainrate(i) = 0.
            snowrate(i) = totrate(i)
         end if

      end do
      !----------------------------------------------------------------
      return
   end subroutine surf_precip_consun

   
   !/@*
   subroutine surf_precip3_consun(totrate, fneige, fip, rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! totrate  total precipitation rate
      ! fneige   snow fraction from Bourge1 subroutine
      ! fip      Fraction of liquid precipitation in the form of ice pellets.
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: totrate, fneige, fip
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Notes
      !   The phase of precipitation at the surface
      !   is determined using the bourge snow fraction
      !   fneige and ice pellet fraction fip.
      !*@/
      integer :: i
      real :: f_solid
      !----------------------------------------------------------------
         do i=1,ni
            f_solid     = fneige(i) + (1.-fneige(i))* fip(i)
            snowrate(i) = totrate(i) * f_solid
            rainrate(i) = totrate(i) * (1.-f_solid)
         end do
      !----------------------------------------------------------------
      return
   end subroutine surf_precip3_consun
   

   !/@*
   subroutine surf_precip_other(tt, tls, tss, totrate, rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! tt       low-level air temperature
      ! tls      total liquid "stratiform" precipitation rate
      ! tss      total solid  "stratiform" precipitation rate
      ! totrate  total precipitation rate
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: tt, tls, tss, totrate
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Notes
      !   The phase of precipitation at the surface
      !   is determined by results from the "explicit"
      !   condensation scheme (MiXED-PHASE, EXMO, KONG-YAU).
      !   it was noted in the LAM2.5 that KONG-YAU
      !   could produce a non zero tss over Quebec in summer
      !   and therefore OPTiON 1 was no longer appropriate
      !*@/
      integer :: i
      real, dimension(ni) :: expli
      !----------------------------------------------------------------
      do i=1,ni
         expli(i) = tls(i) + tss(i)

         ! if "stratiform" (stable condensation) is greater
         ! than 0, then the explicit scheme determines if
         ! the precipitation is rain or snow

         if (expli(i) >= EPSiLON_4) then
            rainrate(i) = tls(i)/expli(i) * totrate(i)
            snowrate(i) = tss(i)/expli(i) * totrate(i)

            ! if "stratiform" precipitation is null, then
            ! we use low-level air temperature to specify
            ! the phase of precip.

         else if (tt(i) >= TCDK) then
            rainrate(i) = totrate(i)
            snowrate(i) = 0.0
         else if (tt(i) < TCDK) then
            rainrate(i) = 0.
            snowrate(i) = totrate(i)
         end if

      end do
      !----------------------------------------------------------------
      return
   end subroutine surf_precip_other


   !/@*
   subroutine surf_precip_sps_w19(tt, hu, pp, totrate, rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! tt       low-level air temperature
      ! hu       low-level air specific humidity
      ! pp       low-level air pressure
      ! totrate  total precipitation rate
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: tt, hu, pp, totrate
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Notes
      !   Rain-snow partitionong scheme based on wet-bulb temperature 
      !   developped by Wang et al. (GRL, 2019)
      !   To be used only for SPS
      !*@/
      integer :: i
      real :: za, zb, zc
      real, dimension(ni) :: tw, prain, ta
      !----------------------------------------------------------------
      ! Constant for the SPS_W19 options
      za = 6.99e-5  ! [-]
      zb = 2.0  ! [K-1]
      zc = 3.97 ! [K]

      ! Temperature in C
      ta(:) = tt(:) - TCDK

      ! Compute wet-bulb temperature (in deg C) using an iterative method 

      tw(:) = wetbulbt(pp(:), ta(:), hu(:)) 

      do i=1,ni           
         ! Compute rain fraction
         !
         ! Note that the formulation for the snow fraction differs from Eq (1) in the W19 paper 
         ! We use here:
         !                           1.
         !    fs = -----------------------------
         !          1.+ ZA*EXP{ ZB * (TW + ZC)} 
         ! with TW in deg C
         ! This formulation gives the same solution as the one shown on Fig 1b 
         ! of Wang et al. (2019)
         if (tw(i) > 5.) then
            prain(i) =  1.
         else if (tw(i) < -5.) then
            prain(i) = 0.
         else
            prain(i) = 1. - 1./(1.+za*exp(zb*(tw(i) + zc)))
         endif

         if (prain(i) < 0.01) prain(i) = 0.   ! remove rain fraction less that 1 percent
         if (prain(i) > 0.99) prain(i) = 1.   ! remove snow fraction less that 1 percent

         ! Compute rain and snow rates
         rainrate(i) = prain(i)*totrate(i)
         snowrate(i) = max(0.,(1. - prain(i))*totrate(i))
      end do
      !----------------------------------------------------------------
      return
   end subroutine surf_precip_sps_w19
  
   subroutine surf_precip_sps_h13(tt, hu, pp, totrate, rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! tt       low-level air temperature
      ! hu       low-level air specific humidity
      ! pp       low-level air pressure
      ! totrate  total precipitation rate
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: tt, hu, pp, totrate
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Notes
      !   Rain-snow partitionong scheme based on hydrometeor temperature 
      !   developped by Harder and Pomeroy (HP, 2013)
      !   To be used only for SPS
      !*@/
      integer :: i
      real :: zd, ze 
      real, dimension(ni) :: ti, prain, ta
      !----------------------------------------------------------------
      !      Constant for the SPS_H13 options
      zd = 2.630006 ! [-]
      ze = 0.09336 ! [-]

      ! Temperature in C
      ta(:) = tt(:) - TCDK

      ! Compute hydrometeor temperature (in deg C) using an iterative method 
      ti(:) =  hydromet_temp(pp(:), ta(:), hu(:))       


      do i=1,ni         

         ! Compute rain fraction 
         !               1.
         !    fr = -----------------
         !          1.+ ZD* ZE^TI  
         ! with TI in deg C in the range [-4, 4]

         if(ti(i) > 4) then
             prain(i) =  1.
         else if(ti(i)<-4) then
             prain(i) =  0.
         else
             prain(i) =  1./(1.+zd* ze**(ti(i)))
         endif

         if (prain(i) < 0.01) prain(i) = 0.   ! remove rain fraction less that 1 percent
         if (prain(i) > 0.99) prain(i) = 1.   ! remove snow fraction less that 1 percent

         ! Compute rain and snow rates
         rainrate(i) = prain(i)*totrate(i)
         snowrate(i) = max(0.,(1. - prain(i))*totrate(i))
      end do
      !----------------------------------------------------------------
      return
   end subroutine surf_precip_sps_h13
   !/@*
   subroutine surf_precip_sps_frc(tt, hu, pp, totrate, &
        rainfrac, snowfrac, frfrac, pefrac, &
        rainrate, snowrate, ni)
      implicit none
!!!#include <arch_specific.hf>
      !@Object Determine the phase of precipitation reaching the surface
      !@Arguments
      !          - Input -
      ! tt       low-level air temperature
      ! hu       low-level air specific humidity
      ! pp       low-level air pressure
      ! totrate  total precipitation rate
      ! rainfrac fraction of precip. falling as rain (only used for SPS_FRC)
      ! snowfrac fraction of precip. falling as snow (only used for SPS_FRC)
      ! frfrac   fraction of precip. falling as freezing rain (only used for SPS_FRC)
      ! pefrac   fraction of precip. falling as ice pellets (only used for SPS_FRC
      !          - Output -
      ! rainrate rate of precipitation at the surface (liquid)
      ! snowrate rate of precipitation at the surface (solid)

      integer, intent(in) :: ni
      real, dimension(ni), intent(in) :: tt, hu, pp, totrate
      real, dimension(:), pointer :: rainfrac, snowfrac, frfrac, pefrac  !# intent(in)
      real, dimension(ni), intent(out) :: rainrate, snowrate

      !@Notes
      !   Rain/snow partitioning is provided directly in the meteorological forcing used to drive SPS.
      !   Option useful to use the phase separation from an atmospheric model to drive SPS 
      !   See Vionnet et al (2022, WRR) for examples
      !   The partioning of Wang et al (2019) is used as a backup if the fraction is not defined
      !*@/
      integer :: i
      real :: za, zb, zc
      real, dimension(ni) :: tw, prain, ta
      !----------------------------------------------------------------
      ! Constant for the SPS_W19 options
      za = 6.99e-5  ! [-]
      zb = 2.0      ! [K-1]
      zc = 3.97     ! [K]

      ! Temperature in C
      ta(:) = tt(:) - TCDK

      ! Compute wet-bulb temperature for backup
      tw(:) = wetbulbt(pp(:), ta(:), hu(:))

      do i=1,ni   

         rainrate(i) = rainfrac(i)*totrate(i)

         ! So far freezing rain and ice pellets are considerd as solid precipitation since the land surface schame is not 
         ! equipped to deal specifically with freezing rain and ice pellets. 
         snowrate(i) = (snowfrac(i)+frfrac(i)+pefrac(i))*totrate(i)

         if (rainfrac(i)+snowfrac(i)+frfrac(i)+pefrac(i) < 0.99 .and. totrate(i) > 0.) then
            !
            ! the fractions are not defined properly for this grid point. in this case, use w19 
            ! to compute the rain/snow partitioning. 

            prain(i) = 1. - 1./(1.+za*exp(zb*(tw(i) + zc)))

            if (prain(i) < 0.01) prain(i) = 0.   ! remove rain rate less that 1 percent
            if (prain(i) > 0.99) prain(i) = 1.   ! remove snow rate less that 1 percent

            ! compute rain and snow rates
            rainrate(i) = prain(i)*totrate(i)
            snowrate(i) = max(0.,(1. - prain(i))*totrate(i))
         endif

      end do
      !----------------------------------------------------------------
      return
   end subroutine surf_precip_sps_frc

end module surf_precip
