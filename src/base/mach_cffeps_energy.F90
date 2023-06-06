!============================================================================!
!
! Projet / Project : GEM-MACH
! Fichier / File   : mach_cffeps_energy.ftn90
! Creation         : Kerry Anderson
!                  : (C To Fortran)
!                    A. Akingunola, J. Chen, and P. Makar - Fall 2018
! Description      : Energy Balance of a Fire required for Plume Rise
!                    Anderson, K.R.; Pankratz, A; Mooney, C. 2011.
!                    A thermodynamic approach to estimating smoke plume heights.
!                    In 9th Symp. on Fire and Forest Meteorology,
!                    Oct 18-20, 2011.  Palm Springs, CA.
!                    Am. Meteorol. Soc., Boston, MS.
!
!                    Bulk densities for FBP fuel types are based on;
!
!                  Anderson, K.R.  2000.
!                  Incorporating smoldering into fire growth modelling.
!                  Pages 31-36 in 3rd Symp. on Fire and Forest Meteorology,
!                  Jan. 9-14, 2000, Long Beach, CA. Am. Meteorol. Soc., Boston, MS.
!============================================================================
!!if_on
 subroutine mach_cffeps_energy(fbp, ifuel, qplume, tsurf, area, perimeter, &
                               growth, fmc, mcL, mcF, mcH)
   use mach_cffeps_mod,           only: fbp_type
!!if_off
   use mach_cffeps_mod,           only: fuel_size, cffeps_sinks,           &
                                        cffeps_timestep, cffeps_radiation, &
                                        cffeps_fueltypes
   use mach_cffeps_emissions_mod, only: rhoBs, rho2_consts, &
                                        rho4_consts, rho6_consts, rho8_consts
   implicit none
!!if_on
   type(fbp_type),  intent(in)  :: fbp
   real(kind=4),    intent(out) :: qplume
   real(kind=4),    intent(in)  :: tsurf
   real(kind=4),    intent(in)  :: mcL, mcF, mcH
   real(kind=4),    intent(in)  :: area, perimeter, growth, fmc
   integer(kind=4), intent(in)  :: ifuel
!!if_off
   real(kind=4), parameter  :: cpw = 4.1855 * 1.0e3 ! /* convert J g^-1 K^-1 to J kg-1 K-1 */,
   real(kind=4), parameter  :: cpwood = 1.7 * 1.0e3 !  /* J kg-1 K-1 */,
   real(kind=4), parameter  :: lv = 2.501E6 !/* J kg^-1 */,
                                    ! heat of combustion [J kg^-1] */
   real(kind=4), parameter  :: H = 18000.0 * 1.0e3 !/* convert kJ kg^-1 to J kg^-1 */
   real(kind=4), parameter  :: pi = 3.1415926, stefan = 5.67E-8, tcdk = 273.16

   real(kind=4)  :: As, Aw, Ag, rs, mw, &
                    qfire, qw, qf, qs, qr, qinc, fi, dT, tfire, h2, &
                    rho2, rho4, rho6, rho8, dob, work
   character(len=FUEL_SIZE) :: fueltype
!
   fueltype = cffeps_fueltypes(ifuel)
!
!/* Units are converted to SI units. */
   As = 100.0 * 100.0 * area         !/* convert A to m^2 */
   Ag = 100.0 * 100.0 * growth       !/* convert area growth to m^2 */

!/* use the fire area in the case where growth is not used */
!   if (Ag <= 0.) Ag = As
   if (Ag < 0.0) Ag = As

!/********************************************************************/

!/* Calculate Qfire, Qw, Qf, Qs, Qr, Qinc and from these Qplume  */

   qfire = H * fbp%tfc * Ag       ! /* eqn 2 - based on area growth */

   if (cffeps_sinks .and. cffeps_radiation <= 3) then

!/* Assign default duff characteristics by fuel type based on the literature. */
!/* Use values from Anderson 2000 in g/cm^3 */
      select case (trim(fueltype))
         case('C1', 'C2', 'C3', 'C4', 'C6', 'C7', 'D1', &
              'M1', 'M2', 'M3', 'M4', 'S1', 'S2', 'S3')
            rho2 = rho2_consts(ifuel)
            rho4 = rho4_consts(ifuel)
            rho6 = rho6_consts(ifuel)
            rho8 = rho8_consts(ifuel)
         case default
            rho2 = rhoBs(ifuel)
            rho4 = rhoBs(ifuel)
            rho6 = rhoBs(ifuel)
            rho8 = rhoBs(ifuel)
      end select

!/*    Calculate depth of burn and the duff consumption from the surface fuel
!      consumption (SFC) and densities */
      if (trim(fueltype) == "O1a" .or. trim(fueltype) == "O1b") then
         mw = mcL * fbp%sfc   ! fbp%sfc for the fuel type is grass fuel load, GFL
                              ! as defined in cffeps_FBPCalc
         dob = 0.0
      else
         work = rho2 * 20.0 + rho4 * 20.0
         if (fbp%sfc < (2.0 * rho2 * 10.0)) then
              !/* rho*10 give mass per square metre per cm of depth;
              !   rho*20. gives mass per 2cm layer */
            dob = fbp%sfc / (rho2 * 10.0)
         else if (fbp%sfc < work) then
            dob = 2.0 + (fbp%sfc - rho2 * 20.0) / (rho4 * 10.0)
         else if (fbp%sfc < (work + rho6 * 20.0)) then
            dob = 4.0 + (fbp%sfc - work) / (rho6 * 10.0)
         else
            dob = 6.0 + (fbp%sfc - work - rho6 * 20.0) / (rho8 * 10.0)
         end if

!/*    Water mass is calculated from the moisture content as described by
!      the FFMC, (top 2 cm), DMC (2-5 cm), DC (5+ cm) and the depth of burn.
!      Add the water mass from the crown fuel consumption (TFC-SFC) and the
!      foliar moisture content (FMC). */
!
         if (dob < 2.0) then
            mw = mcL * dob * 0.5 * (rho2 * 20.0) ! Mw = mass of water
         else
            mw = mcL * rho2 * 20.0

            if (dob < 4.0) then
               mw = mw + mcF * (dob - 2.0) * (rho4 * 10.0)
            else
               mw = mw + mcF * (rho4 * 20.0)
               if (dob < 6.0) then
                  mw = mw + mcF * (dob - 4.0) * (rho6 * 10.0)
               else
                  mw = mw + mcF * (rho6 * 20.0)
                  mw = mw + mcH * (dob - 6.0) * (rho8 * 10.0)
               end if
            end if
         end if
      end if

      if (fbp%tfc > fbp%sfc) mw = mw + fmc * 1.0e-2 * (fbp%tfc - fbp%sfc)

      dT = 373.16 - tsurf
      ! /* energy that goes into heating fuel moisture and evaporating it */
      qw = mw * (lv + cpw * dT) * Ag

      !/* new value for heating temperature for fuel (~600 K) based on
      !   Albini et al 1995... but this is just the ignition temperature! */
!//    dT = 600.0 - tsurf
! /* newer value for heating temperature for fuel (mid point between
!    600 K and 777oC) also based on Albini et al 1995 */
      dT = (551.92 + tcdk) - tsurf

      !/* energy that goes into heating fuel to 500oC */
      qf = fbp%tfc * cpwood * dT * Ag

    !/* new approach assuming that 50% amount of heat equal to surface fuel
    !   consumption is injected into surface */
      qs = 0.5 * H * fbp%sfc * Ag

   ! /* assumption that the percent incomplete combustion is half of the
   !    Crown Fraction Burned */
      fi = fbp%cfb * 0.5

!/* Radiation method 1 (default scheme) =- from Byram 1959, 1973: page 161
!   - assumption that 14% (17%?) lost due to radiation */
      qr  = qfire * 1200.0 / 8600.0

!/* Radiation method 2 - Flame height to size of fire */
      if ((cffeps_radiation == 0) .or. (cffeps_radiation == 2)) then ! Method 2
      !  /* radius - eqn 18 (recalculated for case when Qplume was known */
         rs = sqrt(As / pi)

         h2 = sqrt(fbp%hfi / 300.0) !/* Byram (1959, 1973 pg 175) */
         Aw = 2.0 * pi * rs * h2     !/* area of the wall surrounding the fire */

         if (Aw < As .and. As > 0.0) qr = qfire * Aw / As * 0.5
      end if

      if (cffeps_radiation == 3) then !// Method 3
      ! // new thinking... back radiation equal the forward radiation (assuming
      !    a steady state fire)
!         Qr  = Qf + Qw
         h2 = sqrt(fbp%hfi / 300.0) !/* Byram (1959, 1973 pg 175) */
 ! /* http://wildfiretoday.com/2011/02/26/at-what-temperature-does-a-forest-fire-burn/ */
         tfire = tcdk + 800.0 + 400.0 * fbp%cfb  !/* [K] */
        !  /* assumption Qr is the backward radiation, which is based on the
        !     flame temperature */
         qr = h2 * perimeter * stefan * tfire**4 * cffeps_timestep * 3600.0
      end if
!
     !/* relevant only for methods 1 and 2 */
     !/* Qinc is a weak link in the energy budget model */
     !// new approach for Qinc - 10/08/2015 (Stull's eqn occasionally produced
     !   negative values)
!      Qinc = Fi * (Qf + Qw)
     ! newer approach for Qinc - only energy of fuel burned is involved
     !                           in incomplete combustion - 04-06-2018
      qinc = fi * qf

!      Qplume = (1.-Fi)*(Qfire - Qw - Qf - Qr - Qs) ! old equation
       ! corrected equation based on discussions with Roland Stull and
       ! Rosie Howard - 05/28/2015
!      Qplume = Qfire - (1 + Fi) * (Qw + Qf) - Qr - Qs
      qplume = qfire - qw - qf - qr - qs - qinc
!
!/* Radiation method 4 - use the value as a percentage of energy of the fire
!   that enters the plume */
   else if (cffeps_radiation > 3) then
      qplume = real(cffeps_radiation) * 1.0e-2 * qfire
!
   else
      qplume = qfire
   end if

   qplume = max(0.0, qplume)

   return

 end subroutine mach_cffeps_energy
