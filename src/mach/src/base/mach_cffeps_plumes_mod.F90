!============================================================================!
!
! Projet / Project : GEM-MACH
! Fichier / File   : mach_cffeps_plume_mod.ftn90
! Creation         : Kerry Anderson 
!                  : (C To Fortran)
!                    A. Akingunola, J. Chen, and P. Makar - Fall 2018
! Description      : Plume model based on;
!                    Anderson, K.R.; Pankratz, A; Mooney, C. 2011.   
!                    A thermodynamic approach to estimating smoke plume heights.
!                    In 9th Symp. on Fire and Forest Meteorology, 
!                    Oct 18-20, 2011.  Palm Springs, CA.
!                    Am. Meteorol. Soc., Boston, MS.
!
!============================================================================
module mach_cffeps_plume_mod

   use mach_cffeps_mod,  only: MET_LEVELS, fuel_size, met_type,          &
                               cffeps_fire_shape, cffeps_timestep,       &
                               cffeps_alpha, cffeps_sinks, msl_pressure, &
                               cffeps_fire_type
   implicit none
   public  :: cffeps_plumecalc, cffeps_plumecalc2, rh_calc
   private :: r_calc
   
  contains
  
   subroutine cffeps_plumecalc(area, Qt, tsurf, lapse_rate, perimeter, &
                               dz, mplume, Qbb)
      implicit none
       
      real(kind=4),    intent   (in) :: area      
      real(kind=4),    intent   (in) :: Qt    
      real(kind=4),    intent   (in) :: tsurf      
      real(kind=4),    intent   (in) :: lapse_rate      
      real(kind=4),    intent   (in) :: perimeter      
      real(kind=4),    intent  (out) :: dz      
      real(kind=4),    intent  (out) :: mplume      
      real(kind=4),    intent  (out) :: Qbb      
      
      integer(kind=4), parameter :: imax = 10000
      real(kind=4), parameter    :: grav = 9.80025, cpd = 1005.0,       &
                                    rgasd = 287.05, Ld = -9.8 / 1000.0, &
                                    pi = 3.1415926, stefan = 5.67E-8
                                   
      real(kind=4)  :: As, At, d, rs, rt, rho      
      real(kind=4)  :: Qatm, q, ddz, ps, Mo, Le, work, Tm, Tt
!      real(kind=4)  :: pt, rhot, rhos, Mt, Ms      
      integer(kind=4) :: ii
!         
!/* 1. Units are converted to SI units. */
      !/* lapse rate limited to -9.0 oC/km */
      Le = max(-0.009, lapse_rate)  !/* Environmental lapse rate [oC m^-1] */
      ps = msl_pressure * 100.0     !/* convert ps to [pa] */
      As = 100.0 * 100.0 * area     !/* convert A to [m^2] */

!/* I don't use Area growth (per timestep) for anything in this routine */
      if ((trim(cffeps_fire_shape) == "line" .or.     &
           trim(cffeps_fire_shape) == "cresent") .and. (perimeter > 0.0)) then
         d = As / perimeter  ! /* forward spread distance [m] */
      else
         d = 0.0
      end if

!/*****  Testing this section *****/

!/*  This approach turns out to be wrong
!    as the plume height varies with the time step.  */

!//        As = 100.0 * 100.0 * growth !/* convert area growth to [m^2] */

!/*  By keeping the plume volume based on the cumulative fire size,
!    the answer is stable regadless of time step.
!
!    Note that energy is still calculated by area growth over time */
!
! /* radius - eqn 18 (recalculated for case when Qplume was known */
      rs = sqrt(As / pi)
!
!/********************************************************************/
!/* 2.
      dz = 1000.0    ! // first guest at height
      ddz = 1000.0   ! // step used in the iteractive process

      ii = 0
      do while (ii < imax .and. abs(ddz) >= 1.0)
         rt = rs + dz * tan(cffeps_alpha)  !          /* eqn 19 */
         At = pi * rt * rt

!/* 3.a. Calculate the energy per unit mass needed to heat plume to dry adiabat*/
         q = -.5 * cpd * Ld * dz * log(1.0 + dz * (Le - Ld) / tsurf)  
                                                           !/* eqn 8 */

!/* 3.b. Calculate the column mass of the plume */
         if (Le /= 0.0) then
            work = (1.0 + Le * dz / tsurf)**(-grav / Le / rgasd)
         else !   // in the isothermal case, take the average using Le+.5 and Le-.5
            work = ((1.0 + Le * dz / tsurf)**(-grav / (Le + 0.5) / rgasd) +  &
                    (1.0 + Le * dz / tsurf)**(-grav / (Le - 0.5) / rgasd)) * 0.5
         end if
         mplume = ps / grav * As * (1.0 - work) !/* eqn 16 */
               
         Mo = mplume   !// Mo is mass without entrainment

!/* 3.c. Calculate air density in plume */
!            !// new mass weighted scheme to determine rho
!         pt = ps * work                         !/* eqn 14 */
!         Mt = (ps - pt) / grav * At             !/* eqn 15 */
!         Ms = (ps - pt) / grav * As             !/* eqn 15 */
!         if ((Mt + Ms) > 0.0) then
!            Tt = tsurf + Le * dz        !// plume top temperature
!            rhot = pt / (Tt * rgasd)    !// ideal gas law
!            rhos = ps / (tsurf * rgasd) !// ideal gas law
!
!            rho = (Mt * rhot + Ms * rhos) / (Mt + Ms) 
!         else
!            rho = 0.0
!         end if

        !   // original scheme (need to use this to agree with
        !      paper exercise in FFM9 paper)
         if ((dz > 0.0) .and. (As > 0.0)) then
            rho = Mo / (As * dz) 
         else
            rho = 0.0 !   // and later on M will be zero
         end if

!/* 3.d. Adjust mass for entrainment */
         if (cffeps_alpha > 0.0) then ! /* include entrainment */
                    ! /* eqn 20 */
            mplume = 1.0 / 3.0 * pi * rho * dz * (rt**2 + rs * rt + rs**2) 
         end if

         if (trim(cffeps_fire_shape) == "line" .or. &
             trim(cffeps_fire_shape) == "cresent") then
            !        /* this calculation is for the volume of a wedge -
            !           no allowance for curvature */
            !        /* For "cresent" shape, use a wedge following the 
            !            burning perimeter (estimated in CFFEPS) */
            mplume = rho * (d * dz * perimeter +  &
                       perimeter * dz * dz * tan(cffeps_alpha))
         end if

!/* 3.e. Calculate total energy needed to heat plume to dry adiabat */
         Qatm = q * mplume                  !  /* eqn 17 */

!/* 3.f. If energy required to heat the air mass (Qatm) is greater than 
!        plume energy (Qt), halve the height-step size and reduce the height, 
!        otherwise increase the height */
         if (Qatm > Qt) then
            ! // overshot the top, reduce the height-step size
            ddz = 0.5 * ddz
            dz = dz - ddz
         else
            dz = dz + ddz
         end if

         ii = ii + 1
!/* 3.g. Repeat calculations until the height-step size is less than a metre */
            !/* expect this to converge to within less than a metre */
      end do

!/* 4. Calculate the black body radiation loss Qbb. */
      if (cffeps_sinks) then
         Tt = tsurf + Le * dz        ! // plume top temperature
         Tm = Tsurf + (Le - Ld) * dz ! // modified surface temperature
         Qbb = cffeps_timestep * 3600.0 * dz * perimeter * stefan * &
               (((Tm + Tt) * 0.5)**4 - ((tsurf + Tt) * 0.5)**4)
      else
         Qbb = 0.0
      end if
         
      return
   end subroutine cffeps_plumecalc
!
!============================================================================
!
   subroutine cffeps_plumecalc2(met, rh0, area, Qt, mw0, perimeter, dz, mplume, Qbb)
!/*     This routine allows for a piecewise integration of an upper air profile
!       to calculate plume rise.
!    It has been tested under the ICAO standard atmosphere and is within 1% of 
!    the PlumeCalc routine.
!
!    With that said, there appears to be an error in the calculations that I 
!    cannot trace -- apparently in energy calculations (q) and specifically in 
!    the temperature calculations (T1, T2, T3, T4) as absolute temperature 
!    ratio ln(theta1/theta2) is handled consistently.  Mass (M) is correct.
!    The calculation does not precisely measure the area in the trapezoid and,
!    as a result, the addition of two adjacent trapazoids does not equal the 
!    area of the larger trapezoid.  The error tends to propagate through the 
!    calculations resulting in higher plumes with more piecewise calculations.
!    The error is small (tens of meters) but is annoying -- KRA 2018-05-16 */
!
      implicit none
       
      type(met_type),  intent(in)    :: met
      real(kind=4),    intent(in)    :: rh0      
      real(kind=4),    intent(in)    :: area      
      real(kind=4),    intent(in)    :: Qt     
      real(kind=4),    intent(in)    :: mw0     
      real(kind=4),    intent(in)    :: perimeter      
      real(kind=4),    intent(out)   :: dz      
      real(kind=4),    intent(out)   :: mplume      
      real(kind=4),    intent(out)   :: Qbb      
      
      real(kind=4), parameter  :: grav = 9.80665, cpd = 1005.0,       &
                                  rgasd = 287.05, Ld = -9.8 / 1000.0, &
                                  pi = 3.1415926, stefan = 5.67E-8,   &
                                  eps = 0.622, l_v = 2.5008E6
                                   
      real(kind=4)  :: As, At, rs, rt, rho      
      real(kind=4)  :: Qatm, q, qa, qb, qc, qd, Mo
      real(kind=4)  :: Le, T0, Tm, Tt, d, theta1, theta2, theta5, & 
                       Mc, T1, T2, T3, T4, T5, Td, fQ, pQatm, &
                       P0_mb, P1, P2, Z0, Z1, Z2, ZZ, ddz, ddP, ddT
      real(kind=4)  :: Qw, M1, M2, Td0, Tw1, Tw2, thetaP,     &
                       LCL, hu0, hu_sat, hu_new, rh, moist_ad_lr
      integer(kind=4) :: i, iConverge
      
!/* 1. Units are converted to SI units. */
      As = 100.0 * 100.0 * area      !/* convert A to [m^2] */
      ZZ = 1.0                       ! /* the vertical step used in converging [m] */

!/* I don't use Area growth (per timestep) for anything in this routine */
      if ((trim(cffeps_fire_shape) == "line" .or.     &
           trim(cffeps_fire_shape) == "cresent") .and. (perimeter > 0.0)) then
         d = As / perimeter  ! /* forward spread distance [m] */
      else
         d = 0.0
      end if

!/*****  Testing this section *****/

!/*  This approach turns out to be wrong
!    as the plume height varies with the time step.  */

!//        As = 100.0 * 100.0 * growth !/* convert area growth to [m^2] */

!/*  By keeping the plume volume based on the cumulative fire size,
!    the answer is stable regadless of time step.
!
!    Note that energy is still calculated by area growth over time */
!
! /* radius - eqn 18 (recalculated for case when Qplume was known */
      rs = sqrt(As / pi)
!
!/********************************************************************/
!/* 2.   Initialize plume height, plume mass blackbody radiation loss 
      dz = 0.0
      mplume = 0.0
      Qbb = 0.0
      pQatm = 0.0
         
!/* 3. Set array pointers to first upper air data point. */
      T0 = met%T(1)   !  // temperature at/near surface [K]
      Z0 = met%Z(1)   !  // height at/near surface [m] -- constant
      P0_mb = met%ps  !  // pressure at/near surface [mb]

!/* 4. Set cumulative energy per unit mass, cumulative mass and 
!      plume height to zero. */
      qc = 0.0; q = 0.0  ! // cumulative energy
      Mc = 0.0; Mo = 0.0 ! // cumulative mass
      qw = 0.0
      iConverge = 0      ! // track whether top has been reached 
                         !    (Qatm>Qplume) and convergence is required
!
      if (trim(cffeps_fire_type) == "wet") then
         hu0 = met%hus   ! // mixing ratio (specific humidity) at surface [kg/kg]
      ! Estimate the dew point temperature [K] from RH, using CFFEPS formulation
      ! taking from Irabarne and Godson (1973) [VII-6]
         Td0 = 1.0 / (-4.25E-4 * log10(rh0 * 1.0e-2) + 1.0 / T0)
      end if
!
!/* 5. Begin converging on plume height. */
      i = 1
      do 
!/* these calculations are done in [deg K] and [Pa] */
         if (iConverge == 0) then
!/* 5.a.i. Level 2 is the upper air data for the current level */
            T2 = met%T(i)
            Z2 = met%Z(i)
            P2 = met%P(i)

!/* 5.a.ii Level 1 is the upper air data for the next level */
            T1 = met%T(i + 1)
            Z1 = met%Z(i + 1)
            P1 = met%P(i + 1)

            ddZ = Z1 - Z2
            ddP = P1 - P2
            ddT = T1 - T2
!/* correct for geopotential height (plume rise is calculated in gpm) */
!//          Z1 = Z1 * Re*Re/(Z1 + Re)/(Z1+ Re);
!//          Z2 = Z2 * Re*Re/(Z2 + Re)/(Z2+ Re);
!
!/* 5.a.iv. If stepping upward, qc and Mc are updated (when stepping down,
!           qc and MC remain constant) */
             qc = q  !  // old energy
             Mc = Mo !  // Mo is mass without entrainment
         end if
!
         ! // environmental lapse rate in upper layer
         Le = (T1 - T2) / (Z1 - Z2)
         dz = Z1 - Z0
                
         rt = rs + dz * tan(cffeps_alpha)  !   /* eqn 19 */
         At = pi * rt * rt

!/*5.c. Calculate energy required to heat atmosphere to dry adiabat.*/
!   Solution is found by calculating energy contained in trapezoid defined by 
!   (T1, theta1), (T2, theta2), (T3, theta3), (T4, theta4)
!   where the temperatures are at the next level (T1), the current level (T2),
!   while T3 and T4 are those when lowered adiabatically to the 
!   surface (T1 to T4, T2 to T3) */

!/* alternate approach to be consistent with PlumeCalc (note that theta is wrt
!   the surface layer, not 1000 mb) * /
!         theta1 = (T1+273.16) - Ld*(Z1-Z0);
!         theta2 = (T2+273.16) - Ld*(Z2-Z0); */

!/* 5.c.i. Calculate T3 and T4 using dry adiabat */
         T3 = T2 - Ld * (Z2 - Z0)
         T4 = T1 - Ld * (Z1 - Z0)
         T5 = T2 - Le * (Z2 - Z0)

!/* 5.c.ii. Calculate potential temperatures */
         theta1 = T4
         theta2 = T3
         theta5 = T5

!// old strategy, calculate area in trapezoid -- don't know why this is wrong
         ! // use the baseline qc as T1, P1, Z1 is converged on 0
!//       q = qc - .5 * cp * log((theta1)/(theta2)) * (T1+T2-T3-T4)  

!/* 5.c.iv. Calculate area in trapezoid
!       (new strategy, subtracting triangles to calculate area in trapezoid) */
         qd = -0.5 * cpd * log(theta1 / theta5) * (T1 - T4)
         qa = -0.5 * cpd * log(theta2 / theta5) * (T2 - T3)
         qb = qd - qa ! // resulting area of trapezoid

!/* 5.c.v. Add new energy and mass to previous values */
         q = qc + qb
         !/* eqn 15 - Mo is column mass*/
         Mo = Mc + (P2 - P1) / grav * As 

!/* 5.d. Calculate air density in plume. */
         if ((dz > 0.0) .and. (As > 0.0)) then
            !   // original scheme (need to use this to agree with
            !      paper exercise in FFM9 paper)
            rho = Mo / (As * dz) 
         else
            rho = 0.0 !   // and later on M will be zero
         end if

!/* 5.e. Adjust mass for entrainment */
         if (dz == 0.0) then
            mplume = 0.0
         else
            if (cffeps_alpha > 0.0) then ! /* include entrainment */
                               ! /* eqn 20 */
               mplume = 1.0 / 3.0 * pi * rho * dz * (rt**2 + rs * rt + rs**2)
            else
               mplume = Mo
            end if

            if (trim(cffeps_fire_shape) == "line" .or. &
                trim(cffeps_fire_shape) == "cresent") then
               !   /* this calculation is for the volume of a wedge -
               !      no allowance for curvature */
               !   /* For "cresent" shape, use a wedge following the 
               !      burning perimeter (estimated in CFFEPS) */
               mplume = rho * (d * dz * perimeter + &
                          perimeter * dz * dz * tan(cffeps_alpha))
            end if
         end if

         if (trim(cffeps_fire_type) == "wet") then
         
            m2 = mplume  ! // new mass of plume
            
!            mw0 = 0.0    ! setting mw0 to 0 turns off PlumeCalcWetter (returning to PlumeCalcWet)

!/* Adjust mixing ration for moisture released due to fuel moisture and 
!   for chemical decomposition                       
            T4 = max(T4, T0)   ! avoid superadiabatic situations
            !  LCL = 120.0 * (T0 - Td0) => LCL of the base of the cumulus layer [AGL]
            !  LCL = 120.0 * (T4 - Td0) => LCL of the base of the pyrocumulus [AGL]
            LCL = 120.0 * (T4 - Td0)

            if (mw0 > 0.0 .and. m2 > 0.0) then
               hu_new = hu0 + mw0 / m2   ! new mixing ratio [kg/kg] based moisture added from fire
               hu_sat = r_Calc(P0_mb, T4)  ! saturation mixing ratio 
               if (hu_new > hu_sat) hu_new = hu_sat
               rh = rh_calc(hu_new, P0_mb, T4)  ! recalculate surface rh
               rh = min(rh, 100.0)
                !Dew point calculated using CFFEPS' T_d_Calc(T4, rh)
               Td = 1.0 / (-4.25E-4 * log10(rh * 1.0e-2) + 1.0 / T4)
               LCL = 120.0 * (T4 - Td)  !LCL of the base of the pyrocumulus [AGL]
            else
               hu_new = hu0 ! r remains unadjusted
            end if

!/* 5.f. an approximation of condensation */
            if (dz > LCL) then
         ! saturation is reached (note that r0 from surface is originally used)
               Tw2 = T4 + Ld * LCL ! lift parcel from Z0 to LCL using 
                                   ! dry adiabat
                   !  // lift parcel from LCL to Z1 using pseudo-adiabat
!/*** Moist adiabatic lapse rate [K/m], based on; ***/
!/* https://en.wikipedia.org/wiki/Lapse_rate#Moist_adiabatic_lapse_rate */
               moist_ad_lr = (-grav * (rgasd * Tw2**2 + l_v * hu0 * Tw2) / &
                             (cpd * rgasd * Tw2**2 + l_v * l_v * hu0 * eps))
               Tw1 = Tw2 - moist_ad_lr * (Z1 - Z2)
               thetaP = Tw1 - Ld * (Z1 - Z0)
               qw = 0.5 * cpd * log(thetaP / theta1) * (Tw2 - T1)
               rt = rs + LCL * tan(cffeps_alpha)  !  /* eqn 19 */
               m1 = 1.0 / 3.0 * pi * LCL * (rt**2 + rs * rt + rs**2) * rho ! /* 20 */
               qw = qw * (m2 - m1) ! the energy resulting from pseudoadiabatic rise
            else
               qw = 0.0   ! no energy added as plume height is the below the LCL
               m1 = m2
            end if
         end if
         
         Qatm = q * mplume                        !/* eqn 17 */

!/* 5.g. If energy required to heat the air mass (Qatm) is greater than ,
!        plume energy (Qt) begin calculations downward, reducing the top level
!        by 1 metre until Qatm < Qt again, in which case the soluiton is reached. */
         if (iConverge == 0 .and. Qatm < (Qt + Qw) .and. i < MET_LEVELS) then
            i = i + 1
            pQatm = Qatm
            if (trim(cffeps_fire_type) == "wet") then
               m1 = m2  !// previous cumulative mass of plume
            end if
         else if (Qatm <= 0.0) then 
            iConverge = -1

         else if (i >= MET_LEVELS) then  ! maxed out the levels 
            iConverge = -1

! passed equilibrium level
         else if (zz < 0.0) then ! this is the old method that steps down by 1 m 
            if (Z1 + ZZ <= Z2 .or. Qatm < (Qt + Qw)) then
               iConverge = -1
            else   !  ZZ is negative hence +ZZ to step down
               iConverge = 1     ! converge on the answer
               P1 = P2 + (P1 - P2) * (Z1 - Z2 + ZZ) / (Z1 - Z2)
               T1 = T2 + (T1 - T2) * (Z1 - Z2 + ZZ) / (Z1 - Z2)
               Z1 = Z1 + ZZ
            end if
         else if (zz == 0.0) then ! this method simply estimates the value 
                                  ! based on the fractional difference between 
                                  ! previous calculation and current (no 
                                  ! attempt at convergence)
            iConverge = -1 ! // approximate the answer
            fQ = (Qt - pQatm) / (Qatm - pQatm) !// weighted by fractional
            Z1 = Z2 + (Z1 - Z2) * fQ
            P1 = P2 + (P1 - P2) * fQ
            T1 = T2 + (T1 - T2) * fQ
            dz = Z1 - Z0
         else  ! (ZZ > 0.) bisection (interval halving) method
            iConverge = 1 !  // converge on the answer
            ddZ = ddZ * 0.5
            ddP = ddP * 0.5
            ddT = ddT * 0.5
                    
            if (Qatm > (Qt + Qw)) then
               Z1 = Z1 - ddZ
               P1 = P1 - ddP
               T1 = T1 - ddT
            else
               Z1 = Z1 + ddZ
               P1 = P1 + ddP
               T1 = T1 + ddT
            end if
  
            if (ddZ < ZZ) then
               dz = Z1 - Z0                    
               iConverge = -1 !   // exit
            end if
         end if

         if (iConverge < 0) exit 
      end do

!/* 6. From the previous plume energy Qo, calculate the black body radiation
!      loss Qbb. */
      if (cffeps_sinks) then
         Tt = T0 + Le * dz        !// plume top temperature
         Tm = T0 + (Le - Ld) * dz ! // modified surface temperature
         Qbb = cffeps_timestep * 3600.0 * dz * perimeter * stefan * &
               (((Tm + Tt) * 0.5)**4 - ((T0 + Tt) * 0.5)**4)
      end if
      
      return
   end subroutine cffeps_plumecalc2
!
!============================================================================
!
    real(kind=4) function e_w_calc(Td)   ![K]
       implicit none
       
       real(kind=4), intent(in) :: Td
!/*** vapour pressure over water [mb] ***********************/
!/* Goff-Gratch formula (List 1951) */
       e_w_calc = 10** (-7.90298 * (373.16 / Td - 1.0) + &
                        5.02808 * log10(373.16 / Td) -   &
                        1.3816E-07 * (10**(11.344 * (1.0 - (Td /  &
                        373.16))) - 1.0) + 8.1328E-03 *           &
                        (10**(-3.49149 * ((373.16 / Td) - 1.0)) - &
                        1.0) + log10(1013.246))
     
       return
    end function e_w_calc
   
!/* The following equations are taken from Irabarne and Godson (1973)
!   References are by chapter-equation */
!
!/*** mixing ratio [] ******
   real(kind=4) function r_calc(p, &   ! [mb]
                                Td)    ! [K]
       implicit none
       
       real(kind=4), intent(in) :: p, Td
       real(kind=4), parameter  :: eps = 0.622
       
       r_calc = eps * e_w_calc(Td) / p    ! [Eqn. IV - 77]
       
       return
    end function r_calc
 
!/*** Relative humidity [] ******
    real(kind=4) function rh_calc(q, &   ! [kg/kg]
                                  p, &   ! [mb]
                                  Td)    ! [K]
       implicit none
       
       real(kind=4), intent(in) :: q, p, Td
       real(kind=4)             :: r, e
       real(kind=4), parameter  :: eps = 0.622
       
       r = q / (1.0 - q)                      ! [Eqn. IV - 74]
       e = p * r / (eps + r)                  ! [Eqn. IV - 76]
       rh_calc = 100.0 * e / e_w_calc(Td)     ! [Eqn. IV - 83]
       
       return
    end function rh_calc
 
end module mach_cffeps_plume_mod
