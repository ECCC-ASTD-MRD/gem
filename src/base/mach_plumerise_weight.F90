!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http: /  / www.ec.gc.ca                             !
!============================================================================!
!
! Projet / Project : GEM-MACH
! Fichier / File   : mach_plumerise_weight.ftn90
! Creation       : Stephane Gaudreault, Alexander Kallaur, mars 2008 adaptation for GEM-MACH
!                  Wanmin Gong, for AURAMS 2004
!                  Original version from Janusz Pudykiewicz for CHRONOS 1995
!
! Description    : Assembly of major-point-source emission data for calculation
!                  with Eulerian transport model simulating atmospheric
!                  oxidants. Plume rise formula based on Briggs (1984).
!
! Modification   : Deji Akingunola, Paul Makar, Fall 2017:
!                   - Add a new (optional) new layered plumerise algorithm based on
!                     Briggs (1984) irregular profile.
!
! Arguments:
!            IN
!                cur_source          -> ID of the current source
!                z_magl              -> height of the model levels (m, agl)
!                temperature         -> Model temperature (K)
!                uu, vv              -> Wind components
!                boundlayer_height   -> Height of the boundary layer (m, agl)
!                friction_velocity   -> Friction Speed U* (in m / s)!
!                index_above_pbl     -> Index above boundary layer height
!                safe_inv_mo_length  -> Inverse of Monin obukhov_length
!                stb_func            -> Atmospheric stability parameter
!
!             IN / OUT
!                weight              -> Weight of each cell in the column 
!                                       over the point source
!
! Extra info     :
!
!============================================================================
!
!!if_on
subroutine mach_plumerise_weight(cur_source, z_magl, index_above_pbl, &
                                 safe_inv_mo_length, weight, pbl_hgt, &
                                 ustar, z_temp, wnd_spd, stb_func)
   use chm_ptopo_grid_mod,      only: chm_nk
!!if_off
   use chm_utils_mod,           only: global_debug, chm_lun_out, chm_error_l
   use chm_nml_mod,             only: chm_mj_treatment_s
   use chm_mjrpts_sortinfo_mod, only: lstack_info, i_gilc, i_gjlc, i_hgt, &
                                      i_tem, i_vel, i_dia
   use chm_consphychm_mod,      only: grav, pi, cpd

   implicit none
!!if_on
   integer(kind=4), intent (in) :: cur_source
   real(kind=4),    intent (in) :: z_magl (chm_nk + 1)
   real(kind=4),    intent(out) :: weight (chm_nk)
   integer(kind=4), intent (in) :: index_above_pbl
   real(kind=4),    intent (in) :: pbl_hgt
   real(kind=4),    intent (in) :: ustar
   real(kind=4),    intent (in) :: safe_inv_mo_length
   real(kind=4),    intent (in) :: z_temp (chm_nk + 1)
   real(kind=4),    intent (in) :: wnd_spd(chm_nk + 1)
   real(kind=4),    intent (in) :: stb_func(chm_nk)
!!if_off
!
!  Local variables
!
!  automatic arrays
!
   real(kind=4)    :: layer_height(chm_nk)
!
!  scalars
!
   integer(kind=4) :: k                           ! loop index
   integer(kind=4) :: stack_k                     ! Level above the stack height
   real(kind=4)    :: stack_diameter
   real(kind=4)    :: buoyancy_flux
   real(kind=4)    :: air_temp                    ! Air temperature at position of the stack
   real(kind=4)    :: distance                    ! distance used in interpolation
   real(kind=4)    :: delh1, delh2, plumerise
   real(kind=4)    :: dtdz
   real(kind=4)    :: distance_ratio              ! distance from the top of the boundary layer_height to the stack divided by the height of the plume
   real(kind=4)    :: convective_velocity         ! The convective velocity scale is used as an indicator to identify the particular stability regime
   real(kind=4)    :: wind_speed                  ! Horizontal wind speed for the layer_height
   real(kind=4)    :: stability_parameter         ! Stability parameter used to distinguish the different stability regime (Briggs, 1984)
   real(kind=4)    :: penetration                 ! Plume penetration to free troposphere
   real(kind=4)    :: plume_top                   ! height of the plume top in meter above ground
   real(kind=4)    :: plume_bottom                ! height of the plume bottom in meter above ground
   integer(kind=4) :: index_top                   ! height of the plume top in grid coordinate
   integer(kind=4) :: index_bottom                ! height of the plume bottom in grid coordinate
   real(kind=4)    :: ratio
   real(kind=4)    :: stdev                       ! standard deviation used in the gaussian case
   real(kind=4)    :: column_sum
   real(kind=4)    :: stack_temperature           ! Stack temperature in Kelvin
   real(kind=4)    :: stack_magl                  ! Stack height in meter above ground
   real(kind=4)    :: stack_flow_rate             ! Stack Volume Flow rates ( in meters^3 / sec)
   real(kind=4)    :: stack_radius
   real(kind=4)    :: work
   logical(kind=4) :: local_dbg
   integer(kind=4) :: isvolcano
!
   real(kind=4), dimension(chm_nk) :: zprime, buoy_fl1, buoy_fl2
   real(kind=4)    :: lapse_rate, zprime0, zprime_max, buoy_inv
   integer(kind=4) :: inv_levl, buoyancy_level
!
! Constants
!
   real(kind=4), parameter :: THRESHOLD = 1.e-02
   real(kind=4), parameter :: third8 = 8.0 / 3.0
   real(kind=4), parameter :: third1 = 1.0 / 3.0

   character(len=8) :: plume_type

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   stack_magl        = lstack_info(i_hgt, cur_source)
   stack_temperature = lstack_info(i_tem, cur_source)
   stack_flow_rate   = lstack_info(i_vel, cur_source)
   stack_diameter    = lstack_info(i_dia, cur_source)

   stack_radius      =  stack_diameter * 0.5
   
   weight = 0.0

   isVolcano = 0

   if (local_dbg) then
       write (chm_lun_out, *) 'stack diameter: ', stack_diameter
   end if

   if (stack_diameter < 0.0) then 
      !   stack is a volcano!

      isVolcano = 1

   end if

   IF_VOLCANO: if (isVolcano == 1) then  
 
      !   stack_magl contains height of volcano above sea level
      !   stack_flow_rate contains height of column cloud above sea level 
 
      plume_top = stack_flow_rate - stack_magl 
      plume_top = max(1.5 * THRESHOLD, plume_top)
!     plume_top = max(pbl_hgt, plume_top)
      plume_bottom = 0.0
      plumerise = max(1.5 * THRESHOLD, plume_top)
      penetration = 0.

   else  ! Not volcano, but industrial stack

      do k = chm_nk, 1, -1
         if (stack_magl <= z_magl(k)) then
            stack_k = k
            exit
         end if
      end do

      distance = (z_magl(stack_k) - stack_magl) / &
                 (z_magl(stack_k) - z_magl(stack_k + 1))
      air_temp = Linear_Interpolation(z_temp(stack_k + 1), z_temp(stack_k), &
                                      distance)
      wind_speed = Linear_Interpolation(wnd_spd(stack_k + 1), &
                                        wnd_spd(stack_k), distance)
      if (wind_speed <= 0.0) wind_speed = 0.01

      if (stack_temperature > air_temp .and. stack_flow_rate > 0.0) then
         buoyancy_flux = (grav / pi) * stack_flow_rate *                  &
                         (stack_temperature - air_temp) / stack_temperature
      else
         buoyancy_flux = 0.0
      end if
!
      if (chm_mj_treatment_s == 'PLUMERISE' .and. buoyancy_flux > 0.0) then
! Evaluate stability parameter s = (g / theta)*dtheta / dz
! Note that s is defined as (g / Ta)*dtheta / dz in several references, but the 
! former is defined in Briggs (1984) and is consistent with the definition of 
! N^2 (N: Brunt-Vaisala frequency).
! Note also dtheta / dz = (theta / T)(dT / dz + 0.0098)
         dtdz = (air_temp - z_temp(stack_k + 1)) / &
                (stack_magl -  z_magl(stack_k + 1))
   ! Minimum potential temp gradient for stable condition
         dtdz = max(dtdz, -0.005)
         stability_parameter = (grav / air_temp) * (dtdz + 0.0098)
!
   ! Plume rise under different stability conditions

         if (stack_magl < pbl_hgt) then
            if ((safe_inv_mo_length < 0.0) .and. &
                (-stack_magl * safe_inv_mo_length) > 4.0) then
         ! Unstable condition

                work = (buoyancy_flux / wind_speed)**(3.0 / 5.0)
                convective_velocity = -2.5 * (ustar**3) * safe_inv_mo_length

         ! The plume rise formula proposed by Briggs (1984) is given by
                delh1 = 3.0 * work / (convective_velocity**(2.0 / 5.0))

         ! However, Briggs (1983) suggested a reasonable approximation for the 
         ! convective velocity scale. The rational for this simplification is due
         ! to the lack of data for evaluation.
                delh2 = 30.0 * work

         ! We take the minimum from these two formulations
                plumerise = min(delh1, delh2)

            else if ((safe_inv_mo_length >= 0.0) .and. &
                     (stack_magl * safe_inv_mo_length) > 0.5) then
         ! Stable condition

                plumerise = 2.6 * (buoyancy_flux /                        &
                            (wind_speed * stability_parameter))**(third1)
            else
         ! Neutral condition

                work = buoyancy_flux / (ustar**2 * wind_speed)
                delh1 = 39.0 * buoyancy_flux**(3.0 / 5.0) / wind_speed
                delh2 = 1.2 * work**(3.0 / 5.0) *                         &
                        (stack_magl + 1.3 * work)**(2.0 / 5.0)
                plumerise = min(delh1, delh2)
            end if
         else
      ! Outside the boundary layer_height
            plumerise = 2.6 * (buoyancy_flux /                            &
                       (wind_speed * stability_parameter))**(1.0 / 3.0)
         end if
!         
! Calculate plume penetration
         if ((stack_magl < pbl_hgt) .and. (plumerise > 1.0e-02)) then
            distance_ratio = (pbl_hgt - stack_magl) / plumerise
            if (distance_ratio >= 1.5) then
               penetration = 0.0
            else if (distance_ratio <= 0.5) then
               penetration = 1.0
            else
               penetration = 1.5 - distance_ratio
            end if
            plumerise = min((0.62 + 0.38 * penetration) *               &
                           (pbl_hgt - stack_magl), plumerise)
         else
            penetration = 0.0
         end if
!      
      else if (chm_mj_treatment_s == 'PLUMERISE2' .and. &
               buoyancy_flux > 0.0) then
   ! Fix Briggs' definition of buoyancy flux and multiply by pi for initial flux  
         buoyancy_flux = grav * stack_flow_rate *       &
                         (stack_temperature - air_temp) / air_temp
! NEW PLUMERISE
         dtdz = (z_temp(stack_k) - z_temp(stack_k + 1)) / &
                 (z_magl(stack_k) -  z_magl(stack_k + 1))
   ! Minimum potential temp gradient for stable condition
         dtdz = max(dtdz, -0.005)  !TODO:  Not needed for this method
         lapse_rate = - grav / cpd
   ! Determine air parcel stability based on temperature gradient across &
   ! the layer enclosing the stack height
         k = stack_k
!
   ! Set a limit on the plumerise for the initial buoyant rise
         zprime_max = stack_radius * (sqrt(99.0) - 1.0) / 0.125
         
         zprime0 = z_magl(k) - stack_magl
         do while ((lapse_rate > dtdz) .and. (zprime0 <= zprime_max))
         ! Assume air is unstable, therefore determine the inversion layer,
         ! and work out the buoyancy at the inversion level
         !
         ! Find the inversion layer
            dtdz = (z_temp(k - 1) - z_temp(k)) / (z_magl(k - 1) -  z_magl(k))
            zprime0 = z_magl(k) - stack_magl
            k = k - 1
         end do
!         write(*,*) 'plumerise2 ', cur_source, k, stack_k
!
         inv_levl = k
         zprime(k)   = (z_magl(k) - stack_magl) + 0.001
         if (k < stack_k) then
            ! Evaluate flux reduction (buoyant flux left at inversion)
            buoy_inv = stack_radius ** 2 / ((zprime0 + stack_radius)**2) * &
                        buoyancy_flux
         else
            ! K'th level is at stack_k
            buoy_inv = buoyancy_flux
         end if
         buoy_fl1(k) = buoy_inv - (0.015 * stb_func(k) * buoy_inv**third1 * &
                       zprime(k)**third8)
         buoy_fl2(k) = buoy_inv - (0.053 * stb_func(k) * wnd_spd(k + 1) *   &
                       zprime(k)**3)
         if (buoy_fl1(k) <= 0.0 .and. buoy_fl2(k) <= 0.0) then
            plumerise = buoy_inv * zprime(k) / (buoy_inv - buoy_fl1(k))
         else if (buoy_fl1(k) <= 0.0) then
            plumerise = buoy_inv * zprime(k) / (buoy_inv - buoy_fl1(k))
         else if (buoy_fl2(k) <= 0.0) then
            plumerise = buoy_inv * zprime(k) / (buoy_inv - buoy_fl2(k))
         else
!
            do k = inv_levl - 1, 1, -1
               zprime(k)   = (z_magl(k) - stack_magl) + 0.001
!               dtdz = (z_temp(k) - z_temp(k + 1)) / zprime(k)
!               stb_func(k) = (grav / z_temp(k)) * (dtdz + 0.0098)
               buoy_fl1(k) = buoy_fl1(k + 1) - (0.015 * stb_func(k) *       &
                             buoy_fl1(k + 1)**third1 *                      &
                             (zprime(k)**third8 - zprime(k + 1)**third8))
               buoy_fl2(k) = buoy_fl2(k + 1) - (0.053 * stb_func(k) *       &
                             wnd_spd(k + 1) * (zprime(k)**3 - zprime(k + 1)**3))
! Check whether the plume buoyancy is now 0 or -ve; and determine if it is
! Vertical plume or bent-over plume?
               buoyancy_level = k
               if (buoy_fl1(k) <= 0.0 .and. buoy_fl2(k) <= 0.0) then
                  if ((buoy_fl1(k + 1) - buoy_fl1(k)) > (buoy_fl2(k + 1) - buoy_fl2(k))) then
                     plume_type = 'vertical'
                  else 
                     plume_type = 'bentover'
                  end if
                  exit
               else if (buoy_fl1(k) <= 0.0) then
                  plume_type = 'vertical'
                  exit
               else if (buoy_fl2(k) <= 0.0) then
                  plume_type = 'bentover'
                  exit
               else
                  cycle
               end if
            end do
            k = buoyancy_level
            select case (plume_type)
               case ('vertical')
!               write(*,*) 'plume_custom ', k, inv_levl, buoy_fl1(k + 1), buoy_fl2(k + 1), stb_func(k), zprime(k + 1)
                  work = buoy_fl1(k + 1)**(2.0/3.0) / (0.015 * stb_func(k)) + &
                         zprime(k + 1)**third8
                  plumerise = work ** (3.0 / 8.0)
               case ('bentover')
                  work = buoy_fl2(k + 1) / (0.053 * stb_func(k) * wnd_spd(k + 1)) + &
                         zprime(k + 1)**3
                  plumerise = work ** third1
               case default
                   write(*,*) 'plume level unrealistically high ', buoyancy_level
                   chm_error_l = .true.
                   return
            end select
!            
         end if
         penetration = 0.0
! End of new plumerise
      else
         buoyancy_flux = 0.0
         plumerise     = 0.0
         penetration   = 0.0
      end if

      if (plumerise > THRESHOLD) then
         plume_top     = stack_magl + 1.5 * plumerise
         plume_bottom  = stack_magl + 0.5 * plumerise

         if (penetration > 0.0) then
            plume_top = min(pbl_hgt, plume_top )
            plumerise = plume_top - plume_bottom
         end if

         if (chm_mj_treatment_s == 'PLUMERISE') then
      ! Assuming uniformly mixed plume within BL under unstable condition
      ! (mixed between surface and plume top level at the moment)
            if (safe_inv_mo_length < 0.0 .and. &
                (-stack_magl * safe_inv_mo_length) > 4.0) then
               plume_bottom = 0.0
               plumerise    = plume_top - plume_bottom
            end if
         end if
      end if
     
   end if IF_VOLCANO

   if (plumerise > THRESHOLD) then

      ! layer_height containing plume top and bottom
      plume_top = max(plume_top, 0.0)
      plume_top = min(plume_top, z_magl(1))

      do k = chm_nk, 1, -1
         index_top = k
         if (z_magl(k) > plume_top) then
            exit
         end if
      end do

      plume_bottom = max(plume_bottom, 0.0 )
      plume_bottom = min(plume_bottom, z_magl(1))
      do k = chm_nk, 1, -1
         index_bottom = k
         if (z_magl(k) > plume_bottom) then
            exit
         end if
      end do

      ! layer_height emission amongst layer_heights

      layer_height = 0

      ! Check if plumetop is below boundary layer height
      if (penetration > 0.0 .and. index_above_pbl <= index_top) then
         layer_height(index_above_pbl) = penetration
         ratio = (1.0 - penetration) / plumerise
      else
         ratio = 1.0 / plumerise
      end if

      if (index_bottom == index_top) then
!        print *, 'INDEX TOP AND BOTTOM , chm_nk' , index_top , chm_nk
         layer_height(index_top) = layer_height(index_top) + 1.0 - penetration

      else
         layer_height(index_bottom) = (z_magl(index_bottom) - plume_bottom) * ratio
         layer_height(index_top)    = layer_height(index_top) + &
                                      (plume_top - z_magl(index_top + 1)) * ratio
         if ((index_bottom - index_top) > 1) then
            do k = index_bottom-1, index_top + 1, -1
               layer_height(k) = (z_magl(k) - z_magl(k + 1)) * ratio
            end do
         end if
      end if

      do k = 1, chm_nk -1
         weight(k) = layer_height(k) / (z_magl(k) - z_magl(k + 1))
      end do
         weight(chm_nk) = layer_height(chm_nk) / z_magl(chm_nk)

      if (local_dbg) then
         write (chm_lun_out, *) 'in mach_plumerise_weight'
         write (chm_lun_out, *) 'Buoyant plume'
         if (isVolcano == 1) then
            write (chm_lun_out, *) ' Volcano with...'
         else
            write (chm_lun_out, *) 'stack height ', stack_magl
            write (chm_lun_out, *) 'in model layer ', stack_k, ' with upper bound ', z_magl(stack_k)
            if (stack_k == chm_nk )  then 
               write (chm_lun_out, *) ' and surface as lower bound'
            else
               write (chm_lun_out, *) ' and lower bound ', z_magl(stack_k+1)
            end if
         end if
         write (chm_lun_out, *) 'plume top and plume bottom are  ', plume_top, plume_bottom
         write (chm_lun_out, *) 'with respective upper bounds at levels ', index_top, index_bottom
         write (chm_lun_out, *) 'weights are ', weight
      end if

   else

      ! Non buoyant plumes are treated as Gaussian plume
      if (stack_magl > pbl_hgt) then
         stdev = 50.0
      else
         stdev = 500.0
      endif
      column_sum = 0.0
      do k = 1, chm_nk-1
         column_sum = column_sum + exp( -(z_magl(k) - stack_magl)**2 / &
                      (2.0 * stdev * stdev)) * (z_magl(k) - z_magl(k+1))
      end do
      column_sum = column_sum + exp( -(z_magl(chm_nk) - stack_magl)**2 / &
                      (2.0 * stdev * stdev)) * z_magl(chm_nk)

      do k = 1, chm_nk
         weight(k) = exp(-(z_magl(k) - stack_magl)**2 / &
                          (2.0 * stdev * stdev)) / column_sum
      end do
      if (local_dbg) then
         write (chm_lun_out, *) 'in mach_plumerise_weight'
         write (chm_lun_out, *) 'Non-buoyant plume'
         write (chm_lun_out, *) 'stack height ', stack_magl
         write (chm_lun_out, *) 'weights are ', weight
      end if
   end if

   contains
      real function Linear_Interpolation(val1, val2, distance)
         implicit none
         real(kind=4), intent(in) :: val1, val2, distance
         Linear_Interpolation = val1 + ((val2 - val1) * distance)
      end function

end subroutine mach_plumerise_weight
