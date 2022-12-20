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

module lhn_mod
  ! apply Latent Heat Nudging at beginning of model integration
   implicit none
   private
#include <rmn/msg.h>
   public :: lhn2


contains

  !/@*
  subroutine lhn2(dbus, fbus, vbus, &
                  dt, ni, nk, kount)
     use phy_options
     use phybus
     use tendency, only: apply_tendencies
     use debug_mod, only: init2nan
     
     implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf> 
#include "phymkptr.hf" 

     !@Object Apply latent heat nudging during model integration
     !@Arguments
     !
     !          - Input -
     ! dt       timestep (sec.)
     ! ni       horizontal running length
     ! nk       vertical dimension
     ! kount    timestep number
     integer, intent(in) :: kount, ni, nk
     real, intent(in) :: dt

     !          - Input/Output -
     ! dbus     dynamics input field
     ! fbus     historic variables for the physics
     ! vbus     physics tendencies and other output fields from the physics
     real, dimension(:), pointer, contiguous :: dbus, fbus, vbus

     !@Author D.Jacques 

     !@description
     ! latent heat nudging from radar-inferred precipitation rates
     ! This code computes increments of temperature and moisture for a better
     ! agreement between the simulated and observed precipitation.
     !
     ! The link between the incremental form and classic LHN equation is as follows:
     !
     ! With classic LHN we have:
     !     final_profile = (radar_precip / model_precip) * initial_profile
     !                   =           pr_ratio            * initial_profile
     ! 
     ! In incremental form we want:
     !     finalProfile = initialProfile + increment
     !  
     ! So classic LHN in incremental form is:
     !     finalProfile = initialProfile + (pr_ratio * initial_profile - initial_profile)
     !                  = initialProfile + [(pr_ratio  - 1.) * initial_profile]
     !                  = initialProfile + increment
     ! 
     ! Additionnally, modulation factors are added to 
     !   1- Give the user the possibility to control intensity at which LHN is applied (through 'lhn_weight' defined in gem_settings)
     !   2- Slowly start applying LHN at the beginning of model integrations       (time_modulation)
     !   3- Deactivate LHN when temperatures at an altitude of 1km are below 0degC (temperature_modulation)
     !   4- Use the quality of observations to modulate the intensity of LHN being applied (radar_qi)
     !
     ! In the code, LHN increments are computed as:
     ! increment        = modulation * [pr_ratio  - 1.] * initial_profile
     !
     ! These, increments are transformed into tendencies before their application
     !*@/
     
     !parameters
     !coefficient A for saturation vapor pressure
     !    [Pa]   from Rogers and Yau eq. 2.12  (the one in the book is for kPa)
     real, parameter :: YAUA = 2.53e11
     !coefficient B for saturation vapor pressure [K]   from Rogers and Yau eq. 2.12
     real, parameter :: YAUB = 5.42e3
     !ratio of dry/moist gas constant: RGAS/Rv [unitless]
     real, parameter :: EPSIL = 0.622
     !lhn control parameters
     real, parameter :: MIN_PR      = 0.1    !mm/h instantaneous precip rate
     real, parameter :: MAX_RATIO   = 1.5            !max scaling factor
     real, parameter :: MIN_RATIO   = 1./MAX_RATIO   !min scaling factor
     real, parameter :: OBSV_HEIGHT = 1000.  !altitude of radar observations [m] AGL
     
     !internal variables
     integer :: i, k
     real    :: time_modulation, temperature_modulation, modulation_factor
     real    :: this_radar_pr, this_radar_qi, this_model_pr
     real    :: pr_ratio, scale_fact, r, temp_at_obs_height
     real, dimension(ni,nk)        :: pres_pa, es, qvs_old, qvs_new, rh
     real, pointer, dimension(:,:), contiguous :: tinc, ttend, temp_k, heights_agl, sigma, qv, zlhm
     real, pointer, dimension(:), contiguous   :: zlhnr, ztree, psp, zlhs
     real, pointer, dimension(:), contiguous   :: radar_pr, radar_qi, model_rt
     
     !do nothing if LHN not in use
     if (lhn /= 'IRPCP') return

     call init2nan(pres_pa, es, qvs_old, qvs_new, rh)

     !Output diagnostic:  
     ! 
     !code for decision tree
     MKPTR1D(ztree, tree, vbus)
     ztree = 0.

     !lhn_weight must be specified in gem_settings.nml
     !no LHN applied if lhn_weight <= 0.
     if (lhn_weight <= 0.) then
        ztree = 6.1
        return
     endif
  
     !LHN time modulation as a function of model timestep  (kount)
     if (kount < lhn_start) then
           !LHN is not applied during first few timesteps
           ztree = 6.2
           return
     else if (kount  <  lhn_start+lhn_ramp) then
           !ramp up LHN modulation from 0 to 1 for a few timesteps
           time_modulation = real(kount - lhn_start + 1.)/lhn_ramp
     else if (kount <= lhn_stop) then
           !no time modulation during 'normal' application of LHN 
           time_modulation = 1.       
     else 
           !timestep is > lhn_stop ; LHN is not applied
           ztree = 6.2
           return
     end if 

     !Output diagnostics:  
     ! 
     !scaling factor used in LHN 
     MKPTR1D(zlhnr, lhnr, vbus)
     zlhnr = 0.
     !moisture tendency due to LHN
     MKPTR2D(zlhm, tlhm, vbus)
     zlhm  = 0.
     !total moisture added/removed in a column
     MKPTR1D(zlhs, tlhs, vbus)
     zlhs  = 0.
     !temperature tendencies due to Latent Heat Nudging
     MKPTR2D(tinc, tlhn, vbus)
     tinc  = 0.

     !Model variables
     !
     !temperature
     MKPTR2D(temp_k, tplus, dbus)
     !AGL heights (m) on thermo levels
     MKPTR2D(heights_agl, gztherm, vbus)
     !vapor mixing ratio
     MKPTR2D(qv, huplus, dbus)
     !model surface pressure
     MKPTR1D(psp, pmoins, fbus)
     !model level
     MKPTR2D(sigma, sigt, dbus)

     !Horizontally smoothed quantities
     !
     !temperature tendencies due to latent heat release
     MKPTR2D(ttend, tcond_smt, fbus)
     !modeled precip rate
     MKPTR1D(model_rt, rt_smt, fbus)
     !observed precip rate
     MKPTR1D(radar_pr, rdpr_smt, fbus)
     !observation quality index
     MKPTR1D(radar_qi, rdqi_smt, fbus)

  
     !record relative humidity and temperature before LHN
     do k = 1,nk
        do i = 1,ni
           pres_pa(i,k) = sigma(i,k)*psp(i)             !air pressure [Pa]     formulation from cnv_main.F90  
           es(i,k)   = YAUA*exp(-1.*YAUB/temp_k(i,k))   !saturation vapor pressure [Pa]   Eq. 2.12 of R&Y
           qvs_old(i,k)  = EPSIL*es(i,k)/pres_pa(i,k)   !satturation mixing ratio [kg/kg] Eq. 2.18 of R&Y
           rh(i,k)   = qv(i,k)/qvs_old(i,k)
        enddo
     enddo


     !apply LHN
     do i=1,ni

        !
        !check if observation quality is good enough to apply LHN
        this_radar_qi = radar_qi(i)
        if (this_radar_qi <= 0.) then
           !radar quality is poor -> do nothing and move on to next point
           ztree(i) = 5.
           cycle
        endif

        !
        !modulation of LHN as a function of temperature
        !
        !first interpolate temperature at an altitude of 1km AGL
        do k=nk,1,-1
            if (heights_agl(i,k) > OBSV_HEIGHT) exit
        enddo
        !   at this point,
        !   k is index of first level (starting from the ground) above obsv_height 
        !   k+1 is index of level just below obsv_height
        !   (   height(k) - obsv_height    ) / (               delta h               )
        r = (heights_agl(i,k) - OBSV_HEIGHT) / (heights_agl(i,k) - heights_agl(i,k+1))
        temp_at_obs_height =  (1.-r)*temp_k(i,k) + r*temp_k(i,k+1)
        !
        !second, set modulation from 0. to 1 in the interval between 0 and 5 degC (273-278 degK)
        if (temp_at_obs_height <= 273.) then
           !no LHN below freezing 
           ! temperature_modulation = 0. -> do nothing and move on to next point
           ztree(i) = 6.3
           cycle
        else if (temp_at_obs_height >= 278.) then
           !no influence of temperature above 5 degC
           temperature_modulation = 1.
        else
           ! modulate LHN linearly with temperature between 0 and 5 degC
           temperature_modulation = (temp_at_obs_height - 273.) / 5.
        endif

        !if code gets here, we knon that 
        !     radar observation is valid and temperature at 1km AGL > 0 degC

        !group all modulation factors into one 
        modulation_factor = lhn_weight*temperature_modulation*time_modulation*this_radar_qi
  
        !convert model precip rate in m/s to mm/h
        this_model_pr = model_rt(i)*3.6e6
        this_radar_pr = radar_pr(i) !already in mm/h
  
        !
        !walk down the LHN tree
        if (this_radar_pr > MIN_PR) then
            !radar has precip
  
            if (this_model_pr > MIN_PR) then
               !model has precip
  
               pr_ratio = this_radar_pr / this_model_pr
               ztree(i) = 1.
               !insure ratio is within bounds, adjust status accordingly
               if (pr_ratio < MIN_RATIO) then
                   pr_ratio = MIN_RATIO
                   ztree(i) = .9
               else
                   if (pr_ratio > MAX_RATIO) then
                       pr_ratio = MAX_RATIO
                       ztree(i) = 1.1
                   endif
               endif
               !temperature increment due to LHN
               scale_fact = modulation_factor*(pr_ratio - 1.)
               do k=1,nk
                  !only modify LHR profiles where +ve
                  if (ttend(i,k) > 0.) then
                      !temperature increment due to LHN
                      tinc(i,k) = scale_fact*ttend(i,k)
                  endif
               enddo
               zlhnr(i) = scale_fact
               
            else
               !model has NO precip
               !in this case, use predefined typical profile
               call use_avg_profile(this_radar_pr, pres_pa, modulation_factor, i, ni, nk, &
                                    tinc, ztree(i) )
            endif
        else
            !radar has NO precip
            if (this_model_pr > MIN_PR) then
               !model has precip
               scale_fact = modulation_factor*(MIN_RATIO - 1.)
               do k=1,nk
                  if (ttend(i,k) > 0.) then
                      !cool existing profile for less precip
                      tinc(i,k) = scale_fact*ttend(i,k)  
                  endif
               enddo
               zlhnr(i) = scale_fact
               ztree(i) = 3.
            else
              !model has NO precip -> do nothing
               ztree(i) = 4.
            endif
        endif
     enddo
  
     call apply_tendencies(dbus, vbus, fbus, tplus, tlhn, ni, nk)
  
     !compute increments to humidity to conserve RH
     do k = 1,nk
        do i = 1,ni
  
           !adjust moisture only where temperature has changed
           if (abs(tinc(i,k)*dt) >= 1e-3) then
               
               es(i,k)   = YAUA*exp(-1.*YAUB/temp_k(i,k))     !saturation vapor pressure  [Pa] Eq. 2.12 of R&Y
               qvs_new(i,k)  = EPSIL*es(i,k)/pres_pa(i,k)     !saturation mixing ratio [kg/kg] Eq. 2.18 of R&Y
  
               !increment to moisture      
               !     want = have + increment -> increment = want - have
               zlhm(i,k) = rh(i,k)*qvs_new(i,k) - qv(i,k)
  
               !output diagnostic
               zlhs(i) =  zlhs(i) + zlhm(i,k)
           endif

        enddo
     enddo
  
     !change increment into a tendency
     zlhm = zlhm/dt
     call apply_tendencies(dbus, vbus, fbus, huplus, tlhm, ni, nk)
     
  end subroutine lhn2


  !/@*
  subroutine use_avg_profile(radar_pr, model_pres_pa, modulation_factor, i, ni, nk, &
                             lhn_profile, stat) 
     use phy_options
     implicit none
#include <arch_specific.hf>
     !@Object Use typical heating profiles where model has no precip 
     !@Arguments
     !
     !          - Input -
     ! radar_pr         precip rate observed by radar
     ! model_pres_pa    atmospheric pressure
     ! i                x index
     ! ni               x dimension
     ! nk               z dimension
     integer, intent(in) :: i,ni,nk            !x index ; x dim ;  number of vertical levels
     real,    intent(in) :: radar_pr           !precip rate observed by radar
     real,    intent(in) :: modulation_factor  !modulation factor for LHN
     real,    intent(in), dimension(ni,nk) :: model_pres_pa !model pressure
     
     !          - Output -
     ! lhn_profile      Latent Heating profile interpolated on model pressure
     ! stat             status for decision tree statistics
     real, intent(out), dimension(ni,nk) :: lhn_profile  !Temperature increment due to LHN 
     real, intent(out) :: stat                !status for decision tree statistics

     !*@/

     integer, parameter :: NAVG = 18
     real,    parameter :: delta_p = 5.000e+03  !bin size of average profile [Pa]

     !          - Internal -
     integer               :: k, hp,lp
     real                  :: r 
     real, dimension(NAVG) :: avg_profile, avg_pres_pa !average profile to be interpolated

     !pressure levels for average profiles [Pa]
     avg_pres_pa=(/ 2.000e+04, 2.500e+04, 3.000e+04, 3.500e+04, 4.000e+04, &
                    4.500e+04, 5.000e+04, 5.550e+04, 6.000e+04, 6.500e+04, &
                    7.000e+04, 7.500e+04, 8.000e+04, 8.500e+04, 9.000e+04, &
                    9.500e+04, 1.000e+05, 1.050e+05 /)

     !determine which profile to use based on intensity of radar precipitation
     !
     ! avg_profiles are heating profiles due to latent heat release in K/s
     !
     if (radar_pr < 1.5) then
        avg_profile = 0. !no heating for radar pr < 1.5 mm/h
        stat = 2.0
     else if (radar_pr < 3.0) then
        !< 3.0 millimeters of rain per second
        avg_profile = (/ 0.000e+00, 0.000e+00, 0.000e+00, 5.000e-04, 9.000e-04, &
                         1.200e-03, 1.400e-03, 1.500e-03, 1.600e-03, 1.600e-03, &
                         1.700e-03, 1.800e-03, 1.700e-03, 1.300e-03, 5.000e-04, &
                         0.000e+00, 0.000e+00, 0.000e+00 /)
        stat = 2.1
     else if (radar_pr < 6.12) then
        !between 3. and 6.12 mm/h
        avg_profile = (/ 0.000e+00, 3.000e-04, 7.000e-04, 1.000e-03, 1.300e-03, &
                         1.700e-03, 2.000e-03, 2.100e-03, 2.200e-03, 2.300e-03, &
                         2.500e-03, 2.400e-03, 2.300e-03, 2.100e-03, 1.000e-03, &
                         0.000e+00, 0.000e+00, 0.000e+00 /)
 
        stat = 2.2
     else if (radar_pr < 12.5) then
        !between 6.12 and 12.5 mm/h
        avg_profile = (/ 0.000e+00, 3.000e-04, 7.000e-04, 1.000e-03, 1.300e-03, &
                         1.700e-03, 2.000e-03, 2.100e-03, 2.200e-03, 2.600e-03, &
                         2.800e-03, 3.000e-03, 3.000e-03, 2.500e-03, 1.500e-03, &
                         0.000e+00, 0.000e+00, 0.000e+00 /)
        stat = 2.3
     else if (radar_pr < 25.) then
        !between 12. and 25. mm/h
        avg_profile = (/ 0.000e+00, 8.000e-04, 1.600e-03, 2.300e-03, 3.000e-03, &
                         3.800e-03, 4.300e-03, 5.100e-03, 5.900e-03, 6.200e-03, &
                         5.600e-03, 5.000e-03, 4.100e-03, 2.900e-03, 1.500e-03, &
                         0.000e+00, 0.000e+00, 0.000e+00 /)
        stat = 2.4
     else if (radar_pr < 50.) then
        !between 25. and 50. mm/h
        avg_profile = (/ 0.000e+00, 1.400e-03, 2.600e-03, 4.000e-03, 5.900e-03, &
                         7.300e-03, 8.700e-03, 9.400e-03, 9.700e-03, 9.600e-03, &
                         9.100e-03, 7.400e-03, 5.900e-03, 3.900e-03, 2.000e-03, &
                         0.000e+00, 0.000e+00, 0.000e+00 /)
        stat = 2.5
     else
        !above 50 mm/h
        avg_profile = (/ 0.000e+00, 1.900e-03, 3.800e-03, 7.300e-03, 1.140e-02, &
                         1.390e-02, 1.480e-02, 1.490e-02, 1.300e-02, 1.150e-02, &
                         9.700e-03, 8.000e-03, 5.800e-03, 3.000e-03, 1.500e-03, &
                         0.000e+00, 0.000e+00, 0.000e+00 /)
        stat = 2.6
     endif

     !initialize profile to zero
     lhn_profile(i,:) = 0.

     !skip high altitudes where LHR does not occur
     k=1
     do while (model_pres_pa(i,k) < avg_pres_pa(1)) 
        k = k + 1
        lhn_profile(i,k) = 0.
     enddo
     !linear interp of average profile
     hp = 1
     lp = 2
     do k=k,nk
        if (model_pres_pa(i,k) >= avg_pres_pa(lp)) then
           hp = hp+1
           lp = lp+1
           if (lp > NAVG) exit
        endif
        r = ( model_pres_pa(i,k) - avg_pres_pa(hp) )/delta_p
        lhn_profile(i,k) = modulation_factor * ( (1.-r)*avg_profile(hp)+r*avg_profile(lp) )
     enddo
     
  end subroutine use_avg_profile

end module lhn_mod
