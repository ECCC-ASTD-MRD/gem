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
#include <msg.h>
   public :: lhn2


contains

  !/@*
  subroutine lhn2(d, f, v, dsiz, fsiz, vsiz,&
                  dt, ni, nk, kount)
     use phy_options
     use phybus
     use tendency, only: apply_tendencies

     implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf> 
#include "phymkptr.hf" 

     !@Object Apply latent heat nudging during model integration
     !@Arguments
     !
     !          - Input -
     ! dsiz     dimension of dbus
     ! fsiz     dimension of fbus
     ! vsiz     dimension of vbus
     ! dt       timestep (sec.)
     ! ni       horizontal running length
     ! nk       vertical dimension
     ! kount    timestep number
     integer, intent(in) :: dsiz, fsiz, vsiz, kount, ni, nk
     real, intent(in) :: dt

     !          - Input/Output -
     ! d        dynamics input field
     ! f        historic variables for the physics
     ! v        physics tendencies and other output fields from the physics
     real, intent(inout), target :: d(dsiz), f(fsiz), v(vsiz)

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
     ! Additionnally, weighting factors are added to 
     !   1- Use the quality of observations to modulate the intensity of LHN being applied
     !   2- Give the user the possibility to diminish the influence of LHN
     !
     ! In the code, LHN increments are computed as:
     ! increment        = weights * [pr_ratio  - 1.] * initial_profile
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
     real, parameter :: MIN_PR    = 0.1    !mm/h instantaneous precip rate
     real, parameter :: MAX_RATIO = 1.5            !max scaling factor
     real, parameter :: MIN_RATIO = 1./MAX_RATIO   !min scaling factor
     
     !internal variables
     integer :: i, k
     real    :: master_volume, this_radar_pr, this_radar_qi, this_model_pr
     real    :: pr_ratio, scale_fact
     real, dimension(ni,nk)        :: pres_pa, es, qvs_old, qvs_new, rh, trec
     real, pointer, dimension(:,:) :: tinc, ttend, temp_k, sigma, qv, zlhm, zsm2d
     real, pointer, dimension(:)   :: zlhnr, ztree, psp, zlhs
     real, dimension(ni)           :: radar_pr, radar_qi, model_rt

     !do nothing if LHN not in use
     if (lhn /= 'IRPCP') return
  
     !LHN weight as a function of model timestep during a 6h period
     !lhn_weight must be specified in gem_settings.nml
     !  TODO: this should not be hard coded
     select case (kount)       
        case (360:)                 
           master_volume = 0.      !no lhn after IAU period
        case (20:359)                 
           master_volume = lhn_weight !normal volume when in operation
        case (10:19)                  !ramp up volume linearly for a few timesteps
           master_volume = lhn_weight*(kount - 9.)/10. 
        case (:9)                
           master_volume = 0.      !no lhn during spinup period
     end select

     !Output diagnostic:  
     ! 
     !code for decision tree
     MKPTR1D(ztree, tree, v)
     ztree = 0.

     if (master_volume > 0.) then
  
        !Output diagnostics:  
        ! 
        !scaling factor used in LHN 
        MKPTR1D(zlhnr, lhnr, v)
        zlhnr = 0.
        !moisture tendency due to LHN
        MKPTR2D(zlhm, tlhm, v)
        zlhm  = 0.
        !total moisture added/removed in a column
        MKPTR1D(zlhs, tlhs, v)
        zlhs  = 0.
        !temperature tendencies due to Latent Heat Nudging
        MKPTR2D(tinc, tlhn, v)
        tinc  = 0.

        !Model variables
        !
        !temperature
        MKPTR2D(temp_k, tplus, d)
        !vapor mixing ratio
        MKPTR2D(qv, huplus, d)
        !model surface pressure
        MKPTR1D(psp, pmoins, f)
        !model level
        MKPTR2D(sigma, sigt, d)

        !Horizontally smoothed quantities
        !
        !temperature tendencies due to latent heat release
        MKPTR2D(ttend, smta, f)
        !
        !TODO use 2D smoothing instead of putting them in
        !a 3D array
        !
        !3D container for 2D fields
        MKPTR2D(zsm2d, sm2d, f)
        !model precip rate      
        model_rt = zsm2d(:,1) 
        !radar precip rate
        radar_pr = zsm2d(:,2)
        !radar quality index
        radar_qi = zsm2d(:,3)

  
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

           this_radar_qi = radar_qi(i)

           if (this_radar_qi > 0.) then
              !radar observation is valid
  
              !convert model precip rate in m/s to mm/h
              this_model_pr = model_rt(i)*3.6e6
              this_radar_pr = radar_pr(i) !already in mm/h
  
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
                     scale_fact = this_radar_qi*master_volume*(pr_ratio - 1.)
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
                     call use_avg_profile(this_radar_pr, this_radar_qi, pres_pa, master_volume, i, ni, nk, dt, &
                                          v(tlhn), ztree(i) )
                  endif
              else
                  !radar has NO precip
                  if (this_model_pr > MIN_PR) then
                     !model has precip
                     scale_fact = this_radar_qi*master_volume*(MIN_RATIO - 1.)
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
           else
              !radar quality is poor -> do nothing
              ztree(i) = 5.
           endif 
        enddo
  
        !change temperature increment into a tendency
        tinc = tinc/dt
        call apply_tendencies(d,v,f,tplus,tlhn,ni,nk)
  
        !compute increments to humidity to conserve RH
        do k = 1,nk
           do i = 1,ni
  
              !adjust moisture only where temperature has changed
              if (abs(tinc(i,k)*dt) >= 1e-3) then
                  
                  es(i,k)   = YAUA*exp(-1.*YAUB/temp_k(i,k))     !saturation vapor pressure  [Pa] Eq. 2.12 of R&Y
                  qvs_new(i,k)  = EPSIL*es(i,k)/pres_pa(i,k)     !saturation mixing ratio [kg/kg] Eq. 2.18 of R&Y
  
                  !increment to moisture      
                  !     want = have + increment -> increment = want - have
                  zlhm(i,k) = master_volume*(rh(i,k)*qvs_new(i,k) - qv(i,k))
  
                  !output diagnostic
                  zlhs(i) =  zlhs(i) + zlhm(i,k)
              endif

           enddo
        enddo
  
        !change increment into a tendency
        zlhm = zlhm/dt
        call apply_tendencies(d,v,f,huplus,tlhm,ni,nk)
  
     else
        !master volume = 0 ->  LHN is not applied
        ztree = 6.
     endif
     
  end subroutine lhn2


  !/@*
  subroutine use_avg_profile(radar_pr, radar_qi, model_pres_pa, master_volume, i, ni, nk, dt, &
                             lhn_profile, stat) 
     implicit none
#include <arch_specific.hf>
     !@Object Use typical heating profiles where model has no precip 
     !@Arguments
     !
     !          - Input -
     ! radar_pr         precip rate observed by radar
     ! radar_qi         radar observation quality index
     ! model_pres_pa    atmospheric pressure
     ! i                x index
     ! ni               x dimension
     ! nk               z dimension
     ! dt               model timestep (s) 
     integer, intent(in) :: i,ni,nk            !x index ; x dim ;  number of vertical levels
     real,    intent(in) :: dt                 !model timestep
     real,    intent(in) :: radar_pr           !precip rate observed by radar
     real,    intent(in) :: radar_qi           !radar observation quality index
     real,    intent(in) :: master_volume      !modulation factor for LHN
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

     !transform heating rates (K/s) into increments (K) for consistency with the rest of LHN code
     avg_profile = avg_profile * dt

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
        lhn_profile(i,k) = master_volume * ( (1.-r)*avg_profile(hp)+r*avg_profile(lp) )
     enddo
     
  end subroutine use_avg_profile

end module lhn_mod
