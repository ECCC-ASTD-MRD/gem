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

module sfclayer
  use, intrinsic :: iso_fortran_env, only: REAL64
  use tdpack
  use sfclayer_funcs
  implicit none
  private

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

  ! Public parameters
  integer, parameter, public :: SL_OK    = RMN_OK           !Return code for success
  integer, parameter, public :: SL_ERROR = RMN_ERR          !Return code for error

  ! Private parameters
  integer, parameter :: SL_LONG = 1024                      !Maximum string length
  real,    parameter :: FACTN = 1.2                         !PBL height scaling factor
  real,    parameter :: HMIN = 30.                          !Minimum PBL height
  real,    parameter :: DEFAULT_HGHTT_DIAG = 1.5            !Default screen level
  real,    parameter :: DEFAULT_HGHTM_DIAG = 10.            !Default anemometer level
  real,    parameter :: DEFAULT_L_MIN = -1.                 !Default minimum M-O length

  ! Private variables
  real, save    :: beta = 1.                                !Prandtl number for a neutral profile
  real, save    :: rineutral = 0.                           !Width of neutral Ri regime
  logical, save :: tdiaglim_default = .false.               !Default value for inversion limiter
  logical, save :: z0ref = .true.                           !Use a reference roughness (max of z0m and z0t)

  ! Private procedures
  procedure(stability_function), pointer, save :: &
       sf_stable => sf_stable_delage97, &                   !Stability functions (stable surface layer)
       sf_unstable => sf_unstable_delage92                  !Stability functions (unstable surface layer)

  ! Public subprograms
  public :: sl_put          !Set module variables
  public :: sl_get          !Retrieve module variables/parameters
  public :: sl_prelim       !Preliminary surface layer diagnostics
  public :: sl_sfclayer     !Surface layer parameterization
  public :: sl_adjust       !Adjust diagnostic values

  ! Generic procedures
  interface sl_put
     module procedure sl_put_i4
     module procedure sl_put_r4
     module procedure sl_put_l
     module procedure sl_put_s
  end interface sl_put

  interface sl_get
     module procedure sl_get_r4
     module procedure sl_get_l
  end interface sl_get

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_put_i4(key,val) result(status)
     ! Set a variable in the surface layer module

     ! Input arguments
     character(len=*), intent(in) :: key                  !Name of key to set
     integer, intent(in) :: val                           !Value of the key

     ! Output arguments
     integer :: status                                    !Return status of function
     character(len=256) :: tmp_S

     ! Initialize return value
     status = SL_ERROR

     ! Attempt to set value of requested key
     select case (key)
     case ('stderr','STDERR')
        call msg_toall(MSG_WARNING,'(sl_put_i4) STDERR is ignored by the MSG messaging system')
     case DEFAULT
        write(tmp_S, '(i0)') val
        call msg_toall(MSG_WARNING,'(sl_put_i4) cannot set '//trim(key)//'='//trim(tmp_S))
        return
     end select

     ! Successful completion of subprogram
     status = SL_OK

  end function sl_put_i4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_put_r4(key,val) result(status)
     ! Set a variable in the surface layer module

     ! Input arguments
     character(len=*), intent(in) :: key                  !Name of key to set
     real, intent(in) :: val                              !Value of the key

     ! Output arguments
     integer :: status                                    !Return status of function

     ! Initialize return value
     status = SL_ERROR

     ! Attempt to set value of requested key
     select case (key)
     case ('beta','BETA')
        beta = val
     case ('rineutral','RINEUTRAL')
        rineutral = val
     case DEFAULT
        call msg_toall(MSG_WARNING,'(sl_put_r4) cannot set '//trim(key))
        return
     end select

     ! Successful completion of subprogram
     status = SL_OK

  end function sl_put_r4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_put_l(key,val) result(status)
     ! Set a variable in the surface layer module

     ! Input arguments
     character(len=*), intent(in) :: key                  !Name of key to set
     logical, intent(in) :: val                           !Value of the key

     ! Output arguments
     integer :: status                                    !Return status of function

     ! Initialize return value
     status = SL_ERROR

     ! Attempt to set value of requested key
     select case (key)
     case ('tdiaglim','TDIAGLIM')
        tdiaglim_default = val
     case ('z0ref','Z0REF')
        z0ref = val
     case DEFAULT
        call msg_toall(MSG_WARNING,'(sl_put_l) cannot set '//trim(key))
        return
     end select

     ! Successful completion of subprogram
     status = SL_OK

  end function sl_put_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_put_s(key,val) result(status)
     use clib_itf_mod, only: clib_toupper, CLIB_OK
     ! Set a variable in the surface layer module

     ! Input arguments
     character(len=*), intent(in) :: key                  !Name of key to set
     character(len=*), intent(in) :: val                  !Value of the key

     ! Output arguments
     integer :: status                                    !Return status of function

     ! Local variables
     character(len=len_trim(val)) :: ucval

     ! Initialize return value
     status = SL_ERROR

     ! Convert value to upper case
     ucval = val
     if (clib_toupper(ucval) /= CLIB_OK) then
        call msg_toall(MSG_WARNING,'(sl_put_s) cannot convert input to upper case: '//trim(val))
        return
     endif

     ! Attempt to set value of requested key
     select case (key)
     case ('sl_stabfunc_stab','SL_STABFUNC_STAB')
        if (sf_get(sf_stable, 'stable', ucval) /= SF_OK) then
           call msg_toall(MSG_WARNING,'(sl_put_s) cannot acquire stable SF for '//trim(val))
           return
        endif
     case ('sl_stabfunc_unstab','SL_STABFUNC_UNSTAB')
        if (sf_get(sf_unstable, 'unstable', ucval) /= SF_OK) then
           call msg_toall(MSG_WARNING,'(sl_put_s) cannot acquire unstable SF for '//trim(val))
           return
        endif
     case DEFAULT
        call msg_toall(MSG_WARNING,'(sl_put_s) cannot set '//trim(key))
       return
     end select

     ! Successful completion of subprogram
     status = SL_OK

  end function sl_put_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_get_r4(key,val) result(status)
     use sfclayer_funcs
     ! Retrieve a variable or parameter from the module

     ! Input arguments
     character(len=*), intent(in) :: key                  !Name of key to retrieve

     ! Output arguments
     integer :: status                                    !Return status of function
     real, intent(out) :: val                             !Value of the key

     ! Initialize return value
     status = SL_ERROR
     val = -1.

     ! Attempt to return value of requested key
     select case (key)
     case ('beta','BETA')
        val = beta
     case ('bh91_a','BH91_A')
        val = BH91_A
     case ('bh91_b','BH91_B')
        val = BH91_B
     case ('bh91_c','BH91_C')
        val = BH91_C
     case ('bh91_d','BH91_D')
        val = BH91_D
     case ('d97_as','D97_AS')
        val = D97_AS
     case ('dg92_ci','DG92_CI')
        val = DG92_CI
     case ('l07_ah','L07_AH')
        val = L07_AH
     case ('l07_am','L07_AM')
        val = L07_AM
     case ('rineutral','RINEUTRAL')
        val = rineutral
     case DEFAULT
        call msg_toall(MSG_WARNING,'(sl_get_r4) cannot retrieve '//trim(key))
        return
     end select

     ! Successful completion of subprogram
     status = SL_OK

  end function sl_get_r4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_get_l(key,val) result(status)
     ! Retrieve a variable or parameter from the module

     ! Input arguments
     character(len=*), intent(in) :: key                  !Name of key to retrieve

     ! Output arguments
     integer :: status                                    !Return status of function
     logical, intent(out) :: val                          !Value of the key

     ! Initialize return value
     status = SL_ERROR
     val = .false.

     ! Attempt to return value of requested key
     select case (key)
     case ('tdiaglim','TDIAGLIM')
        val = tdiaglim_default
     case DEFAULT
        call msg_toall(MSG_WARNING,'(sl_get_l) cannot retrieve '//trim(key))
        return
     end select

     ! Successful completion of subprogram
     status = SL_OK

  end function sl_get_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_sfclayer(t_air,q_air,spd_air,dir_air,hghtm_air,hghtt_air,t_sfc,q_sfc,z0m,z0t,lat,fcor, &
       coefm,coeft,rib,flux_t,flux_q,ilmo,ue,h,lzz0m,lzz0t,stabm,stabt,spdlim, &
       t_diag,q_diag,u_diag,v_diag,hghtm_diag,hghtm_diag_row,hghtt_diag, &
       hghtt_diag_row,tdiaglim,L_min,optz0,z0mloc,z0t_optz0) result(status)
    ! Surface layer parameterization

    ! Input arguments
    real, dimension(:), intent(in) :: t_air                     !Lowest level potential temperature (K)
    real, dimension(:), intent(in) :: q_air                     !Lowest level specific humidity (kg/kg)
    real, dimension(:), intent(in) :: spd_air                   !Lowest level wind speed (m/s)
    real, dimension(:), intent(in) :: dir_air                   !Lowest level wind direction (rad)
    real, dimension(:), intent(in) :: hghtm_air                 !Height of the lowest momentum level (m)
    real, dimension(:), intent(in) :: hghtt_air                 !Height of the lowest thermodynamic level (m)
    real, dimension(:), intent(in) :: t_sfc                     !Surface air temperature (K)
    real, dimension(:), intent(in) :: q_sfc                     !Surface specific humidity (kg/kg)
    real, dimension(:), intent(in) :: z0m                       !Momentum roughness length (m)
    real, dimension(:), intent(in) :: z0t                       !Thermal roughness length (m)
    real, dimension(:), intent(in) :: lat                       !Latitude (rad)
    real, dimension(:), intent(in) :: fcor                      !Coriolis factor (/s)
    logical, intent(in), optional :: tdiaglim                   !Limit temperature inversion in lowest layer [get 'tdiaglim']
    real, intent(in), optional :: L_min                         !Minimum stable Obukhov length [none]
    integer, intent(in), optional :: optz0                      !Alternative roughness length adjustment [0]
    real, intent(in), optional :: hghtm_diag                    !Diagnostic level height for momentum (m) [10.]
    real, dimension(:), intent(in), optional :: hghtm_diag_row  !Diagnostic heights for momentum (m; supercedes hghtm_diag)
    real, intent(in), optional :: hghtt_diag                    !Diagnostic level height for thermdynamics (m) [1.5]
    real, dimension(:), intent(in), optional :: hghtt_diag_row  !Diagnostic heights for thermdynamics (m; supercedes hghtt_diag)
    real, dimension(:), intent(in), optional :: z0mloc          !Local Momentum roughness length (no orography) (m)


    ! Output arguments
    integer :: status                                           !Return status of function
    real, dimension(:), intent(out), optional :: ilmo           !Inverse of the Obukov length (/m)
    real, dimension(:), intent(out), optional :: h              !Boundary layer height (m)
    real, dimension(:), intent(out), optional :: ue             !Friction velocity (m/s)
    real, dimension(:), intent(out), optional :: flux_t         !Temperature flux (Km/s)
    real, dimension(:), intent(out), optional :: flux_q         !Moisture flux (kgm/kgs)
    real, dimension(:), intent(out), optional :: coefm          !Momentum exchange coefficient (m/s)
    real, dimension(:), intent(out), optional :: coeft          !Thermal exchange coefficient (m/s)
    real, dimension(:), intent(out), optional :: rib            !Bulk Richardson number
    real, dimension(:), intent(out), optional :: lzz0m          !Log of adjusted momentum roughness
    real, dimension(:), intent(out), optional :: lzz0t          !Log of adjusted thermodynamic roughness
    real, dimension(:), intent(out), optional :: stabm          !Integrated momentum stability function
    real, dimension(:), intent(out), optional :: stabt          !Integrated thermodynamic stability function
    real, dimension(:), intent(out), optional :: spdlim         !Adjusted wind speed at the lowest momentum level (m/s)
    real, dimension(:), intent(out), optional :: t_diag         !Diagnostic level temperature (K)
    real, dimension(:), intent(out), optional :: q_diag         !Diagnostic level moisture (kg/kg)
    real, dimension(:), intent(out), optional :: u_diag         !Diagnostic level u-wind (m/s)
    real, dimension(:), intent(out), optional :: v_diag         !Diagnostic level v-wind (m/s)
    real, dimension(:), intent(out), optional :: z0t_optz0      !Thermodynamic roughness calculated using optz0


    ! Internal variables
    integer :: my_optz0
    real :: my_L_min
    real, dimension(size(t_air)) :: my_coefm,my_coeft,my_rib,my_flux_t,my_flux_q, &
         my_ilmo,my_ue,my_h,my_lzz0m,my_lzz0t,my_stabm,my_stabt,my_t_diag,my_q_diag, &
         my_u_diag,my_v_diag,my_hghtm_diag,my_hghtt_diag,my_spdlim, my_z0mloc, my_z0t_optz0
    logical :: my_tdiaglim

    ! Initialize return value
    status = SL_ERROR

    ! Set default values
    my_optz0 = 0
    if (present(optz0)) my_optz0 = optz0
    my_hghtm_diag(:) = DEFAULT_HGHTM_DIAG
    if (present(hghtm_diag)) my_hghtm_diag(:) = hghtm_diag
    if (present(hghtm_diag_row)) my_hghtm_diag(:) = hghtm_diag_row(:)
    my_hghtt_diag(:) = DEFAULT_HGHTT_DIAG
    if (present(hghtt_diag)) my_hghtt_diag(:) = hghtt_diag
    if (present(hghtt_diag_row)) my_hghtt_diag(:) = hghtt_diag_row(:)
    my_tdiaglim = tdiaglim_default
    if (present(tdiaglim)) my_tdiaglim = tdiaglim
    my_L_min = DEFAULT_L_MIN
    if (present(L_min)) my_L_min = L_min
    my_z0mloc(:) = z0m(:)
    if (present(z0mloc)) my_z0mloc(:) = z0mloc(:)


    ! Compute surface fluxes on request
    call flxsurf(my_coefm,my_coeft,my_rib,my_flux_t,my_flux_q,my_ilmo,my_ue, &
         fcor,t_air,q_air,hghtm_air,hghtt_air,spd_air,t_sfc,q_sfc,my_h,z0m,z0t,my_L_min, &
         my_lzz0m,my_lzz0t,my_stabm,my_stabt,my_spdlim,size(t_air),my_optz0,my_z0mloc,my_z0t_optz0)

    ! Diagnostic level calculations
    if (any((/present(hghtm_diag),present(hghtt_diag),present(hghtm_diag_row),present(hghtt_diag_row)/))) then
       ! Compute diagnostic level quantities
       call diasurf(my_u_diag,my_v_diag,my_t_diag,my_q_diag,size(t_air), &
            dir_air,t_sfc,q_sfc,z0m,z0t,my_ilmo,hghtm_air,my_h,my_ue,my_flux_t, &
            my_flux_q,my_hghtm_diag,my_hghtt_diag,lat)
       ! Apply diagnostic adjustments if requested
       if (sl_adjust(t_air,hghtt_air,my_t_diag,hghtt_diag_row=my_hghtt_diag,tdiaglim=my_tdiaglim, &
            adj_t_diag=my_t_diag) /= SL_OK) then
          call msg_toall(MSG_ERROR,'(sl_sfclayer) unable to adjust diagnostic values')
          return
       endif
    else
       if (any((/present(t_diag),present(q_diag),present(u_diag),present(v_diag)/))) then
          call msg_toall(MSG_ERROR,'(sl_sfclayer) diagnostic heights must be provided for calculation')
          return
       endif
    endif

    ! Fill requested return variables
    if (present(coefm)) coefm = my_coefm
    if (present(coeft)) coeft = my_coeft
    if (present(rib)) rib = my_rib
    if (present(flux_t)) flux_t = my_flux_t
    if (present(flux_q)) flux_q = my_flux_q
    if (present(ilmo)) ilmo = my_ilmo
    if (present(ue)) ue = my_ue
    if (present(h)) h = my_h
    if (present(lzz0m)) lzz0m = my_lzz0m
    if (present(lzz0t)) lzz0t = my_lzz0t
    if (present(stabm)) stabm = my_stabm
    if (present(stabt)) stabt = my_stabt
    if (present(spdlim)) spdlim = my_spdlim
    if (present(t_diag)) t_diag = my_t_diag
    if (present(q_diag)) q_diag = my_q_diag
    if (present(u_diag)) u_diag = my_u_diag
    if (present(v_diag)) v_diag = my_v_diag
    if (present(z0t_optz0)) then
       if (my_optz0 > 0) then
          z0t_optz0 = my_z0t_optz0
       else
          call msg_toall(MSG_ERROR,'(sl_sfclayer) optz0 > 0 must be provided for z0t_optz0 calculation')
          return
       endif
    endif
    ! Successful completion of subprogram
    status = SL_OK

  end function sl_sfclayer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_prelim(t_air,q_air,u_air,v_air,p_sfc,hghtm_air,spd_air,dir_air,tv_air,rho_air, &
                     min_wind_speed,min_wind_reduc) result(status)
    ! Preliminary calculations for derived surface layer fields

    ! Input arguments
    real, dimension(:), intent(in) :: t_air               !Lowest level temperature (K)
    real, dimension(:), intent(in) :: q_air               !Lowest level specific humidity (kg/kg)
    real, dimension(:), intent(in) :: u_air               !Lowest level u-component wind speed (m/s)
    real, dimension(:), intent(in) :: v_air               !Lowest level v-component wind speed (m/s)
    real, dimension(:), intent(in) :: p_sfc               !Surface pressure (Pa)
    real, dimension(:), intent(in) :: hghtm_air           !Height of the lowest momentum level (m)
    real, intent(in), optional :: min_wind_speed          !Minimum wind speed to apply (m/s) [0.]
    character(len=*), intent(in), optional :: min_wind_reduc !Minimum wind speed reduction type (['none'],'linear')

    ! Output arguments
    integer :: status                                     !Return status of function
    real, dimension(:), intent(out), optional :: spd_air  !Lowest level wind speed (m/s)
    real, dimension(:), intent(out), optional :: dir_air  !Lowest level wind direction (rad)
    real, dimension(:), intent(out), optional :: tv_air   !Lowest level virtual temperature (K)
    real, dimension(:), intent(out), optional :: rho_air  !Lowest level air density (kg/m3)

    ! Internal variables
    integer :: i
    real, dimension(size(u_air)) :: my_min_wind_speed
    real, dimension(size(t_air)) :: my_tv_air
    character(len=SL_LONG) :: my_min_wind_reduc

    ! Initialize return value
    status = SL_ERROR

    ! Set default values
    my_min_wind_speed = 0.
    if (present(min_wind_speed)) my_min_wind_speed = min_wind_speed
    my_min_wind_reduc = 'none'
    if (present(min_wind_reduc)) my_min_wind_reduc = min_wind_reduc

    ! Adjust minimum wind speed as function of momentum level height if requested
    if (my_min_wind_reduc == 'linear') my_min_wind_speed = min(my_min_wind_speed,max(my_min_wind_speed*hghtm_air/40.,1.e-6))

    ! Compute requested wind fields
    if (present(spd_air)) spd_air = max(sqrt(u_air*u_air + v_air*v_air),my_min_wind_speed)
    if (present(dir_air)) dir_air = atan2(v_air,sign(max(abs(u_air),epsilon(u_air)),u_air))

    ! Compute requested mass fields
    if (any((/present(tv_air),present(rho_air)/))) then
       do i=1,size(my_tv_air)
          my_tv_air(i) = fotvt(t_air(i),q_air(i))
       enddo
    endif
    if (present(tv_air)) tv_air = my_tv_air
    if (present(rho_air)) rho_air = p_sfc/(RGASD*my_tv_air)

    ! Succesful completion of subprogram
    status = SL_OK

  end function sl_prelim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_adjust(t_air,hghtt_air,t_diag,hghtt_diag,hghtt_diag_row,tdiaglim,adj_t_diag) result(status)
    ! Adjust diagnostic level values according to prescribed operations

    ! Input arguments
    real, dimension(:), intent(in) :: t_air                     !Lowest level temperature (K)
    real, dimension(:), intent(in) :: hghtt_air                 !Height of the lowest thermodynamic level (m)
    real, dimension(:), intent(in) :: t_diag                    !Diagnostic level temperature (K)
    real, intent(in), optional :: hghtt_diag                    !Diagnostic level height for thermdynamics (m) [1.5]
    real, dimension(:), intent(in), optional :: hghtt_diag_row  !Diagnostic height for thermdynamics (m; supercedes hghtt_diag)
    logical, intent(in), optional :: tdiaglim                   !Limit temperature inversion in lowest layer [get 'tdiaglim']

    ! Output arguments
    integer :: status                                           !Return status of function
    real, dimension(:), intent(out), optional :: adj_t_diag     !Adjusted diagnostic level temperature (K)

    ! Internal variables
    real, dimension(size(t_diag)) :: my_adj_t_diag,my_hghtt_diag
    logical :: my_tdiaglim

    ! Initialize return values
    status = SL_ERROR

    ! Set default values
    my_tdiaglim = tdiaglim_default
    if (present(tdiaglim)) my_tdiaglim = tdiaglim
    my_hghtt_diag(:) = DEFAULT_HGHTT_DIAG
    if (present(hghtt_diag)) my_hghtt_diag(:) = hghtt_diag
    if (present(hghtt_diag_row)) my_hghtt_diag(:) = hghtt_diag_row(:)
    my_adj_t_diag(:) = t_diag(:)

    ! Modify diagnostic level temperature on request (~ 8K/40m maximum)
    if (my_tdiaglim) my_adj_t_diag(:) = max(t_diag(:),t_air(:)-0.2*(hghtt_air(:)-my_hghtt_diag(:)))
    if (present(adj_t_diag)) adj_t_diag(:) = my_adj_t_diag(:)

    ! Succesful completion of subprogram
    status = SL_OK

 end function sl_adjust

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine flxsurf(cmu, ctu, rib, ftemp, fvap, ilmo, &
       ue, fcor, ta , qa , zu, zt, vmod, tg , qg , h , z0 , z0t,  &
       lmin, lzz0, lzz0t, fm, fh, va, n, optz0, z0loc, optz0_z0t)
    use sfclayer_funcs, only: stability_function, D97_AS, SF_MOMENTUM, SF_HEAT
!!!#include <arch_specific.hf>
    ! Represent surface layer similarity state, turbulent tranfer coefficients and fluxes

    ! Input arguments
    integer, intent(in) :: n                                    !Grid dimension
    real, dimension(:), intent(in) :: ta                        !Lowest level potential temperature (K)
    real, dimension(:), intent(in) :: qa                        !Lowest level specific humidity (kg/kg)
    real, dimension(:), intent(in) :: vmod                      !Lowest-level wind speed (m/s)
    real, dimension(:), intent(in) :: zu                        !Height of the lowest momentum level (m)
    real, dimension(:), intent(in) :: zt                        !Height of the lowest thermodynamic level (m)
    real, dimension(:), intent(in) :: tg                        !Surface air temperature (K)
    real, dimension(:), intent(in) :: qg                        !Surface specific humidity (kg/kg)
    real, dimension(:), intent(in) :: z0                        !Momentum roughness length (m)
    real, dimension(:), intent(in) :: z0t                       !Thermal roughness length (m)
    real, dimension(:), intent(in) :: fcor                      !Coriolis factor (/s)
    real, intent(in) :: lmin                                    !Minimum stable Obukhov length [none]
    integer, intent(in) :: optz0                                !Alternative roughness length adjustment [0]
    real, dimension(:), intent(in) :: z0loc                     !Local momentum roughness length (m) for compz0() calc.

    ! Output arguments
    real, dimension(:), intent(out) :: ilmo                     !Inverse of the Obukov length (/m)
    real, dimension(:), intent(out) :: h                        !Boundary layer height (m)
    real, dimension(:), intent(out) :: ue                       !Friction velocity (m/s)
    real, dimension(:), intent(out) :: ftemp                    !Temperature flux (Km/s)
    real, dimension(:), intent(out) :: fvap                     !Moisture flux (kgm/kgs)
    real, dimension(:), intent(out) :: cmu                      !Momentum exchange coefficient (m/s)
    real, dimension(:), intent(out) :: ctu                      !Thermal exchange coefficient (m/s)
    real, dimension(:), intent(out) :: rib                      !Bulk Richardson number
    real, dimension(:), intent(out) :: lzz0                     !Log of adjusted momentum roughness
    real, dimension(:), intent(out) :: lzz0t                    !Log of adjusted thermodynamic roughness
    real, dimension(:), intent(out) :: fm                       !Integrated momentum stability function
    real, dimension(:), intent(out) :: fh                       !Integrated thermodynamic stability function
    real, dimension(:), intent(out) :: va                       !Adjusted lowest-level wind speed (m/s)
    real, dimension(:), intent(out) :: optz0_z0t                !Thermal roughness length (m) calculated with compz0()

    ! Internal parameters
    real, parameter :: VMIN=1.e-6                               !Minimum wind speed
    real, parameter :: RIBMIN=1.e-5                             !Minimum absolute bulk Richardson number
    real, parameter :: RIBMAX=1.e5                              !Maximum absolute bulk Richardson number
    real, parameter :: HMAX=1500.                               !Maximum PBL height estimate
    real, parameter :: EPSLN=1e-5                               !Small value

    ! Internal variables
    integer :: j, it, itmax
    real :: cm, ct, zp
    real :: ilmoj, fmj, fhj, hj, z0j, z0tj, zuj, ztj, vaj, fcorj
    real :: z0rmj, z0rtj, lzz0j, lzz0tj, ribj, z0locj
    real(REAL64) :: dthv, tva, tvs
    real :: ilmax, vlmin, hc, ribc, ilmm
    real :: am, ah, dfm, dfh, g, dg
    real, dimension(n) :: z0rt, z0rm  !#, lz0,lz0t
    real, dimension(2) :: ff_zm, ff_zh
    procedure(stability_function), pointer :: sf

    ! Setup for numerical solution
    itmax  = 3
    ilmax = -1.
    if (lmin > 0.) then
       ilmax = 1./lmin
       itmax = 5
    endif

    ! Establish reference roughness length estimates
    if (z0ref) then
       do j=1,n
          z0rm(j) = max(z0(j),z0t(j))
          z0rt(j) = z0rm(j)
       enddo
    else
       z0rm(:) = z0(:)
       z0rt(:) = z0t(:)
    endif

    ! Compute neutral stability functions
    do j=1,n
       lzz0(j) = log((z0rm(j) + zu(j)) / z0(j))
       lzz0t(j) = log((z0rt(j) + zt(j)) / z0t(j))
    enddo

    ! Update estimate of the Obukhov length and stability functions
    ROW: do j=1,n

       ! Copy inputs to scalars as an optimization
       z0j = z0(j)
       z0locj = z0loc(j)
       z0tj = z0t(j)
       z0rmj = z0rm(j)
       z0rtj = z0rt(j)
       fcorj = fcor(j)
       lzz0j = lzz0(j)
       lzz0tj = lzz0t(j)
       zuj = zu(j)
       ztj = zt(j)

       ! Compute bulk Ri and consistent lowest-level wind
       vaj = max(vmod(j), VMIN)
       zp  = zuj**2 / ztj
       tva = (1.d0 + delta * qa(j)) * ta(j)
       tvs = (1.d0 + delta * qg(j)) * tg(j)
       dthv   = tva - tvs
       vlmin  = 0.
       ! Apply Obukhov length limiter under stable conditions
       if (lmin > 0. .and. dthv > 0.) then
          ilmm = sign(ilmax, real(dthv))
          hc = pblheight(zuj, z0j, vaj, ilmm, fcorj, lzz0j)
          sf => sf_stable
          ff_zm = (/zuj+z0rmj, z0j/)
          ff_zh = (/ztj+z0rtj, z0tj/)
          call sf(fmj, fhj, lzz0j, lzz0tj, ilmm, hc, beta, &
               F_zm=ff_zm, F_zh=ff_zh)
          ribc = zp * ilmm * fhj / fmj**2
          ! Wind speed is adjusted to keep Obukhov length above the minimum
          vlmin = sqrt( max(0.d0, 1./ribc * GRAV * zp * dthv/(tvs + 0.5*dthv)) )
       endif
       vaj = max(vaj, vlmin)
       ribj = GRAV / (tvs + 0.5*dthv) * zp * dthv/vaj**2
       ribj = sign(max(abs(ribj), RIBMIN), ribj)
       ribj = sign(min(abs(ribj), RIBMAX), ribj)

       ! Prescribe a smooth-surface neutral regime below a threshold bulk Ri
       if (rineutral > 0. .and. z0j < 0.01 .and. abs(ribj) < rineutral) ribj = 0.

       ! First guess of neutral stability function adjustments and Obukhov length
       ! uses the Delage (1997) stability functions
       if (ribj >= 0.)  then
          am = 2.5*D97_AS * ribj / max(2*z0j, 1.0)
          ah = 2.5*D97_AS * ribj / max(sqrt(z0j*z0tj), 1.0)
       else
          am = -min(0.7 + log(1.-ribj), lzz0j-1.)
          ah = -min(0.7 + log(1.-ribj), lzz0tj-1.)
       endif
       fmj = lzz0j + am
       fhj = beta * (lzz0tj + ah)
       ilmoj = ribj * fmj * fmj / (zp*fhj)

       ! Find the (inverse) Obukhov length (byproducts: PBL height, integrated stability functions)
       ITERATIONS: do it=1,itmax

          ! Update roughness and neutral stability functions on request
          if(optz0 > 0) then
             call compz0_a(optz0, z0j, z0locj, z0tj, fmj, vaj, fcorj, 1)
             lzz0j = log((z0rmj + zuj) / z0j)
             lzz0tj = log((z0rtj + ztj) / z0tj)
          endif

          ! Ensure that the Obukhov length is finite
          ilmoj = sign(max(abs(ilmoj), EPSLN), ilmoj)

          ! Estimate PBL height based on surface layer properties
          hj  = pblheight(zuj, z0j, vaj, ilmoj, fcorj, fmj)

          ! Update stability functions using appropriate stability branch
          if (ilmoj > 0.) then
             sf => sf_stable
          else
             sf => sf_unstable
          endif
          ff_zm = (/zuj+z0rmj, z0j/)
          ff_zh = (/ztj+z0rtj, z0tj/)
          call sf(fmj, fhj, lzz0j, lzz0tj, ilmoj, hj, beta, &
               F_zm=ff_zm, F_zh=ff_zh, &
               F_dfm=dfm, F_dfh=dfh)
          ! Use Newton-Raphson iterations to solve the equation
          !   rib - (fh/fm^2)*zp*ilmo = 0
          ! This is done by defining an auxiliary function
          !   g(ilmo) = rib - (fh/fm^2)*zp*ilmo
          ! and finding its zero using a Newton-Raphson method
          !   ilmo_(n+1) = ilmo_n - g/dg
          ! where
          !   dg  = derivative of g w.r.t. ilmo
          !       = -fh/(fm^2)*zp*(1+dfh/fh-2*dfm/fm)
          if (it < itmax) then
             g = ribj - fhj / (fmj * fmj) * zp * ilmoj
             dg = -fhj / (fmj*fmj) * zp * (1. + dfh/fhj - 2.*dfm/fmj)
             dg = sign(max(abs(dg),EPSLN),dg)
             ilmoj = ilmoj - g/dg
          endif

       enddo ITERATIONS

       ! Diagnose surface layer properties
       fm(j) = fmj                              !Integrated momentum stability function
       fh(j) = fhj                              !Integrated heat stability function
       rib(j) = ribj                            !Bulk Richardson number
       ilmo(j) = ilmoj                          !Inverse Obukhov length
       h(j) = hj                                !Boundary layer height estimate
       va(j) = vaj                              !Adjusted first-level wind speed

       ! Diagnose turbulent surface exhange properties
       cm = karman / fmj                        !Momentum exchange coefficient
       ct = karman / fhj                        !Heat exchange coefficient
       ue(j) = cm * vaj                         !Friction velocity
       cmu(j) = cm * ue(j)                      !Momentum exchange coefficient * ustar
       ctu(j) = ct * ue(j)                      !Heat exchange coefficient * ustar
       ftemp(j) = -ctu(j) * (ta(j) - tg(j))     !Turbulent temperature flux estimate
       fvap(j) = -ctu(j) * (qa(j) - qg(j))      !Turbulent moisture flux estimate
       h(j)     = min(hj, HMAX)                 !Boundary layer height estimate

    enddo ROW

    ! Recompute consistent roughness length estimates on request
    if(optz0.gt.0) then
       call compz0_a(optz0, z0, z0loc, z0t, fm, va, fcor, n)
       optz0_z0t = z0t
!!$       lz0(1:n) = z0(1:n)
!!$       lz0t(1:n) = z0t(1:n)
!!$       call compz0(optz0, lz0, z0loc, lz0t, fm, va, fcor, n)
!!$       optz0_z0t = lz0t
    endif

    ! End of subprogram
    return
  end subroutine flxsurf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine diasurf(uz, vz, tz, qz, ni, angi, tg, qg, z0, z0t, ilmo, za, &
       h, ue, ftemp, fvap, zu, zt, lat)
    use sfclayer_funcs, only: stability_function, D97_AS, SF_MOMENTUM, SF_HEAT
!!!#include <arch_specific.hf>
    ! Compute state variables at the diagnostic level

    ! Input arguments
    integer, intent(in) :: ni                                   !Grid dimension
    real, dimension(:), intent(in) :: angi                      !Lowest level wind direction (rad)
    real, dimension(:), intent(in) :: tg                        !Surface air temperature (K)
    real, dimension(:), intent(in) :: qg                        !Surface specific humidity (kg/kg)
    real, dimension(:), intent(in) :: z0                        !Momentum roughness length (m)
    real, dimension(:), intent(in) :: z0t                       !Thermal roughness length (m)
    real, dimension(:), intent(in) :: ilmo                      !Inverse of the Obukov length (/m)
    real, dimension(:), intent(in) :: za                        !Height of the lowest momentum level (m)
    real, dimension(:), intent(in) :: h                         !Boundary layer height (m)
    real, dimension(:), intent(in) :: ue                        !Friction velocity (m/s)
    real, dimension(:), intent(in) :: ftemp                     !Temperature flux (Km/s)
    real, dimension(:), intent(in) :: fvap                      !Moisture flux (kgm/kgs)
    real, dimension(:), intent(in) :: zu                        !Diagnostic heights for momentum (m)
    real, dimension(:), intent(in) :: zt                        !Diagnostic heights for thermdynamics (m)
    real, dimension(:), intent(in) :: lat                       !Latitude (rad)

    ! Output arguments
    real, dimension(:), intent(out) :: tz                       !Diagnostic level temperature (K)
    real, dimension(:), intent(out) :: qz                       !Diagnostic level moisture (kg/kg)
    real, dimension(:), intent(out) :: uz                       !Diagnostic level u-wind (m/s)
    real, dimension(:), intent(out) :: vz                       !Diagnostic level v-wind (m/s)

    ! Internal parameters
    real, parameter :: ANGMAX = 0.85                            !Maximum frictional deflection (rad)

    ! Internal variables
    integer :: j
    real :: fh, fm, h1, h2, h3, hh, ct, ctu, cm, vits, dang, ang, hi
    real, dimension(ni) :: lzz0t, lzz0
    real, dimension(2) :: ff_zm, ff_zh
    procedure(stability_function), pointer :: sf

    ! Initialize neutral stability functions
    do j=1,ni
       lzz0(j) = log(1 + zu(j) / z0(j))
       lzz0t(j) = log(1 + zt(j) / z0t(j))
    enddo

    ! Compute state variables at the diagnostic levels
    ROW: do j=1,ni

       ! Combine boundary layer height estimates under stable conditions
       if (ilmo(j) > 0.) then
          h1 = (za(j) + 10.*z0(j)) * factn
          h2 = h(j)
          h3 = factn / (4.*D97_AS*beta * ilmo(j))
          hh = max(HMIN, h1, h2, h3)
       else
          hh = h(j)
       endif

       ! Compute integrated stability functions
       if (ilmo(j) > 0.) then
          sf => sf_stable
       else
          sf => sf_unstable
       endif
       ff_zm = (/zu(j)+z0(j), z0(j)/)
       ff_zh = (/zt(j)+z0t(j), z0t(j)/)
       call sf(fm, fh, lzz0(j), lzz0t(j), ilmo(j), hh, beta, &
            F_zm=ff_zm, F_zh=ff_zh)

       ! Compute exchange coefficients
       ct = karman / fh
       cm = karman / fm

       ! Diagnose temperature and moisture by reversing the bluk flux equations
       !   ftemp = -ctu*((tz + g*zt/cp) - tg) ==> tz = tg - ftemp/ctu - g*zt/cp
       !   fvap  = -ctu*( qz            - qg) ==> qz = qg - fvap /ctu
       ctu = ct * ue(j)
       tz(j) = tg(j) - ftemp(j)/ctu - GRAV/CPD * zt(j)
       qz(j) = qg(j) - fvap (j)/ctu

       ! Diagnose wind speed by reversing the friction velocity equation
       !   ue = cm*vits ==> vits = ue/cm
       vits  = ue(j)/cm

       ! Frictional turning only under stable conditions
       if (ilmo(j).gt.0.) then
          hi = 1./MAX(HMIN,hh)
       else
          hi = 0.
       endif

       ! Diagnose wind components by computing frictional turning angle
       dang = (za(j)-zu(j)) * hi * ANGMAX * sin(lat(j))
       ang = angi(j) + dang
       uz(j) = vits*COS(ang)
       vz(j) = vits*SIN(ang)

    enddo ROW

    ! End of subprogram
    return
  end subroutine diasurf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pblheight(F_zu, F_z0, F_u, F_ilmo, F_fcor, F_fm) result(F_nhb)
    use sfclayer_funcs, only: D97_AS
!!!#include <arch_specific.hf>
    ! Compute the planetary boundary layer height.

    ! Input arguments
    real, intent(in) :: F_zu                !height of wind input (m)
    real, intent(in) :: F_z0                !roughness length for momentum (m)
    real, intent(in) :: F_u                 !wind speed at F_zu (m/s)
    real, intent(in) :: F_ilmo              !inverse of Obukhov length (1/m)
    real, intent(in) :: F_fcor              !Coriolis factor (1/s)
    real, intent(in) :: F_fm                !integrated stability function for momemtum

    ! Output arguments
    real :: F_nhb                           !height of the PBL (m)

    ! Local parameters
    real, parameter :: BS = 1.              !stability function coefficient (stable)
    real, parameter :: CORMIN = .7e-4       !minimum Coriolis parameter value (/s)

    ! Local variables
    real :: h1, h2, h3, f

    ! Boundary layer height estimate based on surface layer quantities
    f = max(abs(F_fcor), cormin)
    if (F_ilmo >= 0.) then
       ! Stable branch
       h1 = (F_zu + 10.*F_z0) * factn
       h2 = BS * SQRT(karman * F_u / (F_ilmo*f*F_fm))
       h3 = factn / (4*D97_AS * beta * F_ilmo)
       F_nhb = MAX(HMIN, h1, h2, h3)
    else
       ! Unstable branch
       h1 = 0.3 * (F_u * karman / F_fm) / f
       F_nhb = max(HMIN, h1)
    endif

    ! End of subprogram
    return
  end function pblheight

end module sfclayer
