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

module sfclayer_mod
  use, intrinsic :: iso_fortran_env, only: REAL64
  use tdpack
  implicit none
  private

#include <rmnlib_basics.hf>
#include <msg.h>

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

  ! Private parameters for the Beljaars and Holtslag (1991) stability functions (stable)
  real, parameter :: BH91_A = 1.                            !A coefficient
  real, parameter :: BH91_B = 2./3.                         !B coefficient
  real, parameter :: BH91_C = 5.                            !C coefficient
  real, parameter :: BH91_D = 0.35                          !D coefficient

  ! Private parameters for the Lock (2007) stability functions (stable)
  real, parameter :: L07_AH = 1.                            !AH coefficient [was 4 in Lock (2007)]
  real, parameter :: L07_AM = 1.                            !AM coefficient [was 4 in Lock (2007)]

  ! Private parameters for the Delage (1997) stability functions (stable)
  real, parameter :: D97_AS = 12.                           !AS coefficient

  ! Private parameters for the Delage and Girard (1991) stability functions (unstable)
  real, parameter :: DG92_CI = 40.                          !CI coefficient
  real, parameter :: DG92_RAC3 = 1.732050776                !RAC3 coefficient

  ! Private variables
  real    :: beta = 1.                                      !Prandtl number for a neutral profile
  real    :: rineutral = 0.                                 !Width of neutral Ri regime
  logical :: tdiaglim_default = .false.                     !Default value for inversion limiter
  logical :: z0ref = .true.                                 !Use a reference roughness (max of z0m and z0t)
  character(len=SL_LONG) :: sl_stabfunc_stab = 'DELAGE97'   !Class of stability functions for the stable case
  character(len=SL_LONG) :: sl_stabfunc_unstab = 'DELAGE92' !Class of stability functions for the unstable case

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
        write(tmp_S, *) val
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
        sl_stabfunc_stab = ucval
     case ('sl_stabfunc_unstab','SL_STABFUNC_UNSTAB')
        sl_stabfunc_unstab = ucval
     case DEFAULT
        call msg_toall(MSG_WARNING,'(sl_put_s) cannot set '//trim(key))
       return
     end select

     ! Successful completion of subprogram
     status = SL_OK

  end function sl_put_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sl_get_r4(key,val) result(status)
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
       hghtt_diag_row,tdiaglim,L_min,optz0) result(status)
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
    real, intent(in), optional :: L_min                         !Minimum stable Monin-Obukhov length [none]
    integer, intent(in), optional :: optz0                      !Alternative roughness length adjustment [0]
    real, intent(in), optional :: hghtm_diag                    !Diagnostic level height for momentum (m) [10.]
    real, dimension(:), intent(in), optional :: hghtm_diag_row  !Diagnostic heights for momentum (m; supercedes hghtm_diag)
    real, intent(in), optional :: hghtt_diag                    !Diagnostic level height for thermdynamics (m) [1.5]
    real, dimension(:), intent(in), optional :: hghtt_diag_row  !Diagnostic heights for thermdynamics (m; supercedes hghtt_diag)

    ! Output arguments
    integer :: status                                           !Return status of function
    real, dimension(:), intent(out), optional :: ilmo           !Inverse of the Monin-Obukov length (/m)
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

    ! Internal variables
    integer :: my_optz0
    real :: my_L_min
    real, dimension(size(t_air)) :: my_coefm,my_coeft,my_rib,my_flux_t,my_flux_q, &
         my_ilmo,my_ue,my_h,my_lzz0m,my_lzz0t,my_stabm,my_stabt,my_t_diag,my_q_diag, &
         my_u_diag,my_v_diag,my_hghtm_diag,my_hghtt_diag,my_spdlim
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

    ! Compute surface fluxes on request
    call flxsurf5(my_coefm,my_coeft,my_rib,my_flux_t,my_flux_q,my_ilmo,my_ue, &
         fcor,t_air,q_air,hghtm_air,hghtt_air,spd_air,t_sfc,q_sfc,my_h,z0m,z0t,my_L_min, &
         my_lzz0m,my_lzz0t,my_stabm,my_stabt,my_spdlim,size(t_air),my_optz0)

    ! Diagnostic level calculations
    if (any((/present(hghtm_diag),present(hghtt_diag),present(hghtm_diag_row),present(hghtt_diag_row)/))) then
       ! Compute diagnostic level quantities
       call diasurf5(my_u_diag,my_v_diag,my_t_diag,my_q_diag,size(t_air), &
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
  subroutine flxsurf5(cmu, ctu, rib, ftemp, fvap, ilmo, &
       ue, fcor, ta , qa , zu, zt, vmod, &
       tg , qg , h , z0 , z0t, lmin, &
       lzz0, lzz0t, fm, fh, va, n, optz0 )

    implicit none
!!!#include <arch_specific.hf>
    integer n
    real cmu(n),ctu(n),rib(n),fcor(n),ilmo(n)
    real ftemp(n),fvap(n),ta(n),qa(n),zu(n),zt(n),va(n)
    real tg(n),qg(n),h(n),z0(n),ue(n),vmod(n)
    real z0t(n),lzz0(n),lzz0t(n)
    real fm(n),fh(n)
    real lmin
    integer optz0
    !
    !Author
    !          Y.Delage (Jul 1990)
    !Revision
    ! 001      G. Pellerin (Jun 94) New function for unstable case
    ! 002      G. Pellerin (Jui 94) New formulation for stable case
    ! 003      B. Bilodeau (Nov 95) Replace VK by KARMAN
    ! 004      M. Desgagne (Dec 95) Add safety code in function ff
    !                               and ensures that RIB is non zero
    ! 005      R. Sarrazin (Jan 96) Correction for H
    ! 006      C. Girard (Nov 95) - Diffuse T instead of Tv
    ! 007      G. Pellerin (Feb 96) Revised calculation for H (stable)
    ! 008      G. Pellerin (Feb 96) Remove corrective terms to CTU
    ! 009      Y. Delage and B. Bilodeau (Jul 97) - Cleanup
    ! 010      Y. Delage (Feb 98) - Addition of HMIN
    ! 011      D. Talbot and Y. Delage (Jan 02) -
    !             Correct bug of zero divide by dg in loop 35
    ! 012      Y. Delage (Oct 03) - Set top of surface layer at ZU +Z0
    !                   - Output UE instead of UE**2 and rename subroutine
    !                   - Change iteration scheme for stable case
    !                   - Introduce log-linear profile for near-neutral stable cases
    !                   - set VAMIN inside flxsurf and initialise ILMO and H
    !                   - Put stability functions into local functions via stabfunc.cdk
    ! 013      Y. Delage (Sep 04) - Input of wind and temperature/humidity
    !                                at different levels
    ! 014      R. McTaggart-Cowan and B. Bilodeau (May 2006) -
    !             Clean up stabfunc.cdk
    ! 015      L. Spacek (Dec 07) - Correction of the log-linear profile
    !                               Double precision for rib calculations
    ! 016      A. Zadra (Feb 11) - clean-up and generalization
    ! 017      A. Zadra (Sep 11) - convert to fortran 90
    ! 018      S. Leroyer (Apr 11) - Argument optz0  for options for roughness lenghts computation
    ! 019      M. Abrahmowicz (May 13) - optz0=0 now means no extra roughness calculation (no call to compz0)
    !
    !Object
    !          to calculate surface layer transfer coefficients and fluxes
    !
    !Arguments
    !
    !          - Output -
    ! cmu      transfer coefficient of momentum times UE
    ! ctu      transfer coefficient of temperature times UE
    ! rib      bulk Richardson number
    ! ftemp    temperature flux
    ! fvap     vapor flux
    ! ilmo     (1/length of Monin-Obukov)
    ! ue       friction velocity
    ! h        height of the boundary layer
    ! fm       momentum stability function
    ! fh       heat stability function
    ! lzz0     log ((zu+z0)/z0)
    ! lzz0t    log ((zt+z0t)/z0t)
    !
    !          - Input -
    ! fcor     Coriolis factor
    ! zu       height of wind input (measured from model base at topo height + Z0)
    ! zt       height of temperature and humidity input
    ! ta       potential temperature at ZT
    ! qa       specific humidity at ZT
    ! va       wind speed at ZU
    ! tg       surface temperature
    ! qg       specific humidity at the surface
    ! z0       roughness length for momentum      flux calculations
    ! z0t      roughness length for heat/moisture flux calculations
    ! n        horizontal dimension
    ! optz0    option for z0, z0t formulations -- optz0=0 : NO EXTRA CALCULATIONS.
    !
    !
    EXTERNAL COMPZ0
    !
    integer j,it,itmax
    real cm,ct
    real x1,x0,y1,y0
    real(REAL64) :: dthv,tva,tvs
    real vmin,ribmin,ribmax,hmax,epsln,ilmax,vlmin,hc,ribc,ilmm
    real am,ah,dfm,dfh,g,dg
    real, dimension(n) :: zp,z0rt,z0rm

    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     Initialize internal parameters:
    !     - maximum number of iterations (itmax)
    !     - minimum value of wind speed (vamin)
    !     - minimum absolute value of bulk Richardson number (ribmin)
    !     - maximum absolute value of bulk Richardson number (ribmax)
    !     - maximum value of boundary layer depth (hmax)
    !     - minimum value of inverse Monin-Obukhov length (epsln)
    !
    itmax  = 3
    vmin   = 1.0e-6
    ribmin = 1.0e-05
    ribmax = 1.0e5
    hmax   = 1500.0
    epsln  = 1.0e-05
    ilmax = -1.
    if (lmin > 0.) then
       ilmax = 1./lmin
       itmax = 5
    endif
    !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     STEP 1: Impose minimum value on wind speed (va) and initialize neutral
    !             stability functions (lzz0 and lzz0t)
    !                 lzz0  = log( (zu+z0 )/z0  ) = log( 1+zu/z0  )
    !                 lzz0t = log( (zt+z0t)/z0t ) = log( 1+zt/z0t )
    !             z0, z0t are from the previous time step
    !
    if (z0ref) then
       z0rm(:) = max(z0(:),z0t(:))
       z0rt(:) = z0rm(:)
    else
       z0rm(:) = z0(:)
       z0rt(:) = z0t(:)
    endif
    do j=1,n
       va(j)    = MAX(vmod(j),vmin)
       lzz0 (j) = (z0rm(j)+zu(j))/z0(j)
       lzz0t(j) = (z0rt(j)+zt(j))/z0t(j)
!!$       lzz0 (j) = 1+zu(j)/z0(j)
!!$       lzz0t(j) = 1+zt(j)/z0t(j)
    enddo
    call vslog(lzz0t,lzz0t,n)
    call vslog(lzz0 ,lzz0 ,n)
    !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     STEP 2: Calculate bulk Richardson number (rib)
    !              rib = A/(B^2)
    !              A = (g/tvm)*(tva-tvs)/zt = buoyancy term
    !                 tva = atmospheric virtual temperature
    !                 tvs = surface virtual temperature
    !                 tvm = 0.5*(tva+tvs)
    !              B = va/zu = shear
    !
    do j=1,n
       zp(j)  = zu(j)**2/zt(j)
       tva    = (1.d0+delta*qa(j))*ta(j)
       tvs    = (1.d0+delta*qg(j))*tg(j)
       dthv   = tva-tvs
       vlmin  = 0.
       if (ilmax > 0. .and. dthv > 0.) then
          ilmm   = sign(ilmax,real(dthv))
          hc    = sf_pblheight(zu(j),z0(j),va(j),ilmm,fcor(j),lzz0(j))
          fm(j) =        lzz0 (j) &
               + sf_momentum(zu(j)+z0rm(j),ilmm,hc,x1) &
               - sf_momentum(      z0  (j),ilmm,hc,x0)
          fh(j) = beta*( lzz0t(j) &
               + sf_heat(zt(j)+z0rt(j),ilmm,hc,y1) &
               - sf_heat(      z0t (j),ilmm,hc,y0))
          ribc  = zp(j)*ilmm * fh(j) / (fm(j)*fm(j))
          vlmin = sqrt( max(0d0, 1d0/ribc * grav*zp(j) * dthv/(tvs+0.5*dthv)) )
       endif
       va(j) = max(va(j), vlmin)
       rib(j) = grav/(tvs+0.5*dthv)*zp(j)*dthv/(va(j)*va(j))
       rib(j) = sign(max(abs(rib(j)),ribmin),rib(j))
       rib(j) = sign(min(abs(rib(j)),ribmax),rib(j))
    enddo

    ! Prescribe a smooth-surface neutral regime below a threshold bulk Ri
    if (rineutral > 0.) then
       where (z0 < 0.01 .and. abs(rib) < rineutral)
          rib = 0.
       endwhere
    endif

    !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     STEP 3: Calculate inverse length of Monin-Obukhov (ilmo), integrated
    !             stability functions (fm, fh) and boundary-layer height (h)
    !             using an iterative method to solve the equation:
    !                  rib - (fh/fm^2)*zp*ilmo = 0
    !             This is done by defining an auxiliary function
    !                  g(ilmo) = rib - (fh/fm^2)*zp*ilmo
    !             and finding its zero using a Newton-Raphson method
    !                  ilmo_(n+1) = ilmo_n - g/dg
    !             where
    !                  dg  = derivative of g w.r.t. ilmo
    !                      = -fh/(fm^2)*zp*(1+dfh/fh-2*dfm/fm)
    !                  dfh = derivative of fh w.r.t. ilmo
    !                  dfm = derivative of fm w.r.t. ilmo
    !
    !------ first guess for fm, fh (chosen as small deviations am,ah from the
    !       neutral case) and the corresponding first guess for ilmo.
    !       TODO: replace this with a function-specific first guess.
    !
    do j=1,n
       if(rib(j).ge.0.)  then
          am = 2.5*D97_AS*rib(j)/MAX(      2*z0(j)     ,1.0)
          ah = 2.5*D97_AS*rib(j)/MAX(SQRT(z0(j)*z0t(j)),1.0)
       else
          am = -MIN(0.7+LOG(1-rib(j)),lzz0 (j)-1)
          ah = -MIN(0.7+LOG(1-rib(j)),lzz0t(j)-1)
       endif
       fm(j)   =       lzz0 (j) + am
       fh(j)   = beta*(lzz0t(j) + ah)
       ilmo(j) = rib(j)*fm(j)*fm(j)/(zp(j)*fh(j))
    enddo
    !
    !------ iterative solution using Newton-Raphson method
    !       (no convergence criterion, but fixed number of
    !        iterations = itmax-1)
    !
    do it=1,itmax

       if(optz0.gt.0) then
          !-- new computation of z0 and z0t from the last computed fm
          ! and repeat step 1 for new lzz0 lzz0t
          call compz0(optz0, z0, z0t, fm, va, fcor,n)
          !
          do j=1,n
             lzz0 (j) = (z0rm(j)+zu(j))/z0(j)
             lzz0t(j) = (z0rt(j)+zt(j))/z0t(j)
!!$             lzz0 (j) = 1+zu(j)/z0(j)
!!$             lzz0t(j) = 1+zt(j)/z0t(j)
          enddo
          call vslog(lzz0t,lzz0t,n)
          call vslog(lzz0 ,lzz0 ,n)
       endif
       !
       do j=1,n
          if (rib(j).ge.0.) ilmo(j) = MAX( epsln,ilmo(j))
          if (rib(j).lt.0.) ilmo(j) = MIN(-epsln,ilmo(j))
          h(j)  = sf_pblheight(zu(j),z0(j),va(j),ilmo(j),fcor(j),fm(j))
          fm(j) =        lzz0 (j) &
               + sf_momentum(zu(j)+z0rm(j),ilmo(j),h(j),x1) &
               - sf_momentum(      z0  (j),ilmo(j),h(j),x0)
          fh(j) = beta*( lzz0t(j) &
               + sf_heat(zt(j)+z0rt(j),ilmo(j),h(j),y1) &
               - sf_heat(      z0t (j),ilmo(j),h(j),y0))
          dfm =       x1 - x0
          dfh = beta*(y1 - y0)
          if (it.lt.itmax) then
             g  = rib(j) - fh(j)/(fm(j)*fm(j))*zp(j)*ilmo(j)
             dg = -fh(j)/(fm(j)*fm(j))*zp(j)*(1+dfh/fh(j)-2*dfm/fm(j))
             dg = sign(max(abs(dg),epsln),dg)
             ilmo(j) = ilmo(j) - g/dg
          endif
       enddo
    enddo
    !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     STEP 4: Calculate friction velocity (ue), transfer coefficients (cmu,ctu),
    !             temperature and vapor fluxes (ftemp, fvap):
    !               ue    = (k/fm)*va      ;  cmu = (k/fm)*ue
    !               ftemp = -ctu*(ta-tg)   ;  ctu = (k/fh)*ue
    !               fvap  = -ctu*(qa-qg)
    !             Note: the momemtum flux is not output separately because it is
    !             simply given by
    !               fmom = cmu*va = (k/fm)*ue*va = ue^2
    !             Finaly, we impose a maximum on the BL height (h).
    !
    do j=1,n
       cm       = karman/fm(j)
       ct       = karman/fh(j)
       ue(j)    = cm*va(j)
       cmu(j)   = cm*ue(j)
       ctu(j)   = ct*ue(j)
       ftemp(j) = -ctu(j)*(ta(j)-tg(j))
       fvap(j)  = -ctu(j)*(qa(j)-qg(j))
       !
       h(j)     = MIN(h(j),hmax)
    enddo
    !
    !
    if(optz0.gt.0) then
       !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !     STEP 5: Calculate aerodynamic (z0) and thermal (z0t) roughness lenghts
       !              z0 and z0t computed with the last fm for use of the next time step
       !
       call compz0(optz0, z0, z0t, fm, va, fcor, n)
       !
       !
    endif
    !
    !
    return
  end subroutine flxsurf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine diasurf5(uz,vz,tz,qz,ni,angi,tg,qg,z0,z0t,ilmo,za, &
       h,ue,ftemp,fvap,zu,zt,lat)
    implicit none
!!!#include <arch_specific.hf>
    integer ni
    real zt(ni),zu(ni)
    real uz(ni),vz(ni),tz(ni),qz(ni),za(ni),angi(ni)
    real tg(ni),qg(ni),ue(ni),ftemp(ni),fvap(ni)
    real lat(ni),ilmo(ni),z0t(ni),z0(ni),h(ni)
    !
    !Author
    !          Yves Delage  (Aug1990)
    !
    !Revision
    ! 001      G. Pellerin(JUN94)
    !          Adaptation to new surface formulation
    ! 002      B. Bilodeau (Nov 95) - Replace VK by KARMAN
    ! 003      R. Sarrazin (Jan 96) - Prevent problems if zu < za
    ! 004      G. Pellerin (Feb 96) - Rewrite stable formulation
    ! 005      Y. Delage and B. Bilodeau (Jul 97) - Cleanup
    ! 006      Y. Delage (Feb 98) - Addition of HMIN
    ! 007      G. Pellerin (Mai 03) - Conversion IBM
    !               - calls to vslog routine (from massvp4 library)
    ! 008      Y. Delage (Oct 03) - Change UE2 by UE and rename subroutine
    !             - Introduce log-linear profile for near-neutral cases
    !             - Put stability functions into local functions via stabfunc.cdk
    ! 009      R. McTaggart-Cowan and B. Bilodeau (May 2006)
    !             - Clean up stabfunc.cdk
    ! 010      A. Zadra (Apr 11) - clean up, simplify and add comments
    ! 011      A. Zadra (Sep 11) - convert to fortran 90
    !
    !Object
    !          to calculate the diagnostic values of U, V, T, Q
    !          near the surface (ZU and ZT)
    !
    !Arguments
    !
    !          - Output -
    ! uz       U component of the wind at z=zu
    ! vz       V component of the wind at z=zu
    ! tz       temperature in kelvins at z=zt
    ! qz       specific humidity at z=zt
    !
    !          - Input -
    ! ni       number of points to process
    ! angi     wind direction at z=za (rad)
    ! tg       temperature at the surface (z=0) in Kelvins
    ! qg       specific humidity
    ! ps       surface pressure at the surface
    ! ilmo     inverse of Monin-Obukhov lenth
    ! h        height of boundary layer
    ! ue       friction velocity
    ! z0       roughness lenth for winds
    ! z0t      roughness lenth for temperature and moisture
    ! ftemp    temperature flux at surface
    ! fvap     vapour flux at surface
    ! za       height of the lowest momentum level (m)
    ! zu       height for computation of wind components
    ! zt       height for computation of temperature and moisture
    ! lat      latitude
    !
    real, parameter :: ANGMAX = 0.85
    !
    integer j
    real fh,fm,x0,x1,y0,y1,h1,h2,h3,hh
    real ct,ctu,cm,vits
    real dang,ang,hi
    !
    REAL, dimension(ni) :: lzz0t
    REAL, dimension(ni) :: lzz0

    !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     STEP 1: Initialize neutral stability functions (lzz0 and lzz0t)
    !
    do j=1,ni
       lzz0 (j) = 1+zu(j)/z0(j)
       lzz0t(j) = 1+zt(j)/z0t(j)
    enddo
    call vslog(lzz0t,lzz0t,ni)
    call vslog(lzz0 ,lzz0 ,ni)
    !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     STEP 2: Compute diagnostic fields, given the fluxes (ftemp,fvap,ue):
    !
    !     (a) for tz and qz, this is done by reversing the flux equations
    !       ftemp = -ctu*((tz + g*zt/cp) - tg) ==> tz = tg - ftemp/ctu - g*zt/cp
    !       fvap  = -ctu*( qz            - qg) ==> qz = qg - fvap /ctu
    !
    !     (b) for vits, this is done reversing the friction velocity equation
    !       ue = cm*vits ==> vits = ue/cm
    !
    do j=1,ni
       !
       !---- integrated stability functions and transfer coefficients
       !     for the diagnostic layers
       if ( ilmo(j).gt.0.) then
          h1    = (za(j)+10.*z0(j))*factn
          h2    = h(j)
          h3    = factn/(4*D97_AS*beta*ilmo(j))
          hh    = MAX(hmin,h1,h2,h3)
       else
          hh    = h(j)
       endif
       !
       fm    =        lzz0 (j) &
            + sf_momentum(zu(j)+z0(j),ilmo(j),hh,x1) &
            - sf_momentum(      z0(j),ilmo(j),hh,x0)
       fh    = beta*( lzz0t(j) &
            + sf_heat(zt(j)    +z0t(j),ilmo(j),hh,y1) &
            - sf_heat(          z0t(j),ilmo(j),hh,y0))
       ct    = karman/fh
       cm    = karman/fm
       !
       !---- diagnostic temperature tz and specific humidity qz
       !
       ctu   = ct*ue(j)
       tz(j) = tg(j) - ftemp(j)/ctu - grav/cpd*zt(j)
       qz(j) = qg(j) - fvap (j)/ctu
       !
       !---- diagnostic wind speed
       !
       vits  = ue(j)/cm
       !
       !---- diagnostic wind components
       !     Note: We assume that the wind direction changes in the
       !           stable layer only, according to the formula
       !             (ang2-ang1)/ang_max = - (z2-z1)/h * sin(lat)
       !
       !
       !     (1) angle difference between levels zu and za (dang)
       if (ilmo(j).gt.0.) then
          hi = 1./MAX(hmin,hh)
       else
          hi = 0.
       endif
       dang = (za(j)-zu(j)) * hi * ANGMAX * SIN(lat(j))
       !
       !     (2) angle at diagnostic level zu
       ang = angi(j) + dang
       !
       !     (3) wind components at diagnostic level zu
       uz(j) = vits*COS(ang)
       vz(j) = vits*SIN(ang)
       !
    enddo
    !
    return
  end subroutine diasurf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sf_heat(F_z,F_ilmo,F_h,F_deriv) result(F_nfh)
    implicit none
!!!#include <arch_specific.hf>
    !@objective Compute stability function for heat/moisture and its derivative.
    !@arguments
    real, intent(in) :: F_z                 !height (m)
    real, intent(in) :: F_ilmo              !inverse of Monin-Obukhov length (1/m)
    real, intent(in) :: F_h                 !height of the PBL (m)
    real, intent(out) :: F_deriv            !derivative of stability function w.r.t. log(F_ilmo)
    real :: F_nfh                           !primitive of stability function (excluding neutral term)
    !@author  A. Zadra, 2011-10
    !@revisions
    !  2011-10, A. Zadra; original code
    !@description
    !   Fusion of stability functions for heat based on PSI and FHI,
    !   also outputing its derivative w.r.t. log(ILMO).
    !   Stable branch: default based on Delage 1997, BLM 82, 23-48
    !   Unstable branch: default based on Delage & Girard 1992, BLM 58, 19-31
    !*@/

    REAL hi,a,b,c,d,ah
    !
    if (F_ilmo.ge.0.) then
       !------ stable branch
       select case (sl_stabfunc_stab)
       case ('LOCK07')
          ah = L07_AH
          F_nfh = (1.+ 0.5*ah*F_z*F_ilmo)**2
          F_deriv = ah*F_z*F_ilmo*(1.+ 0.5*ah*F_z*F_ilmo)
       case ('BELJAARS91')
          a = BH91_A
          b = BH91_B
          c = BH91_C
          d = BH91_D
          F_nfh = (1.+b*a*F_z*F_ilmo)**1.5 - 1 + b*(F_z*F_ilmo-c/d)*EXP(-d*F_z*F_ilmo) + b*c/d
          F_deriv = a*F_z*F_ilmo*(1.+b*a*F_z*F_ilmo)**0.5 + b*F_z*F_ilmo*EXP(-d*F_z*F_ilmo) &
               - d*F_z*F_ilmo*b*(F_z*F_ilmo-c/d)*EXP(-d*F_z*F_ilmo)
       case ('DELAGE97')
          hi = 1./F_h
          d  = 4*D97_AS*beta*F_ilmo
          c  = d*hi - hi**2
          b  = d - 2*hi
          a  = SQRT(1 + b*F_z - c*F_z**2)
          F_nfh = 0.5*( a - F_z*hi - LOG(1+b*F_z*0.5+a) &
               - b/(2*SQRT(c))*ASIN((b-2*c*F_z)/d) )
          F_deriv = 0.5*( 1/(2*a)*((d*F_z-b*hi/c)*(1-F_z*hi) &
               - d*F_z*(1-F_z*hi+a)/(1+0.5*b*F_z+a)) &
               -(d**2*hi/(4*c**1.5))*ASIN((b-2*c*F_z)/d) )
       case DEFAULT
          call msg_toall(MSG_WARNING, '(sf_heat) unknown SL stable function '// &
               trim(sl_stabfunc_stab))
          F_nfh = 1.; F_deriv = 0.
          return
       end select

    else
       !------ unstable branch
       a = (1-DG92_CI*F_z*beta*F_ilmo)**(0.33333333)
       select case (sl_stabfunc_unstab)
       case ('DYER74')
          a = (1-16*F_z*F_ilmo)**(0.5)
          F_nfh = - 2.*LOG(a+1)
          F_deriv = (1./a) - 1
       case ('DELAGE92')
          F_nfh = -1.5*LOG(a**2+a+1)+DG92_RAC3*ATAN((2*a+1)/DG92_RAC3)
          F_deriv = (1./a) - 1
       case DEFAULT
          call msg_toall(MSG_WARNING, '(sf_heat) unknown SL unstable function '// &
               trim(sl_stabfunc_unstab))
          F_nfh = 1.; F_deriv = 0.
          return
       end select

    endif

    return
  end function sf_heat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sf_momentum(F_z,F_ilmo,F_h,F_deriv) result(F_nfm)
    implicit none
!!!#include <arch_specific.hf>
    !@objective Compute stability function for momentum and its derivative.
    !@arguments
    real, intent(in) :: F_z                 !height (m)
    real, intent(in) :: F_ilmo              !inverse of Monin-Obukhov length (1/m)
    real, intent(in) :: F_h                 !height of the PBL (m)
    real, intent(out) :: F_deriv            !derivative of stability function w.r.t. log(F_ilmo)
    real :: F_nfm                           !primitive of stability function (excluding neutral term)
    !@author  A. Zadra, 2011-10
    !@revisions
    !  2011-10, A. Zadra; original code
    !@description
    !   Fusion of stability functions for momemtum based on PSI and FMI,
    !   also outputing the derivative w.r.t. the log(ILMO).
    !   Stable branch: default based on Delage 1997, BLM 82, 23-48
    !   Unstable branch: default based on Delage & Girard 1992, BLM 58, 19-31
    !*@/

    REAL hi,a,b,c,d,am
    !
    if (F_ilmo.ge.0.) then
       !------ stable branch
       select case (sl_stabfunc_stab)
       case ('LOCK07')
          am = L07_AM
          F_nfm = am*F_z*F_ilmo
          F_deriv = am*F_z*F_ilmo
       case ('BELJAARS91')
          a = BH91_A
          b = BH91_B
          c = BH91_C
          d = BH91_D
          F_nfm = a*F_z*F_ilmo + b*(F_z*F_ilmo-c/d)*EXP(-d*F_z*F_ilmo) + b*c/d
          F_deriv = a*F_z*F_ilmo + b*F_z*F_ilmo*EXP(-d*F_z*F_ilmo) &
               - d*F_z*F_ilmo*b*(F_z*F_ilmo-c/d)*EXP(-d*F_z*F_ilmo)
       case ('DELAGE97')
          hi = 1./F_h
          d  = 4*D97_AS*beta*F_ilmo
          c  = d*hi - hi**2
          b  = d - 2*hi
          a  = SQRT(1 + b*F_z - c*F_z**2)
          F_nfm = 0.5*( a - F_z*hi - LOG(1+b*F_z*0.5+a) &
               - b/(2*SQRT(c))*ASIN((b-2*c*F_z)/d) )
          F_deriv = 0.5*( 1/(2*a)*((d*F_z-b*hi/c)*(1-F_z*hi) &
               - d*F_z*(1-F_z*hi+a)/(1+0.5*b*F_z+a)) &
               -(d**2*hi/(4*c**1.5))*ASIN((b-2*c*F_z)/d) )
       case DEFAULT
          call msg_toall(MSG_WARNING, '(sf_momentum) unknown SL stable function '// &
               trim(sl_stabfunc_stab))
          F_nfm = 1.; F_deriv = 0.
          return
       end select

    else
       !------ unstable branch
       a =(1-DG92_CI*F_z*beta*F_ilmo)**(0.16666666)
       select case (sl_stabfunc_unstab)
       case ('DYER74')
          a = (1-16*F_z*F_ilmo)**(0.25)
          F_nfm = - 2.*LOG(a+1) - LOG(a**2+1) + 2.*ATAN(a)
          F_deriv = (1./a) - 1
       case ('DELAGE92')
          F_nfm = - LOG( (a+1)**2*SQRT(a**2-a+1)*(a**2+a+1)**1.5 ) &
               +DG92_RAC3*ATAN((a**2-1)/(DG92_RAC3*a))
          F_deriv = (1./a) - 1
       case DEFAULT
          call msg_toall(MSG_WARNING, '(sf_momentum) unknown SL unstable function '// &
               trim(sl_stabfunc_unstab))
          F_nfm = 1.; F_deriv = 0.
          return
       end select

    endif

    return
  end function sf_momentum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sf_pblheight(F_zu,F_z0,F_u,F_ilmo,F_fcor,F_fm) result(F_nhb)
    implicit none
!!!#include <arch_specific.hf>
    !@objective Compute the planetary boundary layer height.
    !@arguments
    real, intent(in) :: F_zu                !height of wind input (m)
    real, intent(in) :: F_z0                !roughness length for momentum (m)
    real, intent(in) :: F_u                 !wind speed at F_zu (m/s)
    real, intent(in) :: F_ilmo              !inverse of Monin-Obukhov length (1/m)
    real, intent(in) :: F_fcor              !Coriolis factor (1/s)
    real, intent(in) :: F_fm                !integrated stability function for momemtum
    real :: F_nhb                           !height of the PBL (m)
    !@author  A. Zadra, 2011-10
    !@revisions
    !  2011-10, A. Zadra; original code
    !@description
    !   Compute the planetary boundary layer height for
    !   both stable and unstable cases.
    !*@/
    real,    parameter :: BS = 1.
    real h1,h2,h3,cormin,f
    !
    cormin = 0.7e-4
    !
    f = MAX(ABS(F_fcor),cormin)
    !
    if (F_ilmo.ge.0.) then
       !------ stable branch
       h1 = (F_zu+10.*F_z0)*factn
       h2 = BS*SQRT(karman*F_u/(F_ilmo*f*F_fm))
       h3 = factn/(4*D97_AS*beta*F_ilmo)
       !
       F_nhb = MAX(hmin,h1,h2,h3)
    else
       !------ unstable branch
       h1 = 0.3*(F_u*karman/F_fm)/f
       !
       F_nhb = MAX(hmin,h1)
    endif
    !
    return
  end function sf_pblheight

end module sfclayer_mod
