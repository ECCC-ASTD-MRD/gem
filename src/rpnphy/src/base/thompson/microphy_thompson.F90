module microphy_thompson
   use, intrinsic :: iso_fortran_env
   use phymem, only: phyvar
   use tdpack, only: GRAV, RGASD, EPS1, RAUW
   use phybusidx, only : mg,me_moins,dxdy,glsea
   use phymem, only: phyvar
   use machine, only : kind_phys
   use mpi_f08, only : MPI_Comm, MPI_COMM_WORLD
   use mp_thompson, only : mp_thompson_run, mp_thompson_init
   use module_mp_thompson, only: thompson_init, modelinput
   use phy_status, only: PHY_OK, PHY_ERROR, physeterror
#include "phymkptr.hf"
   implicit none
   private
   ! Public constants
   integer, parameter, public :: THOMPSON_ERROR  = -1
   integer, parameter, public :: THOMPSON_OK     = 0
   ! Public procedures
   public :: thompson_wrapper_gem !Thompson condensation scheme wrapper
   public :: thompson_wrapper_init
   public :: thompson_phybusinit  !Define bus requirements
   public :: thompson_lwc         !Compute liquid water content
   public :: thompson_iwc         !Compute ice water content  

   ! Private constants
   integer, parameter :: IMP_PHYSICS = 1
   logical, parameter :: NO_AEROSOL  = .false.
   logical, parameter :: CONVERT_DRY_RHO = .false.
   type(MPI_Comm), parameter :: mpicomm =   MPI_COMM_WORLD

   !==========
   ! MPI stuff
   !----------
   ! mpicomm, mpirank, mpiroot should not be used if MPI is undef.
   ! However, there is a test on MPI parameter blkno and it must be 1
   ! for the scemne to run I think.
   ! In order suppress message like the following, mpiroot must differ from mpirank
   ! Thompson MP is using N inner loops per time step with an effective time step of  xx.xx seconds
   ! If you want to see the message, make mpiroot and mpirank equal
   integer, parameter :: mpirank = 0
   integer, parameter :: mpiroot = -1
   integer, parameter :: blkno   = 1
   integer, parameter :: threads = 1

    ! Private vars
    logical, save :: is_initialized_L = .false.

contains
   
   subroutine thompson_wrapper_init(F_basepath)
      implicit none
      character(len=*), intent(in) :: F_basepath
      integer :: errflg
      character(len=256) :: errmsg_S

      if (is_initialized_L) return

      modelinput = F_basepath
      call thompson_init(is_aerosol_aware_in=NO_AEROSOL,              &
           merra2_aerosol_aware_in=NO_AEROSOL,      &
           mpicomm=mpicomm, mpirank=mpirank, mpiroot=mpiroot, &
           threads=threads, errmsg=errmsg_S, errflg=errflg)
     
      if(errflg /= 0)then
         call physeterror('thompson_wrapper_init', 'call to thompson_init error: '//trim(errmsg_S))
         return
      end if
      return
   end subroutine thompson_wrapper_init
   

  integer function thompson_wrapper_gem(F_pvars,&
       F_dt_inner,F_sedi_semilag_L,F_decfl,F_cldfrac_S,&
       F_ni,F_nk,F_qv,F_qc,F_qr,F_qi,F_qn,F_qg,F_nice,F_nr,&
       F_tt,F_sigma,F_ps,F_zm,F_ww,F_dt,F_kount,&
       F_liquid_rt,F_solid_rt,F_tttnd,F_qvtnd,F_qctnd,F_qrtnd,F_qitnd,F_cldfrac) result(status)
    implicit none    
#include <arch_specific.hf>
#include <rmnlib_basics.hf>
    ! Subroutine Input Units
    ! All units are in MKS (meter-kilogram-second) system.
    ! Arguments
    type(phyvar), pointer, contiguous, intent(in) :: F_pvars(:)
    real, intent(in) ::    F_dt_inner                     !Thompson dt_inner substep length (-1 for model time step)
    logical, intent(in) :: F_sedi_semilag_L               !Thompson Semi-Lagrangian Sedimentation switch
    integer, intent(in) :: F_decfl                        !Thompson deformation CFL
    character(len=*), intent(in) :: F_cldfrac_S           !Type of Ad Hoc cloud fraction
    integer, intent(in) :: F_ni                           !Number horizontal grid points
    integer, intent(in) :: F_nk                           !Number of levels
    real,target, dimension(:,:), intent(in)    :: F_qv    !(mass of vapor         ) / (mass of air + vapor + condensates) [kg/kg]
    real,target, dimension(:,:), intent(in)    :: F_qc    !(mass of cloud water   ) / (mass of air + vapor + condensates) [kg/kg]
    real,target, dimension(:,:), intent(in)    :: F_qr    !(mass of rainwater     ) / (mass of air + vapor + condensates) [kg/kg]
    real,target, dimension(:,:), intent(inout) :: F_qi    !(mass of cloud ice     ) / (mass of air + vapor + condensates) [kg/kg]
    real,target, dimension(:,:), intent(inout) :: F_qn    !(mass of snow          ) / (mass of air + vapor + condensates) [kg/kg]
    real,target, dimension(:,:), intent(inout) :: F_qg    !(mass of graupel       ) / (mass of air + vapor + condensates) [kg/kg]
    real,target, dimension(:,:), intent(inout) :: F_nice  !(number of ice crystals) / (mass of air + vapor + condensates) [#/kg]
    real,target, dimension(:,:), intent(inout) :: F_nr    !(number of raindrops   ) / (mass of air + vapor + condensates) [#/kg]
    real,target, dimension(:,:), intent(in)    :: F_tt    !Temperature grid values [K]
    real, dimension(:,:), intent(in)    :: F_sigma        !Sigma p/ps [dimensionless]
    real, dimension(:), intent(in)      :: F_ps           !Surface pressure [pa]
    real, dimension(:,:), intent(in)    :: F_zm           !Height above ground of the model momentum levels [m]
    real, dimension(:,:), intent(in)    :: F_ww           !Wind speed along-z [m/s]
    real, intent(in)    :: F_dt                           !Time step [s]
    integer, intent(in) :: F_kount                        !Step number
    real, dimension(:), intent(inout) :: F_liquid_rt      !Liquid precipitation rate [kg m-2 s-1] [m s-1]
    real, dimension(:), intent(inout) :: F_solid_rt       !Solid  precipitation rate [kg m-2 s-1] [m s-1]
    real, dimension(:,:), intent(out) :: F_tttnd          !Temperature tendency (K/s)
    real, dimension(:,:), intent(out) :: F_qvtnd          !Specific humidity tendency (kg/kg/s)
    real, dimension(:,:), intent(out) :: F_qctnd          !Cloud water content mass fraction tendency (kg/kg/s)
    real, dimension(:,:), intent(out) :: F_qrtnd          !Rainwater content   mass fraction tendency (kg/kg/s)
    real, dimension(:,:), intent(out) :: F_qitnd          !Ice content         mass fraction tendency (kg/kg/s)
    real, dimension(:,:), intent(out) :: F_cldfrac        !Cloud Fraction Ad Hoc simple calculation not from Thompson scheme.
    
    ! Local variables
    ! Note kind_phys comes from file Thompson code machine.F90 
    integer :: i,k,spp_mp,n_var_spp,errflg
    integer, dimension(F_ni) :: islmsk
    real :: dt_inner
    real, dimension(F_ni,F_nk) :: pres,omega,factor,refl_10cm,qv,qc,qr,qi,qn,qg,nice,nr,tt,ww
    real, dimension(F_ni,F_nk+1) :: phi
    real, pointer, dimension(:), contiguous :: zmg,zme,zglsea
    real, dimension(F_ni) :: prcp,rain,graupel,ice,snow,sr,max_hail_diam_sfc,area
    real, pointer, contiguous:: aerfld(:,:,:)
    real, pointer, contiguous :: zdxdy(:)
    real,target :: dummy(1,1,1)
    real :: qsmall,switch
    logical :: fullradar_diag_L, &
         do_radar_ref_L,ext_diag_L,reset_diag3d_L,cplchm_L,restart_L
    character(len=256) :: errmsg_S

    status = THOMPSON_ERROR
    
    nullify(zmg,zme,aerfld,zdxdy)

    !==========================================================================
    ! Variables F_tt, F_qv, F_qc and F_qr are not changed in this routine.
    ! Their tendenties, F_tttnd, F_qvtnd, F_qctnd and F_qrtnd respectively, are
    ! computed in this routine and passed trough the interface to be
    ! applies in condensation.F90
    ! Variables F_qi, F_qn and F_qg are updated in this routine. A total
    ! tendency for these are computed in F_qitnd. This is only use
    ! for conservation
    !--------------------------------------------------------------------------

    !================================
    ! Copy Temperature in local array
    !--------------------------------
    call my_invert2(tt,F_tt)
    
    !=======================
    ! Compute pressure in Pa
    !-----------------------
    
    do k=1,F_nk
       pres(:,F_nk-k+1)=F_sigma(:,k)*F_ps(:)
    end do
    
    !==================
    ! convert_dry_rho_L
    !------------------
    ! RPN phy uses the mass ratio
    !
    ! Mass ratio of substance i
    !       mi
    ! qi = ----, mT = ma + mv + sum(mj)
    !       mT
    !
    ! where mi is the mass of substance i
    !       ma is the mass of dry air
    !       mv is the mass of vapor
    !       sum(mj) the sum of all condensates:
    !          mc   mass of cloud water
    !          mice mass of cloud ice (pristine)
    !          mr   mass of rain
    !          mn   mass of snow
    !          mg   mass of gropels
    !
    ! Thompson scheme uses mixing ratio (except for vapor see below)
    !       mi 
    ! ri = ----
    !       ma 
    !
    ! Passing from qi to ri
    !       mi        qi * mT                qi                      qi
    ! ri = ---- = ----------------- = ---------------------- = ----------------
    !       ma    mT - mv - sum(mj)   1 - mv/mT - sum(mj/mT)   1 - qv - sum(qj)
    !
    ! Passing from ri to qi
    !       mi        ri * ma                ri                      ri
    ! qi = ---- = ----------------- = ---------------------- = ----------------
    !       mT    ma + mv + sum(mj)   1 + mv/ma + sum(mj/ma)   1 + rv + sum(rj)
    !
    ! If convert_dry_rho_L=.true., thompson_run will do the following conversion:
    !    convert qv and the other condentates from [(mass of substance x) / (dry air + vapor)] to mixing ratio [(mass of substance x) / (dry air)]
    !    Therefore, to use convert_dry_rho_L=.true. all must be converted in units of [(mass of substance x) / (dry air + vapor)] in this routine
    !
    ! If convert_dry_rho_L=.false. thompson_run will do the following conversion:
    !    convert qv from [(mass of vapor) / (dry air + vapor)] to mixing ratio [(mass of vapor) / (dry air)]
    !    The rest must aldeaddy be in mixing ratio [(mass of substance x) / (dry air)] and no conversion will be done in thompson_run for these.
    !    Therefore, to use convert_dry_rho_L=.false. qv must be converted to [(mass of vapor) / (mass of dry air + vapor)] in this routine
    !    and the rest to mixing ratio [(mass of substance x) / (dry air)].
    !
    ! Vertical velocity
    !------------------
    ! Compute omega in Pa s-1 from vertical velocity in m s-1.
    ! Here we inverse the formula used in mp_thompson_run to get ww from omega.
    ! To compute rho like in the scheme, humidity must be in the form of mixing ratio
    ! omega = -w * g * rho
    
    call my_invert2(qv,F_qv)
    call my_invert2(qc,F_qc)
    call my_invert2(qr,F_qr)
    call my_invert2(qi,F_qi)
    call my_invert2(qn,F_qn)
    call my_invert2(qg,F_qg)
    call my_invert2(nice,F_nice)
    call my_invert2(nr,F_nr)
    ! First we put all in mixing ratio
    factor = 1./(1.-(qv+qc+qr+qi+qn+qg))
    qv=qv*factor
    qc=qc*factor
    qr=qr*factor
    qi=qi*factor
    qn=qn*factor
    qg=qg*factor
    nice=nice*factor
    nr=nr*factor

    call my_invert2(ww,F_ww)
    omega = -ww*GRAV*EPS1*pres/(RGASD*tt*(qv+EPS1))
    ! Second, Thompson expects qv to be in units of [ (mass of vapor) / (dry air + vapor) ]
    ! even if flag convert_dry_rho_L is set to false!
    qv=qv/(1.0+qv)

    !========
    ! Aerosol
    !--------
    if (NO_AEROSOL) then
       call physeterror('thompson_wrapper_gem','allocate array nc,nwfa,nifa,nwfa2d,nifa2d,aerfld')
       return
    else
       ! Not allocating aerfld makes the model crash in debug mode so allocating to shut it up
       aerfld => dummy
    endif

    !====================
    ! Land sea mask in mg
    !--------------------
    ! islmsk,  long_name = sea/land/ice mask (=0/1/2)
    ! I think the land sea mask if only used for is_aerosol_aware_L
    MKPTR1D(zmg, mg, F_pvars)
    MKPTR1D(zglsea, glsea, F_pvars)
    do i=1,F_ni
       if(zmg(i)>=.5)then
          islmsk(i)=1
       else
          islmsk(i)=0
       end if
       if(zglsea(i)>=.5)islmsk(i)=2
    end do

    !==============
    ! Level heights
    !--------------
    ! In the Thompson scheme, phi must be defined at the layer interfaces, and there must be F_nk+1 of them
    ! in order to compute the thickness of the layer. The thermodynamic variables are
    ! located at the center of the layers, while the interfaces correspond to the momentum levels
    ! `me_moins`. These are heights above the surface, but tests in the Thompson scheme use actual heights
    ! for the nucleation code. Therefore, we add the surface heights to the profiles.
    ! Phi must be in units of m²/s².
    ! Note that `me_moins`, (the height of the surface times grav) is already in units of m²/s².
    ! Compute phi = g * z with respect to see level.
    MKPTR1D(zme, me_moins, F_pvars)
    phi(:,1)=zme(:)
    do k=1,F_nk
       phi(:,k+1)=F_zm(:,F_nk-k+1)*GRAV + zme(:)
    end do

    !============
    ! Substepping
    !------------
    ! There are two possible substepping methods that are equivalent.
    ! * First, using dt_inner
    !      When using dt_inner, nsteps must be set to 1.
    !      The scheme will iterate multiple times to complete the model time step.
    ! * Second, using nsteps
    !      When using nsteps, dt_inner must be set to the model time step.
    !      The scheme must be called nsteps times with istep being the iteration index.
    !      first_time_step_L must be set to true when istep == 1 and false afterward.
    ! Tests have shown that these two options lead to the same results.
    ! Therefore, only the dt_inner option is available to the user via the physics_cfgs namelist
    ! under the name thompson_dt_inner.
    !
    ! Note that there is an internal automatic substepping for hydrometeor sedimentation that
    ! is independent of the above parameters. If sedi_semilag_L is set to true, the number of
    ! iteration can be controller via decfl. sedi_semilag_L and decfl are part of
    ! physics_cfgs namelist under the names of thompson_sedi_semilag_L and thompson_decfl
    ! A larger CFL will require fewer iterations but the solution may lose precision.
    ! Parameter decfl is not used if sedi_semilag_L is set to .false.
    
    if(F_dt_inner < 0.)then
       dt_inner=F_dt
    else       
       if(F_dt_inner > F_dt)then
          call physeterror('thompson_wrapper_gem','thompson_dt_inner from physics_cfgs must be smaler equal to the model time step')
          return
       endif
       ! Check if F_dt_inner is a divisor of F_dt
       if( F_dt/real(max(nint(F_dt/F_dt_inner),1)) /= F_dt_inner)then
          call physeterror('thompson_wrapper_gem','The value of thompson_dt_inner from physics_cfgs must evenly divide the model time step, but it does not.')
          return
       endif
       dt_inner=F_dt_inner
    end if

    !=========================
    !Set precipitation to zero
    !-------------------------
    prcp=0.D0; rain=0.D0; graupel=0.D0; ice=0.D0; snow=0.D0; sr=0.D0  !#TODO: check every step????

    !======
    ! Radar
    !------
    fullradar_diag_L=.false.; do_radar_ref_L=.false.
    
    !=====
    ! Diag
    !-----
    ext_diag_L=.false.; reset_diag3d_L=.false.
    if(ext_diag_L)then
       call physeterror('thompson_wrapper_gem','Creat and allocate array diag3d and add it to call')
       return
    endif

    !=============================
    ! Optional Random perturbation
    !-----------------------------
    ! If the optional parametres are not used, spp_mp, n_var_spp should not be used either
    spp_mp=0; n_var_spp=0

    !=============
    ! Coupled Chem
    !--------------
    cplchm_L=.false.
    if(ext_diag_L)then
       call physeterror('thompson_wrapper_gem','Creat and allocate arrays pfi_lsan and pfl_lsan an add it to call')
       return
    endif

    !===============
    ! Initialization
    !---------------
    ! Do it once at beging of integration
    ! imp_physics and imp_physics_thompson are integer that are compared and
    ! must be equal, so put 1 for both
    ! Parameter area is not used but not optional, so we put it
    MKPTR1D(zdxdy, dxdy, F_pvars)
    area=zdxdy
    
    if (.not.is_initialized_L) then

       restart_L = (F_kount == 0)
       call mp_thompson_init(&
            ncol=F_ni, nlev=F_nk, con_g=GRAV, con_rd=RGASD, con_eps=EPS1,      &
            restart=restart_L, imp_physics=IMP_PHYSICS,                        &
            imp_physics_thompson=IMP_PHYSICS, convert_dry_rho=CONVERT_DRY_RHO, &
            spechum=qv, qc=qc, qr=qr, qi=qi, qs=qn, qg=qg, ni=nice, nr=nr,     &
            is_aerosol_aware=NO_AEROSOL,                                       &
            merra2_aerosol_aware=NO_AEROSOL,                                   &
            tgrs=tt, prsl=pres, phil=phi, area=area,                           &
            aerfld=aerfld, mpicomm=mpicomm, mpirank=mpirank, mpiroot=mpiroot,  &
            threads=threads, ext_diag=ext_diag_L,                                    &
            is_initialized=is_initialized_L, errmsg=errmsg_S, errflg=errflg)
       if(errflg /= 0)then
          call physeterror('thompson_wrapper_gem', 'call to mp_thompson_init error: '//trim(errmsg_S))
          return
       end if
    endif
   
    !print*,'F_dt=',F_dt,', dt_inner=',dt_inner,'sedi_semilag=',F_sedi_semilag_L,',F_decfl=',F_decfl,', nsteps=',nsteps,','
    
    call mp_thompson_run(&
         ncol=F_ni, nlev=F_nk, con_g=GRAV, con_rd=RGASD,                          &
         con_eps=EPS1, convert_dry_rho=CONVERT_DRY_RHO,                           &
         spechum=qv, qc=qc, qr=qr, qi=qi, qs=qn, qg=qg, ni=nice, nr=nr,           &
         is_aerosol_aware=NO_AEROSOL,                                             &
         merra2_aerosol_aware=NO_AEROSOL,                                         &
         aero_ind_fdb=NO_AEROSOL,                                                 &
         tgrs=tt, prsl=pres, phii=phi, omega=omega,                               &
         sedi_semi=F_sedi_semilag_L, decfl=F_decfl, islmsk=islmsk, dtp=F_dt,      &
         dt_inner=dt_inner,                                                       &
         first_time_step=.true., istep=1, nsteps=1,                               &
         prcp=prcp, rain=rain, graupel=graupel, ice=ice, snow=snow, sr=sr,        &
         refl_10cm=refl_10cm, fullradar_diag=fullradar_diag_L,                    &
         max_hail_diam_sfc=max_hail_diam_sfc,                                     &
         do_radar_ref=do_radar_ref_L, aerfld=aerfld,                              &
         mpicomm=mpicomm, mpirank=mpirank, mpiroot=mpiroot, blkno=blkno,          &
         ext_diag=ext_diag_L, reset_diag3d=reset_diag3d_L,                        &
         spp_mp=spp_mp, n_var_spp=n_var_spp,                                      &
         cplchm=cplchm_L,                                                         &
         is_initialized=is_initialized_L, errmsg=errmsg_S, errflg=errflg)
    
    if(errflg /= 0)then
       call physeterror('thompson_wrapper_gem','call to mp_thompson_run error: '//trim(errmsg_S))
       return
    end if
    
!!$      subroutine mp_thompson_run(ncol, nlev, con_g, con_rd,        &
!!$                              con_eps, convert_dry_rho,            &
!!$                              spechum, qc, qr, qi, qs, qg, ni, nr, &
!!$                              is_aerosol_aware,                    &
!!$                              merra2_aerosol_aware, nc, nwfa, nifa,&
!!$                              nwfa2d, nifa2d, aero_ind_fdb,        &
!!$                              tgrs, prsl, phii, omega,             &
!!$                              sedi_semi, decfl, islmsk, dtp,       &
!!$                              dt_inner,                            &
!!$                              first_time_step, istep, nsteps,      &
!!$                              prcp, rain, graupel, ice, snow, sr,  &
!!$                              refl_10cm, fullradar_diag,           &
!!$                              max_hail_diam_sfc,                   &
!!$                              do_radar_ref, aerfld,                &
!!$                              mpicomm, mpirank, mpiroot, blkno,    &
!!$                              ext_diag, diag3d, reset_diag3d,      &
!!$                              spp_wts_mp, spp_mp, n_var_spp,       &
!!$                              spp_prt_list, spp_var_list,          &
!!$                              spp_stddev_cutoff,                   &
!!$                              cplchm, pfi_lsan, pfl_lsan,          &
!!$                              is_initialized, errmsg, errflg)

    ! Thompson is returning qv in units of  [ (mass of vapor) / (dry air + vapor) ]
    ! even it flag convert_dry_rho_L is set to false!
    ! We convert qv back to mixing ration

    call my_invert1(tt)
    call my_invert1(qv)
    call my_invert1(qc)
    call my_invert1(qr)
    call my_invert1(qi)
    call my_invert1(qn)
    call my_invert1(qg)
    call my_invert1(nice)
    call my_invert1(nr)
    
    qv = qv/(1.0-qv)
    ! Convert all to mass ratio
    
    factor = 1./(1.+(qv+qc+qr+qi+qn+qg))
    qv=qv*factor
    qc=qc*factor       !MPQC !CQ
    qr=qr*factor       !MPQR
    qi=qi*factor       !MPQI
    qn=qn*factor       !MPQS
    qg=qg*factor       !MPQG
    F_nice=nice*factor !MPNI
    F_nr=nr*factor     !MPNR
    
    ! Ad hoc cloud fraction
    qsmall = 1.e-14
    if(trim(F_cldfrac_S) == 'liq_ice_snow')then
       switch=1.
    else if(trim(F_cldfrac_S) == 'liq_ice')then
       switch=0.
    else
       call physeterror('thompson_wrapper_gem','Wrong choice for THOMPSON_cldfrac')
       return
    endif
    do k=1,F_nk
       do i=1,F_ni
          if(qc(i,k)+qi(i,k)+switch*qn(i,k) > qsmall)then
             F_cldfrac(i,k)=1.
          else
             F_cldfrac(i,k)=0.
          endif
       end do
    end do

    ! Precipitation
    ! In the scheme, the precipitation units are [m] accumulation during the model step.
    ! In the RPN physics, the precipitation rate comming out of the condensation schemes
    ! like Thonpson, must be mass flux [kg m-2 s-1].
    ! Therefore we must do the following conversion:
    !
    ! Mutiply the value out of the scheme by the water density 1000 [kg m-3] (tdpack_const RAUW)
    ! [m] [kg m-3] = [kg m-2]
    !
    ! Then we must devide by the time step to get the flux
    ! [kg m-2] [s] = [kg m-2 s-1]

    F_liquid_rt = rain*RAUW/F_dt                  !TLS P2
    F_solid_rt = (graupel + ice + snow)*RAUW/F_dt !TSS P4

    ! Compute tendencies
    F_tttnd=(tt-F_tt)/F_dt !STE
    F_qvtnd=(qv-F_qv)/F_dt !SQE
    F_qctnd=(qc-F_qc)/F_dt !SQCE
    F_qrtnd=(qr-F_qr)/F_dt !SQRE
    ! Ice tendencies is only used to compute conservation.
    F_qitnd=( (qi-F_qi) + (qn-F_qn) + (qg-F_qg) )/F_dt !MPQI
    ! Apply ice snow and graupel tendencies
    F_qi=qi; F_qn=qn; F_qg=qg;
    
    status = THOMPSON_OK
    
  end function thompson_wrapper_gem


  subroutine my_invert1(F_in)
     implicit none
     real, intent(inout) :: F_in(:,:)
     real :: rtmp(size(F_in,1))
     integer :: k, nk
     nk = size(F_in,2)
     do k = 1, nk/2
        rtmp = F_in(:,k)
        F_in(:,k) = F_in(:,nk-k+1)
        F_in(:,nk-k+1) = rtmp
     enddo
     return
  end subroutine my_invert1

  
  subroutine my_invert2(F_out, F_in)
     implicit none
     real, intent(out) :: F_out(:,:)
     real, intent(in)  :: F_in(:,:)
     integer :: nk
     nk = size(F_in,2)
     F_out = F_in(:,nk:1:-1)
     return
  end subroutine my_invert2

 
  ! Define bus requirements
  function thompson_phybusinit() result(F_istat)
    use bus_builder, only: bb_request
    implicit none
    integer :: F_istat                          !Function return status
    F_istat = PHY_ERROR
    if (bb_request((/ &
         'CLOUD_WATER_MASS ', &
         'CLOUD_WATER_NUM  ', &
         'RAIN_MASS        ', &
         'RAIN_NUM         ', &
         'CLOUD_ICE_MASS   ', &
         'CLOUD_ICE_NUM    ', &
         'SNOW_MASS        ', &
         'SNOW_NUM         ', &
         'GRAUPEL_MASS     ', &
         'GRAUPEL_NUM      ', &
         'RATE_PRECIP_TYPES', &
         'ICE_MASS_TEND    ', &
         'REFLECTIVITY     '&
       /)) /= PHY_OK) then
       call physeterror('microphy_thompson::thompson_phybusinit', &
            'Cannot construct bus request list')
       return
    endif
    F_istat = PHY_OK
  end function thompson_phybusinit

  ! Compute total water mass
  function thompson_lwc(F_qltot, F_pvars, F_tminus) result(F_istat)
    use phybusidx
    implicit none
    real, dimension(:,:), intent(out) :: F_qltot        !Total water mass (kg/kg)
    type(phyvar), pointer, contiguous :: F_pvars(:)     !All phy vars (meta + slab data)
    logical, intent(in), optional :: F_tminus           !Compute fields at time-minus [false]
    integer :: F_istat                                  !Return status
    integer :: ni, nkm1
    real, dimension(:,:), pointer, contiguous :: zqc, zqr
    logical :: my_tminus
    F_istat = PHY_ERROR
    my_tminus = .false.
    if (present(F_tminus)) my_tminus = F_tminus
    ni = size(F_qltot, dim=1); nkm1 = size(F_qltot, dim=2)
    if (my_tminus) then
       MKPTR2Dm1(zqc, qcmoins, F_pvars)
       MKPTR2Dm1(zqr, qrmoins, F_pvars)
    else
       MKPTR2Dm1(zqc, qcplus, F_pvars)
       MKPTR2Dm1(zqr, qrplus, F_pvars)
    endif
    F_qltot(:,:) = max(zqc(:,:), 0.) + max(zqr(:,:), 0.)
    F_istat = PHY_OK
    return
  end function thompson_lwc

  ! Compute total ice mass
  function thompson_iwc(F_qitot, F_pvars, F_tminus) result(F_istat)
    use phybusidx
    implicit none
    real, dimension(:,:), intent(out) :: F_qitot        !Total ice mass (kg/kg)
    type(phyvar), pointer, contiguous :: F_pvars(:)     !All phy vars (meta + slab data)
    logical, intent(in), optional :: F_tminus           !Compute fields at time-minus [false]
    integer :: F_istat                                  !Return status
    integer :: ni, nkm1
    real, dimension(:,:), pointer, contiguous :: zqi, zqn, zqg
    logical :: my_tminus
    F_istat = PHY_ERROR
    my_tminus = .false.
    if (present(F_tminus)) my_tminus = F_tminus
    ni = size(F_qitot, dim=1); nkm1 = size(F_qitot, dim=2)
    if (my_tminus) then
       MKPTR2Dm1(zqi, qimoins, F_pvars)
       MKPTR2Dm1(zqn, qnmoins, F_pvars)
       MKPTR2Dm1(zqg, qgmoins, F_pvars)
    else
       MKPTR2Dm1(zqi, qiplus, F_pvars)
       MKPTR2Dm1(zqn, qnplus, F_pvars)
       MKPTR2Dm1(zqg, qgplus, F_pvars)
    endif
    F_qitot(:,:) = max(zqi(:,:), 0.) + max(zqn(:,:), 0.) + &
         max(zqg(:,:), 0.)
    F_istat = PHY_OK
    return
  end function thompson_iwc

end MODULE microphy_thompson
