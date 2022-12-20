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

module mixing_length
   ! Container for calculation and manipulation of mixing and dissipation length scales in the PBL.
   use, intrinsic :: iso_fortran_env, only: INT64
   use phy_status, only: PHY_OK, PHY_ERROR
   use tdpack, only: GRAV, KARMAN, PI, DELTA, CAPPA
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   private

   ! External symbols
   integer, external :: neark
  
   ! Private parameters
   integer, parameter :: ML_LONG=1024              !Long character string
   integer, parameter :: ML_BLAC62=1               !Key for 'blac62' mixing length
   integer, parameter :: ML_BOUJO=2                !Key for 'boujo' mixing length
   integer, parameter :: ML_TURBOUJO=3             !Key for 'turboujo' mixing length
   integer, parameter :: ML_LH=4                   !Key for 'lh' mixing length
   real, parameter :: ML_MIN=1.                    !Minimum mixing length (m)
   
   ! API types
   type closureMap
      sequence
      character(len=ML_LONG) :: name
      integer :: key
   end type closureMap
   
   ! API parameters
   real, parameter, public :: ML_LMDA=200.         !Maximum stable asymptotic mixing length for Blackadar (1962) estimate
   real, parameter, public :: ML_LMDAS=40.         !Minimum stable asymptotic mixing length for Lock (2007) estimate
   real, parameter, public :: ML_LMDA_UNSTAB=5000. !Maximum unstable mixing length for Blackadar (1962) estimate
   type(closureMap), dimension(4), parameter, public :: ML_CLOSURES = (/ & !Available closures and mappings
        closureMap('BLAC62', ML_BLAC62), &
        closureMap('BOUJO', ML_BOUJO), &
        closureMap('TURBOUJO', ML_TURBOUJO), &
        closureMap('LH', ML_LH)/)  !Remember to update the documentation strings in phy_options() if changes are made here
   
   ! Private variables
   logical :: initialized=.false.                  !Package initialization status
   character(len=ML_LONG) :: mlblac_max = 'BLAC62' !Asymptotic formulation for Blacadar (1962) mixing length with clipping

   ! API subprograms
   public :: ml_init
   public :: ml_put
   public :: ml_key
   public :: ml_tfilt
   public :: ml_blend
   public :: ml_calc_blac
   public :: ml_calc_boujo
   public :: ml_calc_lh
   public :: ml_compute

   ! Generic procedures
   interface ml_put
      module procedure ml_put_s
   end interface ml_put

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ml_init(mlen, imlen) result(status)
    ! Initialize the mixing length module and return the integer key

    ! Input arguments
    character(len=*), intent(in) :: mlen                   !Mixing length closure

    ! Output arguments
    integer, intent(out) :: imlen                          !Integer key for mixing length closure
    integer :: status                                      !Return status of function

    ! Initialize return value
    status = PHY_ERROR

    ! Set initialization status
    if (initialized) then
       call msg_toall(MSG_ERROR,'(ml_init) Attempt to reinitialized mixing length package')
       return
    endif
    initialized = .true.
    
    ! Map closure name to key
    imlen = -1
    if (ml_key(mlen, imlen) /= PHY_OK) return

    ! Successful completion of subprogram
    status = PHY_OK
    
  end function ml_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ml_key(mlen, imlen) result(status)
    use clib_itf_mod, only: clib_toupper, CLIB_OK
    ! Map the mixing length closure to an integer key

    ! Input arguments
    character(len=*), intent(in) :: mlen                   !Mixing length closure

    ! Output arguments
    integer, intent(out) :: imlen                          !Integer key for mixing length closure
    integer :: status                                      !Return status of function

    ! Internal variables
    integer :: i
    character(len=ML_LONG) :: ucmlen
    
    ! Initialize return value
    status = PHY_ERROR

    ! Check initialization status
    if (.not.initialized) then
       call msg_toall(MSG_ERROR,'(ml_key) mixing length package not initialized')
       return
    endif
    
    ! Convert value to upper case
    ucmlen = mlen
    if (clib_toupper(ucmlen) /= CLIB_OK) then
       call msg_toall(MSG_ERROR,'(ml_key) cannot convert mixing length to upper case: '//trim(mlen))
       return
    endif

    ! Map closure name to key
    imlen = -1
    i = 1
    do while (imlen < 0 .and. i <= size(ML_CLOSURES))
       if (ucmlen == ML_CLOSURES(i)%name) then
          imlen = ML_CLOSURES(i)%key
       else
          i = i+1
       endif
    enddo
    if (imlen < 0) then
       call msg_toall(MSG_ERROR,'(ml_key) invalid mixing length request: '//trim(mlen))
       return
    endif
       
    ! Successful completion of subprogram
    status = PHY_OK

  end function ml_key
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function ml_put_s(key,val) result(status)
      use clib_itf_mod, only: clib_toupper, CLIB_OK
      ! Set a variable in the mixing length module

      ! Input arguments
      character(len=*), intent(in) :: key                  !Name of key to set
      character(len=*), intent(in) :: val                  !Value of the key

      ! Output arguments
      integer :: status                                    !Return status of function

      ! Local variables
      character(len=len_trim(val)) :: ucval

      ! Initialize return value
      status = PHY_ERROR

      ! Check initialization status
      if (.not.initialized) then
         call msg_toall(MSG_ERROR,'(ml_put_s) mixing length package not initialized')
         return
      endif
      
      ! Convert value to upper case
      ucval = val
      if (clib_toupper(ucval) /= CLIB_OK) then
         call msg_toall(MSG_ERROR,'(ml_put_s) cannot convert input to upper case: '//trim(val))
         return
      endif

      ! Attempt to set value of requested key
      select case (key)
      case ('mlblac_max','MLBLAC_MAX')
         mlblac_max = ucval
      case DEFAULT
         call msg_toall(MSG_ERROR,'(ml_put_s) cannot set value for '//trim(key)//' to '//trim(ucval))
         return
      end select

      ! Successful completion of subprogram
      status = PHY_OK

   end function ml_put_s
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ml_compute(zn, zd, pri, mlen, t, qe, qc, z, gzmom, s, se, ps, &
        enold, buoy, rig, w_cld, f_cs, fm, turbreg, z0, &
        hpbl, lh, hpar, mrk2, dxdy, tau, kount) result(stat)
     use phy_options, only: ilongmel, pbl_diss, pbl_mlturb_diss, timings_L
     use ens_perturb, only: ens_spp_get
     ! High-level function to compute and blend mixing lengths

     ! Argument declaration
     integer, dimension(:), intent(in) :: mlen             !mixing length closure key
     real, dimension(:,:), intent(in) :: t                 !dry air temperature (K)
     real, dimension(:,:), intent(in) :: qe                !specific humidity on e-levels (kg/kg)
     real, dimension(:,:), intent(in) :: qc                !PBL cloud water content (kg/kg)
     real, dimension(:,:), intent(in) :: z                 !height AGL of thermodynamic levels (m)
     real, dimension(:,:), intent(in) :: gzmom             !height AGL of momentum levels (m)
     real, dimension(:,:), intent(in) :: s                 !sigma values for full momentum levels
     real, dimension(:,:), intent(in) :: se                !sigma values for half (energy) levels
     real, dimension(:), intent(in) :: ps                  !surface pressure (Pa)
     real, dimension(:,:), intent(in) :: enold             !TKE from previous time step (m2/s2)
     real, dimension(:,:), intent(in) :: buoy              !buoyancy flux (m2/s2)
     real, dimension(:,:), intent(in) :: rig               !gradient Richardson number
     real, dimension(:,:,:), intent(in) :: w_cld           !cloud-layer velocity scales (m/s)
     real, dimension(:,:), intent(in) :: f_cs              !factor for l_s coeff c_s (nonlocal)
     real, dimension(:,:), intent(in) :: fm                !PBL momentum stability function
     real, dimension(:,:), intent(in) :: turbreg           !turbulence regime flag
     real, dimension(:), intent(in) :: z0                  !momentum roughness length (m)
     real, dimension(:), intent(in) :: hpbl                !depth of the PBL (m)
     real, dimension(:), intent(in) :: lh                  !launching height for gravity waves (m)
     real, dimension(:), intent(in) :: hpar                !height of parcel ascent (m)
     real, dimension(:,:), intent(in) :: mrk2              !Markov chains for stochastic parameters
     real, dimension(:), intent(in) :: dxdy                !horizontal grid area (m2)
     real, intent(in) :: tau                               !time step (s)
     integer, intent(in) :: kount                          !step number
     real, dimension(:,:), intent(inout) :: pri            !inverse Prandtl number
     real, dimension(:,:), intent(inout) :: zn             !mixing length (m)
     real, dimension(:,:), intent(out) :: zd               !dissipation length (m)

     ! Local variables and common blocks
     include "phyinput.inc"
#include "clefcon.cdk"
     integer :: n, nk, j, k
     real, dimension(size(t,dim=1)) :: mlemod, mlemodt, mlmult
     real, dimension(size(t,dim=1),size(t,dim=2)) :: zn_blac, zd_blac, pri_blac, &
          zn_boujo, zd_boujo, pri_boujo, zn_turboujo, zd_turboujo, pri_turboujo, &
          zn_lh, zd_lh, pri_lh, blend_hght, przn, te, tv, qce, weight, &
          rif, znold
     logical :: one_ml_form, mlemod_calc, mlemodt_calc
     logical, dimension(size(t,dim=1),size(t,dim=2)) :: boujo_valid
     
     ! Set return status
     stat = PHY_ERROR

     ! Check initialization status
     if (.not.initialized) then
        call msg_toall(MSG_ERROR,'(ml_compute) mixing length package not initialized')
        return
     endif

     ! Timings for mixing length calculation
     if (timings_L) call timing_start_omp(500, 'mixing_length', 430)
     
     ! Set array bounds
     n = size(t, dim=1)
     nk = size(t, dim=2)

     ! Retain mixing length information from previous step
     znold(:,:) = zn(:,:)
     
     ! Obtain stochastic parameter information
     mlmult(:) = ens_spp_get('ml_mult', mrk2, 1.)
     mlemod(:) = ens_spp_get('ml_emod', mrk2, 0.)
     mlemodt(:) = ens_spp_get('ml_emodt', mrk2, 0.)
     mlemod_calc = any(mlemod(:) /= 0.)     !model based on 'blac62' and 'boujo'
     mlemodt_calc = any(mlemodt(:) /= 0.)   !model based on 'blac62' and 'boujo' below eta=0.7
     one_ml_form = all(mlen(:) == ilongmel) .and. .not.(mlemod_calc .or. mlemodt_calc)

     ! Precompute state variables if required
     if (any(mlen(:) == ML_BOUJO) .or. any(mlen(:) == ML_TURBOUJO) .or. &
          mlemod_calc .or. mlemodt_calc) then
        te(1:n,1:nk)  = t(1:n,1:nk)
        qce(1:n,1:nk) = qc(1:n,1:nk)
        tv = te*(1.0+DELTA*qe-qce)*(se**(-CAPPA))
     endif
     boujo_valid(:,:) = .false.
     zn_boujo(:,:) = -1.

     ! Mixing length estimate based on regime-dependent Bougeault and Lacarrere (1989)
     if (any(mlen(:) == ML_TURBOUJO)) then
        where (nint(turbreg(:,:)) /= LAMINAR)
           boujo_valid(:,:) = .true.
        endwhere
        if (ml_calc_boujo(zn_boujo, tv, enold, w_cld, z, s, ps, &
             mask=boujo_valid, init=.false.) /= PHY_OK) then
           call physeterror('moistke', 'error returned by B-L mixing length estimate')
           return
        endif
        if (ml_calc_blac(zn_blac, rig, w_cld, z, z0, fm, hpbl, lh, przn=przn) /= PHY_OK) then
           call physeterror('moistke', 'error returned by Blackadar mixing length estimate')
           return
        endif
        blend_hght(:,:nk-1) = gzmom(:,2:)
        blend_hght(:,nk) = 0.
        if (ml_blend(zn_turboujo, zn_blac, zn_boujo, blend_hght, s, ps) /= PHY_OK) then
           call physeterror('moistke','error returned by mixing length blending')
           return
        endif
        where (boujo_valid(:,:))
           przn(:,:) = 1.                  ! Nonlocal (Bougeault) mixing length for turbulent flows
        elsewhere
           zn_turboujo(:,:) = zn_blac(:,:) ! Local (Blackadar) mixing length for laminar flows
        endwhere
        do k=1,nk
           zn_turboujo(:,k) = zn_turboujo(:,k) * mlmult(:)
        enddo        
        if (ml_tfilt(zn_turboujo, znold, tau, kount) /= PHY_OK) then
           call physeterror('moistke', 'error returned by mixing length time filtering')
           return
        endif
        pri_turboujo = pri/przn
        rif = pri_turboujo*rig
        zd_turboujo = zn_turboujo * (1.-min(rif,0.4)) / (1.-2.*min(rif,0.4))
        zd_turboujo = max(zd_turboujo,1.e-6)
        if (pbl_mlturb_diss) then
           where (.not.boujo_valid(:,:)) zd_turboujo = max(zn_turboujo,1.E-6)
        endif
        if (one_ml_form) then
           zn = zn_turboujo
           zd = zd_turboujo
           pri = pri_turboujo
        endif
     endif

     ! Mixing length estimate based on Bougeault and Lacarrere (1989)
     if (any(mlen(:) == ML_BOUJO) .or. mlemod_calc .or. mlemodt_calc) then
        ! Do not recompute at points already done by 'turboujo'
        do k=1,nk
           do j=1,n
              if (boujo_valid(j,k)) then
                 boujo_valid(j,k) = .false.
              else
                 if (mlemodt_calc .and. se(j,k) > 0.7) boujo_valid(j,k) = .true.
                 if (mlen(j) == ML_BOUJO .or. mlemod_calc) boujo_valid(j,k) = .true.
              endif
           enddo
        enddo
        if (ml_calc_boujo(zn_boujo, tv, enold, w_cld, z, s, ps, &
             mask=boujo_valid, init=.false.) /= PHY_OK) then
           call physeterror('moistke', 'error returned by B-L mixing length estimate')
           return
        endif
        if (ml_calc_blac(zn_blac, rig, w_cld, z, z0, fm, hpbl, lh, przn=przn) /= PHY_OK) then
           call physeterror('moistke', 'error returned by Blackadar mixing length estimate')
           return
        endif
        where (zn_boujo(:,:) < 0.)
           zn_boujo(:,:) = zn_blac(:,:)
        elsewhere
           boujo_valid(:,:) = .true.
        endwhere
        blend_hght(:,:nk-1) = gzmom(:,2:)
        blend_hght(:,nk) = 0.
        if (ml_blend(zn_boujo, zn_blac, zn_boujo, blend_hght, s, ps) /= PHY_OK) then
           call physeterror('moistke','error returned by mixing length blending')
           return
        endif
        do k=1,nk
           zn_boujo(:,k) = zn_boujo(:,k) * mlmult(:)
        enddo
        przn = 1.
        if (ml_tfilt(zn_boujo, znold, tau, kount) /= PHY_OK) then
           call physeterror('moistke', 'error returned by mixing length time filtering')
           return
        endif
        pri_boujo = pri/przn
        rif = pri_boujo*rig
        zd_boujo = max(zn_boujo * (1.-min(rif,0.4)) / (1.-2.*min(rif,0.4)), 1.e-6)
        if (one_ml_form) then
           zn = zn_boujo
           zd = zd_boujo
           pri = pri_boujo
        endif
     endif

     ! Mixing length estimate based on Lenderink and Holtslag (2004) modified for no-local Cu
     if (any(mlen(:) == ML_LH)) then
        przn = 1.
        pri_lh = pri/przn
        if (ml_calc_lh(zd_lh, zn_lh, pri_lh, buoy, enold, w_cld, rig, z, f_cs, hpar) /= PHY_OK) then
           call physeterror('moistke', 'error returned by L-H mixing length estimate')
           return
        endif
        do k=1,nk
           zd_lh(:,k) = zd_lh(:,k) * mlmult(:)
        enddo
        if (kount > 0) then ! Blend towards Blackadar length scale at upper levels
           if (ml_calc_blac(zn_blac, rig, w_cld, z, z0, fm, hpbl, lh) /= PHY_OK) then
              call physeterror('moistke', 'error returned by Blackadar mixing length estimate')
              return
           endif
           do k=1,nk
              zn_blac(:,k) = zn_blac(:,k) * mlmult(:)
           enddo
           call blweight2(weight, s, ps, n, nk)
           zn_lh(:,:) = zn_blac(:,:) + (zn_lh(:,:)-zn_blac(:,:))*weight(:,:)
           zd_lh(:,:) = zn_blac(:,:) + (zd_lh(:,:)-zn_blac(:,:))*weight(:,:)
           if (ml_tfilt(zn_lh, znold, tau, kount) /= PHY_OK) then
              call physeterror('moistke', 'error returned by mixing length time filtering')
              return
           endif
        endif
        if (one_ml_form) then
           zn = zn_lh
           zd = zd_lh
           pri = pri_lh
        endif
     endif

     ! Mixing length estimate based on Blackadar (1962)
     if (any(mlen(:) == ML_BLAC62) .or. mlemod_calc .or. mlemodt_calc) then
        if (ml_calc_blac(zn_blac, rig, w_cld, z, z0, fm, hpbl, lh, dxdy=dxdy, przn=przn) /= PHY_OK) then
           call physeterror('moistke', 'error returned by Blackadar mixing length estimate')
           return
        endif
        do k=1,nk
           zn_blac(:,k) = zn_blac(:,k) * mlmult(:)
        enddo
        if (ml_tfilt(zn_blac, znold, tau, kount) /= PHY_OK) then
           call physeterror('moistke', 'error returned by mixing length time filtering')
           return
        endif
        zd_blac = max(zn_blac,1.E-6)
        pri_blac = pri/przn
        if (one_ml_form) then
           zn = zn_blac
           zd = zd_blac
           pri = pri_blac
        endif
     endif
     
     ! Merge mixing length estimates if required
     if (.not.one_ml_form) then
        do j=1,n
           if (mlen(j) == ML_TURBOUJO) then
              zn(j,:) = zn_turboujo(j,:)
              zd(j,:) = zd_turboujo(j,:)
              pri(j,:) = pri_turboujo(j,:)
           elseif (mlen(j) == ML_BOUJO) then
              zn(j,:) = zn_boujo(j,:)
              zd(j,:) = zd_boujo(j,:)
              pri(j,:) = pri_boujo(j,:)
           elseif (mlen(j) == ML_LH) then
              zn(j,:) = zn_lh(j,:)
              zd(j,:) = zd_lh(j,:)
              pri(j,:) = pri_lh(j,:)
           elseif (mlen(j) == ML_BLAC62) then
              zn(j,:) = zn_blac(j,:)
              zd(j,:) = zd_blac(j,:)
              pri(j,:) = pri_blac(j,:)
           else
              call physeterror('moistke', 'invalid mixing length ')
              return
           endif
        enddo
     endif

     ! Apply stochastic error models on request
     if (mlemod_calc) then
        do k=1,nk
           zn(:,k) = min( max( zn(:,k) + mlemod(:) * (zn_boujo(:,k) - zn_blac(:,k)), &
                ML_MIN), ML_LMDA_UNSTAB)
        enddo
     endif
     if (mlemodt_calc) then
        do k=1,nk
           where (boujo_valid(:,k))
              zn(:,k) = min( max( zn(:,k) + mlemodt(:) * (zn_boujo(:,k) - zn_blac(:,k)), &
                   ML_MIN), ML_LMDA_UNSTAB)
           endwhere
        enddo
     endif

     ! Adjust dissipation length scale on user request
     if (pbl_diss == 'LIM50') zd = min(zd,50.)

     ! Recycling of length scales
     if (any('zn'==phyinread_list_s(1:phyinread_n))) zn(:,:) = znold(:,:)
     
     ! Completed timings
     if (timings_L) call timing_stop_omp(500)
     
     ! Successful end of subprogram
     stat = PHY_OK
     return
   end function ml_compute
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ml_tfilt(zn,znold,tau,kount) result(stat)
      ! Time filtering of mixing length scale
      use phy_options, only: pbl_zntau

      ! Argument declaration
      integer, intent(in) :: kount                      !time step number
      real, intent(in) :: tau                           !time step (s)
      real, dimension(:,:), intent(in) :: znold         !mixing length from previous step (m)
      real, dimension(:,:), intent(inout) :: zn         !mixing length (m)

      ! Local variables and common blocks
#include "phyinput.inc"

      ! Set return status
      stat = PHY_ERROR

      ! Check initialization status
      if (.not.initialized) then
         call msg_toall(MSG_ERROR,'(ml_tfilt) mixing length package not initialized')
         return
      endif
      
      ! Do not filter if field has been read on this step
      READ_INPUT: if (any('zn'==phyinread_list_s(1:phyinread_n))) then
         zn = znold
         stat = PHY_OK
         return
      endif READ_INPUT

      ! Apply time filter during integration
      if (kount > 0 .and. pbl_zntau > 0.) zn = zn + (znold-zn)*exp(-tau/pbl_zntau)

      ! Successful end of subprogram
      stat = PHY_OK
      return
   end function ml_tfilt

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ml_blend(znblend,zn1,zn2,zz,sigma,ps) result(stat)
      ! Blend two mixing lengths using a prescribed vertical profile

      ! Argument declaration
      real, dimension(:,:), intent(in) :: zn1           !mixing length used near surface above mid-troposphere (m)
      real, dimension(:,:), intent(in) :: zn2           !mixing length used in lower to mid-troposphere (m)
      real, dimension(:,:), intent(in) :: zz            !height of e-levs (m)
      real, dimension(:,:), intent(in) :: sigma         !coordinate values for e-levs
      real, dimension(:), intent(in) :: ps              !surface pressure (Pa)
      real, dimension(:,:), intent(out) :: znblend      !blended mixing length (m)

      ! Local parameters
      real, parameter :: ZMIX=500.,PLOW=550E2,PHIGH=450E2

      ! Local variables and common blocks
      integer :: j,k,nj,nk
      real :: pres

      ! Set return status
      stat = PHY_ERROR
      
      ! Check initialization status
      if (.not.initialized) then
         call msg_toall(MSG_ERROR,'(ml_blend) mixing length package not initialized')
         return
      endif
    
      ! Initialization
      nj = size(zn1,dim=1)
      nk = size(zn1,dim=2)

      ! Blend mixing length estimates
      znblend = zn2
      do j=1,nj
         do k=1,nk
            ! Near-surface blending
            if ( zz(j,k) < ZMIX ) znblend(j,k) = zn1(j,k)+zz(j,k)/ZMIX*(znblend(j,k)-zn1(j,k))
            ! Mid-tropospheric blending
            pres   = sigma(j,k) * ps(j)
            znblend(j,k) = zn1(j,k) + ( min( max( pres, PHIGH ) , PLOW ) - PHIGH ) &
                 * ( znblend(j,k) - zn1(j,k) ) / ( PLOW - PHIGH )
            znblend(j,k) = max(znblend(j,k),1.e-6)
         enddo
      enddo

      ! Successful end of subprogram
      stat = PHY_OK
      return
   end function ml_blend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ml_calc_blac(zn,rig,w_cld,zz,z0,fm,hpbl,lh,dxdy,przn) result(stat)
      ! Compute the Blackadar (1962) mixing length scale

      ! Argument declaration
      real, dimension(:,:), intent(in) :: rig             !gradient Richardson number
      real, dimension(:,:,:), intent(in) :: w_cld         !Cu cloud velocity scale profile (1 L_heat; 2 L_momentum; 3 L_dissipation)
      real, dimension(:,:), intent(in) :: zz              !height of e-levs (m)
      real, dimension(:,:), intent(in) :: fm              !momentum stability function
      real, dimension(:), intent(in) :: z0                !roughness length for momentum (m)
      real, dimension(:), intent(in) :: hpbl              !height of the PBL (m)
      real, dimension(:), intent(in) :: lh                !launching height (m)
      real, dimension(:), intent(in), optional :: dxdy    !horizontal grid cell area (m^2)
      real, dimension(:,:), intent(out) :: zn             !mixing length (m)
      real, dimension(:,:), intent(out), optional :: przn !mixing length-based Prandtl number adjustment

      ! Local variables and common blocks
      integer :: k,nk
      real, dimension(size(rig,dim=1)) :: lmda_t
      real, dimension(size(rig,dim=1),size(rig,dim=2)) :: lmda,my_przn

      ! Set return status
      stat = PHY_ERROR

      ! Check initialization status
      if (.not.initialized) then
         call msg_toall(MSG_ERROR,'(ml_calc_blac) mixing length package not initialized')
         return
      endif
      
      ! Initialization
      nk = size(zz,dim=2)

      ! Compute asymptotic mixing length
      select case (mlblac_max)
      case ('LOCK07')
         do k=1,nk
            where (rig(:,k) >= 0.)
               lmda_t(:) = min(max(ML_LMDAS,0.3*hpbl(:)),ML_LMDA)
               lmda(:,k) = min(max(lmda_t(:),sqrt(2.)*lh(:)),ML_LMDA)
               my_przn(:,k) = ( 1. + KARMAN*(zz(:,k)+z0(:))/lmda_t(:) ) &
                    / ( 1. + KARMAN*(zz(:,k)+z0(:))/lmda(:,k) )
            elsewhere
               lmda(:,k) = ML_LMDA
               my_przn(:,k) = 1.
            endwhere
         enddo
      case DEFAULT
         lmda(:,:) = ML_LMDA
         my_przn(:,:) = 1.
      end select

      ! Adjust mixing length limit for high resolution (<850m) grids (Mason and Brown 1999)
      if (present(dxdy)) then
         do k=1,nk
            lmda(:,k) = min(lmda(:,k),0.23*sqrt(dxdy(:)))
         enddo
      endif

      ! Compute Blackadar mixing length
      select case (mlblac_max)
      case ('BLAC62')
         do k=1,nk
            zn(:,k) = min(KARMAN*(zz(:,k)+z0(:)),lmda(:,k))
         enddo
      case DEFAULT
         do k=1,nk
            zn(:,k) = KARMAN*(zz(:,k)+z0(:))/(1.+KARMAN*(zz(:,k)+z0(:))/lmda(:,k))
         enddo
      end select
      zn = zn / fm
      where (rig < 0.) zn = min(zn,ML_LMDA_UNSTAB)

      ! Adjust mixing length for non-local cloud effects (Lock and Mailhot 2006)
      zn = (zn**3 + w_cld(:,:,2)**3)**(1./3.)

      ! Fill optional return requests
      if (present(przn)) przn(:,:) = my_przn(:,:)

      ! Successful end of subprogram
      stat = PHY_OK
      return
   end function ml_calc_blac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ml_calc_boujo(zn,th,en,w_cld,zz,sigma,ps,mask,init) result(stat)
      ! Compute the Bougeault and Lacarrere (1989) length scales
      use integrals, only: int_profile, int_solve, INT_ERR

      ! Argument declaration
      real, dimension(:,:), intent(in) :: th                !virtual potential temperature on e-levs (K)
      real, dimension(:,:), intent(in) :: en                !turbulent kinetic energy (m2/s2)
      real, dimension(:,:,:), intent(in) :: w_cld           !Cu cloud velocity scale profile (1 L_heat; 2 L_momentum; 3 L_dissipation)
      real, dimension(:,:), intent(in) :: zz                !height of e-levs (m)
      real, dimension(:,:), intent(in) :: sigma             !coordinate values for e-levs
      real, dimension(:), intent(in) :: ps                  !surface pressure (Pa)
      logical, dimension(:,:), intent(in), optional :: mask !calculation mask [.true.]
      logical, intent(in), optional :: init                 !initialize the mixing length [.true.]
      real, dimension(:,:), intent(inout) :: zn             !mixing length (m)

      !Author
      !          S. Belair (November 1996)
      !
      !Revision
      ! 001      J. Mailhot (July 1999) - version using the virtual potential
      !                                   temperature; change name to MIXLEN1
      ! 002      J. Mailhot (Sept 1999) - clipping of RIF maximum for computation of ZE
      ! 003      S. Belair (Oct 1999)   - staggerred values for the virtual
      !                                   potential temperature and the heights
      ! 004      J. Mailhot (July 2000) - correct limits for solution of quadratic eqn.
      ! 005      J. Mailhot (Aug 2000) - add relaxation option (RELAX = .T. or .F.)
      ! 006      S. Belair, J. Mailhot (March 2001)
      !                                  blend between local (i.e.,
      !                                  Bougeault-Lacarrere) and
      !                                  background (i.e., input) mixing and
      !                                  dissipation lengths
      ! 007      A-M Leduc  (Oct 2001)   - Automatic arrays
      ! 008      J. Mailhot (May 2002) - restrict local mixing to convective case
      ! 009      S. Belair, J. Mailhot (June 2002) - use fixed heights for blend
      !                                              and remove stability considerations
      ! 010      S. Belair (Jan 2003)   -reorganization and modification of bougeault
      !                                   mixlen1--->mixlen2
      ! 011      B. Bilodeau (Aug 2003) - IBM conversion (scalar version)
      ! 012      S. Belair (March 2004) - Relax the mixing length towards the
      !                                   Blackadar value in the upper troposphere
      !                                   (i.e., between plow and phigh)
      ! 002      L. Spacek (Dec 2007)   - all calculations on energy levels
      !
      !          C. Girard              - memory optimization
      !
      !Object
      !           Calculates the mixing length ZN and the dissipation
      !           length ZE based on the Bougeault and Lacarrere method.
      !

      ! Local variables and common blocks
      integer :: n,nk,j,ki,istat,nc
      integer, dimension(size(th,dim=1)) :: slk,indx
      real :: gravinv
      real, dimension(size(th,dim=1)) :: a,zdep
      real, dimension(size(th,dim=1),size(th,dim=2)) :: y,lup,ldown,zcoord
      logical :: myInit
      logical, dimension(size(th,dim=1)) :: intok
      logical, dimension(size(th,dim=1),size(th,dim=2)) :: myMask

      ! Set return status
      stat = PHY_ERROR

      ! Check initialization status
      if (.not.initialized) then
         call msg_toall(MSG_ERROR,'(ml_calc_boujo) mixing length package not initialized')
         return
      endif
      
      ! Initialization
      n = size(th,dim=1)
      nk = size(th,dim=2)
      gravinv = 1./GRAV
      istat = neark(sigma,ps,1000.,n,nk,slk) !determine "surface layer" vertical index for buoyancy coefficient
      myMask = .true.
      if (present(mask)) myMask = mask
      myInit = .true.
      if (present(init)) myInit = init

      ! Expend TKE through vertical displacements
      if (myInit) zn = -1.
      do ki=2,nk-1

         ! Prepare integral equation components
         nc = 0
         indx = -1
         do j=1,n
            if (.not. myMask(j,ki)) cycle
            nc = nc+1
            indx(j) = nc
            a(nc) = min(en(j,ki),4.)*th(j,slk(j))*gravinv
            y(nc,:) = th(j,:) - th(j,ki)
            zcoord(nc,:) = zz(j,:)
            zdep(nc) = zz(j,ki)
         end do

         ! Solve integral equation in both vertical directions
         if (int_solve(lup(1:nc,ki),y(1:nc,:),zcoord(1:nc,:),zdep(1:nc),a(1:nc),'up') == INT_ERR) then
            call msg(MSG_ERROR,'(ml_calc_boujo) error returned by integral solver for upwards displacement')
            return
         endif
         if (any(lup(1:nc,ki) > ML_MIN)) then
         if (int_solve(ldown(1:nc,ki),-y(1:nc,:),zcoord(1:nc,:),zdep(1:nc),a(1:nc),'down',found=intok(1:nc)) == INT_ERR) then
               call msg(MSG_ERROR,'(ml_calc_boujo) error returned by integral solver for downwards displacement')
               return
            endif
            where (.not.intok(1:nc))
               ldown(1:nc,ki) = zdep(1:nc)
            elsewhere
               ldown(1:nc,ki) = min(ldown(1:nc,ki),zdep(1:nc))
            endwhere
         else        
            ldown(1:nc,ki) = ML_MIN
         endif
         
         ! Select final mixing length
         do j=1,n
            nc = indx(j)
            if (nc < 0) cycle
            ! Compute "average" mixing length
            zn(j,ki) = min(lup(nc,ki),ldown(nc,ki))
            zn(j,ki) = min(zn(j,ki),zz(j,ki))
            zn(j,ki) = max(zn(j,ki),ML_MIN)
            ! Adjust mixing length for non-local cloud effects (Lock and Mailhot 2006)
            zn(j,ki) = (zn(j,ki)**3 + w_cld(j,ki,2)**3)**(1./3.)
         enddo

      enddo

      ! Copy values into extremes
      zn(:,1) = zn(:,2)
      zn(:,nk) = zn(:,nk-1)

      ! Successful end of subprogram
      stat = PHY_OK
      return
   end function ml_calc_boujo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function ml_calc_lh(l_e,l_m,pri,dthv,en,w_cld,ri,zs,f_cs,hpar) result(stat)
      ! Compute the Lenderink and Hostlag (2004) length scales.

      ! Arguments
      real, dimension(:,:), intent(in) :: dthv                  !buoyancy flux (m2/s2)
      real, dimension(:,:), intent(in) :: en                    !turbulent kinetic energy (m2/s2)
      real, dimension(:,:,:), intent(in) :: w_cld               !Cu cloud velocity scale profile (1 L_heat; 2 L_momentum; 3 L_dissipation)
      real, dimension(:,:), intent(in) :: ri                    !gradient Richardson number
      real, dimension(:,:), intent(in) :: zs                    !height of energy levels (m)
      real, dimension(:,:), intent(in) :: f_cs                  !interpolation factor for length scale coefficients
      real, dimension(:), intent(in) :: hpar                    !height of parcel ascent (unstable) or pbl depth (stable)
      real, dimension(:,:), intent(out) :: l_e                  !dissipation length scale (m)
      real, dimension(:,:), intent(out) :: l_m                  !momentum mixing length scale (m)
      real, dimension(:,:), intent(out) :: pri                  !inverse Prandtl number (ratio of l_h/l_m)

      !Author
      !          A. Lock (May 2004)

      !Revision

      !Object
      !           Calculates the mixing lengths
      !           based on the method of Lenderink and Holtslag (2004)
      !           includes non-local velocity scale contribution in cumulus

      !Notes
      !           Based on Mailhot and Lock (2004; LM04)

      ! Local parameters
      real, parameter :: B_LH=4.0     ! L&H's constants b
      real, parameter :: C_N=0.516    ! L&H's constants c_n (same as RPN "AA")

      ! Local declarations
      integer :: j,k,n,nk
      integer :: imh       ! counter to allow same code to be used for momentum and heat
      real :: alpha_n,alpha_c,alpha_r ! L&H's constants
      real :: alpha_c_fac             ! ratio of alpha_c to apha_n
      real :: c_stab_d                ! "c_m" for dissipn in stable conditions
      real :: c_stab_h,c_stab_m       ! c_h, c_m for stable conditions
      real :: l_min                   ! minimum for "integral" length scale
      real :: recip_l_stable          ! 1/l_stable
      real :: wf                      ! 
      real, dimension(size(en,dim=1)) :: &
           ri_int         ! vertically integrated function of Ri
      real, dimension(size(en,dim=1)) :: &
           recip_inf      ! 1/l_inf, where l_inf is length scale away from surface
      real, dimension(size(en,dim=1),size(en,dim=2)) :: &
           lup          ,&! mixing length from upward integration
           ldown        ,&! mixing length from downward integration
           f_ri         ,&! function of Ri
           l_h            ! mixing length for heat

      ! Set return status
      stat = PHY_ERROR

      ! Check initialization status
      if (.not.initialized) then
         call msg_toall(MSG_ERROR,'(ml_key) mixing length package not initialized')
         return
      endif
   
      ! Initializations and constants
      n = size(en,dim=1)
      nk = size(en,dim=2)
      alpha_n = C_N * KARMAN

      ! Length scale estimates for momentum (imh=1) and heat (imh=2) sequentially
      MH_LOOP: do imh=1,2

         ! Coefficients depend on momentum or heat context
         if (imh == 1) then
            alpha_c_fac = 3.0
         else
            alpha_c_fac = 5.0
         endif
         alpha_c = alpha_c_fac * alpha_n
         alpha_r = PI*alpha_n*B_LH/(alpha_c-alpha_n)  ! Less conservative estimate for alpha_r (Eq. 13 in L&H)

         ! Compute Richarson number function (simplified version of L&H formula)
         do k=1,nk
            do j=1,n
               if (ri(j,k) > 0.0) then
                  f_ri(j,k) =  alpha_n - 2.0*alpha_n*b_lh*ri(j,k)
               else
                  f_ri(j,k) =  alpha_n - (2.0/pi)*(alpha_c-alpha_n) &
                       *atan(alpha_r*ri(j,k))
               endif
            enddo
         enddo

         ! Integrate upwards for lup using trapezoidal integration
         ri_int = 0.
         do k=nk-1,2,-1
            do j=1,n
               ri_int(j) = ri_int(j) + 0.5 * (zs(j,k)-zs(j,k+1)) &
                    * (f_ri(j,k)+f_ri(j,k+1))
               ri_int(j) = max(   0.0, ri_int(j) )
               lup(j,k)  = max( 1.e-6, ri_int(j) )
            enddo
         enddo

         ! Integrate downwards for ldown using trapezoidal integration
         ri_int = 0.
         do k=2,nk-1
            do j=1,n
               ri_int(j) = ri_int(j) + 0.5 * (zs(j,k-1)-zs(j,k)) &
                    * (f_ri(j,k-1)+f_ri(j,k))
               ri_int(j) = max(   0.0, ri_int(j) )
               ! Near-surface increase to reduce its influence in stable case (L&H Section 6)
               ldown(j,k)= max( 75.*exp(-zs(j,k)/500.) , ri_int(j) )
            enddo
         enddo

         ! Combine length estimates into a single integral length scale
         select case (imh)
         case (1) ! momentum length scale (l_m)
            do k=2,nk-1
               do j=1,n
                  if (lup(j,k) > 0.1 .and. ldown(j,k) > 0.1) then
                     l_m(j,k) = 1.0/( 1.0/lup(j,k) + 1.0/ldown(j,k) )
                  else
                     l_m(j,k) = 0.
                  endif
               enddo
            enddo
            do j=1,n
               l_m(j,1)  = 0.0
            enddo
         case (2) ! heat length scale (l_h)
            do k=2,nk-1
               do j=1,n
                  if (lup(j,k) > 0.00001 .and. ldown(j,k) > 0.00001) then
                     l_h(j,k) = 1.0/( 1.0/lup(j,k) + 1.0/ldown(j,k) )
                  else
                     l_h(j,k) = 0.
                  endif
               enddo
            enddo
            do j=1,n
               l_h(j,1)  = 0.0
            enddo
         end select

      enddo MH_LOOP

      ! Asymptotic mixing length (ensure larger than L&H 75m in cumulus layers with 0.1*hpar)
      recip_inf = 1./max(75.,min(200.,0.1*hpar))

      ! Blend integral length scales with stable and minimum scales (Appendix B of L&H 2004)
      do k=1,nk-1
         do j=1,n

            ! Average inverse square of integral and minimum length scales
            l_min = 1.0/(recip_inf(j)+(1.0/(0.5*C_N*KARMAN*zs(j,k))))
            l_h(j,k) = l_h(j,k)**2. + l_min**2.
            l_m(j,k) = l_m(j,k)**2. + l_min**2.
            l_e(j,k) = l_m(j,k)

            ! Add stable length scale estimate adjusted by cloud properties
            STABLE: if ( dthv(j,k) > 0. .and. en(j,k) > 0. ) then
               !         ! lh04 standard coeffs:
               !          c_stab_h = 0.2
               !          c_stab_m = 0.2 * min( 3., 1.+2.*ri(j,k) )
               !          c_stab_d = c_stab_m
               !         ! comparison with gabls les suggests following coeffs (lm04 eq. 14)
               c_stab_h  = 0.072                             ! 0.36*cs_lh
               c_stab_m  = 0.05 * min( 3., 1.+2.*ri(j,k) )   ! 0.25*cs_lh
               c_stab_d  = 0.1  * min( 3., 1.+5.*ri(j,k) )   !~0.5 *cs_lh
               !         ! -------------------------------------------------------------
               !         ! f_cs is a factor for the coefficient c_s that interpolates
               !         ! between the larger original values (f_cs=1) and the smaller
               !         ! values proposed by mailhot and lock (2004) from les of gabls.
               !         ! set in nlcalc.
               !         ! -------------------------------------------------------------               
               wf = max(0.,min(1.,f_cs(j,k)))
               c_stab_h  = (1.-wf) * c_stab_h + wf * 0.2 
               c_stab_m  = (1.-wf) * c_stab_m + wf * 0.2 * min( 3., 1.+2.*ri(j,k) )
               c_stab_d  = (1.-wf) * c_stab_d + wf * 0.2 * min( 3., 1.+2.*ri(j,k) )
               
               recip_l_stable = sqrt( dthv(j,k) ) /  &
                    ( (c_stab_h*sqrt(en(j,k)))**3. + w_cld(j,k,1)**3. )**(1./3.)
               l_h(j,k) = 1.0/sqrt( (1.0/l_h(j,k)) + recip_l_stable**2. )
               recip_l_stable = sqrt( dthv(j,k) ) /  &
                    ( (c_stab_m*sqrt(en(j,k)))**3. + w_cld(j,k,2)**3. )**(1./3.)
               l_m(j,k) = 1.0/sqrt( (1.0/l_m(j,k)) + recip_l_stable**2. )
               recip_l_stable = sqrt( dthv(j,k) ) /  &
                    ( (c_stab_d*sqrt(en(j,k)))**3. + w_cld(j,k,3)**3. )**(1./3.)
               l_e(j,k) = 1.0/sqrt( (1.0/l_e(j,k)) + recip_l_stable**2. )
            else
               l_h(j,k) = sqrt(l_h(j,k))
               l_m(j,k) = sqrt(l_m(j,k))
               l_e(j,k) = l_m(j,k)
            endif STABLE
         enddo
      enddo

      ! Conversions for output
      pri(:,1:nk-1) = l_h(:,1:nk-1) / l_m(:,1:nk-1)
      l_m(:,1:nk-1) = l_m(:,1:nk-1) / 0.516         ! allow for AA factor
      l_e(:,1:nk-1) = l_e(:,1:nk-1) * 0.14/0.071    ! allow for different constants in LH04 to RPN
      pri(:,nk) = 0.
      l_m(:,nk) = 0.
      l_e(:,nk) = 0.

      ! Successful completion of subprogram
      stat = PHY_OK
      
   end function ml_calc_lh

end module mixing_length

