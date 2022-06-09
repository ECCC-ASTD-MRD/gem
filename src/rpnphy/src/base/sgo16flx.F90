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
!-------------------------------------- LICENCE END --------------------------

!/@*
subroutine sgo16flx3(uu, vv, utend, vtend, ttend, utendgwd4, vtendgwd4, &
     tth, ttf, ss, ssh,  &
     ilev, ilg, il1, il2,  &
     tau, taufac,  &
     gc, height, slope, xcent, mtdir,  &
     cdmin, fc, phic, stabfactor, nldirfactor, mrk2, &
     gwdrag, blocking, split_tend_L)
   use integrals, only: int_profile, int_solve, INT_OK
   use phy_options, only: sgo_windfac
   use tdpack, only: CAPPA, CPD, GRAV, RGASD
   use ens_perturb, only: ens_nc2d, ens_spp_get
   implicit none
!!!#include <arch_specific.hf>
   !@Arguments
   logical, intent(in) :: gwdrag                   ! .true. for gravity wave drag
   logical, intent(in) :: blocking                 ! .true. for low-level blocking
   logical, intent(in) :: stabfactor               ! .true. for stability enhancement factor
   logical, intent(in) :: nldirfactor              ! .true. for directional enhacement factor
   logical, intent(in) :: split_tend_L             ! .true. for keeping GWD uvtend seprated from total
   integer, intent(in) :: ilev                     ! number of levels
   integer, intent(in) :: ilg                      ! total horizontal dimension
   integer, intent(in) :: il1                      ! horiz. index to start calculation
   integer, intent(in) :: il2                      ! horiz. index to stop calculation 
   real, intent(in) :: tau                         ! timestep
   real, intent(in) :: taufac                      ! 1/(length scale) for GWD scheme 
   real, intent(in) :: cdmin                       ! minimum bulk drag coeff for blocking
   real, intent(in) :: fc                          ! tuning factor for blocking height 
   real, intent(in) :: phic                        ! critical phase value for blocking height
   real, dimension(ilg), intent(in) :: gc          ! land-cover fraction
   real, dimension(ilg), intent(in) :: height      ! mountain (launching) height
   real, dimension(ilg), intent(in) :: slope       ! mountain slope
   real, dimension(ilg), intent(in) :: xcent       ! mountain eccentricity
   real, dimension(ilg), intent(in) :: mtdir       ! mountain direction
   real, dimension(ilg,ilev), intent(in) :: uu     ! u-component of wind
   real, dimension(ilg,ilev), intent(in) :: vv     ! v-component of wind
   real, dimension(ilg,ilev), intent(in) :: tth    ! virtual temperature at thermo levels
   real, dimension(ilg,ilev), intent(in) :: ttf    ! virtual temperature at mom levels 
   real, dimension(ilg,ilev), intent(in) :: ss     ! sigma at thermo levels
   real, dimension(ilg,ilev), intent(in) :: ssh    ! sigma at mom levels
   real, dimension(ilg,ens_nc2d), intent(in) :: mrk2 !Markov chains for perturbations
   real, dimension(ilg,ilev), intent(out) :: utend ! total tendency for u
   real, dimension(ilg,ilev), intent(out) :: vtend ! total tendency for v
   real, dimension(ilg,ilev), intent(out) :: ttend ! total tendency for t
   real, dimension(ilg,ilev), intent(out) :: utendgwd4 ! GWD tendency for u
   real, dimension(ilg,ilev), intent(out) :: vtendgwd4 ! GWD tendency for v
    
   !@Author A. Zadra - May 2002 - From lin_sgoflx1
   !@Revisions
   ! 001    B. Bilodeau and A. Zadra (Apr 2003) - Add smoothing of BV frequency.
   ! 002    J.-P. Toviessi (May 2003) - IBM conversion
   !              - calls to vsqrt routine (from massvp4 library)
   !              - calls to vdiv  routine (from massvp4 library)
   !              - calls to vrec  routine (from massvp4 library)
   !              - calls to vlog  routine (from massvp4 library)
   ! 003    A. Zadra and B. Bilodeau (Oct 2003) - Modifications to the 
   !                new code in order to validate with oper. code
   ! 004    A. Zadra (Apr 2004) - Implicit formulation for utendllb and 
   !                vtendllb
   ! 005    G. Pellerin (March 2006) - calculation of piotwo
   ! 006    A. Zadra (Dec 2008) - Changes in the calculation of directional
   !                factor for blocking term
   ! 007    A. Zadra (Jan 2011) - In the GWD scheme, use reference level at a 
   !                a fixed height (href = 130m)
   ! 008    A. Zadra (Apr 2011) - For blocking term:
   !                - new enhancement factor due to stability
   !                - new nonlinear amplification factor
   ! 009    A. Zadra & R. McTaggart-Cowan (June 2011) - top level index of tendency
   !                calculation is now provided through the interface (TOPLEV)
   ! 010    A. Zadra - add logical switches (stabfactor,nldirfactor) for optional
   !                amplification factors, and minimum value (cdmin) of 
   !                effective drag coefficient for blocking term
   ! 011    A. Plante - Removed toplev since notop is not a option anymore
   !
   ! 012    A. Zadra - Critical value phic (used in calculation of blocking height)
   !                is adjusted according to the resolution of the lowest layer
   ! 013    A. Zadra (Nov 2016) - Various corrections and simplifications:
   !                - from double to single precision
   !                - use of integral module
   !                - simplified calculations in blocking scheme
   !                - correction to eccentricity factor
   !                - 2 new input tuning parameters: fc and phic
   !                - removed orographic lift and lee-wave breaking
   !                - removed apply tendency section
   ! 014    A. Zadra (Apr 2017) - Add tendency for temperature
   !@Object
   !        Compute total wind tendencies from 2 processes:
   !          - gravity-wave drag (McFarlane 87)
   !          - low-level blocking (Lott & Miller 97)
   !        with 2 optional enhancements for drag coefficient in blocking scheme
   !          - enhancement factor due to stability (Vosper et al. 2009)
   !          - nonlinear directional amplification (Wells et al. 2008)
   
   !  Local parameter definitions
   real, parameter :: vmin1  = 2.     ! minimum wind speed to produce GWD tendencies
   real, parameter :: vmin2  = 0.01   ! minimum wind speed to produce blocking tendencies
   real, parameter :: v0    = 1.e-6   ! epsilon for wind speed
   real, parameter :: hmin1  = 10.    ! minimum height to launch gravity waves
   real, parameter :: hmin2  = 3.     ! minimum height to allow low-level blocking
   real, parameter :: cdblk = 1.      ! bulk drag coefficient of blocking (approx. 1)
   real, parameter :: href  = 130.    ! reference height for GWD
   real, parameter :: nmin  = 1.e-8   ! epsilon for buoyancy freqency
   
   !  Local variable declarations
   integer :: i, l, ii, len, lref, lrefm, jyes, jno, istat
   integer, dimension(ilg) :: drag, izt0

   real :: tmp, aux, eta, piotwo, psi, cpsi, spsi, ratio, fvert, fwind, amp, gmsi, &
        B, C, ks2, vent
   real, dimension(ilg) :: ub, vb, vmodb, env, slp2, slpf, gamma, theta, &
        hblk, uav, vav, velav, nav, fdir, dscale, fdscale, cdf, &
        ampd, uavd, vavd, velavd, navd, vphic, phicni
   real, dimension(ilg, ilev) :: u, v, tf, th, s, sh, sexpk, shexpk, bvf, &
        bvfreq, veln, asq, asqi, asqs, dfac, depfac, denfac, grad, zb, &
        utendgwd, vtendgwd, utendllb, vtendllb, utendtot, vtendtot

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     INITIAL STEPS
   !--------------------------------------------------------------------
   !
   utend(il1:il2,:) = 0.
   vtend(il1:il2,:) = 0.
   ttend(il1:il2,:) = 0.
   
   !-------------------------------------------------------------------
   !     Find and gather active grid columns

   len   = il2 - il1 + 1
   jyes = 0
   jno  = len + 1

   do i=il1,il2
      if ( gc(i) == -1. ) then
         jyes       = jyes + 1
         drag(jyes) = i
      else
         jno         = jno - 1
         drag(jno)  = i
      endif
   enddo

   !     Check if there is AT LEAST ONE active column

   if (jyes <= 0) return

   !     Retrieve stochastic perturbations on request
   phicni = ens_spp_get('sgo_phic', mrk2, default=phic)
   
   do i=1,len
      ii = drag(i) + il1 - 1
      env(i)   = height(ii) 
      slp2(i)  = slope(ii)
      gamma(i) = xcent(ii)
      theta(i) = mtdir(ii)
      vphic(i) = phicni(ii)
   enddo

   do i=1,len
      slpf(i) = slp2(i)/max(env(i),hmin2) 
   enddo

   do l=1,ilev
      do i=1,len
         ii      = drag(i) + il1 - 1
         u(i,l)  = uu(ii,l)
         v(i,l)  = vv(ii,l)
         tf(i,l) = ttf(ii,l)
         th(i,l) = tth(ii,l)
         s(i,l)  = ss(ii,l)
         sh(i,l) = ssh(ii,l)
         sexpk(i,l)  = s(i,l)**CAPPA
         shexpk(i,l) = sh(i,l)**CAPPA
      enddo
   enddo

   !     Compute auxiliary profiles

   IF_GWD_BLO: if (gwdrag .or. blocking) then
      !--------------------------------------------------------------------
      !     Build elevation field:

      do i=1,len
         zb(i,ilev)  = -(RGASD/GRAV)*th(i,ilev)*alog(s(i,ilev))
      enddo
      do l=ilev-1,1,-1
         do i=1,len
            zb(i,l)  = zb(i,l+1) + &
                 (RGASD/GRAV)*th(i,l)*alog(s(i,l+1)/s(i,l))          
         enddo
      enddo

      !-------------------------------------------------------------------
      !     Calculate BF frequency:

      bvf = 0.
      do l=2,ilev
         do i=1,len
            aux = (th(i,l)/shexpk(i,l)-th(i,l-1)/shexpk(i,l-1)) &
                 /(sh(i,l)-sh(i,l-1))
            if ( aux >= (-5./s(i,l)) ) aux = -5./s(i,l)
            aux = -aux*s(i,l)*GRAV/(RGASD*tf(i,l))
            bvf(i,l) = GRAV*aux*sexpk(i,l)/tf(i,l)
            bvf(i,l) = sqrt(bvf(i,l))
         enddo
      enddo
      do i=1,len
         bvf(i,1) = bvf(i,2)
      enddo

   endif IF_GWD_BLO

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     GRAVITY-WAVE DRAG
   !--------------------------------------------------------------------
   !     Initialize arrays

   utendgwd(1:len,:) = 0.
   vtendgwd(1:len,:) = 0.
   !     
   izt0(1:len)     = ilev-1
   vmodb(1:len)    = v0
   ub(1:len)       = v0
   vb(1:len)       = v0
   veln(1:len,:)     = v0
   asq(1:len,:)    = hmin1*hmin1
   asqs(1:len,:)   = hmin1*hmin1
   asqi(1:len,:)   = hmin1*hmin1
   dfac(1:len,:)   = 0.
   depfac(1:len,:) = 0.
   denfac(1:len,:) = 0.
   grad(1:len,:)   = 0.     

   IF_GWD: if (gwdrag) then
      !--------------------------------------------------------------------
      !     Smooth B V frequency
      bvfreq = bvf
      do l=2,ilev
         do i=1,jyes
            ratio = 5.*alog(s(i,l)/s(i,l-1))            
            bvfreq(i,l) = (bvfreq(i,l-1) + ratio*bvfreq(i,l)) &
                 /(1.+ratio)
         enddo
      enddo

      !     Find reference level:
      do l=ilev-1,1,-1
         do i=1,jyes
            if (zb(i,l) <= href)    izt0(i)  = l
         enddo
      enddo
      !--------------------------------------------------------------------
      !     Wind and unit vector at level LREFM:

      do i=1,jyes
         lrefm = izt0(i)
         vmodb(i) = u(i,lrefm)*u(i,lrefm) + v(i,lrefm)*v(i,lrefm)
         vmodb(i) = sqrt(vmodb(i))
      enddo
      
      do i=1,jyes
         lrefm = izt0(i)
         if (vmodb(i) <= 1.)  vmodb(i) = 1.
         ub(i) = u(i,lrefm)/vmodb(i)
         vb(i) = v(i,lrefm)/vmodb(i)
      enddo

      !--------------------------------------------------------------------
      !     Project wind field on reference wind:

      do l=1,ilev
         do i=1,jyes 
            veln(i,l) = u(i,l)*ub(i)+v(i,l)*vb(i)
            if (veln(i,l) <= v0)  veln(i,l) = v0
         enddo
      enddo

      !--------------------------------------------------------------------
      !     Stress field
      !     Compute stress at reference level:

      do i=1,jyes
         lrefm = izt0(i)
         lref  = lrefm + 1
         asq(i,lref)  = env(i)*env(i)
         asqs(i,lref) = asq(i,lref)
         asqi(i,lref) = asq(i,lref)
         dfac(i,lref) = bvfreq(i,lrefm)*sh(i,lrefm)*vmodb(i) &
              /tf(i,lrefm)
      enddo

      !     Compute stress at other levels (bottom-up):

      depfac = 0.
      do i=1,jyes
         lrefm = izt0(i)
         do l=lrefm,1,-1
            dfac(i,l) = bvfreq(i,l)*s(i,l)*veln(i,l)/tf(i,l)
            asqi(i,l) = asq(i,l+1)*dfac(i,l+1)/dfac(i,l)
            if (veln(i,l) >= 1.) then
               tmp = veln(i,l)/bvfreq(i,l)
               asqs(i,l) = 0.5*tmp*tmp
            else
               asqs(i,l) = 0.
            endif
            if (asqi(i,l) <= asqs(i,l)) then
               asq(i,l)    = asqi(i,l)
            else
               asq(i,l)    = asqs(i,l)
            endif
            depfac(i,l) = (taufac*GRAV/RGASD)*asq(i,l)*dfac(i,l)
         enddo
      enddo
      do i=1,jyes
         lrefm = izt0(i)
         lref  = lrefm + 1
         do l=ilev,lref,-1
            depfac(i,l) = depfac(i,lrefm)
         enddo
      enddo

      !--------------------------------------------------------------------
      !     Compute gwd tendencies:

      do i=1,jyes
         if ((depfac(i,2) - depfac(i,1)).gt.1.e-10) then
            eta = 1.
         else
            eta = 0.
         endif

         tmp = 2.*sh(i,1) + eta*3.*tau*depfac(i,1)/veln(i,1) 
         tmp = 1.0/tmp
         grad(i,1) = 2.*eta*depfac(i,1) &
              *tmp
         utendgwd(i,1) = -ub(i)*grad(i,1)
         vtendgwd(i,1) = -vb(i)*grad(i,1)
         denfac(i,1) = grad(i,1)*3.*tau*depfac(i,1)/veln(i,1)

         lrefm = izt0(i)
         lref  = lrefm + 1
         do l=ilev,lref,-1
            utendgwd(i,l) = 0.
            vtendgwd(i,l) = 0.
         enddo
      enddo

      do i=1,jyes
         lrefm = izt0(i)
         do l=2,lrefm
            if ((depfac(i,l) - depfac(i,l-1)).gt.1.e-10) then
               eta = 1.
            else
               eta = 0.
            endif
            tmp =  2.*(sh(i,l)-sh(i,l-1)) + &
                 eta*3.*tau*depfac(i,l)/veln(i,l)
            tmp = 1.0/tmp
            grad(i,l) = ( 2.*depfac(i,l)-2.*depfac(i,l-1) + &
                 eta*denfac(i,l-1) )*tmp
            utendgwd(i,l) = -ub(i)*grad(i,l)
            vtendgwd(i,l) = -vb(i)*grad(i,l)
            denfac(i,l) = grad(i,l)*3.*tau*depfac(i,l)/veln(i,l)  
         enddo
      enddo

      do l=1,ilev
         do i=1,jyes
            if ( vmodb(i).lt.vmin1 .or. env(i).lt.hmin1 ) then
               utendgwd(i,l) = 0.
               vtendgwd(i,l) = 0.
            endif
         enddo
      enddo

   endif IF_GWD


   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     LOW-LEVEL BLOCKING
   !--------------------------------------------------------------------
   !     Initialize arrays

   utendllb(1:len,:) = 0.
   vtendllb(1:len,:) = 0.

   fdir(1:len)  = 1.
   ampd(1:len)  = 1.
   cdf(1:len)   = cdmin
   hblk(1:len)  = 0.              
   uav(1:len)   = v0
   vav(1:len)   = v0
   velav(1:len) = v0
   nav(1:len)   = nmin        
   dscale(1:len) = hmin2        
   uavd(1:len)   = v0
   vavd(1:len)   = v0
   velavd(1:len) = v0
   navd(1:len) = nmin

   IF_BLO: if (blocking) then
      !--------------------------------------------------------------------
      !     Take the original (non-smoothed) B-F frequency
      bvfreq = bvf

      !     Compute average wind, wind speed and buoyancy

      uav(1:jyes) = v0
      istat = int_profile(uav(1:jyes),u(1:jyes,:),zb(1:jyes,:), &
           0.5*env(1:jyes),env(1:jyes))
      uav(1:jyes) = 2.*uav(1:jyes)/max(env(1:jyes),hmin2)

      vav(1:jyes) = v0
      istat = int_profile(vav(1:jyes),v(1:jyes,:),zb(1:jyes,:), &
           0.5*env(1:jyes),env(1:jyes))
      vav(1:jyes) = 2.*vav(1:jyes)/max(env(1:jyes),hmin2)

      velav(1:jyes) = sqrt(uav(1:jyes)*uav(1:jyes) &
           + vav(1:jyes)*vav(1:jyes))
      velav(1:jyes) = max(velav(1:jyes),v0)

      nav(1:jyes) = nmin
      istat = int_profile(nav(1:jyes),bvfreq(1:jyes,:),zb(1:jyes,:), &
           0.*env(1:jyes),env(1:jyes))
      nav(1:jyes) = nav(1:jyes)/max(env(1:jyes),hmin2)
      nav(1:jyes) = max(nav(1:jyes),nmin)
      
      !     Compute blocking height
      hblk(1:jyes) = max(fc*env(1:jyes) & 
           - vphic(1:jyes)*velav(1:jyes)/nav(1:jyes), 0.)

      !--------------------------------------------------------------------
      !     Compute directional factor 
      !     (from Lott & Miller 97, modified by A. Zadra)
      !     and optional nonlinear amplification
      !     (from Wells et al. 2008, QJRMS)

      piotwo = .5*acos(-1.0)

      do i=1,jyes

         !     Angle between mean wind and topography:
         if ( abs(vav(i)) .lt. v0 .and. abs(uav(i)) .lt. v0) then 
            psi = 0.
         else
            psi = theta(i) - atan2(vav(i),uav(i))
            if (psi .gt.   piotwo )  psi = psi - 2.*piotwo 
            if (psi .lt. (-piotwo))  psi = psi + 2.*piotwo
         endif
         cpsi = (cos(psi))*(cos(psi))
         spsi = (sin(psi))*(sin(psi))
         gmsi = (gamma(i))*(gamma(i))

         !     Directional factor:
         tmp = cpsi + gmsi*spsi 
         if (tmp .lt. 1.e-6) then 
            amp = 0.
         else 
            tmp = 1.0/(cpsi + gmsi*spsi)
            tmp = max( ( gmsi*cpsi + spsi )*tmp , 0. )
            amp = max( 2. - sqrt( tmp ) , 0. )
         endif
         
         !     Eccentricity factor
         B = 1. - 0.18*gamma(i) - 0.04*gmsi
         C = 0.48*gamma(i) + 0.30*gmsi
         aux = max(B*cpsi + C*spsi,0.)

         !     Combined directional/eccentricity factor

         fdir(i) = amp*aux

         !     Direction-dependent amplification factor due to nonlinearities
         !     (based on Fig.12b of Wells et al. 2008)

         if (nldirfactor) then
            ampd(i) = 1.0 + 0.5*abs(cos(psi))
         endif

      enddo
      !--------------------------------------------------------------------
      !     Compute the enhancement factor due to stability:
      !     (based on Vosper et al. 2009, QJRMS)

      IF_STABFAC: if (stabfactor) then

         !     estimate bulk depth

         dscale = hmin2
         dscale(1:jyes) = env(1:jyes) + velav(1:jyes)/nav(1:jyes)
         dscale(1:jyes) = max(dscale(1:jyes),hmin2)
         dscale(1:jyes) = min(dscale(1:jyes),zb(1:jyes,1))

         !     estimate means in bulk depth

         uavd(1:jyes) = v0
         istat = int_profile(uavd(1:jyes),u(1:jyes,:),zb(1:jyes,:), &
              0.*dscale(1:jyes),dscale(1:jyes))
         uavd(1:jyes) = uavd(1:jyes)/dscale(1:jyes)

         vavd(1:jyes) = v0
         istat = int_profile(vavd(1:jyes),v(1:jyes,:),zb(1:jyes,:), &
              0.*dscale(1:jyes),dscale(1:jyes))
         vavd(1:jyes) = vavd(1:jyes)/dscale(1:jyes)     

         velavd(1:jyes) = sqrt(uavd(1:jyes)*uavd(1:jyes) &
              + vavd(1:jyes)*vavd(1:jyes))
         velavd(1:jyes) = max(velavd(1:jyes),v0)

         navd(1:jyes) = nmin
         istat = int_profile(navd(1:jyes),bvfreq(1:jyes,:),zb(1:jyes,:), &
              0.*dscale(1:jyes),dscale(1:jyes))
         navd(1:jyes) = navd(1:jyes)/dscale(1:jyes)
         navd(1:jyes) = max(navd(1:jyes),nmin)


         !     5- compute Froude number Fd and enhancement factor
         !        (empirical relation based on fig.4 of Vosper et al. 2009)

         do i=1,jyes
            cdf(i) = cdmin
            fdscale(i) = velavd(i) / max( env(i)*navd(i), 1.e-6 )
            if ( fdscale(i) <= 3.0 ) then
               cdf(i) = (6.0-cdmin)*exp(-2.8*fdscale(i)*fdscale(i)) + cdmin
            endif
         enddo

      endif IF_STABFAC

      !--------------------------------------------------------------------
      !     Compute llb tendencies:

      do l=ilev,1,-1
         do i=1,jyes
            if ( velav(i) >= vmin2 .and. env(i) >= hmin2 .and. &
                 hblk(i) >= hmin2 ) then

               !           Height factor:
               fvert = sqrt( max((hblk(i) - zb(i,l)),0.) &
                    /(0.5*env(i)   + zb(i,l)) )
                    
               !           Mountain-base reduction factor
!!$               fbas = tanh( zb(i,l)/(0.025*env(i)) )
                    
               !           Wind-speed factor:
               fwind = 0.5*( 1. + tanh((velav(i) - sgo_windfac(1))/ &
                    max(sgo_windfac(2),0.01)) )

               !           Implicit calculation of llb tendencies:
               ks2 = 0.5*cdblk*cdf(i)*slpf(i)*ampd(i)*fdir(i)*fvert*fwind!*fbas
               vent = sqrt(u(i,l)*u(i,l) + v(i,l)*v(i,l))
               utendllb(i,l) = -ks2*u(i,l)*vent/(1.0+ks2*vent*tau) 
               vtendllb(i,l) = -ks2*v(i,l)*vent/(1.0+ks2*vent*tau)               

            endif
         enddo
      enddo

   endif IF_BLO
   
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     TOTAL DRAG AND RESULTING WIND FIELD
   !--------------------------------------------------------------------

   utendtot(1:len,:) = 0.
   vtendtot(1:len,:) = 0.

   !     Add and scatter tendencies
   utendtot(1:len,:) = utendgwd(1:len,:) + utendllb(1:len,:)
   vtendtot(1:len,:) = vtendgwd(1:len,:) + vtendllb(1:len,:)

   do l=1,ilev
      do i=1,len
         ii = drag(i) + il1 - 1
         utend(ii,l) = utendtot(i,l)
         vtend(ii,l) = vtendtot(i,l)
         ttend(ii,l) = (-1./CPD) * (u(i,l)*utendtot(i,l) + v(i,l)*vtendtot(i,l))
      enddo
   enddo
 
   if (split_tend_L) then
      utendgwd4 = 0.
      vtendgwd4 = 0.
      do l=1,ilev
         do i=1,len
            ii = drag(i) + il1 - 1
            utendgwd4(ii,l) = real(utendgwd(i,l))
            vtendgwd4(ii,l) = real(vtendgwd(i,l))
         enddo
      enddo
   endif

   !--------------------------------------------------------------------
   return
end subroutine sgo16flx3
