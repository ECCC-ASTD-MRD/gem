module sgo16flx
   use integrals, only: int_profile
   use phy_options, only: sgo_windfac
   use tdpack, only: CPD, GRAV, RGASD
   use ens_perturb, only: ens_nc2d, ens_spp_get
   implicit none
   private
   public :: sgo16flx3

contains


!/@*
subroutine sgo16flx3(uu, vv, utend, vtend, ttend, utendgwd4, vtendgwd4, &
     tth, ttf, ss, ssh, ssexpk, sshexpk,  &
     ilev, ni,  &
     tau, taufac,  &
     gc, height, slope, xcent, mtdir,  &
     cdmin, fc, phic, stabfactor, nldirfactor, mrk2, &
     gwdrag, blocking, split_tend_L)
   implicit none
!!!#include <arch_specific.hf>
   !@Arguments
   logical, intent(in) :: gwdrag                   ! .true. for gravity wave drag
   logical, intent(in) :: blocking                 ! .true. for low-level blocking
   logical, intent(in) :: stabfactor               ! .true. for stability enhancement factor
   logical, intent(in) :: nldirfactor              ! .true. for directional enhacement factor
   logical, intent(in) :: split_tend_L             ! .true. for keeping GWD uvtend seprated from total
   integer, intent(in) :: ilev                     ! number of levels
   integer, intent(in) :: ni                       ! total horizontal dimension
   real, intent(in) :: tau                         ! timestep
   real, intent(in) :: taufac                      ! 1/(length scale) for GWD scheme 
   real, intent(in) :: cdmin                       ! minimum bulk drag coeff for blocking
   real, intent(in) :: fc                          ! tuning factor for blocking height 
   real, intent(in) :: phic                        ! critical phase value for blocking height
   real, dimension(ni), intent(in) :: gc          ! land-cover fraction
   real, dimension(ni), intent(in) :: height      ! mountain (launching) height
   real, dimension(ni), intent(in) :: slope       ! mountain slope
   real, dimension(ni), intent(in) :: xcent       ! mountain eccentricity
   real, dimension(ni), intent(in) :: mtdir       ! mountain direction
   real, dimension(ni,ilev), intent(in) :: uu     ! u-component of wind
   real, dimension(ni,ilev), intent(in) :: vv     ! v-component of wind
   real, dimension(ni,ilev), intent(in) :: tth    ! virtual temperature at thermo levels
   real, dimension(ni,ilev), intent(in) :: ttf    ! virtual temperature at mom levels 
   real, dimension(ni,ilev), intent(in) :: ss     ! sigma at thermo levels
   real, dimension(ni,ilev), intent(in) :: ssh    ! sigma at mom levels
   real, dimension(ni,ilev), intent(in) :: ssexpk
   real, dimension(ni,ilev), intent(in) :: sshexpk
   real, dimension(ni,ens_nc2d), intent(in) :: mrk2 !Markov chains for perturbations
   real, dimension(ni,ilev), intent(out) :: utend ! total tendency for u
   real, dimension(ni,ilev), intent(out) :: vtend ! total tendency for v
   real, dimension(ni,ilev), intent(out) :: ttend ! total tendency for t
   real, dimension(ni,ilev), intent(out) :: utendgwd4 ! GWD tendency for u
   real, dimension(ni,ilev), intent(out) :: vtendgwd4 ! GWD tendency for v
    
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
   
   !  Local variable declarations
   integer :: i, l, jyes, jno
   integer, dimension(ni) :: drag
   
   !-------------------------------------------------------------------
   !     Init output fields

   do l=1,ilev
      do i=1,ni
         utend(i,l) = 0.
         vtend(i,l) = 0.
         ttend(i,l) = 0.
         if (split_tend_L) then
            utendgwd4(i,l) = 0.
            vtendgwd4(i,l) = 0.
         endif
      enddo
   enddo

   !     Find and gather active grid columns

   jyes = 0
   jno  = ni + 1
   do i=1,ni
      if (abs(nint(gc(i))) == 1.) then
         jyes       = jyes + 1
         drag(jyes) = i
      else
         jno        = jno - 1
         drag(jno)  = i
      endif
   enddo

   !     Check if there is AT LEAST ONE active column
   if (jyes == 0) return

   call priv_sgo16flx(uu, vv, utend, vtend, ttend, utendgwd4, vtendgwd4, &
        tth, ttf, ss, ssh, ssexpk, sshexpk, &
        ilev, ni,  jyes, &
        tau, taufac,  &
        drag, height, slope, xcent, mtdir,  &
        cdmin, fc, phic, stabfactor, nldirfactor, mrk2, &
        gwdrag, blocking, split_tend_L)

   !-------------------------------------------------------------------
   return
end subroutine sgo16flx3


!/@*
subroutine priv_sgo16flx(uu, vv, utend, vtend, ttend, utendgwd4, vtendgwd4, &
     tth, ttf, ss, ssh, ssexpk, sshexpk, &
     ilev, ni,  jyes, &
     tau, taufac,  &
     drag, height, slope, xcent, mtdir,  &
     cdmin, fc, phic, stabfactor, nldirfactor, mrk2, &
     gwdrag, blocking, split_tend_L)
   implicit none
!!!#include <arch_specific.hf>
   !@Arguments
   logical, intent(in) :: gwdrag                   ! .true. for gravity wave drag
   logical, intent(in) :: blocking                 ! .true. for low-level blocking
   logical, intent(in) :: stabfactor               ! .true. for stability enhancement factor
   logical, intent(in) :: nldirfactor              ! .true. for directional enhacement factor
   logical, intent(in) :: split_tend_L             ! .true. for keeping GWD uvtend seprated from total
   integer, intent(in) :: ilev                     ! number of levels
   integer, intent(in) :: ni                       ! total horizontal dimension
   integer, intent(in) :: jyes                     ! nb of gathered points
   real, intent(in) :: tau                         ! timestep
   real, intent(in) :: taufac                      ! 1/(length scale) for GWD scheme 
   real, intent(in) :: cdmin                       ! minimum bulk drag coeff for blocking
   real, intent(in) :: fc                          ! tuning factor for blocking height 
   real, intent(in) :: phic                        ! critical phase value for blocking height
   integer, dimension(ni), intent(in) :: drag
   real, dimension(ni), intent(in) :: height      ! mountain (launching) height
   real, dimension(ni), intent(in) :: slope       ! mountain slope
   real, dimension(ni), intent(in) :: xcent       ! mountain eccentricity
   real, dimension(ni), intent(in) :: mtdir       ! mountain direction
   real, dimension(ni,ilev), intent(in) :: uu     ! u-component of wind
   real, dimension(ni,ilev), intent(in) :: vv     ! v-component of wind
   real, dimension(ni,ilev), intent(in) :: tth    ! virtual temperature at thermo levels
   real, dimension(ni,ilev), intent(in) :: ttf    ! virtual temperature at mom levels 
   real, dimension(ni,ilev), intent(in) :: ss     ! sigma at thermo levels
   real, dimension(ni,ilev), intent(in) :: ssh    ! sigma at mom levels
   real, dimension(ni,ilev), intent(in) :: ssexpk
   real, dimension(ni,ilev), intent(in) :: sshexpk
   real, dimension(ni,ens_nc2d), intent(in) :: mrk2 !Markov chains for perturbations
   real, dimension(ni,ilev), intent(out) :: utend ! total tendency for u
   real, dimension(ni,ilev), intent(out) :: vtend ! total tendency for v
   real, dimension(ni,ilev), intent(out) :: ttend ! total tendency for t
   real, dimension(ni,ilev), intent(out) :: utendgwd4 ! GWD tendency for u
   real, dimension(ni,ilev), intent(out) :: vtendgwd4 ! GWD tendency for v

   !  Local parameter definitions
   real, parameter :: vmin1  = 2.     ! minimum wind speed to produce GWD tendencies
   real, parameter :: vmin2  = 0.01   ! minimum wind speed to produce blocking tendencies
   real, parameter :: v0    = 1.e-6   ! epsilon for wind speed
   real, parameter :: hmin1  = 10.    ! minimum height to launch gravity waves
   real, parameter :: hmin2  = 3.     ! minimum height to allow low-level blocking
   real, parameter :: cdblk = 1.      ! bulk drag coefficient of blocking (approx. 1)
   real, parameter :: href  = 130.    ! reference height for GWD
   real, parameter :: nmin  = 1.e-8   ! epsilon for buoyancy freqency
   real, parameter :: pi = acos(-1.0)

   !  Local variable declarations
   integer :: i, l, ii, lref, lrefm, istat
   integer, dimension(ni) :: izt0

   real :: tmp, aux, eta, psi, cpsi, spsi, ratio, fvert, fwind, amp, gmsi, &
        B, C, ks2, vent, fac, utendtot, vtendtot, rnldirfactor, &
        asqi, asqs, velavd, fdscale
   real :: phicni(ni)
   real, dimension(jyes) :: ub, vb, vmodb, env, slp2, slpf, gamma, theta, &
        hblk, uav, vav, velav, nav, fdir, dscale, cdf, &
        ampd, uavd, vavd, navd, vphic, work1, zeros
   real, dimension(jyes, ilev) :: u, v, tf, th, s, sh, sexpk, shexpk, bvf, &
        bvfreq, veln, asq, dfac, depfac, denfac, grad, zb, &
        utendgwd, vtendgwd, utendllb, vtendllb

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     INITIAL STEPS
   !--------------------------------------------------------------------
   rnldirfactor = 0.
   if (nldirfactor) rnldirfactor = 1.
   
   !     Retrieve stochastic perturbations on request
   phicni = ens_spp_get('sgo_phic', mrk2, default=phic)
   
   ! Gather active columns
   do i=1,jyes
      ii = drag(i)
      env(i)   = height(ii) 
      slp2(i)  = slope(ii)
      gamma(i) = xcent(ii)
      theta(i) = mtdir(ii)
      vphic(i) = phicni(ii)
      slpf(i) = slp2(i)/max(env(i),hmin2) 
      zeros(i) = 0.
   enddo

   do l=1,ilev
      do i=1,jyes
         ii      = drag(i)
         u(i,l)  = uu(ii,l)
         v(i,l)  = vv(ii,l)
         tf(i,l) = ttf(ii,l)
         th(i,l) = tth(ii,l)
         s(i,l)  = ss(ii,l)
         sh(i,l) = ssh(ii,l)
         sexpk(i,l)  = ssexpk(ii,l)
         shexpk(i,l) = sshexpk(ii,l)
      enddo
   enddo
   
   !     Compute auxiliary profiles

   IF_GWD_BLO: if (gwdrag .or. blocking) then
      !--------------------------------------------------------------------
      !     Build elevation field:

      do i=1,jyes
         zb(i,ilev)  = -(RGASD/GRAV)*th(i,ilev)*alog(s(i,ilev))
      enddo
      do l=ilev-1,1,-1
         do i=1,jyes
            zb(i,l)  = zb(i,l+1) + &
                 (RGASD/GRAV)*th(i,l)*alog(s(i,l+1)/s(i,l))          
         enddo
      enddo

      !-------------------------------------------------------------------
      !     Calculate BF frequency:

      !#TODO: compute in gwd, shared in gwspectrum?
      do l=2,ilev
         do i=1,jyes
            aux = (th(i,l)/shexpk(i,l)-th(i,l-1)/shexpk(i,l-1)) &
                 /(sh(i,l)-sh(i,l-1))
            aux = min(aux, -5./s(i,l))
            aux = -aux*s(i,l)*GRAV/(RGASD*tf(i,l))
            bvf(i,l) = GRAV*aux*sexpk(i,l)/tf(i,l)
            bvf(i,l) = sqrt(bvf(i,l))
         enddo
      enddo
      do i=1,jyes
         bvf(i,1) = bvf(i,2)
      enddo

   endif IF_GWD_BLO

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     GRAVITY-WAVE DRAG
   !--------------------------------------------------------------------
   !     Initialize arrays

   do i=1,jyes
      izt0(i)  = ilev-1
   enddo
   do l=1,ilev
      do i=1,jyes
         utendgwd(i,l) = 0.
         vtendgwd(i,l) = 0.
      enddo
   enddo

   IF_GWD: if (gwdrag) then
      !--------------------------------------------------------------------
      !     Smooth B V frequency
      !#TODO: this smoothing is also done in gwspectrum
      do l=1,ilev
         bvfreq(1:jyes,l) = bvf(1:jyes,l)
      enddo
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
         vmodb(i) = max(1., sqrt(vmodb(i)))
         ub(i) = u(i,lrefm)/vmodb(i)
         vb(i) = v(i,lrefm)/vmodb(i)
      enddo

      !--------------------------------------------------------------------
      !     Project wind field on reference wind:

      do l=1,ilev
         do i=1,jyes 
            veln(i,l) = u(i,l)*ub(i)+v(i,l)*vb(i)
            veln(i,l) = max(v0, veln(i,l))
         enddo
      enddo

      !--------------------------------------------------------------------
      !     Stress field
      !     Compute stress at reference level:

      do i=1,jyes
         lrefm = izt0(i)
         lref  = lrefm + 1
         asq(i,lref)  = env(i)*env(i)
         dfac(i,lref) = bvfreq(i,lrefm)*sh(i,lrefm)*vmodb(i) &
              /tf(i,lrefm)
      enddo

      !     Compute stress at other levels (bottom-up):
      
      do l=ilev-1,1,-1
         do i=1,jyes
            lrefm = izt0(i)
            if (l <= lrefm) then
               dfac(i,l) = bvfreq(i,l)*s(i,l)*veln(i,l)/tf(i,l)
               asqi = asq(i,l+1)*dfac(i,l+1)/dfac(i,l)
               if (veln(i,l) >= 1.) then
                  tmp = veln(i,l)/bvfreq(i,l)
                  asqs = 0.5*tmp*tmp
               else
                  asqs = 0.
               endif
               asq(i,l)    = min(asqi, asqs)
               depfac(i,l) = (taufac*GRAV/RGASD)*asq(i,l)*dfac(i,l)
            endif
         enddo
      enddo
      do l=ilev,1,-1
         do i=1,jyes
            lrefm = izt0(i)
            if (l > lrefm) then
               depfac(i,l) = depfac(i,lrefm)
            endif
         enddo
      enddo

      !--------------------------------------------------------------------
      !     Compute gwd tendencies:

      do i=1,jyes
         work1(i) = 1.
         if (vmodb(i).lt.vmin1 .or. env(i).lt.hmin1) work1(i) = 0.
         if ((depfac(i,2) - depfac(i,1)).gt.1.e-10) then
            eta = 1.
         else
            eta = 0.
         endif

         tmp = 2.*sh(i,1) + eta*3.*tau*depfac(i,1)/veln(i,1) 
         tmp = 1.0/tmp
         grad(i,1) = 2.*eta*depfac(i,1) &
              *tmp
         utendgwd(i,1) = work1(i) * (-ub(i)*grad(i,1))
         vtendgwd(i,1) = work1(i) * (-vb(i)*grad(i,1))
         denfac(i,1) = grad(i,1)*3.*tau*depfac(i,1)/veln(i,1)
      enddo
      
      do l=2,ilev-1
         do i=1,jyes
            if (l <= izt0(i)) then
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
               utendgwd(i,l) = work1(i) * (-ub(i)*grad(i,l))
               vtendgwd(i,l) = work1(i) * (-vb(i)*grad(i,l))
               denfac(i,l) = grad(i,l)*3.*tau*depfac(i,l)/veln(i,l)
            endif
         enddo
      enddo

   endif IF_GWD


   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     LOW-LEVEL BLOCKING
   !--------------------------------------------------------------------
   !     Initialize arrays

   do l=1,ilev
      do i=1,jyes
         utendllb(i,l) = 0.
         vtendllb(i,l) = 0.
      enddo
   enddo

   IF_BLO: if (blocking) then
      !--------------------------------------------------------------------
      !     Compute average wind, wind speed and buoyancy

      do i=1,jyes
         work1(i) = 0.5 * env(i)
      enddo

      istat = int_profile(uav, u, zb, work1, env)
      istat = int_profile(vav, v, zb, work1, env)
      istat = int_profile(nav, bvf, zb, zeros, env)

      do i=1,jyes
         uav(i) = 2.*uav(i)/max(env(i),hmin2)
         vav(i) = 2.*vav(i)/max(env(i),hmin2)
         velav(i) = max(sqrt(uav(i)*uav(i) + vav(i)*vav(i)), v0)
         nav(i) = max(nav(i)/max(env(i), hmin2), nmin)

         !     Compute blocking height
         hblk(i) = max(fc*env(i) - vphic(i)*velav(i)/nav(i), 0.)
      enddo

      !--------------------------------------------------------------------
      !     Compute directional factor 
      !     (from Lott & Miller 97, modified by A. Zadra)
      !     and optional nonlinear amplification
      !     (from Wells et al. 2008, QJRMS)

      do i=1,jyes

         !     Angle between mean wind and topography:
         !#TODO: can we avoid this protection? with the scaling factor for psi, we always come back to the x>0 branch of atan2
         if ( abs(vav(i)) .lt. v0 .and. abs(uav(i)) .lt. v0) then 
            psi = 0.
         else
            psi = theta(i) - atan2(vav(i),uav(i))
            fac = psi/(0.5*pi)
            fac = sign(float(floor(abs(fac))),fac)
            psi = psi - pi * max(-1.,min(fac,1.))
         endif
         cpsi = (cos(psi))*(cos(psi))
         spsi = (sin(psi))*(sin(psi))
         gmsi = (gamma(i))*(gamma(i))

         !     Directional factor:
         tmp = cpsi + gmsi*spsi
         !#TODO: can we avoid this protection?
         if (tmp .lt. 1.e-6) then 
            amp = 0.
         else 
            tmp = 1.0/tmp
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

         ampd(i) = 1.0 + rnldirfactor*0.5*abs(cos(psi))

      enddo
      !--------------------------------------------------------------------
      !     Compute the enhancement factor due to stability:
      !     (based on Vosper et al. 2009, QJRMS)

      IF_STABFAC: if (stabfactor) then

         !     estimate bulk depth

         do i=1,jyes
            dscale(i) = env(i) + velav(i)/nav(i)
            dscale(i) = max(dscale(i),hmin2)
            dscale(i) = min(dscale(i),zb(i,1))
         enddo
         
         !     estimate means in bulk depth

         istat = int_profile(uavd, u, zb, zeros, dscale)
         istat = int_profile(vavd, v, zb, zeros, dscale)
         istat = int_profile(navd, bvf, zb, zeros, dscale)

         do i=1,jyes
            uavd(i) = uavd(i)/dscale(i)

            vavd(i) = vavd(i)/dscale(i)     

            velavd = sqrt(uavd(i)*uavd(i) &
                 + vavd(i)*vavd(i))
            velavd = max(velavd,v0)

            navd(i) = navd(i)/dscale(i)
            navd(i) = max(navd(i),nmin)

            !     5- compute Froude number Fd and enhancement factor
            !        (empirical relation based on fig.4 of Vosper et al. 2009)

            cdf(i) = cdmin
            fdscale = velavd / max( env(i)*navd(i), 1.e-6 )
            if ( fdscale <= 3.0 ) then
               cdf(i) = (6.0-cdmin)*exp(-2.8*fdscale*fdscale) + cdmin
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

   !     Add and scatter tendencies

   do l=1,ilev
      do i=1,jyes
         ii = drag(i)
         utendtot = utendgwd(i,l) + utendllb(i,l)
         vtendtot = vtendgwd(i,l) + vtendllb(i,l)
         utend(ii,l) = utendtot
         vtend(ii,l) = vtendtot
         ttend(ii,l) = (-1./CPD) * (u(i,l)*utendtot + v(i,l)*vtendtot)
      enddo
   enddo
 
   if (split_tend_L) then
      do l=1,ilev
         do i=1,jyes
            ii = drag(i)
            utendgwd4(ii,l) = utendgwd(i,l)
            vtendgwd4(ii,l) = vtendgwd(i,l)
         enddo
      enddo
   endif

   !--------------------------------------------------------------------
   return
end subroutine priv_sgo16flx

end module sgo16flx
