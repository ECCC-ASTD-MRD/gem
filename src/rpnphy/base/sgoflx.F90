!-------------------------------------- LICENCE BEGIN ------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html

!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENSEE END --------------------------

!/@*
subroutine sgoflx8(uu,vv,utend,vtend,ttend,utendgwd4,vtendgwd4,&
     tth,ttf,ss,ssh, &
     ilev,ilg,il1,il2, &
     tau,taufac, &
     gc,height,slope,xcent,mtdir, &
     gwdrag,blocking, &
     stabfactor,nldirfactor,cdmin, split_tend_L)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use  tdpack, only: CAPPA, CPD, GRAV, RGASD
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   logical gwdrag,blocking, &
        stabfactor,nldirfactor, split_tend_L

   integer ilev,ilg,il1,il2
   real tau,taufac,cdmin
   real(REAL64) :: uu(ilg,ilev),      vv(ilg,ilev),     utend(ilg,ilev), &
        vtend(ilg,ilev),   ttend(ilg,ilev),  tth(ilg,ilev),    ttf(ilg,ilev), &
        ss(ilg,ilev),      ssh(ilg,ilev),    gc(ilg), &
        height(ilg),       slope(ilg),       xcent(ilg), &
        mtdir(ilg)
   real, intent(out) :: utendgwd4(ilg,ilev), vtendgwd4(ilg,ilev)
   !@Author
   !        A. Zadra - May 2002 - From lin_sgoflx1
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

   ! 012    A. Zadra - Critical value phic (used in calculation of blocking height)
   !                is adjusted according to the resolution of the lowest layer
   !@Object
   !        Simplified version of subgrid orographic drag (sgoflx2) scheme:
   !        - reduced, non-smoothed buoyancy frequency
   !        - shortened gravity-wave drag (McFarlane 87)
   !        - shortened low-level blocking (Lott & Miller 97)
   !        - orographic lift (not yet included)
   !        - lee-wave breaking (not yet included)
   !@Arguments
   !          - Input/Output -
   ! UU       U component of wind as input
   ! VV       V component of wind as input
   !          - Input -
   ! TTH      virtual temperature level means
   ! TTF      virtual temperature
   ! SS       sigma levels
   ! SSH      intermediate levels
   !          - Output -
   ! UTEND    total tendency on U (gwd + blocking)
   ! VTEND    total tendency on V (gwd + blocking)
   ! TTEND    total tendency on T (gwd + blocking)
   ! UTENDGWD4  GWD tendency on U
   ! VTENDGWD4  GWD tendency on V
   !          - Input -
   ! ILEV     number of levels
   ! ILG      total horizontal dimension
   ! IL1      1st dimension of U,V,T to start calculation
   ! IL2      index on horizontal dimension to stop calculation
   ! TAU      timestep times a factor: 1 for two time-level models
   !                                   2 for three time-level models
   ! TAUFAC   1/(length scale).
   ! GC       ground cover (land-sea mask) : between 0(sea) and 1(land)
   ! HEIGHT   launching height (variance associated with unresolved orography)
   ! SLOPE    mountain slope
   ! XCENT    mountain eccentricity
   ! MTDIR    mountain direction
   ! GWDRAG   .true. for gravity wave drag
   ! BLOCKING .true. for blocking
   !*@/

!!$      real(REAL64), dimension(ILG     ) :: VMOD
!!$      real(REAL64), dimension(ILG     ) :: UUB
!!$      real(REAL64), dimension(ILG     ) :: VVB
   INTEGER, dimension(ILG     ) :: DRAG
   real(REAL64), dimension(ILG     ) :: UB
   real(REAL64), dimension(ILG     ) :: VB
   real(REAL64), dimension(ILG     ) :: VMODB

   real(REAL64), dimension(ILG     ) :: ENV

   real(REAL64), dimension(ILG     ) :: SLP2
   real(REAL64), dimension(ILG     ) :: SLPF
   real(REAL64), dimension(ILG     ) :: GAMMA
   real(REAL64), dimension(ILG     ) :: THETA
   INTEGER, dimension(ILG     ) :: IZT1
   INTEGER, dimension(ILG     ) :: IZT2
   INTEGER, dimension(ILG     ) :: IZT3
   INTEGER, dimension(ILG     ) :: IZB
   real(REAL64), dimension(ILG     ) :: HBLK
   real(REAL64), dimension(ILG     ) :: UAV
   real(REAL64), dimension(ILG     ) :: VAV
   real(REAL64), dimension(ILG     ) :: VELAV
   real(REAL64), dimension(ILG     ) :: DELZ
   real(REAL64), dimension(ILG     ) :: FDIR
   real(REAL64), dimension(ILG     ) :: BLOFF
   INTEGER, dimension(ILG     ) :: IZT0
   INTEGER, dimension(ILG     ) :: IZTD
   real(REAL64), dimension(ILG     ) :: DSCALE
   real(REAL64), dimension(ILG     ) :: NDSCALE
   real(REAL64), dimension(ILG     ) :: FDSCALE
   real(REAL64), dimension(ILG     ) :: CDF
   real(REAL64), dimension(ILG     ) :: NAV
   INTEGER, dimension(ILG     ) :: LZI
   INTEGER, dimension(ILG     ) :: ZOFF
   real(REAL64), dimension(ILG     ) :: AMPD
   real(REAL64), dimension(ILG     ) :: HZN
   real(REAL64), dimension(ILG     ) :: CONV

   real(REAL64), dimension(ILG,ILEV) :: U
   real(REAL64), dimension(ILG,ILEV) :: V
   real(REAL64), dimension(ILG,ILEV) :: TF
   real(REAL64), dimension(ILG,ILEV) :: TH
   real(REAL64), dimension(ILG,ILEV) :: S
   real(REAL64), dimension(ILG,ILEV) :: SH
   real(REAL64), dimension(ILG,ILEV) :: SEXPK
   real(REAL64), dimension(ILG,ILEV) :: SHEXPK
   real(REAL64), dimension(ILG,ILEV) :: BVFREQ
   real(REAL64), dimension(ILG,ILEV) :: UTENDGWD
   real(REAL64), dimension(ILG,ILEV) :: VTENDGWD
   real(REAL64), dimension(ILG,ILEV) :: VELN
   real(REAL64), dimension(ILG,ILEV) :: ASQ
   real(REAL64), dimension(ILG,ILEV) :: ASQI
   real(REAL64), dimension(ILG,ILEV) :: ASQS
   real(REAL64), dimension(ILG,ILEV) :: DFAC
   real(REAL64), dimension(ILG,ILEV) :: DEPFAC
   real(REAL64), dimension(ILG,ILEV) :: GRAD
   real(REAL64), dimension(ILG,ILEV) :: DENFAC
   real(REAL64), dimension(ILG,ILEV) :: ZB
   real(REAL64), dimension(ILG,ILEV) :: UTENDLLB
   real(REAL64), dimension(ILG,ILEV) :: VTENDLLB
   real(REAL64), dimension(ILG,ILEV) :: UTENDTOT
   real(REAL64), dimension(ILG,ILEV) :: VTENDTOT

   real(REAL64), dimension(ILG     ) :: QDvtmp1
   real(REAL64), dimension(ILG     ) :: QDvtmp2
!!$      real(REAL64), dimension(ILG     ) :: QDvtmp3

   !***********************************************************************

   integer i,l,ii,len,lref,lrefm,jyes,jno,it
   real(REAL64) :: aux,eta,dz,uparl,piotwo,psi,cpsi,spsi,ratio, &
        fvert,amp,vmin,v0,zero,unit,cdblk, &
        hmin1,hmin2,gmsi
   real(REAL64) :: QDrgas, QDgrav,QDhmin2
   real(REAL64) :: QDtmp
   real(REAL64) :: ks2,vent
   real(REAL64) :: href,nd,fc,phic


   vmin  = 2.
   v0    = 1.e-30
   hmin1  = 10.
   hmin2  = 3.
   QDhmin2= 1.0/hmin2
   QDrgas = 1.0/dble(RGASD)
   QDgrav = 1.0/dble(GRAV)
   zero  = 0.
   unit  = 1.
   cdblk = 1.
   href  = 130.
   fc    = 1.5

   !     Note: The values of some parameters are not uniquely defined
   !           and are subject to tuning, such as:
   !              hmin1 = minimum height to launch gravity waves
   !              hmin2 = minimum height to trigger low-level blocking
   !              cdblk = bulk drag coefficient of blocking (approx. 1)

   len   = il2 - il1 + 1

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     INITIAL STEPS

   !--------------------------------------------------------------------
   !     Initialize total tendency

   do l=1,ilev
      do i=il1,il2
         utend(i,l) = zero
         vtend(i,l) = zero
         ttend(i,l) = zero
      enddo
   enddo

   !-------------------------------------------------------------------
   !     Find and gather active grid columns


   jyes = 0
   jno  = len + 1

   do i=il1,il2
      if ( gc(i).eq.-1. ) then
         jyes       = jyes + 1
         drag(jyes) = i
      else
         jno         = jno - 1
         drag(jno)  = i
      endif
   enddo

   !     Check if there is AT LEAST ONE active column

   if (jyes.le.0) then
      goto 600
   endif

   do i=1,len
      ii = drag(i) + il1 - 1
      env(i)   = height(ii)
      slp2(i)  = slope(ii)
      gamma(i) = xcent(ii)
      theta(i) = mtdir(ii)
   enddo

   do i=1,len
      if (env(i) .lt. hmin2) then
         slpf(i) = slp2(i)*Qdhmin2
      else
         slpf(i) = slp2(i)/env(i)
      endif
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
      enddo
   enddo

   call vpown1 (sexpk ,s ,dble(CAPPA),len*ilev)
   call vpown1 (shexpk,sh,dble(CAPPA),len*ilev)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     GRAVITY-WAVE DRAG

   !--------------------------------------------------------------------
   !     Initialize arrays
   do l=1,ilev
      do i=1,len
         utendgwd(i,l) = zero
         vtendgwd(i,l) = zero
      enddo
   enddo

   do i=1,len
      izt0(i)  = ilev-1
   enddo

   if (gwdrag) then
      !--------------------------------------------------------------------
      !     Build elevation field:

      call VLOG(QDvtmp1(1),sh(1,ilev),jyes)
      do i=1,jyes
         zb(i,ilev)  = -(RGASD*QDgrav)*tf(i,ilev)*QDvtmp1(i)
      enddo
      do l=ilev-1,1,-1
         call VLOG(QDvtmp1(1),sh(1,l+1),jyes)
         call VLOG(QDvtmp2(1),sh(1,l  ),jyes)
         do i=1,jyes
            zb(i,l)  = zb(i,l+1) + &
                 (RGASD*QDgrav)*tf(i,l)*(QDvtmp1(i)-QDvtmp2(i))
         enddo
      enddo
      !-------------------------------------------------------------------
      !     Calculate BF frequency:

      bvfreq = 0.
      do l=2,ilev
         do i=1,jyes
            aux = (th(i,l)/shexpk(i,l)-th(i,l-1)/shexpk(i,l-1)) &
                 /(sh(i,l)-sh(i,l-1))
            if ( aux.ge.(-5./s(i,l)) ) aux = -5./s(i,l)
            aux = -aux*s(i,l)*GRAV/(RGASD*tf(i,l))
            bvfreq(i,l) = GRAV*aux*sexpk(i,l)/tf(i,l)
         enddo
         call vsqrt(bvfreq(1,l),bvfreq(1,l),jyes)
      enddo
      do i=1,jyes
         bvfreq(i,1) = bvfreq(i,2)
      enddo

      !     * Smooth B V frequency

      do l=2,ilev
         call VLOG(QDvtmp2(1),s(1,L),  len)
         call VLOG(QDvtmp1(1),s(1,L-1),len)
         do i=1,jyes
            ratio = 5.*(QDvtmp2(i) - QDvtmp1(i))
            bvfreq(i,l) = (bvfreq(i,l-1) + ratio*bvfreq(i,l)) &
                 /(1.+ratio)
         enddo
      enddo
      !--------------------------------------------------------------------
      !     Find reference level:
      do l=ilev-1,1,-1
         do i=1,jyes
            if (zb(i,l).le.href)    izt0(i)  = l
         enddo
      enddo
      !--------------------------------------------------------------------
      !     Wind and unit vector at level LREFM:

      do i=1,jyes
         lrefm = izt0(i)
         vmodb(i) = u(i,lrefm)*u(i,lrefm) + v(i,lrefm)*v(i,lrefm)
      enddo
      call VSQRT(vmodb,vmodb,jyes)
      do i=1,jyes
         lrefm = izt0(i)
         if (vmodb(i).le.unit)  vmodb(i) = unit
         QDtmp = 1.d0/vmodb(i)
         ub(i) = u(i,lrefm)*QDtmp
         vb(i) = v(i,lrefm)*QDtmp
      enddo
      !--------------------------------------------------------------------
      !     Project wind field on reference wind:

      do l=1,ilev
         do i=1,jyes
            veln(i,l) = u(i,l)*ub(i)+v(i,l)*vb(i)
            if (veln(i,l).le.v0)  veln(i,l) = v0
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
            if (veln(i,l).ge.1.) then
               QDtmp = veln(i,l)/bvfreq(i,l)
               asqs(i,l) = 0.5*Qdtmp*Qdtmp
            else
               asqs(i,l) = zero
            endif
            if (asqi(i,l).le.asqs(i,l)) then
               asq(i,l)    = asqi(i,l)
            else
               asq(i,l)    = asqs(i,l)
            endif
            depfac(i,l) = (taufac*GRAV*QDrgas)*asq(i,l)*dfac(i,l)
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

      utendgwd = 0.
      vtendgwd = 0.
      call VREC(QDvtmp1(1),veln(1,1),jyes)
      do i=1,jyes
         if ((depfac(i,2) - depfac(i,1)).gt.1.e-10) then
            eta = 1.
         else
            eta = 0.
         endif

         QDtmp = 2.*sh(i,1) + eta*3.*tau*depfac(i,1)*QDvtmp1(i)
         QDtmp = 1.0/QDtmp
         grad(i,1) = 2.*eta*depfac(i,1) &
              *QDtmp
         utendgwd(i,1) = -ub(i)*grad(i,1)
         vtendgwd(i,1) = -vb(i)*grad(i,1)
         denfac(i,1) = grad(i,1)*3.*tau*depfac(i,1)*QDvtmp1(i)

         lrefm = izt0(i)
         lref  = lrefm + 1
         do l=ilev,lref,-1
            utendgwd(i,l) = zero
            vtendgwd(i,l) = zero
         enddo
      enddo

      do i=1,jyes
         lrefm = izt0(i)
         do l=2,lrefm
            call VREC(QDvtmp1(1),veln(1,l),jyes)
            if ((depfac(i,l) - depfac(i,l-1)).gt.1.e-10) then
               eta = 1.
            else
               eta = 0.
            endif
            QDtmp =  2.*(sh(i,l)-sh(i,l-1)) + &
                 eta*3.*tau*depfac(i,l)*QDvtmp1(i)
            QDtmp = 1.0/QDtmp
            grad(i,l) = ( 2.*depfac(i,l)-2.*depfac(i,l-1) + &
                 eta*denfac(i,l-1) )*QDtmp
            utendgwd(i,l) = -ub(i)*grad(i,l)
            vtendgwd(i,l) = -vb(i)*grad(i,l)
            denfac(i,l) = grad(i,l)*3.*tau*depfac(i,l)*QDvtmp1(i)
         enddo
      enddo

      do l=1,ilev
         do i=1,jyes
            if ( vmodb(i).lt.vmin .or. env(i).lt.hmin1 ) then
               utendgwd(i,l) = zero
               vtendgwd(i,l) = zero
            endif
         enddo
      enddo

   endif

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     LOW-LEVEL BLOCKING

   !--------------------------------------------------------------------
   !     Initialize arrays
   do l=1,ilev
      do i=1,len
         zb(i,l)       = zero
         utendllb(i,l) = zero
         vtendllb(i,l) = zero
      enddo
   enddo

   do i=1,len
      izt1(i)  = ilev-1
      izt2(i)  = ilev
      izt3(i)  = ilev
      hblk(i)  = zero
      izb(i)   = ilev
      uav(i)   = v0
      vav(i)   = v0
      velav(i) = v0
      delz(i)  = zero
      fdir(i)  = unit
      bloff(i) = 0
      iztd(i)  = ilev-1
      nav(i)   = zero
      lzi(i)   = ilev-1
      zoff(i)  = 0
      ampd(i)  = zero
      hzn(i)   = zero
      conv(i)  = unit
   enddo

   if (blocking) then
      !--------------------------------------------------------------------
      !     Recalculate non-smoothed buoyancy frequency


      bvfreq = 0.
      do l=2,ilev
         do i=1,jyes
            aux = (th(i,l)/shexpk(i,l)-th(i,l-1)/shexpk(i,l-1)) &
                 /(sh(i,l)-sh(i,l-1))
            if ( aux.ge.(-5./s(i,l)) ) aux = -5./s(i,l)
            aux = -aux*s(i,l)*GRAV/(RGASD*tf(i,l))
            bvfreq(i,l) = GRAV*aux*sexpk(i,l)/tf(i,l)
         enddo
         call vsqrt(bvfreq(1,l),bvfreq(1,l),jyes)
      enddo
      do i=1,jyes
         bvfreq(i,1) = bvfreq(i,2)
      enddo

      !     Build elevation field:

      call VLOG(QDvtmp1(1),sh(1,ilev),jyes)
      do i=1,jyes
         zb(i,ilev)  = -(RGASD*QDgrav)*tf(i,ilev)*QDvtmp1(i)
      enddo
      do l=ilev-1,1,-1
         call VLOG(QDvtmp1(1),sh(1,l+1),jyes)
         call VLOG(QDvtmp2(1),sh(1,l  ),jyes)
         do i=1,jyes
            zb(i,l)  = zb(i,l+1) + &
                 (RGASD*QDgrav)*tf(i,l)*(QDvtmp1(i)-QDvtmp2(i))
         enddo
      enddo

      !--------------------------------------------------------------------
      !     Blocking height

      !     Find maximum blocking level, upper level for averaging and
      !     lower level for averaging:
      do l=ilev-2,1,-1
         do i=1,jyes
            if (zb(i,l).lt.( fc*env(i)))    izt3(i) = l
            if (zb(i,l).lt.     env(i) )    izt1(i)  = l
         enddo
      enddo

      do l=ilev-1,1,-1
         do i=1,jyes
            if (zb(i,l).lt.(0.5*env(i)))    izt2(i) = l
         enddo
      enddo

      !     Compute averages:
      do l=ilev,2,-1
         do i=1,jyes
            if (l.le.izt2(i) .and. l.ge.izt1(i))   then
               dz          = zb(i,l-1) - zb(i,l)
               delz(i)     = delz(i) + dz
               uav(i)      = uav(i)  + dz*u(i,l)
               vav(i)      = vav(i)  + dz*v(i,l)
            endif
         enddo
      enddo
      call VREC(QDvtmp1(1),delz(1),jyes)
      do i=1,jyes
         uav(i)      = uav(i)*QDvtmp1(i)
         vav(i)      = vav(i)*QDvtmp1(i)
         if (abs(vav(i)).lt.v0 .and. abs(uav(i)).lt.v0) then
            velav(i)    = v0
         else
            velav(i)    = uav(i)*uav(i) + vav(i)*vav(i)
         endif
      enddo

      call vsqrt(velav,velav,jyes)
      !
      !     Compute blocking height and blocking level:

      do l=2,ilev
         do i=1,jyes
            if (l.ge.izt3(i) .and. bloff(i).eq.0) then
               dz    = zb(i,l-1) - zb(i,l)
               uparl = (u(i,l)*uav(i) + v(i,l)*vav(i))/velav(i)
               if (uparl .lt. v0) then
                  izb(i)   = l-1
                  bloff(i) = 1
               else
                  hblk(i) = hblk(i) + dz*bvfreq(i,l)/uparl
                  phic = 0.44 + 0.06*tanh((zb(i,ilev)-15.)/4.);
                  if (hblk(i) .gt. phic) then
                     izb(i)   = l-1
                     bloff(i) = 1
                  endif
               endif
            endif
         enddo
      enddo

      !--------------------------------------------------------------------
      !     Compute directional factor
      !     (from Lott & Miller 97, modifed by A. Zadra)
      !     and include nonlinear amplification
      !     (from Wells et al. 2008, QJRMS)

      piotwo = .5*acos(-1.d0)
      do i=1,jyes

         !     Angle between mean wind and topography:
         if ( abs(vav(i)) .lt. v0 .and. abs(uav(i)) .lt. v0) then
            psi = zero
         else
            psi = theta(i) - atan2(vav(i),uav(i))
            if (psi .gt.   piotwo )  psi = psi - 2.*piotwo
            if (psi .lt. (-piotwo))  psi = psi + 2.*piotwo
         endif
         cpsi = (cos(psi))*(cos(psi))
         spsi = (sin(psi))*(sin(psi))
         gmsi = (gamma(i))*(gamma(i))

         !     Directional factor:
         amp = cpsi + gmsi*spsi
         if (amp .lt. 1.e-10) then
            amp = zero
         else
            QDtmp = cpsi + gmsi*spsi
            QDtmp = 1.0/QDtmp
            QDtmp = ( gmsi*cpsi + spsi )*QDtmp
            if (QDtmp.lt.zero)   QDtmp = zero
            amp = 2. - sqrt( QDtmp )
            if (amp.lt.zero)     amp = zero
         endif

         QDtmp = cpsi + gmsi*spsi
         if (QDtmp.lt.zero)   QDtmp = zero
         QDtmp = sqrt( QDtmp )
         fdir(i) = amp*QDtmp

         !     Direction-dependent amplification factor due to nonlinearities
         !     (based on Fig.12b of Wells et al. 2008)

         ampd(i) = 1.0
         if (nldirfactor) then
            ampd(i) = 1.0 + 0.5*abs(cos(psi))
         endif

      enddo
      !--------------------------------------------------------------------
      !     Compute the enhancement factor due to stability:
      !     (based on Vosper et al. 2009, QJRMS)

      !     1- find level corresponding to the depth of the
      !        neutral layer or weakly stratified BL, if it extends
      !        above the mountain height
      !       [Note: here that level is defined from bottom up,
      !        as the first level where the buyoancy frequency
      !        differs by more than 5% w.r.t the lowest level]

      do i=1,jyes
         lzi(i)  = ilev-1
         if ( bvfreq(i,ilev).le.0.005 ) then
            zoff(i) = 0
            do l=ilev-1,2,-1
               if ( zoff(i).eq.0 ) then
                  aux = abs( bvfreq(i,l)/bvfreq(i,ilev) - 1.0 )
                  if ( aux.le.0.05 ) then
                     lzi(i)  = l-1
                  else
                     zoff(i) = 1
                  endif
               endif
            enddo
         endif
         hzn(i) = zb(i,lzi(i))
         if ( env(i) .ge. hzn(i) )   hzn(i) = env(i)
      enddo

      !     2- estimate depth scale D to compute Froude number
      !        - minimum value: 1m
      !        - maximum value: 20000m
      !        and the corresponding nearest level iztd

      do i=1,jyes
         dscale(i) = hzn(i) + velav(i)/bvfreq(i,lzi(i)-1)
         if ( dscale(i).lt.1.    )   dscale(i) = 1.
         if ( dscale(i).ge.2.e+4 )   dscale(i) = 2.e+4
      enddo

      do l=ilev-1,1,-1
         do i=1,jyes
            if (zb(i,l).lt.(dscale(i)))    iztd(i) = l
         enddo
      enddo

      !     3- compute first guess of bulk buoyancy frequency
      !        between surface and depth scale D
      !        - minimum value: 4.e-3 s^-1

      do i=1,jyes
         ndscale(i) = ( GRAV/dscale(i) ) &
              * ( ( th(i,iztd(i))*shexpk(i,ilev   ) ) &
              /( th(i,ilev   )*shexpk(i,iztd(i)) ) - 1.0 )
         if ( ndscale(i).gt.0.0   )  ndscale(i) = sqrt(ndscale(i))
         if ( ndscale(i).lt.4.e-3 )  ndscale(i) = 4.e-3
      enddo


      !     4- interative estimate of buoyancy frequency Nd
      !        Note: convergence criteria are
      !        - 5% difference between consecutive values of Nd
      !        - or maximum 5 iterations

      do i=1,jyes
         conv(i) = 1.
         do it=1,5
            if ( conv(i).ge.0.05 ) then
               dscale(i) = hzn(i) + velav(i)/ndscale(i)
               if ( dscale(i).lt.1.    )   dscale(i) = 1.
               if ( dscale(i).ge.2.e+4 )   dscale(i) = 2.e+4
               iztd(i)  = ilev-1
               do l=ilev-1,1,-1
                  if (zb(i,l).lt.(dscale(i)))    iztd(i) = l
               enddo
               nd = ndscale(i)
               ndscale(i) = ( GRAV/dscale(i) ) &
                    * ( ( th(i,iztd(i))*shexpk(i,ilev   ) ) &
                    /( th(i,ilev   )*shexpk(i,iztd(i)) ) - 1.0 )
               if ( ndscale(i).gt.0.0 )  ndscale(i) = sqrt(ndscale(i))
               if ( ndscale(i).lt.nd  )  ndscale(i) = nd
               conv(i) = abs(ndscale(i)/nd - 1.0)
            endif
         enddo
      enddo

      !     5- compute Froude number Fd and enhancement factor
      !        (empirical relation based on fig.4 of Vosper et al. 2009)

      do i=1,jyes
         aux = env(i)*ndscale(i)
         if ( aux.lt.1.e-6 )      aux = 1.e-6
         fdscale(i) = velav(i) / aux
         cdf(i) = 1.0
         if (stabfactor) then
            if ( fdscale(i).le.3.0 ) then
               cdf(i) = (6.0-cdmin)*exp(-2.8*fdscale(i)*fdscale(i)) + cdmin
            else
               cdf(i) = cdmin
            endif
         endif
      enddo
      !--------------------------------------------------------------------
      !     Compute llb tendencies:

      do l=ilev,1,-1
         do i=1,jyes
            if ( velav(i).ge.vmin .and. &
                 l.gt.izb(i)      .and. zb(i,izb(i)).ge.hmin2 ) then

               !           Vertical factor:
               fvert = sqrt( (zb(i,izb(i)) - zb(i,l)) &
                    /(0.5*env(i)   + zb(i,l)) )

               !           Implicit calculation of llb tendencies:

               ks2 = 0.5*cdblk*cdf(i)*slpf(i)*ampd(i)*fdir(i)*fvert
               vent = sqrt(u(i,l)*u(i,l) + v(i,l)*v(i,l))
               utendllb(i,l) = -ks2*u(i,l)*vent/(1.0+ks2*vent*tau)
               vtendllb(i,l) = -ks2*v(i,l)*vent/(1.0+ks2*vent*tau)

            endif
         enddo
      enddo

   endif

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !     TOTAL DRAG AND RESULTING WIND FIELD
   !--------------------------------------------------------------------
   do l=1,ilev
      do i=1,len
         utendtot(i,l) = zero
         vtendtot(i,l) = zero
      enddo
   enddo

   !     Add and scatter tendencies
   do l=1,ilev
      do i=1,len
         utendtot(i,l) = utendgwd(i,l) + &
              utendllb(i,l)
         vtendtot(i,l) = vtendgwd(i,l) + &
              vtendllb(i,l)
      enddo
   enddo

   do l=1,ilev
      do i=1,len
         ii = drag(i) + il1 - 1
         utend(ii,l) = utendtot(i,l)
         vtend(ii,l) = vtendtot(i,l)
         ttend(ii,l) = (-1./CPD) * (u(i,l)*utendtot(i,l) + v(i,l)*vtendtot(i,l))
         utendgwd4(ii,l) = real(utendgwd(i,l))
         vtendgwd4(ii,l) = real(vtendgwd(i,l))
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
600 continue
   !--------------------------------------------------------------------
   return
end subroutine sgoflx8
