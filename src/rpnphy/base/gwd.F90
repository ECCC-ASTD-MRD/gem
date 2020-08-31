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

module gwd
   implicit none
   private
   public :: gwd9

contains

   !/@*
   subroutine gwd9(d, f, vb, sized, sizef, sizev, &
        std_p_prof, tau, kount, trnch, ni, nk, nkm1)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use tdpack_const, only: CAPPA, GRAV, RGASD
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use series_mod, only: series_xst
      use tendency, only: apply_tendencies
      use ens_perturb, only: ens_nc2d, ens_spp_get
      implicit none
!!!#include <arch_specific.hf>
      !@author J.Mailhot RPN(May1990)
      !@Object model the gravity wave drag
      !@Arguments

      integer, intent(in) :: sized, sizef, sizev, trnch, ni, nk, nkm1, kount
      real, target, intent(inout) :: d(sized), f(sizef), vb(sizev)
      real, intent(in) :: std_p_prof(nk)
      real, intent(in) :: tau

      !          - Input/Output -
      ! f        field of permanent physics variables
      ! sizef    dimension of F
      ! u        u component of wind as input
      !          u component of wind modified by the gravity wave
      !          drag as output
      ! v        v component of wind as input
      !          v component of wind modified by the gravity wave
      !          drag as output
      !          - Input -
      ! t        virtual temperature
      ! s        local sigma levels
      ! std_p_prof  standard pressure profil
      !          - Output -
      ! rug      gravity wave drag tendency on the u component of real
      !          wind
      ! rvg      gravity wave drag tendency on the v component of real
      !          wind
      ! run      non-oro. gravity wave drag tendency on the u component of real
      !          wind
      ! rvn      non-oro. gravity wave drag tendency on the v component of real
      !          wind
      ! rtn      non-oro. gravity wave drag tendency on temperature
      !          - Input -
      ! tau      timestep times a factor: 1 for two time-level models
      !                                   2 for three time-level models
      ! trnch    index of the vertical plane(ni*nk) for which the
      !          calculations are done.
      ! ni       horizontal dimension
      ! nk       vertical dimension
      ! nkm1     vertical scope of the operator

      !@notes
      !          this routine needs at least:
      !          ( 12*nk + 12 )*n + 3*nk words in dynamic allocation
      !            +3*nk       -"-         for local sigma
      !            +2*nk       -"-         for gather on s, sh
      !                           - 3*nk   s1,s2,s3 change from 1d to 2d
      !*@/
#include <msg.h>
#include "phymkptr.hf"
#include <rmnlib_basics.hf>

      logical, parameter :: DO_GWDRAG   = .true.
      logical, parameter :: DO_BLOCKING = .true.

      real, dimension(ni) :: rland, rmscons
      real, dimension(ni,nkm1) :: work, rte, rs1, rtgm, ttendgw_r4
      real, dimension(ni,nk) :: tvirt
      real(REAL64), dimension(ni) :: fcorio, land, launch, mtdir8, pp, slope8, &
           sxx8, syy8, sxy8, xcent8
      real(REAL64), dimension(ni,nkm1) :: tt, te, uu, vv, sigma, s1, s2, s3, &
           utendgw, vtendgw, ttendgw, utendno , vtendno, ttendno

      integer :: j, k
      logical :: split_tend_L

      real, pointer, dimension(:)   :: p, zdhdxdy, zdhdx, zdhdy, zfcor, &
           zlhtg, zmg, zmtdir, zslope, zxcent, ztdmask
      real, pointer, dimension(:,:) :: u, v, ztplus, s, rug, rvg, rtg, run, rvn, &
           rtn, zhuplus, zugwdtd1, zvgwdtd1, zmrk2
      real, pointer, dimension(:,:,:) :: zvcoef
      !--------------------------------------------------------------------

      if (gwdrag == 'NIL' .and. .not.non_oro) return
      call msg_toall(MSG_DEBUG, 'gwd [BEGIN]')
      if (timings_L) call timing_start_omp(417, 'gwd', 46)

      MKPTR1D(p, pmoins, f)
      MKPTR1D(zdhdx, dhdx, f)
      MKPTR1D(zdhdxdy, dhdxdy, f)
      MKPTR1D(zdhdy, dhdy, f)
      MKPTR1D(zfcor, fcor, vb)
      MKPTR1D(zlhtg, lhtg, f)
      MKPTR1D(zmg, mg, f)
      MKPTR1D(zmtdir, mtdir, f)
      MKPTR1D(zslope, slope, f)
      MKPTR1D(ztdmask, tdmask, f)
      MKPTR1D(zxcent, xcent, f)

      MKPTR2Dm1(rtg, tgwd, vb)
      MKPTR2Dm1(rtn, tgno, vb)
      MKPTR2Dm1(rug, ugwd, vb)
      MKPTR2Dm1(run, ugno, vb)
      MKPTR2Dm1(rvg, vgwd, vb)
      MKPTR2Dm1(rvn, vgno, vb)
      MKPTR2Dm1(s, sigm, d)
      MKPTR2Dm1(u, uplus, d)
      MKPTR2Dm1(v, vplus, d)

      MKPTR2Dm1(zugwdtd1, ugwd_td1, vb)
      MKPTR2Dm1(zvgwdtd1, vgwd_td1, vb)

      MKPTR2D(zhuplus, huplus, d)
      MKPTR2D(ztplus, tplus, d)

      MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, f)
      
      MKPTR3D(zvcoef, vcoef, 2, vb)

      split_tend_L = (associated(zugwdtd1) .and. associated(zvgwdtd1))
      if (.not.split_tend_L) then
         zugwdtd1 => rug
         zvgwdtd1 => rvg
      endif

      ! Retrieve stochastic parameter information on request
      rmscons(:) = ens_spp_get('rmscon', zmrk2, default=rmscon)

      call init2nan(rland)
      call init2nan(work, rte, rs1, rtgm, tvirt, ttendgw_r4)
      call init2nan(fcorio, land, launch, mtdir8, pp, slope8, sxx8, syy8, sxy8, xcent8)
      call init2nan(tt, te, uu, vv, sigma, s1, s2, s3, utendgw, vtendgw)
      call init2nan(ttendgw, utendno , vtendno, ttendno)

      IF_INIT: if (kount == 0) then
         rland = -abs(nint(zmg))
         call equivmount1(rland, zdhdx, zdhdy, zdhdxdy, &
              ni, 1, ni, zslope, zxcent, zmtdir)
      endif IF_INIT

      ! tt   - virtual temperature on momentum levels

      call mfotvt(tvirt, ztplus, zhuplus, ni, nk, ni)
      call vint_thermo2mom(work, tvirt, zvcoef, ni, nkm1)
      tt    = work
      uu    = u
      vv    = v
      pp    = p
      sigma = s

      ! te   - virtual temperature on thermo/energy levels

      work(1:ni,1:nkm1) = tvirt(1:ni,1:nkm1)
      te = work

      ! Auxiliary sigma and sigma-related fields:
      ! s1    - sigma on half (thermo) levels
      ! s2    - sigma**cappa on full (mom) levels
      ! s3    - sigma**cappa on half (thermo) levels

      do k=1,nkm1-1
         do j=1,ni
            s1(j,k) = 0.5*(s(j,k)+s(j,k+1))
            s2(j,k) = dble(s(j,k))
         enddo
         call vlog(s3(1,k), s1(1,k), ni)
         call vlog(s2(1,k), s2(1,k), ni)
         do j=1,ni
            s3(j,k) = CAPPA*s3(j,k)
            s2(j,k) = CAPPA*s2(j,k)
         enddo
         call vexp(s3(1,k), s3(1,k), ni)
         call vexp(s2(1,k), s2(1,k), ni)
      enddo

      do j=1,ni
         s1(j,nkm1) = 0.5*(s(j,nkm1)+1.)
         s2(j,nkm1) = dble(s(j,nkm1))
      enddo
      call vlog(s3(1,nkm1), s1(1,nkm1), ni)
      call vlog(s2(1,nkm1), s2(1,nkm1), ni)
      do j=1,ni
         s3(j,nkm1) = CAPPA*s3(j,nkm1)
         s2(j,nkm1) = CAPPA*s2(j,nkm1)
      enddo
      call vexp(s3(1,nkm1), s3(1,nkm1), ni)
      call vexp(s2(1,nkm1), s2(1,nkm1), ni)

      IF_GWD86: if (gwdrag == 'GWD86') then

         launch  = zlhtg
         sxx8    = zdhdx
         syy8    = zdhdy
         sxy8    = zdhdxdy
         mtdir8  = zmtdir
         slope8  = zslope
         xcent8  = zxcent
         fcorio  = zfcor
         land = - abs(nint(zmg))  !# land-water index: -1=land, 0=water

         call sgoflx8(uu, vv, utendgw, vtendgw, ttendgw   , &
              zugwdtd1, zvgwdtd1, te, tt, sigma, s1       , &
              nkm1, ni, 1, ni                             , &
              tau, taufac                                 , &
              land, launch, slope8, xcent8, mtdir8        , &
              DO_GWDRAG, DO_BLOCKING,           &
              sgo_stabfac, sgo_nldirfac, sgo_cdmin, split_tend_L)

         rug = real(utendgw)
         rvg = real(vtendgw)
         ttendgw_r4 = real(ttendgw)
         call vint_mom2thermo(rtg, ttendgw_r4, zvcoef, ni, nkm1)

         call apply_tendencies(u, rug, ztdmask, ni, nk, nkm1)
         call apply_tendencies(v, rvg, ztdmask, ni, nk, nkm1)
         call apply_tendencies(ztplus, rtg, ztdmask, ni, nk, nkm1)

         call series_xst(rug, 'GU', trnch)
         call series_xst(rvg, 'GV', trnch)

      endif IF_GWD86

      IF_SGO16: if (gwdrag == 'SGO16') then

         rte(:,:) = real(te(:,:))
         rs1(:,:) = real(s1(:,:))
         rland = - abs(nint(zmg))  !# land-water index: -1=land, 0=water
         call sgo16flx3(u, v, rug, rvg, rtgm, zugwdtd1, zvgwdtd1, &
              rte, ztplus, s, rs1, &
              nkm1, ni, 1, ni,                                   &
              tau, taufac,               &
              rland, zlhtg, zslope, zxcent, zmtdir,          &
              sgo_cdmin, sgo_bhfac, sgo_phic,                &
              sgo_stabfac, sgo_nldirfac, zmrk2,              &
              DO_GWDRAG, DO_BLOCKING, split_tend_L)

         call vint_mom2thermo(rtg, rtgm, zvcoef, ni, nkm1)

         call apply_tendencies(u, rug, ztdmask, ni, nk, nkm1)
         call apply_tendencies(v, rvg, ztdmask, ni, nk, nkm1)
         call apply_tendencies(ztplus, rtg, ztdmask, ni, nk, nkm1)

         call series_xst(rug, 'GU', trnch)
         call series_xst(rvg, 'GV', trnch)

      endif IF_SGO16

      IF_NON_ORO: if (non_oro) then

         call gwspectrum6(sigma, s1, s2, s3, pp, te, &
              tt, uu, vv, ttendno, vtendno, utendno, hines_flux_filter, &
              GRAV, RGASD, rmscons, iheatcal, std_p_prof, non_oro_pbot, &
              ni, nkm1)
         if (phy_error_L) return

         run = real(utendno)
         rvn = real(vtendno)
         rtn = real(ttendno)

         call apply_tendencies(u, run, ztdmask, ni, nk, nkm1)
         call apply_tendencies(v, rvn, ztdmask, ni, nk, nkm1)

      endif IF_NON_ORO

      if (timings_L) call timing_stop_omp(417)
      call msg_toall(MSG_DEBUG, 'gwd [END]')
      !--------------------------------------------------------------------
      return
   end subroutine gwd9

end module gwd
