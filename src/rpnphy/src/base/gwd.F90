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
   subroutine gwd9(dbus, fbus, vbus, &
        std_p_prof, tau, kount, trnch, ni, nk, nkm1)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use tdpack_const, only: CAPPA_8, GRAV, RGASD
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

      integer, intent(in) :: trnch, ni, nk, nkm1, kount
      real, dimension(:), pointer, contiguous :: dbus, fbus, vbus
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
#include <rmn/msg.h>
#include "phymkptr.hf"
#include <rmnlib_basics.hf>

      logical, parameter :: DO_GWDRAG   = .true.
      logical, parameter :: DO_BLOCKING = .true.

      real, dimension(ni) :: rland, rmscons
      real, dimension(ni,nkm1) :: work, rtgm
      real, dimension(ni,nk) :: tvirt
      real(REAL64), dimension(ni) :: pp
      real(REAL64), dimension(ni,nkm1) :: tt, te, uu, vv, sigma, s1, s2, s3, &
           utendno , vtendno, ttendno

      integer :: j, k
      logical :: split_tend_L

      real, pointer, dimension(:), contiguous   :: p, zdhdxdy, zdhdx, zdhdy, &
           zlhtg, zmg, zmtdir, zslope, zxcent, ztdmask
      real, pointer, dimension(:,:), contiguous :: u, v, ztplus, s, rug, rvg, rtg, run, rvn, &
           rtn, zhuplus, zugwdtd1, zvgwdtd1, zmrk2
      real, pointer, dimension(:,:,:), contiguous :: zvcoef
      !--------------------------------------------------------------------

      if (gwdrag == 'NIL' .and. .not.non_oro) return
      
      call msg_toall(MSG_DEBUG, 'gwd [BEGIN]')
      if (timings_L) call timing_start_omp(417, 'gwd', 46)

      call init2nan(rland, rmscons)
      call init2nan(work, rtgm)
      call init2nan(tvirt)
      call init2nan(pp)
      call init2nan(tt, te, uu, vv, sigma, s1, s2, s3)
      call init2nan(utendno , vtendno, ttendno)
      
      MKPTR1D(p, pmoins, fbus)
      MKPTR1D(zdhdx, dhdx, fbus)
      MKPTR1D(zdhdxdy, dhdxdy, fbus)
      MKPTR1D(zdhdy, dhdy, fbus)
      MKPTR1D(zlhtg, lhtg, fbus)
      MKPTR1D(zmg, mg, fbus)
      MKPTR1D(zmtdir, mtdir, fbus)
      MKPTR1D(zslope, slope, fbus)
      MKPTR1D(ztdmask, tdmask, fbus)
      MKPTR1D(zxcent, xcent, fbus)

      MKPTR2Dm1(rtg, tgwd, vbus)
      MKPTR2Dm1(rtn, tgno, vbus)
      MKPTR2Dm1(rug, ugwd, vbus)
      MKPTR2Dm1(run, ugno, vbus)
      MKPTR2Dm1(rvg, vgwd, vbus)
      MKPTR2Dm1(rvn, vgno, vbus)
      MKPTR2Dm1(s, sigm, dbus)
      MKPTR2Dm1(u, uplus, dbus)
      MKPTR2Dm1(v, vplus, dbus)

      MKPTR2Dm1(zugwdtd1, ugwd_td1, vbus)
      MKPTR2Dm1(zvgwdtd1, vgwd_td1, vbus)

      MKPTR2D(zhuplus, huplus, dbus)
      MKPTR2D(ztplus, tplus, dbus)

      MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, fbus)
      
      MKPTR3D(zvcoef, vcoef, 2, vbus)

      split_tend_L = (associated(zugwdtd1) .and. associated(zvgwdtd1))
      if (.not.split_tend_L) then
         zugwdtd1 => rug
         zvgwdtd1 => rvg
      endif
 
      IF_INIT: if (kount == 0) then
         rland = -abs(nint(zmg))
         call equivmount1(rland, zdhdx, zdhdy, zdhdxdy, &
              ni, 1, ni, zslope, zxcent, zmtdir)
      endif IF_INIT

      ! tvirt - virtual temperature on thermo/energy levels
      call mfotvt(tvirt, ztplus, zhuplus, ni, nk, ni)

      ! Auxiliary sigma and sigma-related fields:
      ! s1    - sigma on half (thermo) levels
      ! s2    - sigma**cappa on full (mom) levels
      ! s3    - sigma**cappa on half (thermo) levels

      do k=1,nkm1-1
         do j=1,ni
            s1(j,k) = 0.5*(s(j,k)+s(j,k+1))
            s3(j,k) = exp(CAPPA_8*log(s1(j,k)))
            s2(j,k) = exp(CAPPA_8*log(dble(s(j,k))))
         enddo
      enddo

      do j=1,ni
         s1(j,nkm1) = 0.5*(s(j,nkm1)+1.)
         s3(j,nkm1) = exp(CAPPA_8*log(s1(j,nkm1)))
         s2(j,nkm1) = exp(CAPPA_8*log(dble(s(j,nkm1))))
      enddo

      IF_SGO16: if (gwdrag == 'SGO16') then
         if (timings_L) call timing_start_omp(418, 'gwdsgo', 417)

         work(:,:) = real(s1(:,:))
         rland = - abs(nint(zmg))  !# land-water index: -1=land, 0=water
         call sgo16flx3(u, v, rug, rvg, rtgm, zugwdtd1, zvgwdtd1, &
              tvirt, ztplus, s, work, &
              nkm1, ni, 1, ni,                                   &
              tau, taufac,               &
              rland, zlhtg, zslope, zxcent, zmtdir,          &
              sgo_cdmin, sgo_bhfac, sgo_phic,                &
              sgo_stabfac, sgo_nldirfac, zmrk2,              &
              DO_GWDRAG, DO_BLOCKING, split_tend_L)
         if (phy_error_L) return

         call vint_mom2thermo(rtg, rtgm, zvcoef, ni, nkm1)

         call series_xst(rug, 'GU', trnch)
         call series_xst(rvg, 'GV', trnch)

         if (timings_L) call timing_stop_omp(418)
      endif IF_SGO16

      IF_NON_ORO: if (non_oro) then
         if (timings_L) call timing_start_omp(419, 'gwdnon', 417)

         ! tt   - virtual temperature on momentum levels
         call vint_thermo2mom(work, tvirt, zvcoef, ni, nkm1)
         tt    = work

         !# Copy to Double
         !# TO CHECK: is 64bit really needed for non_oro
         te(1:ni,1:nkm1) = tvirt(1:ni,1:nkm1)
         uu    = u
         vv    = v
         pp    = p
         sigma = s

         ! Retrieve stochastic parameter information on request
         rmscons(:) = ens_spp_get('rmscon', zmrk2, default=rmscon)
         
         call gwspectrum6(sigma, s1, s2, s3, pp, te, &
              tt, uu, vv, ttendno, vtendno, utendno, hines_flux_filter, &
              GRAV, RGASD, rmscons, iheatcal, std_p_prof, non_oro_pbot, &
              ni, nkm1)
         if (phy_error_L) return

         run = real(utendno)
         rvn = real(vtendno)
         rtn = real(ttendno)

         if (timings_L) call timing_stop_omp(419)
      endif IF_NON_ORO

      !# Apply tendencies
      if (gwdrag /= 'NIL') then
         call apply_tendencies(u, rug, ztdmask, ni, nk, nkm1)
         call apply_tendencies(v, rvg, ztdmask, ni, nk, nkm1)
         call apply_tendencies(ztplus, rtg, ztdmask, ni, nk, nkm1)
      endif
      
      if (non_oro) then
         call apply_tendencies(u, run, ztdmask, ni, nk, nkm1)
         call apply_tendencies(v, rvn, ztdmask, ni, nk, nkm1)
      endif
      
      if (timings_L) call timing_stop_omp(417)
      call msg_toall(MSG_DEBUG, 'gwd [END]')
      !--------------------------------------------------------------------
      return
   end subroutine gwd9

end module gwd
