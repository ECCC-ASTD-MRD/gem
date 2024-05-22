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
   subroutine gwd9(pvars, std_p_prof, tau, kount, ni, nk, nkm1, trnch)
      use tdpack, only: CAPPA, PI
      use debug_mod, only: init2nan
      use phy_options
      use phy_status, only: phy_error_L
      use phybusidx
      use phymem, only: phyvar
      use series_mod, only: series_xst
      use equivmount, only: equivmount2
      use sgo16flx,only: sgo16flx3
      use tendency, only: apply_tendencies
      use ens_perturb, only: ens_spp_get
      use vintphy, only: vint_mom2thermo, vint_thermo2mom
      use gwspectrum, only: gwspectrum6
      implicit none
!!!#include <arch_specific.hf>
      !@author J.Mailhot RPN(May1990)
      !@Object model the gravity wave drag
      !@Arguments

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: trnch, ni, nk, nkm1, kount
      real, intent(in) :: std_p_prof(nk)
      real, intent(in) :: tau

      !          - Input/Output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - Input -
      ! std_p_prof  standard pressure profil
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

      real, dimension(ni) :: rmscons
      real, dimension(ni,nkm1) :: tvirtmom, rtgm, sh, sexpk, shexpk
      real, dimension(ni,nk) :: tvirt

      integer :: j, k
      logical :: split_tend_L
      real :: lat1, lat2, val1, val2, lat

      real, pointer, dimension(:), contiguous   :: p, zdhdxdy, zdhdx, zdhdy, &
           zlhtg, zmg, zmtdir, zslope, zxcent, ztdmaskxdt, zdlat
      real, pointer, dimension(:,:), contiguous :: u, v, ztplus, s, rug, rvg, rtg, run, rvn, &
           rtn, zhuplus, zugwdtd1, zvgwdtd1, zmrk2
      real, pointer, dimension(:,:,:), contiguous :: zvcoef
      !--------------------------------------------------------------------

      if (gwdrag == 'NIL' .and. .not.non_oro) return
      
      call msg_toall(MSG_DEBUG, 'gwd [BEGIN]')
      if (timings_L) call timing_start_omp(417, 'gwd', 46)

      call init2nan(rmscons)
      call init2nan(tvirtmom, rtgm, sh, sexpk, shexpk)
      call init2nan(tvirt)
      
      MKPTR1D(p, pmoins, pvars)
      MKPTR1D(zdhdx, dhdx, pvars)
      MKPTR1D(zdhdxdy, dhdxdy, pvars)
      MKPTR1D(zdhdy, dhdy, pvars)
      MKPTR1D(zlhtg, lhtg, pvars)
      MKPTR1D(zmg, mg, pvars)
      MKPTR1D(zmtdir, mtdir, pvars)
      MKPTR1D(zslope, slope, pvars)
      MKPTR1D(ztdmaskxdt, tdmaskxdt, pvars)
      MKPTR1D(zxcent, xcent, pvars)
      MKPTR1D(zdlat, dlat, pvars)

      MKPTR2DN(rtg, tgwd, ni, nkm1, pvars)
      MKPTR2DN(rtn, tgno, ni, nkm1, pvars)
      MKPTR2DN(rug, ugwd, ni, nkm1, pvars)
      MKPTR2DN(run, ugno, ni, nkm1, pvars)
      MKPTR2DN(rvg, vgwd, ni, nkm1, pvars)
      MKPTR2DN(rvn, vgno, ni, nkm1, pvars)
      MKPTR2DN(s, sigm, ni, nkm1, pvars)
      MKPTR2DN(u, uplus, ni, nkm1, pvars)
      MKPTR2DN(v, vplus, ni, nkm1, pvars)

      MKPTR2DN(zugwdtd1, ugwd_td1, ni, nkm1, pvars)
      MKPTR2DN(zvgwdtd1, vgwd_td1, ni, nkm1, pvars)

      MKPTR2DN(zhuplus, huplus, ni, nk, pvars)
      MKPTR2DN(ztplus, tplus, ni, nk, pvars)

      MKPTR2DN(zmrk2, mrk2,  ni, nk, pvars)
      
      MKPTR3D(zvcoef, vcoef, pvars)

      split_tend_L = (associated(zugwdtd1) .and. associated(zvgwdtd1))
      if (.not.split_tend_L) then
         zugwdtd1 => rug
         zvgwdtd1 => rvg
      endif
 
      IF_INIT: if (kount == 0) then
         call equivmount2(zmg, zdhdx, zdhdy, zdhdxdy, &
              zslope, zxcent, zmtdir, ni)
      endif IF_INIT

      ! tvirt - virtual temperature on thermo/energy levels
      call mfotvt(tvirt, ztplus, zhuplus, ni, nk, ni)

      ! Auxiliary sigma and sigma-related fields:
      ! sh    - sigma on half (thermo) levels
      do k=1,nkm1-1
         do j=1,ni
            sh(j,k) = 0.5*(s(j,k)+s(j,k+1))
            sexpk(j,k) = s(j,k)**CAPPA
            shexpk(j,k) = sh(j,k)**CAPPA
         enddo
      enddo
      do j=1,ni
         sh(j,nkm1) = 0.5*(s(j,nkm1)+1.)
         sexpk(j,nkm1) = s(j,nkm1)**CAPPA
         shexpk(j,nkm1) = sh(j,nkm1)**CAPPA
      enddo
      !#TODO: use vect pow fn for sexpk, shexpk?

      IF_SGO16: if (gwdrag == 'SGO16') then
         if (timings_L) call timing_start_omp(418, 'gwdsgo', 417)
         
         call sgo16flx3(u, v, rug, rvg, rtgm, zugwdtd1, zvgwdtd1, &
              tvirt, ztplus, s, sh, sexpk, shexpk, &
              nkm1, ni,                                  &
              tau, taufac,               &
              zmg, zlhtg, zslope, zxcent, zmtdir,          &
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
         call vint_thermo2mom(tvirtmom, tvirt, zvcoef, ni, nkm1)

         ! Retrieve stochastic parameter information on request
         rmscons = ens_spp_get('rmscon', zmrk2, default=rmscon)

         ! Modulate rmscons with latitude
         IF_RMSCONLAT: if (all(rmscon_lat_weights >= 0.)) then
            lat1 = rmscon_lat_weights(1)
            lat2 = rmscon_lat_weights(2)
            val1 = rmscon_lat_weights(3)
            val2 = rmscon_lat_weights(4)
            do j=1,ni
               lat = abs(zdlat(j)*180./PI)
               if (lat >= lat2)then
                  rmscons(j) = rmscons(j) * val2
               else if (lat <= lat1)then
                  rmscons(j) = rmscons(j) * val1
               else
                  rmscons(j) = rmscons(j) * ( &
                       val2 + (lat2-lat)*(val1-val2)/(lat2-lat1))
               end if
            end do
         endif IF_RMSCONLAT
         
         call gwspectrum6(s, sh, p, tvirt, &
              tvirtmom, u, v, rvn, run, hines_flux_filter, &
              rmscons, std_p_prof, non_oro_pbot, &
              ni, nkm1)
         rtn = 0.
         if (phy_error_L) return

         if (timings_L) call timing_stop_omp(419)
      endif IF_NON_ORO

      !# Apply tendencies
      if (gwdrag /= 'NIL') then
         call apply_tendencies(u,   v,   ztplus, &
              &                rug, rvg, rtg, ztdmaskxdt, ni, nk, nkm1)
      endif
      
      if (non_oro) then
         call apply_tendencies(u,   v, &
              &                run, rvn, ztdmaskxdt, ni, nk, nkm1)
      endif
      
      if (timings_L) call timing_stop_omp(417)
      call msg_toall(MSG_DEBUG, 'gwd [END]')
      !--------------------------------------------------------------------
      return
   end subroutine gwd9

end module gwd
