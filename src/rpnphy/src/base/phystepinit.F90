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

module phystepinit
   implicit none
   private
   public :: phystepinit3

contains

   !/@*
   subroutine phystepinit3(pvars, uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, &
        dt, kount, ni, nk, trnch)
      use, intrinsic :: iso_fortran_env, only: INT64, REAL64
      use debug_mod, only: init2nan, assert_not_naninf
      use tdpack_const, only: CAPPA, GRAV, OMEGA
      use calcz0_mod, only: calcz0
      use phy_options
      use phymem, only: phyvar, phymeta, nphyvars, phymem_find, phymem_busreset, PHY_VBUSIDX, PHY_DBUSIDX
      use phybusidx
      use series_mod, only: series_xst
      use sigmalev, only: sigmalev3
      use surf_precip, only: surf_precip1, surf_precip3
      use phybudget, only: pb_compute, pb_residual
      use phy_status, only: PHY_OK
      use vintphy, only: vint_mom2thermo
      implicit none
!!!#include <arch_specific.hf>
      !@Author L. Spacek (Oct 2011)
      !@Object Setup for physics timestep
      !@Arguments
      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: kount, trnch, ni, nk
      real, dimension(ni,nk), intent(out) :: uplus0, vplus0, wplus0, tplus0, &
           huplus0, qcplus0
      real, intent(in) :: dt
      !          - Input / output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - Input -
      ! dt       timestep (parameter) (sec.)
      ! trnch    slice number
      ! ni       horizontal running length
      ! nk       vertical dimension
      !          - Output -
      ! uplus0   initial value of dbus(uplus)
      ! vplus0   initial value of dbus(vplus)
      ! tplus0   initial value of dbus(tplus)
      ! huplus0  initial value of dbus(huplus)
      ! qcplus0  initial value of dbus(qcplus)
      !*@/
#include <rmn/msg.h>
#include <rmnlib_basics.hf>
#include "phymkptr.hf"
      include "surface.cdk"
      include "phyinput.inc"

      logical,parameter:: SHORTMATCH_L = .true.
      logical,parameter:: QUIET_L = .true.

      character(len=32) :: prefix_S,basename_S,time_S,ext_S
      integer                :: i,k,ivar,nvars,istat, ivarlist(nphyvars),ind_sfc, idxv1(1)
      real                   :: rcdt1
      real, dimension(ni,nk) :: work
      real, target :: tmp1d(ni)
      real, pointer :: tmpptr(:)
      real(kind=REAL64), dimension(ni) :: en0, pw0, en1, pw1
      
      type(phymeta), pointer :: vmeta, var_m

      real, pointer, dimension(:), contiguous :: &
           zdlat, zfcor, zpmoins, ztdiag, &
           zqdiag, zza, zztsl, zzusl, zme, zp0, &
           zudiag, zvdiag, zmg, zz0, zz1, zz2, zz3, &
           zz4, ztls, ztss, zrainrate, zsnowrate, zpplus, &
           zrsc, zrlc, zrainfrac, zsnowfrac, zfrfrac, zpefrac, &
           zcone0, zconq0, zcone1, zconq1, zconedyn, zconqdyn, &
           ztdmask, ztdmaskxdt
      real, pointer, dimension(:,:), contiguous :: &
           zgzmom, zgz_moins, zhumoins, zhuplus, &
           zqadv, zqcmoins, zqcplus, zsigm, zsigt, ztadv, ztmoins, ztplus, &
           zuadv, zumoins, zuplus, zvadv, zvmoins, zvplus, zwplus, zze, &
           zgztherm, zfneige, zfip, &
           zqrp, zqrm, zqti1p, zqti1m, zqti2p, zqti2m, zqti3p, zqti3m, zqti4p, zqti4m, &
           zqnp, zqnm, zqgp, zqgm, zqhp, zqhm, zqip, zqim, zsige
      real, pointer, dimension(:,:), contiguous :: tmp1, tmp2
      real, pointer, dimension(:,:,:), contiguous :: zvcoef
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'phystepinit [BEGIN]')
      if (timings_L) call timing_start_omp(405, 'phystepinit', 46)
          
      MKPTR1D(zcone0, cone0, pvars)
      MKPTR1D(zcone1, cone1, pvars)
      MKPTR1D(zconedyn, conedyn, pvars)
      MKPTR1D(zconq0, conq0, pvars)
      MKPTR1D(zconq1, conq1, pvars)
      MKPTR1D(zconqdyn, conqdyn, pvars)
      MKPTR1D(zdlat, dlat, pvars)
      MKPTR1D(zfcor, fcor, pvars)
      MKPTR1D(zpmoins, pmoins, pvars)
      MKPTR1D(zpplus, pplus, pvars)
      MKPTR1D(zqdiag, qdiag, pvars)
      MKPTR1D(ztdiag, tdiag, pvars)
      MKPTR1D(zudiag, udiag, pvars)
      MKPTR1D(zvdiag, vdiag, pvars)
      MKPTR1D(zza, za, pvars)
      MKPTR1D(zztsl, ztsl, pvars)
      MKPTR1D(zzusl, zusl, pvars)
      MKPTR1D(zme, me_moins, pvars)
      MKPTR1D(zp0, p0_plus, pvars)

      MKPTR1D(zmg, mg, pvars)
      MKPTR1D(zz0, z0, pvars)
      MKPTR1D(zz1, z1, pvars)
      MKPTR1D(zz2, z2, pvars)
      MKPTR1D(zz3, z3, pvars)
      MKPTR1D(zz4, z4, pvars)

      MKPTR1D(zrlc, rlc, pvars)
      MKPTR1D(zrsc, rsc, pvars)
      MKPTR1D(ztls, tls, pvars)
      MKPTR1D(ztss, tss, pvars)
      MKPTR1D(zrainrate, rainrate, pvars)
      MKPTR1D(zsnowrate, snowrate, pvars)

      MKPTR1D(zrainfrac, rainfrac, pvars)
      MKPTR1D(zsnowfrac, snowfrac, pvars)
      MKPTR1D(zfrfrac, frfrac, pvars)
      MKPTR1D(zpefrac, pefrac, pvars)

      MKPTR1D(ztdmask, tdmask, pvars)
      MKPTR1D(ztdmaskxdt, tdmaskxdt, pvars)

      ind_sfc = nk
      if (any(pcptype == (/ &
           'NIL   ', &
           'BOURGE'/))) then
         MKPTR2D(zfneige, fneige, pvars)
         MKPTR2D(zfip, fip, pvars)
         if (mpdiag_for_sfc) ind_sfc = 1
      else
         MKPTR2D(zfneige, fneige3d, pvars)
         MKPTR2D(zfip, fip3d, pvars)
      endif

      MKPTR2D(zgzmom, gzmom, pvars)
      MKPTR2D(zgztherm, gztherm, pvars)
      MKPTR2D(zgz_moins, gz_moins, pvars)
      MKPTR2D(zhumoins, humoins, pvars)
      MKPTR2D(zhuplus, huplus, pvars)
      MKPTR2D(zqadv, qadv, pvars)
      MKPTR2D(zqcmoins, qcmoins, pvars)
      MKPTR2D(zqcplus, qcplus, pvars)
      MKPTR2D(zsige, sige, pvars)
      MKPTR2D(zsigm, sigm, pvars)
      MKPTR2D(zsigt, sigt, pvars)
      MKPTR2D(ztadv, tadv, pvars)
      MKPTR2D(ztmoins, tmoins, pvars)
      MKPTR2D(ztplus, tplus, pvars)
      MKPTR2D(zuadv, uadv, pvars)
      MKPTR2D(zumoins, umoins, pvars)
      MKPTR2D(zuplus, uplus, pvars)
      MKPTR2D(zvadv, vadv, pvars)
      MKPTR2D(zvmoins, vmoins, pvars)
      MKPTR2D(zvplus, vplus, pvars)
      MKPTR2D(zwplus, wplus, pvars)
      MKPTR2D(zze, ze, pvars)

      MKPTR3D(zvcoef, vcoef, pvars)
      
      MKPTR2D(zqrp, qrplus, pvars)
      MKPTR2D(zqrm, qrmoins, pvars)
      MKPTR2D(zqti1p, qti1plus, pvars)
      MKPTR2D(zqti2p, qti2plus, pvars)
      MKPTR2D(zqti3p, qti3plus, pvars)
      MKPTR2D(zqti4p, qti4plus, pvars)
      MKPTR2D(zqti1m, qti1moins, pvars)
      MKPTR2D(zqti2m, qti2moins, pvars)
      MKPTR2D(zqti3m, qti3moins, pvars)
      MKPTR2D(zqti4m, qti4moins, pvars)
      MKPTR2D(zqnp, qnplus, pvars)
      MKPTR2D(zqgp, qgplus, pvars)
      MKPTR2D(zqhp, qhplus, pvars)
      MKPTR2D(zqip, qiplus, pvars)
      MKPTR2D(zqnm, qnmoins, pvars)
      MKPTR2D(zqgm, qgmoins, pvars)
      MKPTR2D(zqhm, qhmoins, pvars)
      MKPTR2D(zqim, qimoins, pvars)

      call init2nan(work)
      call init2nan(tmp1d)
      call init2nan(en0, pw0, en1, pw1)

      if (kount == 0) then
         ztdmaskxdt(:) = ztdmask(:) * delt
      endif
      
      IF_DEBUG: if (debug_mem_L) then
         DO_IVAR: do ivar= 1, nphyvars
            vmeta => pvars(ivar)%meta
            if (any(vmeta%vname == phyinread_list_s(1:phyinread_n))) cycle
            if (vmeta%ibus == PHY_DBUSIDX) then
               !# Init diag level of dyn bus (copy down) if var not read
               !# this is needed in debug mode since the bus is init with NaN
               !# It needs to be done at kount==0 and at every restart
               if (vmeta%nk < nk) cycle
               !#TODO: adapt for 4D vars (ni,nk,fmul,nj)
               MKPTR2D(tmp1, ivar, pvars)
               if (RMN_IS_OK(assert_not_naninf(tmp1(:,nk)))) cycle
               if (any(vmeta%vname == (/"pw_gz:m", "pw_wz:p"/)).or.&
                    vmeta%vname(1:3) == "tr/") then
                  tmp1(:,nk) = 0
               else
                  tmp1(:,nk) = tmp1(:,nk-1)
               endif
            elseif (vmeta%ibus == PHY_VBUSIDX) then
               !# Reset Vol bus var to zero if var not read
               !# this is needed in debug mode since the bus is init with NaN
               pvars(ivar)%data(:) = 0.
            endif
         enddo DO_IVAR
         
         !#TODO: check that memgap is still filled with NANs on all buses
      else  !IF_DEBUG
         !#TODO: allow vol bus var recycling
         istat = phymem_busreset(pvars, PHY_VBUSIDX, F_value=0., F_resetonly_L=.true.)
      endif IF_DEBUG

      rcdt1 = 1. / dt

      do k = 1, nk-1
         do i = 1, ni
            zsigm(i,k) = zsigm(i,k) / zp0(i)
            zsigt(i,k) = zsigt(i,k) / zp0(i)
         end do
      end do
      zsigm(:,nk) = 1.
      zsigt(:,nk) = 1.

      zpmoins(:) = zp0(:)  !this is actually the time-plus surface pressure
      zpplus(:) = zp0(:)

      sigw = sigt
      istat = sigmalev3(zvcoef, zsige, zsigm, zsigt, ni, nk)
      if (.not.RMN_IS_OK(istat)) then
         call physeterror('phystepinit', 'Problem in sigmalev')
         return
      endif

      where(zhuplus(:,1:nk-1)  < 0.) zhuplus(:,1:nk-1)  = 0.
      where(zhumoins(:,1:nk-1) < 0.) zhumoins(:,1:nk-1) = 0.

      !# Precompute heights for the paramterizations
      do k = 1, nk-1
         zgzmom(:,k) = (zgz_moins(:,k) - zme) / GRAV
      enddo
      zgzmom(:,nk) = 0.
      call vint_mom2thermo(zgztherm, zgzmom, zvcoef, ni, nk)

      if (fluvert == 'RPNINT') then
         zze = zgzmom
         zze(:,nk) = 0.
      else
         zze = zgztherm
         zze(:,nk-1:nk) = 0.
      endif
      zza = zgzmom(:,nk-1)
      zzusl = zza
      zztsl = zgztherm(:,nk-1)

      !# z0 directionnel
      ! calcul de z0 avec z1,z2,z3,z4 et umoins,vmoins
      if (z0dir) call calcz0(zmg, zz0, zz1, zz2, zz3, zz4, &
           zumoins(:,nk-1), zvmoins(:,nk-1), ni)

      ! Initialize diagnostic level values in the profile
      if (any('pw_tt:p' == phyinread_list_s(1:phyinread_n))) then
         ztdiag = ztplus(:,nk)
      elseif (kount == 0) then
         ztplus(:,nk) = ztplus(:,nk-1)
         ztdiag = ztplus(:,nk)
      else
         ztplus(:,nk) = ztdiag
      endif
      if (any('tr/hu:p' == phyinread_list_s(1:phyinread_n))) then
         zqdiag = zhuplus(:,nk)
      elseif (kount  ==  0) then
         zhuplus(:,nk) = zhuplus(:,nk-1)
         zqdiag = zhuplus(:,nk)
      else
         zhuplus(:,nk) = zqdiag
      endif
      if (any('pw_uu:p' == phyinread_list_s(1:phyinread_n))) then
         zudiag = zuplus(:,nk)
      elseif (kount  ==  0) then
         zuplus(:,nk) = zuplus(:,nk-1)
         zudiag = zuplus(:,nk)
      else
         zuplus(:,nk) = zudiag
      endif
      if (any('pw_vv:p' == phyinread_list_s(1:phyinread_n))) then
         zvdiag = zvplus(:,nk)
      elseif (kount == 0) then
         zvplus(:,nk) = zvplus(:,nk-1)
         zvdiag = zvplus(:,nk)
      else
         zvplus(:,nk) = zvdiag
      endif
      if (associated(zqcplus)) zqcplus(:,nk) = 0.

      ztmoins(:,nk) = ztdiag
      zhumoins(:,nk) = zqdiag
      zumoins(:,nk) = zudiag
      zvmoins(:,nk) = zvdiag
      if (associated(zqcmoins)) zqcmoins(:,nk) = 0.

      !# Save initial time plus values
      do k = 1, nk
         do i = 1, ni
            huplus0(i,k) = zhuplus(i,k)
            uplus0(i,k)  = zuplus (i,k)
            vplus0(i,k)  = zvplus (i,k)
            tplus0(i,k)  = ztplus (i,k)
         enddo
      enddo
      if (diffuw) wplus0(:,1:nk) = zwplus(:,1:nk)
      if (associated(zqcplus)) then
         qcplus0(:,1:nk) = zqcplus(:,1:nk)
      else
         qcplus0(:,1:nk) = 0.
      endif

      ! Microphysics initializations
      if (associated(zqti1m)) zqti1m(:,nk) = 0.
      if (associated(zqti1p)) zqti1p(:,nk) = 0.
      if (associated(zqti2m)) zqti2m(:,nk) = 0.
      if (associated(zqti2p)) zqti2p(:,nk) = 0.
      if (associated(zqti3m)) zqti3m(:,nk) = 0.
      if (associated(zqti3p)) zqti3p(:,nk) = 0.
      if (associated(zqti4m)) zqti4m(:,nk) = 0.
      if (associated(zqti4p)) zqti4p(:,nk) = 0.
      if (associated(zqnp)) zqnp(:,nk) = 0.
      if (associated(zqhp)) zqhp(:,nk) = 0.
      if (associated(zqip)) zqip(:,nk) = 0.
      if (associated(zqgp)) zqgp(:,nk) = 0.
      if (associated(zqnm)) zqnm(:,nk) = 0.
      if (associated(zqhm)) zqhm(:,nk) = 0.
      if (associated(zqim)) zqim(:,nk) = 0.
      if (associated(zqgm)) zqgm(:,nk) = 0.
      if (associated(zqrp)) zqrp(:,nk) = 0.
      if (associated(zqrm)) zqrm(:,nk) = 0.

      !# calcul des tendances de la dynamique
      do k = 1,nk
         do i = 1,ni
            ztadv(i,k) = (ztplus (i,k) - ztmoins (i,k)) * rcdt1
            zqadv(i,k) = (zhuplus(i,k) - zhumoins(i,k)) * rcdt1
            zuadv(i,k) = (zuplus (i,k) - zumoins (i,k)) * rcdt1
            zvadv(i,k) = (zvplus (i,k) - zvmoins (i,k)) * rcdt1
         end do
      end do

      
      IF_CONS: if (lcons) then
         ! Compute pre-physics budget state
         if (associated(zcone0).or.associated(zconq0)) then
            tmp1d = 0.
            if (.not.associated(zcone0)) zcone0 => tmp1d
            if (.not.associated(zconq0)) zconq0 => tmp1d          
            if (pb_compute(zcone0, zconq0, en0, pw0, &
                 pvars, nk-1) /= PHY_OK) then
               call physeterror('phystepinit', &
                    'Problem computing pre-physics budget state')
               return
            endif
            zcone0(:) = real(en0(:))
            zconq0(:) = real(pw0(:))
            if (kount == 0) then
               if (associated(zcone1)) zcone1(:) = zcone0(:)
               if (associated(zconq1)) zconq1(:) = zconq0(:)
            endif
         endif

         ! Compute non-physics (dynamics) budget residual
         if (associated(zconedyn) .or. associated(zconqdyn)) then
            en1(:) = dble(zcone1(:))
            pw1(:) = dble(zconq1(:))
            if (pb_residual(zconedyn, zconqdyn, en1, pw1, pvars, &
                 delt, nk-1) /= PHY_OK) then
               call physeterror('phybusinit', &
                    'Problem computing dynamics budget')
               return
            endif
         endif
      endif IF_CONS
      
      if (associated(ztadv)) call series_xst(ztadv, 'XT', trnch)
      if (associated(zqadv)) call series_xst(zqadv, 'XQ', trnch)

      if (associated(zqcplus)) then
         do k = 1, nk
            do i = 1, ni
               work(i,k) = (zqcplus(i,k) - zqcmoins(i,k)) * rcdt1
            enddo
         enddo
         call series_xst(work, 'XL', trnch)
      endif

      !# Time swap for tracers at diag level (other levels done in dyn)
      nvars = phymem_find(ivarlist, 'tr/', 'V', 'D', QUIET_L, SHORTMATCH_L)
      do ivar = 1, nvars
         vmeta => pvars(ivarlist(ivar))%meta
         call gmmx_name_parts(vmeta%vname, prefix_S, basename_S, time_S, ext_S)
         if (all(time_S /= (/':M', ':m'/)) .and. &
              .not.any(vmeta%vname == (/'tr/hu:m', 'tr/hu:p'/)) &
              ) then
            istat = phymem_find(idxv1, &
                 trim(prefix_S)//trim(basename_S)//':M', F_npath='V', &
                 F_bpath='D', F_quiet=.true., F_shortmatch=.false.)
            if (istat > 0) then
               !#TODO: adapt for 4D vars (ni,nk,fmul,nj)
               tmp1(1:ni,1:nk) => pvars(idxv1(1))%data(:)
               tmp2(1:ni,1:nk) => pvars(ivarlist(ivar))%data(:)
               tmp1(1:ni,nk) = tmp2(1:ni,nk)
!!$            else
!!$               call msg(MSG_WARNING, '(phystepinit) Var not found: '//trim(prefix_S)//trim(basename_S)//':M')
            endif
         endif
      enddo
      
      !# Tracers clipping
      if (clip_tr_L) then
         do ivar = 1, nvars
            vmeta => pvars(ivarlist(ivar))%meta
            if (.not.any(vmeta%vname == (/'tr/hu:m', 'tr/hu:p'/))) then
               !#TODO: adapt for 4D vars (ni,nk,fmul,nj)
               var_m => pvars(ivarlist(ivar))%meta
               tmp2(1:ni,1:nk-1) => pvars(ivarlist(ivar))%data(:)
               tmp2 = min(max(var_m%vmin, tmp2), var_m%vmax)
            endif
         enddo
      endif

      do i=1,ni
         zfcor(i) = 2.*OMEGA*sin(zdlat(i))
      end do
      if (zua > 0. .and. zta > 0.) then
         do i=1,ni
            zztsl(i) = zta
            zzusl(i) = zua
         enddo
      endif

      ! Diagnostic precipitation types
      if (pcptype == 'NIL' .or. & 
           (pcptype == 'BOURGE' .and. .not.mpdiag_for_sfc)) then
         tmp1d = 0.
         nullify(tmpptr)
         call surf_precip1(ztmoins(:,nk-1), &
              tmp1d, tmp1d, &
              zrlc, ztls, zrsc, ztss, &
              tmpptr, tmpptr, tmpptr, tmpptr, &
              zrainrate, zsnowrate, ni)
      elseif (any(pcptype == (/ &
           'SPS_W19', &
           'SPS_H13', &
           'SPS_FRC'/))) then
         tmp1d = zp0*zsigt(:,nk-1)
         call surf_precip1(ztmoins(:,nk-1), &
              zhumoins(:,nk-1), tmp1d, &
              zrlc, ztls, zrsc, ztss, &
              zrainfrac, zsnowfrac, zfrfrac, zpefrac, &
              zrainrate, zsnowrate, ni)
      elseif (pcptype == 'BOURGE3D' .or. &
           (pcptype == 'BOURGE' .and. mpdiag_for_sfc)) then
         call surf_precip3(ztmoins(:,nk-1), &
              zrlc, ztls, zrsc, ztss, zrainrate, zsnowrate, &
              zfneige(:,ind_sfc), zfip(:,ind_sfc), ni)
      endif
      if (timings_L) call timing_stop_omp(405)
      call msg_toall(MSG_DEBUG, 'phystepinit [END]')
      !-------------------------------------------------------------
      return
   end subroutine phystepinit3

end module phystepinit
