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
   subroutine phystepinit3(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, &
        vbus, dbus, fbus, seloc, dt, &
        vsiz, dsiz, fsiz, kount, trnch, ni, nk)
      use, intrinsic :: iso_fortran_env, only: INT64
      use debug_mod, only: init2nan, assert_not_naninf
      use tdpack_const, only: CAPPA, GRAV, OMEGA
      use calcz0_mod, only: calcz0
      use phy_getmeta_mod, only: phy_getmeta
      use phy_options
      use phy_typedef, only: phymeta
      use phybus
      use phygetmetaplus_mod, only: phymetaplus, phygetmetaplus
      use series_mod, only: series_xst
      use sigmalev, only: sigmalev3
      use surf_precip, only: surf_precip1, surf_precip3
      implicit none
!!!#include <arch_specific.hf>
      !@Author L. Spacek (Oct 2011)
      !@Object Setup for physics timestep
      !@Arguments
      integer                :: vsiz, dsiz, fsiz, kount, trnch, ni, nk
      real, dimension(ni,nk) :: uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0
      real, dimension(ni,nk) :: seloc
      real, target           :: vbus(vsiz), dbus(dsiz), fbus(fsiz)
      real                   :: dt
      !          - Input -
      ! dsiz     dimension of dbus
      ! vsiz     dimension of vbus
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
      ! dt       timestep (sec.)
      !          - input/output -
      ! dbus     dynamics input field
      ! vbus     physics tendencies and other output fields from the physics
      !*@/
#include <msg.h>
#include <rmnlib_basics.hf>
#include "phymkptr.hf"
      include "surface.cdk"
      include "phyinput.inc"

      logical,parameter:: SHORTMATCH_L = .true.
      integer,parameter:: MYMAX = 2048

      character(len=1)  :: bus_S
      character(len=32) :: prefix_S,basename_S,time_S,ext_S
      integer                :: i,k,ivar,nvars,istat
      real                   :: sc,rcdt1
      real, dimension(ni,nk) :: work
      real, dimension(ni,nk) :: qe

      type(phymetaplus) :: meta_m, meta_p
      type(phymeta), pointer :: metalist(:)

      real, pointer, dimension(:)   :: zdlat, zfcor, zpmoins, ztdiag, zthetaa, &
           zqdiag, zza, zztsl, zzusl, zme, zp0, &
           zudiag, zvdiag, zmg, zz0, zz1, zz2, zz3, &
           zz4, ztls, ztss, zrainrate, zsnowrate, zthetaap, zpplus, &
           zrsc, zrlc
      real, pointer, dimension(:,:) :: zgzmom, zgz_moins, zhumoins, zhuplus, &
           zqadv, zqcmoins, zqcplus, zsigm, zsigt, ztadv, ztmoins, ztplus, &
           zuadv, zumoins, zuplus, zvadv, zvmoins, zvplus, zwplus, zze, &
           zgztherm, ztve, zfneige, zfip
      real, pointer, dimension(:,:) :: tmp1,tmp2
      real, pointer, dimension(:,:,:) :: zvcoef
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'phystepinit [BEGIN]')
      if (timings_L) call timing_start_omp(405, 'phystepinit', 46)

      MKPTR1D(zdlat, dlat, fbus)
      MKPTR1D(zfcor, fcor, vbus)
      MKPTR1D(zpmoins, pmoins, fbus)
      MKPTR1D(zpplus, pplus, vbus)
      MKPTR1D(zqdiag, qdiag, fbus)
      MKPTR1D(ztdiag, tdiag, fbus)
      MKPTR1D(zudiag, udiag, fbus)
      MKPTR1D(zvdiag, vdiag, fbus)
      MKPTR1D(zthetaa, thetaa, vbus)
      MKPTR1D(zthetaap, thetaap, vbus)
      MKPTR1D(zza, za, vbus)
      MKPTR1D(zztsl, ztsl, vbus)
      MKPTR1D(zzusl, zusl, vbus)
      MKPTR1D(zme, me_moins, dbus)
      MKPTR1D(zp0, p0_plus, dbus)

      MKPTR1D(zmg, mg, fbus)
      MKPTR1D(zz0, z0, fbus)
      MKPTR1D(zz1, z1, fbus)
      MKPTR1D(zz2, z2, fbus)
      MKPTR1D(zz3, z3, fbus)
      MKPTR1D(zz4, z4, fbus)

      MKPTR1D(zrlc, rlc, fbus)
      MKPTR1D(zrsc, rsc, fbus)
      MKPTR1D(ztls, tls, fbus)
      MKPTR1D(ztss, tss, fbus)
      MKPTR1D(zrainrate, rainrate, vbus)
      MKPTR1D(zsnowrate, snowrate, vbus)

      if (any(pcptype == (/ &
           'NIL   ', &
           'BOURGE'/))) then
         MKPTR2DN(zfneige, fneige, ni, 1, fbus)
         MKPTR2DN(zfip, fip, ni, 1, fbus)
      else
         MKPTR2D(zfneige, fneige3d, fbus)
         MKPTR2D(zfip, fip3d, fbus)
      endif

      MKPTR2D(zgzmom, gzmom, vbus)
      MKPTR2D(zgztherm, gztherm, vbus)
      MKPTR2D(zgz_moins, gz_moins, dbus)
      MKPTR2D(zhumoins, humoins, dbus)
      MKPTR2D(zhuplus, huplus, dbus)
      MKPTR2D(zqadv, qadv, vbus)
      MKPTR2D(zqcmoins, qcmoins, dbus)
      MKPTR2D(zqcplus, qcplus, dbus)
      MKPTR2D(zsigm, sigm, dbus)
      MKPTR2D(zsigt, sigt, dbus)
      MKPTR2D(ztadv, tadv, vbus)
      MKPTR2D(ztmoins, tmoins, dbus)
      MKPTR2D(ztplus, tplus, dbus)
      MKPTR2D(ztve, tve, vbus)
      MKPTR2D(zuadv, uadv, vbus)
      MKPTR2D(zumoins, umoins, dbus)
      MKPTR2D(zuplus, uplus, dbus)
      MKPTR2D(zvadv, vadv, vbus)
      MKPTR2D(zvmoins, vmoins, dbus)
      MKPTR2D(zvplus, vplus, dbus)
      MKPTR2D(zwplus, wplus, dbus)
      MKPTR2D(zze, ze, vbus)

      MKPTR3D(zvcoef, vcoef, 2, vbus)

      call init2nan(work, qe)

      IF_DEBUG: if (debug_mem_L) then

         !# Init diag level of dyn bus (copy down) if var not read
         !# this is needed in debug mode since the bus is init with NaN
         !# It needs to be done at kount==0 and at every restart
         bus_S = 'D'
         nullify(metalist)
         nvars = phy_getmeta(metalist, '', F_npath='V', F_bpath=bus_S, &
              F_maxmeta=MYMAX, F_shortmatch=SHORTMATCH_L)
         do ivar = 1, nvars
            istat = phygetmetaplus(meta_m, metalist(ivar)%vname, F_npath='V', &
                 F_bpath=bus_S, F_quiet=.true., F_shortmatch=.false.)
            if (.not.any(meta_m%meta%vname == phyinread_list_s(1:phyinread_n))) then

               if (meta_m%meta%nk >= nk) then
                  !#TODO: adapt for 4D vars (ni,nk,fmul,nj)
                  tmp1(1:ni,1:nk) => meta_m%vptr(:,trnch)
                  if (.not.RMN_IS_OK(assert_not_naninf(tmp1(:,nk)))) then
                     if (any(metalist(ivar)%vname == (/"pw_gz:m", "pw_wz:p"/)).or.&
                          metalist(ivar)%vname(1:3) == "tr/") then
                        tmp1(:,nk) = 0
                     else
                        tmp1(:,nk) = tmp1(:,nk-1)
                     endif
                  endif
               endif
            endif
         enddo
         deallocate(metalist, stat=istat)

         !# Reset Vol bus var to zero if var not read
         !# this is needed in debug mode since the bus is init with NaN
         bus_S = 'V'
         nullify(metalist)
         nvars = phy_getmeta(metalist, '', F_npath='V', F_bpath=bus_S, &
              F_maxmeta=MYMAX, F_shortmatch=SHORTMATCH_L)
         do ivar = 1, nvars
            istat = phygetmetaplus(meta_m, metalist(ivar)%vname, F_npath='V', &
                 F_bpath=bus_S, F_quiet=.true., F_shortmatch=.false.)
            if (.not.any(meta_m%meta%vname == phyinread_list_s(1:phyinread_n))) &
                 meta_m%vptr(:,trnch) = 0.
         enddo
         deallocate(metalist, stat=istat)

         !#TODO: check that memgap is still filled with NANs on all buses
      else
         !#TODO: allow vol bus var recycling
         vbus = 0.
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
      istat = sigmalev3(zvcoef, seloc, zsigm, zsigt, ni, nk)
      if (.not.RMN_IS_OK(istat)) then
         call physeterror('phystepinit', 'Problem in sigmalev')
         return
      endif

      where(zhuplus(:,1:nk-1)  < 0.) zhuplus  = 0.
      where(zhumoins(:,1:nk-1) < 0.) zhumoins = 0.

      !# Precompute heights for the paramterizations
      do k = 1, nk-1
         zgzmom(:,k) = (zgz_moins(:,k) - zme) / GRAV
      enddo
      zgzmom(:,nk) = 0.
      call vint_mom2thermo(zgztherm, zgzmom, zvcoef, ni, nk)

      zze = zgztherm
      zze(:,nk-1:nk) = 0.
      zza = zgzmom(:,nk-1)
      zzusl = zza
      zztsl = zgztherm(:,nk-1)

      !# z0 directionnel
      ! calcul de z0 avec z1,z2,z3,z4 et umoins,vmoins
      if (z0dir) call calcz0(zmg, zz0, zz1, zz2, zz3, zz4, &
           zumoins(:,nk-1), zvmoins(:,nk-1), ni)

      ! calcul du facteur de coriolis (fcor), de la hauteur du
      ! dernier niveau actif (za) et de la temperature potentielle
      ! a ce niveau (thetaa)
      ztve(:,1:nk-1) = ztmoins(:,1:nk-1)
      qe(:,1:nk-1)   = zhumoins(:,1:nk-1)

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
      zqcplus(:,nk) = 0.

      ztmoins(:,nk) = ztdiag
      zhumoins(:,nk) = zqdiag
      zumoins(:,nk) = zudiag
      zvmoins(:,nk) = zvdiag
      zqcmoins(:,nk) = 0.

      !# Save initial time plus values
      do k = 1, nk
         do i = 1, ni
            huplus0(i,k) = zhuplus(i,k)
            qcplus0(i,k) = zqcplus(i,k)
            uplus0 (i,k) = zuplus (i,k)
            vplus0 (i,k) = zvplus (i,k)
            tplus0 (i,k) = ztplus (i,k)
         enddo
      enddo
      if (diffuw) wplus0(:,1:nk) = zwplus(:,1:nk)

      !# calcul des tendances de la dynamique
      do k = 1,nk
         do i = 1,ni
            ztadv(i,k) = (ztplus (i,k) - ztmoins (i,k)) * rcdt1
            zqadv(i,k) = (zhuplus(i,k) - zhumoins(i,k)) * rcdt1
            zuadv(i,k) = (zuplus (i,k) - zumoins (i,k)) * rcdt1
            zvadv(i,k) = (zvplus (i,k) - zvmoins (i,k)) * rcdt1
         end do
      end do

      if (associated(ztadv)) call series_xst(ztadv, 'XT', trnch)
      if (associated(zqadv)) call series_xst(zqadv, 'XQ', trnch)

      if (stcond /= 'NIL') then
         do k = 1, nk
            do i = 1, ni
               work(i,k) = (zqcplus(i,k) - zqcmoins(i,k)) * rcdt1
            enddo
         enddo
         call series_xst(work, 'XL', trnch)
      endif

      nullify(metalist)
      nvars = phy_getmeta(metalist, 'tr/', F_npath='V', F_bpath='D', &
           F_maxmeta=MYMAX, F_shortmatch=SHORTMATCH_L)
      do ivar = 1, nvars
         call gmmx_name_parts(metalist(ivar)%vname, prefix_S, basename_S, &
              time_S, ext_S)
         if (all(time_S /= (/':M', ':m'/)) .and. &
              .not.any(metalist(ivar)%vname == (/'tr/hu:m', 'tr/hu:p'/)) &
              ) then
            istat = phygetmetaplus(meta_m, &
                 trim(prefix_S)//trim(basename_S)//':M', F_npath='V', &
                 F_bpath='D', F_quiet=.true., F_shortmatch=.false.)
            istat = min(phygetmetaplus(meta_p, &
                 trim(prefix_S)//trim(basename_S)//':P', F_npath='V', &
                 F_bpath='D', F_quiet=.true., F_shortmatch=.false.),istat)
            if (RMN_IS_OK(istat)) then
               !#TODO: adapt for 4D vars (ni,nk,fmul,nj)
               tmp1(1:ni,1:nk) => meta_m%vptr(:,trnch)
               tmp2(1:ni,1:nk) => meta_p%vptr(:,trnch)
               tmp1(1:ni,nk) = tmp2(1:ni,nk)
            endif
         endif
      enddo
      if (clip_tr_L) then
         do ivar = 1, nvars
            if (.not.any(metalist(ivar)%vname == (/'tr/hu:m', 'tr/hu:p'/))) then
               istat = phygetmetaplus(meta_m, metalist(ivar)%vname, F_npath='V', &
                    F_bpath='D', F_quiet=.true., F_shortmatch=.false.)
               tmp2(1:ni,1:nk-1) => meta_m%vptr(:,trnch)
               tmp2 = min(max(meta_m%meta%vmin, tmp2), meta_m%meta%vmax)
            endif
         enddo
      endif
      deallocate(metalist, stat=istat)

      call mfotvt(ztve, ztve, qe, ni, nk-1, ni)

      do i=1,ni
         sc = zsigt(i,nk-1)**(-CAPPA)
         zthetaa(i) = sc*ztmoins(i,nk-1)
         zthetaap(i) = sc*ztplus(i,nk-1)
         zfcor  (i)= 2.*OMEGA*sin(zdlat(i))
      end do
      if (zua > 0. .and. zta > 0.) then
         do i=1,ni
            zztsl(i) = zta
            zzusl(i) = zua
         enddo
      endif

      if (any(pcptype == (/ &
           'NIL   ', &
           'BOURGE'/))) then
         call surf_precip1(ztmoins(:,nk-1), &
              zrlc, ztls, zrsc, ztss, zrainrate, zsnowrate, ni)
      elseif (pcptype == 'BOURGE3D') then
         call surf_precip3(ztmoins(:,nk-1), &
              zrlc, ztls, zrsc, ztss, zrainrate, zsnowrate, &
              zfneige(:,nk), zfip(:,nk), ni)
      endif
      if (timings_L) call timing_stop_omp(405)
      call msg_toall(MSG_DEBUG, 'phystepinit [END]')
      !-------------------------------------------------------------
      return
   end subroutine phystepinit3

end module phystepinit
