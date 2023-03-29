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

module turbul
   implicit none
   private
   public :: turbul2

contains

   !/@*
   subroutine turbul2(dbus, fbus, vbus, se, kount, trnch, ni, nk, nkm1)
      use, intrinsic :: iso_fortran_env, only: INT64
      use debug_mod, only: init2nan
      use tdpack_const, only: CPD, DELTA, GRAV, KARMAN
      use series_mod, only: series_xst
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use ens_perturb, only: ens_nc2d
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Object interface for turbulent kinetic energy calculations
      !@Arguments
      !          - Input/Output -
      ! dbus     dynamic             bus
      ! fbus     permanent variables bus
      ! vbus     volatile (output)   bus
      !          - Input -
      ! se       sigma level for turbulent energy
      ! kount    index of timestep
      ! trnch    number of the slice
      ! ni       horizontal dimension
      ! nk       vertical dimension
      ! nkm1     vertical scope of the operator

      integer, intent(in) :: kount, trnch, ni, nk, nkm1
      real, pointer, contiguous :: dbus(:), fbus(:), vbus(:)
      real, intent(in) :: se(ni,nk)

      !@Author J. Mailhot and B. Bilodeau (Jan 1999)
      !*@/
#include "phymkptr.hf"
      include "clefcon.cdk"
      include "surface.cdk"
      include "tables.cdk"
      include "phyinput.inc"

      integer, external :: neark !#TODO make it a module

      integer, parameter :: ICPU = 1 !#TODO: remove this

      real, save :: tfilt = 0.1

      integer :: i, k, stat, ksl(ni)
      real :: cf1, cf2, eturbtau, uet, ilmot, &
           fhz, fim, fit, hst_local, &
           b(ni,nkm1*4), xb(ni), xh(ni), fbsurf(ni)

      real, pointer, dimension(:) :: zbt_ag, zfrv_ag, zftemp_ag, &
           zfvap_ag, zhst_ag, zilmo_ag, zz0_ag, ztsurf  !#TODO: should be contiguous
      real, pointer, dimension(:), contiguous   :: zalfat, zztsl, zh, zlh, &
           zwstar, zudiag, zvdiag, zhpar, zpmoins, zsdtswd, &
           zsdtsws, zwge, zwgmax, zwgmin, zdxdy, ztstar

      real, pointer, dimension(:,:), contiguous :: tke, zenmoins, zkm, zkt, ztmoins, &
           ztve, zze, zzn, zwtng, zwqng, zuwng, zvwng, zqcmoins, zhumoins, &
           zsigm, zsigt, zsigw, zumoins, zvmoins, zfbl, zfblgauss, zfblnonloc, zfnn, &
           zftot, zlwc, zqtbl, zturbreg, zzd, zgq, zgql, zgte, zgzmom, zrif, &
           zrig, zshear2, zvcoef, zc1pbl, zbuoy, zmrk2

      real, dimension(ni,nkm1) :: c, x, wk2d, enold, tmom, qmom, te, qce, qe
      !----------------------------------------------------------------

      call init2nan(xb, xh, fbsurf)
      call init2nan(b, c, x, wk2d, enold, tmom, qmom, te, qce, qe)

      MKPTR1D(zalfat, alfat, vbus)
      MKPTR1DK(zbt_ag, bt, indx_agrege, vbus)
      MKPTR1D(zdxdy, dxdy, fbus)
      MKPTR1DK(zfrv_ag, frv, indx_agrege, fbus)
      MKPTR1DK(zftemp_ag, ftemp, indx_agrege, fbus)
      MKPTR1DK(zfvap_ag, fvap, indx_agrege, fbus)
      MKPTR1D(zh, h, fbus)
      MKPTR1D(zhpar, hpar, fbus)
      MKPTR1DK(zhst_ag, hst, indx_agrege, fbus)
      MKPTR1DK(zilmo_ag, ilmo, indx_agrege, fbus)
      MKPTR1D(zlh, lhtg, fbus)
 
      MKPTR1D(zpmoins, pmoins, fbus)
      MKPTR1D(zsdtswd, sdtswd, vbus)
      MKPTR1D(zsdtsws, sdtsws, vbus)
      MKPTR1D(ztstar, tstar, vbus)
      MKPTR1DK(ztsurf, tsurf, indx_agrege, fbus)
      MKPTR1D(zudiag, udiag, fbus)
      MKPTR1D(zvdiag, vdiag, fbus)
      MKPTR1D(zwge, wge, vbus)
      MKPTR1D(zwgmax, wgmax, vbus)
      MKPTR1D(zwgmin, wgmin, vbus)
      MKPTR1D(zwstar, wstar, vbus)
      MKPTR1DK(zz0_ag, z0, indx_agrege, fbus)
      MKPTR1D(zztsl, ztsl, vbus)

      MKPTR2D(zbuoy, buoy, vbus)
      MKPTR2D(zc1pbl, c1pbl, vbus)
      MKPTR2D(zfbl, fbl, fbus)
      MKPTR2D(zfblgauss, fblgauss, fbus)
      MKPTR2D(zfblnonloc, fblnonloc, vbus)
      MKPTR2D(zfnn, fnn, fbus)
      MKPTR2D(zftot, ftot, fbus)
      MKPTR2D(zgq, gq, vbus)
      MKPTR2D(zgql, gql, vbus)
      MKPTR2D(zgte, gte, vbus)
      MKPTR2D(zgzmom, gzmom, vbus)
      MKPTR2D(zhumoins, humoins, dbus)
      MKPTR2D(zkm, km, vbus)
      MKPTR2D(zkt, kt, vbus)
      MKPTR2D(zlwc, lwc, fbus)
      MKPTR2D(zqcmoins, qcmoins, dbus)
      MKPTR2D(zqtbl, qtbl, fbus)
      MKPTR2D(zrif, rif, vbus)
      MKPTR2D(zrig, rig, vbus)
      MKPTR2D(zshear2, shear2, vbus)
      MKPTR2D(zsigm, sigm, dbus)
      MKPTR2D(zsigt, sigt, dbus)
      MKPTR2D(zsigw, sigw, dbus)
      MKPTR2D(ztmoins, tmoins, dbus)
      MKPTR2D(zturbreg, turbreg, fbus)
      MKPTR2D(ztve, tve, vbus)
      MKPTR2D(zumoins, umoins, dbus)
      MKPTR2D(zuwng, uwng, vbus)
      MKPTR2D(zvcoef, vcoef, vbus)
      MKPTR2D(zvmoins, vmoins, dbus)
      MKPTR2D(zvwng, vwng, vbus)
      MKPTR2D(zwqng, wqng, vbus)
      MKPTR2D(zwtng, wtng, vbus)
      MKPTR2D(zzd, zd, fbus)
      MKPTR2D(zze, ze, vbus)
      MKPTR2D(zzn, zn, fbus)

      MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, fbus)

      if (advectke) then
         MKPTR2D(zenmoins, enmoins, dbus)
         MKPTR2D(tke, enplus, dbus)
      else
         nullify(zenmoins)
         MKPTR2D(tke, en, fbus)
      endif

      eturbtau = delt

      !  filtre temporel

      cf1=1.0-tfilt
      cf2=tfilt

      !  initialiser e avec ea,z et h

      if (kount == 0) then

         INIT_TKE: if (any('en'==phyinread_list_s(1:phyinread_n))) then
            tke = max(tke,0.)
         else
            do k=1,nkm1
!vdir nodep
               do i=1,ni
                  tke(i,k)= max( etrmin, blconst_cu*zfrv_ag(i)**2 * &
                       exp(-(zze(i,k)- zze(i,nkm1))/zh(i)) )
               end do
            end do
         endif INIT_TKE

         if (advectke) then
            do k=1,nkm1
!vdir nodep
               do i=1,ni
                  zenmoins(i,k) = tke(i,k)
               end do
            end do
         endif

      endif

      if (kount.gt.0) then
         call series_xst(zzn, 'LM', trnch)
         call series_xst(tke, 'EN', trnch)
      endif

      do k=1,nkm1
!vdir nodep
         do i=1,ni
            if (advectke) then
               enold(i,k) = max(zenmoins(i,k),etrmin)
               tke(i,k)   = max(tke(i,k),etrmin)
            else
               enold(i,k) = tke(i,k)
            endif
         end do
      end do

      ! convective velocity scale w*
      ! (passed to MOISTKE3 through XH)
      qe(:,1:nkm1) = zhumoins(:,1:nkm1)
      stat = neark(zsigt, zpmoins, 1000., ni, nkm1, ksl)
      do i=1,ni
         xb(i)=1.0+DELTA*qe(i,ksl(i))
         xh(i)=(GRAV/(xb(i)*ztve(i,ksl(i)))) * ( xb(i)*zftemp_ag(i) &
              + DELTA*ztve(i,ksl(i))*zfvap_ag(i) )
         fbsurf(i)=xh(i)
         xh(i)=max(0.0,xh(i))
         xh(i) = (zh(i)*xh(i))**(1./3.)
         zwstar(i) = xh(i)
         ztstar(i) = max(zftemp_ag(i)/max(zwstar(i),WSTAR_MIN),0.)
      enddo

      if (fluvert == 'MOISTKE') then

         ! Initialise non-gradient fluxes to zero
         zwtng = 0.
         zwqng = 0.
         zuwng = 0.
         zvwng = 0.

         call moistke13(tke,enold,zzn,zzd,zrif,zrig,zbuoy,zshear2,zkt,zqtbl,zc1pbl,zfnn, &
              zfblgauss,zfblnonloc,zgte,zgq,zgql,zh,zlh,zhpar,zwtng,zwqng,zuwng,zvwng,&
              zumoins,zvmoins,ztmoins,ztve,zhumoins,qe,zpmoins,zsigm,se,zsigw, &
              zze,zz0_ag,zgzmom,zfrv_ag,zwstar,fbsurf,zturbreg, &
              zmrk2,zvcoef,zdxdy,eturbtau,kount,trnch,ni,nkm1)
         if (phy_error_L) return

      else

         ! Interpolate input profiles to correct set of levels
         call vint_thermo2mom(qmom ,zhumoins, zvcoef, ni, nkm1)
         call vint_thermo2mom(tmom, ztmoins,  zvcoef, ni, nkm1)
         te(:,1:nkm1)  = ztmoins(:,1:nkm1)
         qce(:,1:nkm1) = zqcmoins(:,1:nkm1)

         call eturbl13(tke,enold,zzn,zzd,zrif,zturbreg,zrig,zshear2, &
              zgte,zilmo_ag,zfbl,zgql,zumoins, &
              zvmoins,tmom,te,ztve,qmom,qce,qe,zpmoins, &
              zsigm,se,eturbtau,kount,zgq,zkt,zze,zgzmom, &
              std_p_prof,zfrv_ag,xh,zdxdy,trnch,ni,nkm1, &
              zz0_ag)
         if (phy_error_L) return

         ! implicit diffusion scheme: the diffusion interprets the coefficients of the
         ! surface boundary fluxe condition as those in the alfat+bt*ta expression.
         ! Since the diffusion is made on potential temperature, there is a correction
         ! term that must be included in the non-homogeneous part of the expression -
         ! the alpha coefficient. The correction is of the form za*g/cp. It is
         ! relevant only for the clef option (as opposed to the moistke option) since
         ! in this case, although the diffusion is made on potential temperature,
         ! the diffusion calculation takes an ordinary temperature as the argument.

         if (impflx) then
            do i=1,ni
               zalfat(i) = zalfat(i) + zbt_ag(i)*zztsl(i)*GRAV/CPD
            end do
         endif

      endif

      ! Diagnose variables for turbulent wind (gusts and standard deviations)
      if (associated(zwge)) then
         call twind(zwge, zwgmax, zwgmin, zsdtsws, zsdtswd, &
              ztve, enold, zumoins, zvmoins, zudiag, zvdiag, &
              se, zze, zh, zfrv_ag, &
              zwstar, ni, nkm1)
         call series_xst(zwge   , 'WGE' , trnch)
         call series_xst(zwgmax , 'WGX' , trnch)
         call series_xst(zwgmin , 'WGN' , trnch)
         call series_xst(zsdtsws, 'SDWS', trnch)
         call series_xst(zsdtswd, 'SDWD', trnch)
      endif

      if (kount.eq.0) then
         call series_xst(zzn, 'LM', trnch)
         call series_xst(tke, 'EN', trnch)
      endif

      do k=1,nkm1-1
!vdir nodep
         do i=1,ni
            !# IBM conv. ; pas d'avantage a precalculer sqrt ci-dessous
            zkm(i,k) = blconst_ck*zzn(i,k)*sqrt(tke(i,k))
            zkt(i,k) = zkm(i,k)*zkt(i,k)
         end do
      end do

      ! Set near-surface values (lowest thermo and diagnostic), using the Delage
      ! class of stability functions
      do i=1,ni
         if (fluvert /= 'MOISTKE') tke(i,nkm1) = tke(i,nkm1-1)
         tke(i,nk) = tke(i,nkm1)
         uet = zfrv_ag(i)
         ilmot = zilmo_ag(i)
         if (ilmot > 0.) then
            !# hst_local is used to avoid division by zero
            hst_local = max( zhst_ag(i), zztsl(i)+1.)
            fhz = 1-zztsl(i)/hst_local
            fim = 0.5*(1+sqrt(1+4*d97_as*zztsl(i)*ilmot*beta/fhz))
            fit = beta*fim
         else
            fim=(1-dg92_ci*zztsl(i)*ilmot)**(-.16666666)
            fit=beta*fim**2
            fhz=1
         endif
         zzn(i,nkm1)   = KARMAN*(zztsl(i)+zz0_ag(i))/fim
         zkm(i,nkm1)   = uet*zzn(i,nkm1)*fhz
         zkt(i,nkm1)   = zkm(i,nkm1)*fim/fit
         zzn(i,nk) = KARMAN*zz0_ag(i)
         zkm(i,nk) = uet*zzn(i,nk)
         zkt(i,nk) = zkm(i,nk)/beta
      enddo

      return
   end subroutine turbul2

end module turbul
