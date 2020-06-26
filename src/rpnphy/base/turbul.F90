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
   subroutine turbul2(d, dsiz, f, fsiz, v, vsiz, se, kount, trnch, ni, nk, nkm1)
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
      ! d        dynamic             bus
      ! f        permanent variables bus
      ! v        volatile (output)   bus
      !          - Input -
      ! dsiz     dimension of d
      ! fsiz     dimension of f
      ! vsiz     dimension of v
      ! se       sigma level for turbulent energy
      ! kount    index of timestep
      ! trnch    number of the slice
      ! ni       horizontal dimension
      ! nk       vertical dimension
      ! nkm1     vertical scope of the operator

      integer, intent(in) :: dsiz, fsiz, vsiz
      integer, intent(in) :: kount, trnch, ni, nk, nkm1
      real, target :: d(dsiz), f(fsiz), v(vsiz)
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
           work2(ni), b(ni,nkm1*4), xb(ni), xh(ni)

      real, pointer, dimension(:)   :: zalfat, zbt_ag, zfrv_ag, zftemp_ag, &
           zfvap_ag, zh, zlh, zhst_ag, zilmo_ag, zz0_ag, zztsl, &
           zwstar, zudiag, zvdiag, zhpar, zpmoins, ztsurf, zsdtswd, &
           zsdtsws, zwge, zwgmax, zwgmin, zdxdy, ztstar

      real, pointer, dimension(:,:) :: tke, zenmoins, zkm, zkt, ztmoins, &
           ztve, zze, zzn, zwtng, zwqng, zuwng, zvwng, zqcmoins, zhumoins, &
           zsigm, zsigt, zsigw, zumoins, zvmoins, zfbl, zfblgauss, zfblnonloc, zfnn, &
           zftot, zlwc, zqtbl, zturbreg, zzd, zgq, zgql, zgte, zgzmom, zrif, &
           zrig, zshear2, zvcoef, zqmoins, zc1pbl, zbuoy, zmrk2

      real, dimension(ni,nkm1) :: c, x, wk2d, enold, tmom, qmom, te, qce, qe
      !----------------------------------------------------------------

      call init2nan(work2, xb, xh)
      call init2nan(b, c, x, wk2d, enold, tmom, qmom, te, qce, qe)

      MKPTR1D(zalfat, alfat, v)
      MKPTR1DK(zbt_ag, bt, indx_agrege, v)
      MKPTR1D(zdxdy, dxdy, f)
      MKPTR1DK(zfrv_ag, frv, indx_agrege, f)
      MKPTR1DK(zftemp_ag, ftemp, indx_agrege, f)
      MKPTR1DK(zfvap_ag, fvap, indx_agrege, f)
      MKPTR1D(zh, h, f)
      MKPTR1D(zhpar, hpar, f)
      MKPTR1DK(zhst_ag, hst, indx_agrege, f)
      MKPTR1DK(zilmo_ag, ilmo, indx_agrege, f)
      MKPTR1D(zlh, lhtg, f)
 
      MKPTR1D(zpmoins, pmoins, f)
      MKPTR1D(zsdtswd, sdtswd, v)
      MKPTR1D(zsdtsws, sdtsws, v)
      MKPTR1D(ztstar, tstar, v)
      MKPTR1DK(ztsurf, tsurf, indx_agrege, f)
      MKPTR1D(zudiag, udiag, f)
      MKPTR1D(zvdiag, vdiag, f)
      MKPTR1D(zwge, wge, v)
      MKPTR1D(zwgmax, wgmax, v)
      MKPTR1D(zwgmin, wgmin, v)
      MKPTR1D(zwstar, wstar, v)
      MKPTR1DK(zz0_ag, z0, indx_agrege, f)
      MKPTR1D(zztsl, ztsl, v)

      MKPTR2D(zbuoy, buoy, v)
      MKPTR2D(zc1pbl, c1pbl, v)
      MKPTR2D(zfbl, fbl, f)
      MKPTR2D(zfblgauss, fblgauss, f)
      MKPTR2D(zfblnonloc, fblnonloc, v)
      MKPTR2D(zfnn, fnn, f)
      MKPTR2D(zftot, ftot, f)
      MKPTR2D(zgq, gq, v)
      MKPTR2D(zgql, gql, v)
      MKPTR2D(zgte, gte, v)
      MKPTR2D(zgzmom, gzmom, v)
      MKPTR2D(zhumoins, humoins, d)
      MKPTR2D(zkm, km, v)
      MKPTR2D(zkt, kt, v)
      MKPTR2D(zlwc, lwc, f)
      MKPTR2D(zqcmoins, qcmoins, d)
      MKPTR2D(zqmoins, humoins, d)
      MKPTR2D(zqtbl, qtbl, f)
      MKPTR2D(zrif, rif, v)
      MKPTR2D(zrig, rig, v)
      MKPTR2D(zshear2, shear2, v)
      MKPTR2D(zsigm, sigm, d)
      MKPTR2D(zsigt, sigt, d)
      MKPTR2D(zsigw, sigw, d)
      MKPTR2D(ztmoins, tmoins, d)
      MKPTR2D(zturbreg, turbreg, f)
      MKPTR2D(ztve, tve, v)
      MKPTR2D(zumoins, umoins, d)
      MKPTR2D(zuwng, uwng, v)
      MKPTR2D(zvcoef, vcoef, v)
      MKPTR2D(zvmoins, vmoins, d)
      MKPTR2D(zvwng, vwng, v)
      MKPTR2D(zwqng, wqng, v)
      MKPTR2D(zwtng, wtng, v)
      MKPTR2D(zzd, zd, f)
      MKPTR2D(zze, ze, v)
      MKPTR2D(zzn, zn, f)

      MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, f)

      if (advectke) then
         MKPTR2D(zenmoins, enmoins, d)
         MKPTR2D(tke, enplus, d)
      else
         nullify(zenmoins)
         MKPTR2D(tke, en, f)
      endif

      eturbtau = delt
      qe(:,1:nkm1) = zqmoins(:,1:nkm1)

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
      stat = neark(zsigt, zpmoins, 1000., ni, nkm1, ksl)
      do i=1,ni
         xb(i)=1.0+DELTA*qe(i,ksl(i))
         xh(i)=(GRAV/(xb(i)*ztve(i,ksl(i)))) * ( xb(i)*zftemp_ag(i) &
              + DELTA*ztve(i,ksl(i))*zfvap_ag(i) )
         xh(i)=max(0.0,xh(i))
         work2(i)=zh(i)*xh(i)
      end do
      call vspown1 (xh,work2,1./3.,ni)
      do i=1,ni
         zwstar(i) = xh(i)
         ztstar(i) = max(zftemp_ag(i)/max(zwstar(i),WSTAR_MIN),0.)
      enddo

      if (fluvert == 'MOISTKE') then

         ! Initialise non-gradient fluxes to zero
         zwtng = 0.
         zwqng = 0.
         zuwng = 0.
         zvwng = 0.

         call moistke12(tke,enold,zzn,zzd,zrif,zrig,zbuoy,zshear2,zkt,zqtbl,zc1pbl,zfnn, &
              zfblgauss,zfblnonloc,zgte,zgq,zgql,zh,zlh,zhpar,zwtng,zwqng,zuwng,zvwng,&
              zumoins,zvmoins,ztmoins,ztve,zhumoins,qe,zpmoins,zsigm,se,zsigw, &
              zze,zz0_ag,zgzmom,zfrv_ag,zwstar,zturbreg, &
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

      !    diagnose variables for turbulent wind (gusts and standard deviations)

      if (diag_twind) then
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
