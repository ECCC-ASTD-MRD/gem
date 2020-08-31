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

module ccc1_cccmarad
   implicit none
   private
   public :: cccmarad1

contains

   !/@*
   subroutine cccmarad1(d, dsiz, f, fsiz, v, vsiz , &
        temp, qq, ps, sig, &
        tau, kount , &
        trnch, ni, nkm1, nk, &
        liqwcin, icewcin, liqwpin, icewpin, cldfrac)
      use iso_c_binding
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
      use tdpack_const, only: CAPPA, CONSOL, GRAV, PI, STEFAN
      use cldoppro_MP, only: cldoppro_MP1
      use sfclayer_mod, only: sl_prelim, sl_sfclayer, SL_OK
      use series_mod, only: series_xst, series_isstep
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use ens_perturb, only: ens_nc2d, ens_spp_get
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      !@Arguments

      integer, intent(in) :: dsiz, fsiz, kount, trnch, vsiz, ni, nkm1, nk
      real, intent(inout), target :: d(dsiz), f(fsiz), v(vsiz)
      real, intent(inout) :: temp(ni,nk), qq(ni,nkm1), ps(ni), sig(ni,nkm1+1)
      real, intent(inout) :: liqwcin(ni,nkm1), icewcin(ni,nkm1), cldfrac(ni*nkm1)
      real, intent(inout) :: liqwpin(ni,nkm1), icewpin(ni,nkm1)
      real, intent(in) :: tau

      !          - input/output -
      ! f        field of permanent physics variables
      ! fsiz     dimension of f
      !          - input -
      ! temp     temperature
      ! qq       specific humidity
      ! ps       surface pressure
      ! sig      sigma levels
      ! tau      timestep
      ! kount    number of timesteps
      ! kntrad   frequency of call for infra-red radiation
      ! trnch    index of the vertical plane (ni*nkm1) for which
      !          calculations are to be done.
      ! ni       horizontal dimension
      ! nk       number of flux levels (nkm1+1)
      ! nkm1     number of layers (scope)
      ! liqwcin  in-cloud liquid water content (kg/kg)
      ! icewcin  in-cloud ice    water content (kg/kg)
      ! liqwpin  in-cloud liquid water path (g/m^2)
      ! icewpin  in-cloud ice    water path (g/m^2)
      ! cldfrac  cloud fraction (0.-1.)

      !@Authors p. vaillancourt, d. talbot, j. li, rpn, cmc, cccma; (may 2006)
      !@Object
      !        prepares all inputs for radiative transfer scheme
      !        (cloud optical properties, trace gases, ozone, aerosols..)
      !        executes ccc radiative transfer for infrared and solar radiation

      !@Notes
      !          cccmarad produces:
      !          infra-red rate (ti) of cooling
      !          shortwave rate (t2) of heating
      !          shortwave flux to ground (fdss)
      !          infra-red flux to ground (fdsi)
      !          infra-red flux to the top of the atmosphere (ei)
      !          shortwave flux to the top of the atmosphere (ev)
      !          planetary albedo (ap=ev/incident solar flux)
      !
      ! BEWARE :
      ! Remove comments to the lines at the end of preintp if pressure at model
      ! top is less than .0005 Pa
      ! When pressure a model top is less than 10 hPa then minor bands are used
      ! These variables change values for different topology but do not impact
      ! on validation for different topology : maxc lev1 ncum ncdm in cldifm
      ! mcont               in raddriv
      ! lstart              in qozon3
      !*@/
#include <msg.h>

      include "surface.cdk"
      include "clefcon.cdk"
      include "ozopnt.cdk"
      include "radiation.cdk"
      include "nbsnbl.cdk"
      include "ccc_tracegases.cdk"
      include "tables.cdk"
      include "mcica.cdk"

      logical, parameter :: SLOPE_L = .true.
      real, parameter :: seuil = 1.e-3

      external :: ccc1_ckdlw, ccc1_ckdsw, ccc_dataero, ccc_tracedata

      real :: julien, r0r

      !     pointeurs des variables volatiles de la radiation
      !     determines par une routine de gestion de memoire

      logical,dimension(ni)     :: thold
      integer,dimension(ni)     :: p7,p8
      real,dimension(ni)        :: p1,p3,p4,p5,p6,pbl,albpla,fdl,ful,fslo,&
           rmu0,v1,ws,ws_vs,cosas_vs
      real,dimension(nk)      :: p10,p11
      real,dimension(ni*npcl)   :: p2
      real,dimension(ni,nk)    :: shtj,tfull,s_qrt
      real,dimension(ni,nbs)    :: salb
      real,dimension(ni,nkm1,5)   :: tauae
      real,dimension(ni,nkm1,nbs) :: exta,exoma,exomga,fa,taucs,omcs,gcs
      real,dimension(ni,nkm1,nbl) :: absa,taucl,omcl,gcl

      integer(INT64) :: ncsec_deb, ncsec_now, timestep, csec_in_day, day_reminder
      real(REAL64) :: hz_8
      real :: hz, hzp, ptopoz, alwcap, fwcap, albrmu
!!$      real :: ptop
      integer :: i, k, l, iuv, yy, mo, dd, hh, mn, ss
      logical :: lcsw, lclw, aerosolback
      integer :: il1,il2
      character(len=1) :: niuv

      real, dimension(ni) :: dummy1,dummy2,dummy3,dummy4
      real, dimension(ni) :: vmod2,vdir,th_air,my_tdiag,my_udiag,my_vdiag
      integer :: mpcat

      real :: &
           liqwcin_s(ni,nkm1,nx_loc),  &! subcolumns of cloud liquid water
           icewcin_s(ni,nkm1,nx_loc)    ! subcolumns of cloud ice water

      real :: &
           sigma_qcw(ni,nkm1),         &! std. dev. of cloud water/mean cloud water
           rlc_cf(ni,nkm1),            &! decorelation length for cloud amount (km)
           rlc_cw(ni,nkm1),            &! decorrelation length for cloud condensate(km)
           cldtot(ni)                 ! total cloud fraction as computed using
      ! stochastic cloud generator

      integer :: &
           ncldy(ni),                &! number of cloudy subcolumns
           iseed(ni)                  ! integer pseudo-random number seed

      real :: &
           rseed                     ! real pseudo-random number seed

#define PHYPTRDCL
#include "cccmarad_ptr.hf"

      include "solcons.cdk"

      data lcsw, lclw, aerosolback / .true., .true., .true./

      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'cccmarad1 [BEGIN]')

#undef PHYPTRDCL
#include "cccmarad_ptr.hf"

      if (.not.associated(zo3fk)) zo3fk => zo3s

      call init2nan(p1, p3, p4, p5, p6, pbl, albpla, fdl, ful, fslo)
      call init2nan(rmu0, v1, ws, ws_vs, cosas_vs, p10, p11, p2, cldtot)
      call init2nan(shtj, tfull, s_qrt, salb, sigma_qcw, rlc_cf, rlc_cw)
      call init2nan(tauae, exta, exoma, exomga, fa, taucs, omcs, gcs, absa, taucl)
      call init2nan(omcl, gcl, liqwcin_s, icewcin_s)
      call init2nan(dummy1, dummy2, dummy3, dummy4, vmod2, vdir, th_air, my_tdiag, my_udiag, my_vdiag)

      !  use integer variables instead of actual integers

      il1=1
      il2=ni

      csec_in_day=8640000
      timestep= nint(tau*100.)
      call mu_js2ymdhms(jdateo, yy, mo, dd, hh, mn, ss)
      ncsec_deb= (ss + mn*60 + hh*3600)*100
      ncsec_now= ncsec_deb + kount*timestep
      day_reminder= mod(ncsec_now, csec_in_day)
      hz_8 = day_reminder / 360000.0d0
      hz   = hz_8

      ! redefine co2, ch4, n2o, f11 and f12 concentrations
      ! following corresponding parameters from /OPTIONR/

      co2_ppm = qco2     * 1.e-6
      rmco2   =  co2_ppm * 44d0     / 28.97

      ch4_ppm = qch4     * 1.e-6
      rmch4   =  ch4_ppm * 16.00d0  / 28.97

      n2o_ppm = qn2o     * 1.e-6
      rmn2o   =  n2o_ppm * 44.00d0  / 28.97

      f11_ppm = qcfc11   * 1.e-9
      rmf11   = f11_ppm  * 137.37d0 / 28.97

      f12_ppm = qcfc12   * 1.e-9
      rmf12   = f12_ppm  * 120.91d0 / 28.97

      zt2   = 0.0
      zfdss = 0.0
      zev   = 0.0
      zflusolis = 0.0
      zfsd  = 0.0
      zfsf  = 0.0
      zfsv  = 0.0
      zfsi  = 0.0
      zparr = 0.0
      zfatb = 0.0
      zfadb = 0.0
      zfafb = 0.0
      zfctb = 0.0
      zfcdb = 0.0
      zfcfb = 0.0


      ! calculate the variation of solar constant

      julien = real(jdate_day_of_year(jdateo + kount*int(tau) + MU_JDATE_HALFDAY))
      alf = julien / 365. * 2 * PI
      r0r = solcons(alf)

      ! cosine of solar zenith angle at greenwich hour

      call suncos2(rmu0, dummy1, dummy2, dummy3, dummy4, ni, &
           zdlat, zdlon, hz, julien, .not.SLOPE_L)

      ! calculate cloud optical properties and dependent diagnostic
      ! cloud variables
      ! such as cloud cover, effective and true; cloud top temp and pressure
      ! called every timestep
      if (stcond(1:3)=='MP_') then

         if (stcond(1:5)=='MP_P3')  mpcat = p3_ncat
         if (stcond(1:6)=='MP_MY2') mpcat = 3
         !#TODO: move bus up
         call cldoppro_MP1(d, dsiz, f, fsiz, v, vsiz , &
              taucs, omcs, gcs, taucl, omcl, gcl, &
              liqwcin, icewcin, &
              liqwpin, icewpin, cldfrac, &
              temp, sig, ps, &
              ni, nkm1, nk, mpcat, kount)
         if (phy_error_L) return

      else

         call cldoppro4(taucs, omcs, gcs, taucl, omcl, gcl, &
              ztopthw, ztopthi, zecc,ztcc, &
              zeccl, zeccm, zecch, &
              zctp, zctt, liqwcin, icewcin, &
              liqwpin, icewpin, cldfrac, &
              temp, sig, ps, zmg, zml, zmrk2, ni, &
              ni, nkm1, nk)
      endif

      ! pour les pas de temps radiatifs
      IF_RADIA_ON: if (kount == 0 .or. mod(kount-1, kntrad) .eq. 0) then

         ! Local estimate of screen level temperature and wind

         if (kount == 0 )          then
            my_tdiag = ztdiag
            my_udiag = zudiag
            my_vdiag = zvdiag
         else
            i = sl_prelim(ztmoins(:,nkm1),zhumoins(:,nkm1),zumoins(:,nkm1),zvmoins(:,nkm1), &
                 zp0_plus,zzusl,spd_air=vmod2,dir_air=vdir,min_wind_speed=VAMIN)
            if (i /= SL_OK) then
               call physeterror('cccmarad', 'problem in sl_prelim()')
               return
            endif

            th_air = ztmoins(:,nkm1)*zsigt(:,nkm1)**(-CAPPA)
            i = sl_sfclayer(th_air,zhumoins(:,nkm1),vmod2,vdir,zzusl,zztsl,ztsrad,zqsurf, &
                 zz0_ag,zz0t_ag,zdlat,zfcor,hghtm_diag=10.,hghtt_diag=1.5,t_diag=my_tdiag, &
                 u_diag=my_udiag,v_diag=my_vdiag)
            if (i /= SL_OK) then
               call physeterror('cccmarad', 'problem in sl_sfclayer()')
               return
            endif
         endif

         ! calculte sigma(shtj) and temperature(tfull) at flux levels


         do i = 1, ni
            s_qrt(i,1) = sig(i,1) / sig(i,2)
            s_qrt(i,nk) = 1.0
            ! The following line extrapolates the temperature above model top
            ! for moon layer temperature
            !     tfull(i,1) = 0.5 * (3.0 * temp(i,1) - temp(i,2))
            ! The following line assumes temperature is isothermal above model top
            ! This assumption must also be imposed in raddriv (see calc of a1(i,5)
            ! and planck subroutines
            tfull(i,1) =  temp(i,1)

            !     tfull(i,nk) = 0.5*(my_tdiag(i)+ztsrad(i))
            tfull(i,nk) = my_tdiag(i)

            !     tfull(i,nk) = 0.5*(temp(i,nkm1+1)+ztsrad(i))
            !     tfull(i,nk) = temp(i,nkm1+1)
            !     tfull(i,nk) = ztsrad(i)
            !     tfull(i,nk) : choose either ground temperature or 2m
            !                    temperature (does have an impact)
         enddo
         do k = 2, nkm1
            do i = 1, ni
               s_qrt(i,k) = sig(i,k-1) * sig(i,k)
               tfull(i,k) = 0.5 * (temp(i,k-1) + temp(i,k))
            enddo
         enddo

         call vssqrt (shtj,s_qrt,ni*nk)

         do i = 1, ni
            shtj(i,1)  = sig(i,1) * shtj(i,1)
         enddo

         ! calculate aerosol optical properties

         pbl(:) = ens_spp_get('aero_pbl', zmrk2, default=1500.)
         call ccc_aerooppro2(tauae,exta,exoma,exomga,fa,absa, &
              temp,shtj,sig,ps,zdlat,zmg, zml,pbl,zmrk2, &
              aerosolback,il1, il2, ni, nkm1, nk )

         ! from ozone zonal monthly climatology: interpolate to proper date
         ! and grid, calculate total amount above model top (ptop)

         call radfac4(zo3fk,zoztoit,sig,nk,nkm1,npcl,zdlat,ps,ni,ni, &
              p2, p3, p4, p5, p6, p7, p8, nlacl, &
              goz(fozon), goz(clat), goz(pref))
         if (phy_error_L) return

         ! must modify oztoit to fit the needs of raddriv who expects an average
         ! mixing ratio rather than an integral (convert cm back to kg/kg)
         do i = 1, ni
            !    ptop = sig(i,1)*ps(i)
            ptopoz = -10.0
            !    look for ozone reference pressure level closest to model top
            do k = 0, npcl-1
               if (goz(pref+k) .lt. std_p_prof(1)) then
                  ptopoz = goz(pref+k)
               endif
            enddo
            if (ptopoz.gt.0.0) zoztoit(i)=zoztoit(i)* &
                 GRAV*2.144e-2/ptopoz
         enddo

         ! calculate cosine of solar zenith angle at kount + kntrad - 1

         julien = real(jdate_day_of_year(jdateo + (kount+kntrad-1)*int(tau) + MU_JDATE_HALFDAY))
         ncsec_now= ncsec_deb + (kount + kntrad - 1)*timestep
         day_reminder= mod(ncsec_now, csec_in_day)
         hz_8 = day_reminder / 360000.0d0
         hzp = hz_8

         call suncos2(zcosas, dummy1, dummy2, dummy3, dummy4, ni, &
              zdlat, zdlon, hzp, julien, .not.SLOPE_L)

         do i = 1, ni
            ! albedo (6% to 80%), temporally set the same for all 4 band
            salb(i,1) = amax1 (amin1 (zalvis_ag(i), 0.80), 0.06)
            zsalb6z(i) = salb(i,1)
            !zsalb6z(i) = 0.0
            do l = 2, nbs
               salb(i,l) = salb(i,1)
            enddo

            !...      adjust the cosine of solar zenith angle to radition call time
            zcosas(i) = (rmu0(i) + zcosas(i)) * 0.5

         enddo


         !-----------------------------------------------------------------------
         ! open water albedo adjusted for solar angle and white caps,
         ! fwcap is fraction of white caps, alwcap is albedo of white caps
         ! ws is the 10m wind speed, zcosas is cosine of solar zenith angle,
         ! albrmu is albedo corrected for solar zenith angle
         ! ref for white cap effect is : monahan et al., 1980, jpo, 10,2094-2099
         ! ref for solar angle dependence  : taylor et al., 1996, qjrms,122,839-861
         ! danger: if this code is accepted, it should migrate to where the
         ! agregated albedo is calculated. furthermore, the 10m wind speed
         ! should be recalculated when needed, present ws comes from the beginning
         ! of the previous time step, rather than the end
         !-----------------------------------------------------------------------

         ws_vs=my_udiag*my_udiag+my_vdiag*my_vdiag

         call vspown1(ws, ws_vs, 1.705, ni)
         call vspown1(cosas_vs, zcosas, 1.4, ni)
         alwcap = 0.3
         do i = 1, ni
            ! au pas de temps zero zglsea n est pas defini car
            ! la radiation est faite avant la sfc
            if (zmg(i) .le. 0.01 .and. zglsea(i) .le. 0.01 .and. &
                 zml(i) .le. 0.01 .and. zcosas(i) .gt. seuil ) then
               fwcap      = amin1 (3.84e-06 * ws(i), 1.0)
               albrmu     = 0.037 / (1.1 * cosas_vs(i) + 0.15)
               salb(i,1)  = (1.-fwcap) * albrmu + fwcap * alwcap
               salb(i,1)  = amax1 (amin1 (salb(i,1), 0.80), 0.03)
               do l = 2, nbs
                  salb(i,l) = salb(i,1)
               enddo
            endif
            zsalb6z(i) = salb(i,1)
         enddo

         ! Initialize surface emissivity
         if (.not.rad_esfc) zemisr(:) = 1.0

         IF_SIMISCCP :if (simisccp) then

            ! ISCCP

            ! seed random number generator

            do i = 1, ni

               ! generate the random number based on local latitude,longitude,hour
               ! and julien day.  created so that the size of the seed should not
               ! exceed 2^31-1.  if it does then there will be problems.
               !
               rseed = 1.0e5*((zdlat(i)+(pi/2.0))*2.0*pi+ zdlon(i)) &
                    + hz*1.0e6 &
                    + julien*100.0

               iseed(i) = int(rseed)

            end do

            !  call random_seed(generator=2) ! specific to ibm
            call random_seed(put=iseed)

            ! define the cloud overlap parameters and horizontal variability

            call prep_mcica(rlc_cf, rlc_cw, sigma_qcw, cldfrac, ni, il1, il2, nkm1)

            ! generate sub-olumns of liquid and ice water contents

            call mcica_cld_gen(cldfrac, liqwcin, icewcin, rlc_cf, rlc_cw, &
                 sigma_qcw, temp, sig, ps, ni, il1, il2, nkm1, &
                 ncldy, liqwcin_s, icewcin_s, cldtot)

            ! call the ISCCP simulator

            call isccp_sim_driver( &
                 zitp, zictp, zitau, zicep, zitcf,  &! output
                 zisun, &
                 liqwcin_s, icewcin_s, ps, sig, shtj,         &! input
                 il1, il2, ni, nkm1, nk, &
                 zcosas, ztsrad, temp, qq, zmg, zml)

         endif IF_SIMISCCP

         if (.not.associated(zo3fk, zo3s)) zo3s = zo3fk  !kg /kg air mmr

         ! actual call to the Li & Barker (2005) radiation

         call ccc1_raddriv3 (zfsg,zfsd0,zfsf0,zfsv0,zfsi0, &
              zfatb0,zfadb0,zfafb0,zfctb0,zfcdb0,zfcfb0, &
              albpla,fdl,ful,zt20, zti, &
              zcstt,zcsb,zclt,zclb,zparr0, &
              zfluxds0,zfluxus0,zfluxdl,zfluxul, &
              fslo, zfsamoon, ps, shtj, sig, &
              tfull, temp, ztsrad, zo3s,zoztoit, &
              qq, zcosas, r0r, salb, zemisr, taucs, &
              omcs, gcs, taucl, omcl, gcl, &
              cldfrac, tauae, exta, exoma, exomga, &
              fa, absa, lcsw, lclw, zmrk2, &
              il1, il2, ni, nkm1, nk)

         ! ti (t2): infrared (solar) cooling (heating) rate
         ! fdsi (fdss): infrared (solar) downward flux at surface.
         ! fusi: infrared upward flux at the surface
         ! ei (ev): infrared (solar) upward flux at toa
         ! ap: albedo planetaire.

         thold=(zcosas .gt. seuil .and. rmu0 .gt. seuil)

         do  i = 1, ni
            zfdsi(i)  = fdl(i)
            zfusi(i)  = zemisr(i)*STEFAN*ztsrad(i)**4
            zei(i)    = ful(i)
            zfdss0(i) = zfsg(i)
            zev0(i)   = CONSOL * r0r * zcosas(i) * albpla(i)

            ! moduler les flux et les taux par le cosinus de l'angle solaire.
            ! rapport des cosinus : angle actuel sur angle moyen.

            v1(i) = rmu0(i) / zcosas(i)
            v1(i) = min(v1(i),2.0)
            zvv1(i)= v1(i)
            if (thold(i)) then
               zfdss(i)     = zfdss0(i)         * v1(i)
               zev(i)       = zev0(i)           * v1(i)
               zflusolis(i) = (zfsd0(i)+zfsf0(i))  * v1(i)
               zfsd(i)      = zfsd0(i)          * v1(i)
               zfsf(i)      = zfsf0(i)          * v1(i)
               zfsv(i)      = zfsv0(i)          * v1(i)
               zfsi(i)      = zfsi0(i)          * v1(i)
               zparr(i)     = zparr0(i)         * v1(i)
               zfluxds(i,nk)       = zfluxds0(i,nk)       * v1(i)
               zfluxus(i,nk)       = zfluxus0(i,nk)       * v1(i)
               do iuv=1, RAD_NUVBRANDS
                 zfatb(i,iuv)= zfatb0(i,iuv) * v1(i)
                 zfadb(i,iuv)= zfadb0(i,iuv) * v1(i)
                 zfafb(i,iuv)= zfafb0(i,iuv) * v1(i)
                 zfctb(i,iuv)= zfctb0(i,iuv) * v1(i)
                 zfcdb(i,iuv)= zfcdb0(i,iuv) * v1(i)
                 zfcfb(i,iuv)= zfcfb0(i,iuv) * v1(i)
               enddo
            endif
         enddo

         do  k=1,nkm1
            do  i=1,ni
               if (thold(i)) then
                  zt2(i,k)     = zt20(i,k) * v1(i)
                  zfluxds(i,k) = zfluxds0(i,k) * v1(i)
                  zfluxus(i,k) = zfluxus0(i,k) * v1(i)
               endif
            enddo
         enddo

         !where (zcosas.gt. seuil .and. rmu0 .gt. seuil)

         ! in case mod(kount-1,kntrad) non zero

      else !IF_RADIA_ON

         ! ajustement du solaire aux pas non multiples de kntrad par
         ! modulation avec cosinus de l'angle solaire

         ! moduler par le cosinus de l'angle solaire. mettre a zero les
         ! valeurs appropriees de fdss, ev et t2.

         thold=(zcosas .gt. seuil .and. rmu0 .gt. seuil)
         do i=1,ni
            v1(i) = rmu0(i) / zcosas(i)
            v1(i) = min(v1(i),2.0)
            zvv1(i)= v1(i)
            if (thold(i)) then
               zfdss(i)     = zfdss0(i)         * v1(i)
               zev(i)       = zev0(i)           * v1(i)
               zflusolis(i) = (zfsd0(i)+zfsf0(i))  * v1(i)
               zfsd(i)      = zfsd0(i)          * v1(i)
               zfsf(i)      = zfsf0(i)          * v1(i)
               zfsv(i)      = zfsv0(i)          * v1(i)
               zfsi(i)      = zfsi0(i)          * v1(i)
               zparr(i)     = zparr0(i)         * v1(i)
               do iuv=1, RAD_NUVBRANDS
                 zfatb(i,iuv)= zfatb0(i,iuv) * v1(i)
                 zfadb(i,iuv)= zfadb0(i,iuv) * v1(i)
                 zfafb(i,iuv)= zfafb0(i,iuv) * v1(i)
                 zfctb(i,iuv)= zfctb0(i,iuv) * v1(i)
                 zfcdb(i,iuv)= zfcdb0(i,iuv) * v1(i)
                 zfcfb(i,iuv)= zfcfb0(i,iuv) * v1(i)
               enddo
            endif
         enddo

         do k=1,nkm1
            do i=1,ni
               if (thold(i)) then
                  zt2(i,k) = zt20(i,k) * v1(i)
               endif
            enddo
         enddo

         do k=1,nk
            do i=1,ni
               if (thold(i)) then
                  zfluxds(i,k) = zfluxds0(i,k) * v1(i)
                  zfluxus(i,k) = zfluxus0(i,k) * v1(i)
               endif
            enddo
         enddo

      endif IF_RADIA_ON

      do i=1,ni
         zcang(i) = rmu0(i)

         ! iv represente le flux entrant au sommet de l'atmosphere
         ! if below ensures iv is zero when sun is set

         if (thold(i)) then
            ziv(i) = CONSOL * r0r * rmu0(i)
         else
            ziv(i) = 0.0
         endif

         if (ziv(i).gt. 1.0) then
            zap(i) = zev(i) / ziv(i)
         else
            zap(i) = 0.
         endif

         p1(i) = ziv(i) - zev(i) - zei(i)
      enddo

      ! extraction pour diagnostics

      IF_SERIES: if (series_isstep()) then
         call series_xst(zti    , 'ti', trnch)
         call series_xst(zt2    , 't2', trnch)
         call series_xst(zctp   , 'bp', trnch)
         call series_xst(zctt   , 'be', trnch)
         call series_xst(ztopthw, 'w3', trnch)
         call series_xst(ztopthi, 'w4', trnch)
         call series_xst(ziv    , 'iv', trnch)
         call series_xst(p1     , 'nr', trnch)
         call series_xst(ztcc   , 'tcc', trnch)
         call series_xst(znt    , 'nt', trnch)
         call series_xst(zecc   , 'ecc', trnch)
         call series_xst(zeccl  , 'eccl', trnch)
         call series_xst(zeccm  , 'eccm', trnch)
         call series_xst(zecch  , 'ecch', trnch)
         call series_xst(zev    , 'ev', trnch)
         call series_xst(zei    , 'ei', trnch)
         call series_xst(zap    , 'ap', trnch)
         call series_xst(zfdss  , 'fs', trnch)
         call series_xst(zflusolis, 'fu', trnch)
         call series_xst(zfsd   , 'fsd', trnch)
         call series_xst(zfsf   , 'fsf', trnch)
         call series_xst(zfsv   , 'fsv', trnch)
         call series_xst(zfsi   , 'fsi', trnch)
         call series_xst(zparr  , 'parr', trnch)
         call series_xst(zclb   , 'clb', trnch)
         call series_xst(zclt   , 'clt', trnch)
         call series_xst(zcstt  , 'cst', trnch)
         call series_xst(zcsb   , 'csb', trnch)
         call series_xst(zcosas , 'co', trnch)
         call series_xst(zcang  , 'cx', trnch)
         do iuv=1, RAD_NUVBRANDS
            write(niuv, '(i1)') iuv
            call series_xst(zfatb(:, iuv), 'fat'//trim(niuv), trnch)
            call series_xst(zfadb(:, iuv), 'fad'//trim(niuv), trnch)
            call series_xst(zfafb(:, iuv), 'faf'//trim(niuv), trnch)
            call series_xst(zfctb(:, iuv), 'fct'//trim(niuv), trnch)
            call series_xst(zfcdb(:, iuv), 'fcd'//trim(niuv), trnch)
            call series_xst(zfcfb(:, iuv), 'fcf'//trim(niuv), trnch)
         enddo
         call series_xst(zo3s, 'ozo', trnch)
      endif IF_SERIES

      !  tendances de la radiation
      ztrad = zti + zt2

      call msg_toall(MSG_DEBUG, 'cccmarad1 [END]')
      !----------------------------------------------------------------
      return
   end subroutine cccmarad1

end module ccc1_cccmarad
