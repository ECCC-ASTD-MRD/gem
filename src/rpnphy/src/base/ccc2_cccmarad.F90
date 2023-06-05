!-------------------------------------- LICENCE BEGIN --------------------------
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
!-------------------------------------- LICENCE END ----------------------------

module ccc2_cccmarad
   implicit none
   private
   public :: ccc2_cccmarad2

contains

   !/@*
   subroutine ccc2_cccmarad2(dbus, fbus, vbus, &
        temp, qq, ps, sig, &
        tau, kount , &
        trnch, ni, nkm1, nk, &
        liqwcin, icewcin, liqwpin, icewpin, cldfrac)
      use iso_c_binding
      use, intrinsic :: iso_fortran_env, only: REAL64, INT64
      use debug_mod, only: init2nan
      use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
      use tdpack_const, only: CAPPA, CONSOL2, GRAV, PI, STEFAN
      use cldoppro_MP, only: cldoppro_MP3
      use cldoppro_noMP, only: cldoppro_noMP1
      use diagno_clouds, only: diagno_clouds2
      use prep_cw_rad, only: prep_cw_rad3
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use linoz_param, only: mwt_air, mwt_o3, p_linoz_meso
      use series_mod, only: series_xst, series_isstep
      use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
      use ens_perturb, only: ens_nc2d, ens_spp_get
      use ccc2_uv_raddriv, only: ccc2_uv_raddriv1
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      integer, intent(in) :: kount, trnch, ni, nkm1, nk
      real, pointer, contiguous :: dbus(:), fbus(:), vbus(:)
      real, intent(inout) :: temp(ni,nk), qq(ni,nkm1), ps(ni), sig(ni,nkm1+1)
      real, intent(inout) :: liqwcin(ni,nkm1), icewcin(ni,nkm1), cldfrac(ni*nkm1)
      real, intent(inout) :: liqwpin(ni,nkm1), icewpin(ni,nkm1)
      real, intent(in) :: tau

      !@Authors p. vaillancourt, d. talbot, j. li, rpn, cmc, cccma; (may 2006)
      !@Object Driver routine for radiation
      !        prepares all inputs for radiative transfer scheme
      !        (cloud optical properties, trace gases, ozone, aerosols..)
      !        executes ccc radiative transfer for infrared and solar radiation
      !@Arguments
      !          - input/output -
      ! dbus     field of dyn variables
      ! fbus     field of permanent physics variables
      ! vbus     field of volatile physics variables
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
#include <rmn/msg.h>

      include "surface.cdk"
      include "clefcon.cdk"
      include "ozopnt.cdk"
      include "radiation.cdk"
      include "nbsnbl.cdk"
      include "ccc_tracegases.cdk"
      include "tables.cdk"
      include "phyinput.inc"

      logical, parameter :: SLOPE_L = .true.
      logical, parameter :: DO_UV_ONLY = .true.
      real, parameter :: seuil = 1.e-3
      
      external :: ccc2_ckdlw, ccc2_ckdsw, ccc_dataero, ccc_tracedata

      real :: julien, r0r
      
      logical :: radlinoz_L, radlinghg_L

      logical, dimension(ni) :: thold
      integer, dimension(ni) :: p7, p8
      real, dimension(ni)      :: p1, p3, p4, p5, p6, pbl, albpla, fdl, ful, &
           fslo, rmu0, rmu1, rmu2, avgcos, v1
      real, dimension(nk)      :: p10, p11
      real, dimension(ni*npcl) :: p2
      real, dimension(ni, nk)    :: shtj, tfull
      real, dimension(ni, nkm1)  :: co2, f113, f114, o2
      real, dimension(ni, nbs)   :: salb
      real, dimension(ni, nkm1, 5)   :: tauae
      real, dimension(ni, nkm1, nbs) :: exta, exoma, exomga, fa, taucs, omcs, gcs
      real, dimension(ni, nkm1, nbl) :: absa, taucl, omcl, gcl

      integer(INT64) :: ncsec_deb, ncsec_now, timestep, csec_in_day, day_reminder
      real(REAL64) :: hz_8
      real :: hz, ptopoz, alwcap, fwcap, albrmu, ws
      integer :: i, k, l, iuv, yy, mo, dd, hh, mn, ss, step
      logical :: aerosolback,thisstepisrad,nextstepisrad,thisstepisraduv
      integer :: il1, il2
      character(len=1) :: niuv

      real, dimension(ni) :: dummy1, dummy2, dummy3, dummy4
      real, dimension(ni,nk) :: dum2d, o3uv, o3_vmr, o3_mmr, ch4_vmr, n2o_vmr, cf11_vmr, cf12_vmr
      real, dimension(ni) :: vmod2, vdir, th_air, my_tdiag, my_udiag, my_vdiag
      
      real, target :: dummy1d(ni)
      
#define PHYPTRDCL
#include "cccmarad_ptr.hf"

      include "solcons.cdk"

      data aerosolback /.true./
      ! note: if you use lclw=false; must decomment some lines at end of raddriv.F90

      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'ccc2_cccmarad [BEGIN]')

#undef PHYPTRDCL
#include "cccmarad_ptr.hf"
      
      radlinoz_L = (rad_linoz_L .and. llinoz)
      radlinghg_L = (rad_linoz_L .and. llingh)
      if (.not.associated(zo3fk)) zo3fk => zo3s

      call init2nan(p1, p3, p4, p5, p6, pbl, albpla, fdl, ful, fslo)
      call init2nan(rmu0, v1, p10, p11, p2)
      call init2nan(shtj, tfull, co2, f113, f114, o2, salb)
      call init2nan(tauae, exta, exoma, exomga, fa, taucs, omcs, gcs, absa, taucl)
      call init2nan(omcl, gcl)
      call init2nan(dummy1, dummy2, dummy3, dummy4, vmod2, vdir, th_air, my_tdiag, my_udiag, my_vdiag)
      call init2nan(dum2d, o3uv, o3_vmr, o3_mmr, ch4_vmr, n2o_vmr, cf11_vmr, cf12_vmr)
      
      ! use integer variables instead of actual integers

      il1 = 1
      il2 = ni


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

      do k = 1, nkm1
         do i = 1, ni
            co2(i, k) = rmco2
            zch4(i, k) = rmch4
            zn2o(i, k) = rmn2o
            zcf11(i, k) = rmf11
            zcf12(i, k) = rmf12
            f113(i, k) = rmf113
            f114(i, k) = rmf114
            o2(i, k) = rmo2
         enddo
      enddo

      ! Initialize O3 if not read from a 3D climatology
      if (kount == 0 .and. .not.any(phyinread_list_S(1:phyinread_n) == 'o3ce')) zo3ce = -1.

      if (radghg_L) then
         if (kount == 0) then
            !#TODO: check if var was read from dyn as well
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'ch4c')) zch4c = 1.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'n2oc')) zn2oc = 1.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'cf1c')) zcf1c = 1.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'cf2c')) zcf2c = 1.
         endif
         zch4  = zch4  * zch4c
         zn2o  = zn2o * zn2oc
         zcf11 = zcf11  * zcf1c
         zcf12 = zcf12 * zcf2c

      endif

      ! Initialize O3 tracer for uvindex calculation from climatology, if not found in analysis
      IF_LINOZ0: if (kount == 0 .and. llinoz) then

         if (.not.any(dyninread_list_s == 'o3l') .or.   &
             .not.any(phyinread_list_s(1:phyinread_n) == 'tr/o3l:p')) then

            ! climato_phase2_v2 (Paul V)
            if (minval(zo3ce) >= 0.) zo3lmoins = zo3ce * 1E+9  !micro g /kg air <-- kg /kg air

         endif

      endif IF_LINOZ0

      IF_RDGHG0: if (kount == 0 .and. radghg_L) then

         if (radlinghg_L) then

            if (.not.any(dyninread_list_s == 'ch4l') .or.   &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/ch4l:m')) &
               zch4lmoins = zch4 * 1E+9                      !micro g /kg

            if (.not.any(dyninread_list_s == 'n2ol') .or.   &
                .not.any(phyinread_list_s(1:phyinread_n) == 'tr/n2ol:m')) &
               zn2olmoins = zn2o * 1E+9                      !micro g /kg

            if (.not.any(dyninread_list_s == 'f11l') .or.   &
                .not.any(phyinread_list_s(1:phyinread_n) == 'tr/f11l:m')) &
               zf11lmoins = zcf11 * 1E+9                     !micro g /kg

            if (.not.any(dyninread_list_s == 'f12l') .or.   &
                .not.any(phyinread_list_s(1:phyinread_n) == 'tr/f12l:m')) &
               zf12lmoins = zcf12 * 1E+9                     !micro g /kg

         endif

      endif IF_RDGHG0

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

      ! calculate cloud optical properties and dependent diagnostic

      if (stcond(1:3)=='MP_') then
         call cldoppro_MP3(dbus, fbus, vbus, &
              taucs, omcs, gcs, taucl, omcl, gcl, &
              liqwcin, icewcin, &
              liqwpin, icewpin, cldfrac, &
              temp, sig, ps, ni, nkm1, nk, kount)
         if (phy_error_L) return
      else
         if (.not.associated(ztcsl)) ztcsl => dummy1d
         if (.not.associated(ztcsm)) ztcsm => dummy1d
         if (.not.associated(ztcsh)) ztcsh => dummy1d
         if (.not.associated(ztczl)) ztczl => dummy1d
         if (.not.associated(ztczm)) ztczm => dummy1d
         if (.not.associated(ztczh)) ztczh => dummy1d
         call prep_cw_rad3(fbus, dbus, &
                 temp, qq, ps, sig, &
                 cldfrac, liqwcin, icewcin, liqwpin, icewpin, &
                 kount, ni, nk, nkm1)
         call cldoppro_noMP1(fbus,vbus,taucs, omcs, gcs, taucl, omcl, gcl, &
              liqwcin, icewcin, &
              liqwpin, icewpin, cldfrac, &
              temp, sig, ps, trnch, ni, &
              ni, nkm1, nk)
      endif

      ! calculate diagnostic cloud variables
      ! such as cloud cover, effective and true; cloud top temp and pressure

      call diagno_clouds2(fbus, vbus,taucs, taucl,  &
             zgztherm, cldfrac, &
             temp, sig, ps, trnch, ni, &
             ni, nkm1, nk)

      ! is this or next step a radiation timestep?
      thisstepisraduv = ((kntraduv_S /= '') .and. &
           (kount == 0 .or. mod(kount, kntraduv) == 0))
      thisstepisrad=(kount == 0 .or. mod(kount-1, kntrad) == 0)
      nextstepisrad=(mod(kount, kntrad) == 0)

      csec_in_day  = 8640000
      timestep     =  nint(tau*100.)
      call mu_js2ymdhms(jdateo, yy, mo, dd, hh, mn, ss)
      ncsec_deb    = (ss + mn*60 + hh*3600)*100
      ncsec_now    = ncsec_deb + kount*timestep
      day_reminder =  mod(ncsec_now, csec_in_day)
      hz_8 = day_reminder / 360000.0d0
      hz   = hz_8
      julien = real(jdate_day_of_year(jdateo + kount*int(tau) + MU_JDATE_HALFDAY))

      ! cosine of solar zenith angle at greenwich hour
      call suncos2(rmu0, dummy1, dummy2, dummy3, dummy4, ni, &
           zdlat, zdlon, hz, julien, .not.SLOPE_L)
      if (phy_error_L) return

      ! calculate the variation of solar constant
      alf = julien / 365. * 2 * pi
      r0r = solcons(alf)

      ! pour les pas de temps radiatifs ou pour l'appel diagno a raddriv
      IF_RADIA_RADUV: if (thisstepisrad .or. thisstepisraduv) then

         IF_KOUNT0: if (kount == 0) then
          ! Local estimate of screen level temperature for LW boundary condition
            my_tdiag = ztdiag
            my_udiag = zudiag
            my_vdiag = zvdiag

          !  calculer avgcos(zcosas) et SALB pour kount=0
          !  pour simplifier, on neglige l'effet des vents

            step=kntrad-1
            julien = real(jdate_day_of_year(jdateo + step*int(tau) + MU_JDATE_HALFDAY))
            ncsec_now    = ncsec_deb + step*timestep
            day_reminder = mod(ncsec_now, csec_in_day)
            hz_8 = day_reminder / 360000.0d0
            hz  = hz_8

            call suncos2(rmu2, dummy1, dummy2, dummy3, dummy4, ni, &
                 zdlat, zdlon, hz, julien, .not.SLOPE_L)
            if (phy_error_L) return

            avgcos = (rmu0 + rmu2) * 0.5

            zalwater  = 0.0
            alwcap = 0.3
            do i = 1, ni
               salb(i, 1) = amax1(amin1(zalvis_ag(i), 0.80), 0.03)
               if (zmg(i) <= 0.01 .and. zglsea(i) <= 0.01  &
                    .and. avgcos(i) > seuil) then
                  ws = (my_udiag(i)*my_udiag(i) + my_vdiag(i)*my_vdiag(i))**1.705
                  fwcap      = amin1(3.84e-06 * ws, 1.0)
                  albrmu     = 0.037 / (1.1 * (avgcos(i)**1.4) + 0.15)
                  zalwater(i)  = (1.-fwcap) * albrmu + fwcap * alwcap
                  zalwater(i)  = amax1(amin1(zalwater(i), 0.80), 0.03) ! this max comes from newrad!?!
                  salb(i, 1) = zalwater(i)
               endif
               zsalb6z(i)=salb(i, 1)
               do l = 2, nbs
                  salb(i, l) = salb(i, 1)
               enddo
            enddo

            zcosas=avgcos
         else !IF_KOUNT0
            i = sl_prelim(ztmoins(:, nkm1), zhumoins(:, nkm1), zumoins(:, nkm1), &
                 zvmoins(:, nkm1), zp0_plus, zzusl, spd_air=vmod2, dir_air=vdir, &
                 min_wind_speed=VAMIN)
            if (i /= SL_OK) then
               call physeterror('cccmarad', 'problem in sl_prelim()')
               return
            endif

            th_air = ztmoins(:, nkm1)*zsigt(:, nkm1)**(-cappa)
            i = sl_sfclayer(th_air, zhumoins(:, nkm1), vmod2, vdir, zzusl, zztsl, ztsrad, zqsurf, &
                 zz0_ag, zz0t_ag, zdlat, zfcor, hghtt_diag=1.5, t_diag=my_tdiag)
            if (i /= SL_OK) then
               call physeterror('cccmarad', 'problem in sl_sfclayer()')
               return
            endif
            do i = 1, ni
               ! albedo agregated at previous timestep,set the same for all 4 bands
               salb(i, 1) = amax1(amin1(zalvis_ag(i), 0.80), 0.03)
               zsalb6z(i)=salb(i, 1)
               do l = 2, nbs
                  salb(i, l) = salb(i, 1)
               enddo
            enddo
            zcosas=zcosn  ! use value calculated at previous timestep
         endif IF_KOUNT0

         ! calculte sigma(shtj) and temperature(tfull) at flux levels

         do i = 1, ni
            shtj(i, 1) = sig(i, 1) * sqrt(sig(i, 1) / sig(i, 2))
            shtj(i, nk) = 1.0
            ! The following line extrapolates the temperature above model top
            ! for moon layer temperature
            ! tfull(i, 1) = 0.5 * (3.0 * temp(i, 1) - temp(i, 2))
            ! The following line assumes temperature is isothermal above model top
            ! This assumption must also be imposed in raddriv (see calc of a1(i, 5)
            ! and planck subroutines
            tfull(i, 1) =  temp(i, 1)
            ! Choose boundary condition for LW down flux very near the surface
            ! High vertical resolution SCM tests suggest that average of 2m and ground temperature causes less problems
            tfull(i, nk) = 0.5*(my_tdiag(i)+ztsrad(i))
            ! tfull(i, nk) = my_tdiag(i)
            ! tfull(i, nk) = ztsrad(i)
         enddo
         do k = 2, nkm1
            do i = 1, ni
               shtj(i, k) = sqrt(sig(i, k-1) * sig(i, k))
               tfull(i, k) = 0.5 * (temp(i, k-1) + temp(i, k))
            enddo
         enddo

         ! calculate aerosol optical properties

         pbl(:) = ens_spp_get('aero_pbl', zmrk2, default=1500.)
         call ccc_aerooppro2(tauae, exta, exoma, exomga, fa, absa, &
              temp, shtj, sig, ps, zdlat, zmg, zml, pbl, zmrk2, &
              aerosolback, il1, il2, ni, nkm1, nk)

         ! from ozone zonal monthly climatology: interpolate to proper date
         ! and grid, calculate total amount above model top (ptop)

         call radfac4(zo3fk, zoztoit, sig, nk, nkm1, npcl, zdlat, ps, ni, ni, &
              p2, p3, p4, p5, p6, p7, p8, nlacl, &
              goz(fozon), goz(clat), goz(pref))
         if (phy_error_L) return

         ! must modify oztoit to fit the needs of raddriv who expects an average
         ! mixing ratio rather than an integral (convert cm back to kg/kg)
         do i = 1, ni
            ! ptop = sig(i, 1)*ps(i)
            ptopoz = -10.0
            ! look for ozone reference pressure level closest to model top
            do k = 0, npcl-1
               if (goz(pref+k) < std_p_prof(1)) then
                  ptopoz = goz(pref+k)
               endif
            enddo
            if (ptopoz > 0.) zoztoit(i) = zoztoit(i)* &
                 grav*2.144e-2/ptopoz
         enddo

         ! Initialize surface emissivity
         if (.not.rad_esfc) zemisr(:) = 1.0

         ! Use FK ozone above stratopause 1hPa and era5 (zo3ce) below
         if (.not.associated(zo3fk, zo3s)) zo3s = zo3fk  !kg /kg air mmr
         if (minval(zo3ce) >= 0.) then
            do k = 1, nk
               do i = 1, ni
                  if (sig(i, k) * ps(i) > 100. ) &    !1hPa
                     !    climato+ghg+linoz+era5+ (Irena I) micro g /kg air mmr
                     !zo3s(i, k) = zo3ce(i,k) * 1.E-9      !kg /kg <--  micro g /kg air mmr

                     !    climato_phase2_v2 (Paul V) kg /kg air
                     zo3s(i,k) = zo3ce(i,k)               !kg /kg air
               enddo
            enddo
         endif

         if (llinoz) then
            o3_vmr = zo3lmoins *1.E-9 * mwt_air/ mwt_o3    !mole /mole vmr <-- micro g/kg air
            o3_vmr = max(o3_vmr, 1.E-10)                   !mole /mole vmr
            o3_mmr = o3_vmr * mwt_o3  /mwt_air             !kg /kg air     <-- mole /mole vmr

            ! Use LINOZ ozone (o3_mmr) in radiation
            if (radlinoz_L) zo3s = o3_mmr                                  !kg /kg air

            ! Use LINOZ ozone (o3_mmr) in uvindex calculation
            if (thisstepisraduv) o3uv = o3_mmr                                  !kg /kg air

         else if(thisstepisraduv) then
            ! Ozone climato for index calculation
            o3uv = zo3s 
         endif

         !Use LINOZ ghg (zch4, zn2o, zcf11, zcf12)
         if (radlinghg_L) then
            zch4  = zch4lmoins *1E-9             ! kg /kg air mmr  <- micro g /kg
            zn2o  = zn2olmoins *1E-9             ! kg /kg air mmr  <- micro g /kg
            zcf11 = zf11lmoins *1E-9             ! kg /kg air mmr  <- micro g /kg
            zcf12 = zf12lmoins *1E-9             ! kg /kg air mmr  <- micro g /kg

            zch4  = max(zch4 ,1.E-20)
            zn2o  = max(zn2o ,1.E-20)
            zcf11 = max(zcf11,1.E-20)
            zcf12 = max(zcf12,1.E-20)
         endif

         if (thisstepisraduv) then
            if (timings_L) call timing_start_omp(413, 'rad_uvindex', 410)
            ! appel diagnostique pour calculer fatb,fctb...
            ! utiliser rmu0 (temps courant)
            call ccc2_uv_raddriv1(zfatb, zfadb, zfafb, zfctb, zfcdb, zfcfb, &
                 fslo, zfsamoon, ps, shtj, sig, &
                 temp, o3uv, zoztoit, &
                 qq, co2, zch4,  &
                 o2, rmu0, r0r, salb, taucs, &
                 omcs, gcs, &
                 cldfrac, tauae, exta, exoma, exomga, &
                 fa, zmrk2, &
                 il1, il2, ni, nkm1, nk)

            ! Calcul des l'indices UV
            call ccc2_uvindex2(il1,il2,ni,zfctb,zfatb,salb,ziuvc,ziuva,iuv_method)
            if (timings_L) call timing_stop_omp(413)
         endif


       endif IF_RADIA_RADUV

       IF_RADIA_ON:  if (thisstepisrad) then

         ! call for prognostic outputs
         ! actual call to the Li & Barker (2005) radiation

         call ccc2_raddriv3(zfsg, zfsd0, zfsf0, zfsv0, zfsi0, &
              zfatb0, zfadb0, zfafb0, zfctb0, zfcdb0, zfcfb0, &
              albpla, fdl, ful, zt20, zti, &
              zcstt, zcsb, zclt, zclb, zparr0, &
              zfluxds0, zfluxus0, zfluxdl, zfluxul, &
              fslo, zfsamoon, ps, shtj, sig, &
              tfull, temp, ztsrad, zo3s, zoztoit, &
              qq, co2, zch4, zn2o, zcf11, &
              zcf12, f113, f114, o2, zcosas, r0r, salb, zemisr, taucs, &
              omcs, gcs, taucl, omcl, gcl, &
              cldfrac, tauae, exta, exoma, exomga, &
              fa, absa, rad_sw, rad_lw, zmrk2, .not.DO_UV_ONLY, &
              il1, il2, ni, nkm1, nk)

         ! ti (t2): infrared (solar) cooling (heating) rate
         ! fdsi (fdss): infrared (solar) downward flux at surface.
         ! ei (ev): infrared (solar) upward flux at toa
         ! ap: albedo planetaire.

         thold=(zcosas > seuil .and. rmu0 > seuil)

         do  i = 1, ni
            zfdsi(i)  = fdl(i)
            zei(i)    = ful(i)
            zfusi(i)  = zfluxul(i, nk)
            zfdss0(i) = zfsg(i)
            zev0(i)   = CONSOL2 * r0r * zcosas(i) * albpla(i)

            ! moduler les flux et les taux par le cosinus de l'angle solaire.
            ! rapport des cosinus : angle actuel sur angle moyen.

            v1(i) = rmu0(i) / zcosas(i)
            v1(i) = min(v1(i), 2.0)
            zvv1(i)= v1(i)
            if (thold(i)) then
               zfdss(i)     = zfdss0(i)         * v1(i)
               zev(i)       = zev0(i)           * v1(i)
               zflusolis(i) = (zfsd0(i)+zfsf0(i)) * v1(i)
               zfsd(i)      = zfsd0(i)          * v1(i)
               zfsf(i)      = zfsf0(i)          * v1(i)
               zfsv(i)      = zfsv0(i)          * v1(i)
               zfsi(i)      = zfsi0(i)          * v1(i)
               zparr(i)     = zparr0(i)         * v1(i)
               zfluxds(i, nk) = zfluxds0(i, nk) * v1(i)
               zfluxus(i, nk) = zfluxus0(i, nk) * v1(i)
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

         do k=1, nkm1
            do i=1, ni
               if (thold(i)) then
                  zt2(i, k)     = zt20(i, k) * v1(i)
                  zfluxds(i, k) = zfluxds0(i, k) * v1(i)
                  zfluxus(i, k) = zfluxus0(i, k) * v1(i)
               endif
            enddo
         enddo


      else  !# IF_RADIA_ON

         ! ajustement du solaire aux pas non multiples de kntrad par
         ! modulation avec cosinus de l'angle solaire

         ! moduler par le cosinus de l'angle solaire. mettre a zero les
         ! valeurs appropriees de fdss, ev et t2.

         thold=(zcosas > seuil .and. rmu0 > seuil)
         do i=1, ni
            v1(i) = rmu0(i) / zcosas(i)
            v1(i) = min(v1(i), 2.0)
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

         do k=1, nkm1
            do i=1, ni
               if (thold(i)) then
                  zt2(i, k) = zt20(i, k) * v1(i)
               endif
            enddo
         enddo

         do k=1, nk
            do i=1, ni
               if (thold(i)) then
                  zfluxds(i, k) = zfluxds0(i, k) * v1(i)
                  zfluxus(i, k) = zfluxus0(i, k) * v1(i)
               endif
            enddo
         enddo

         ! end of radiation loop
      endif IF_RADIA_ON

      do i=1, ni

         zcang(i) = rmu0(i)
         ! Need net radiative fluxes of the atmosphere to calculate liquid water static energy residual for radiation scheme,
         ! physics time step and full timestep
         ! NOTE: top of model and NOT top of atmosphere fluxes must be used
         ! znetrad(i) = (SW net TOM - SW net suface) + (LW net TOM - LW net sfc)
         zfnsi(i) =  zfluxdl(i,nk)-zfluxul(i,nk)
         znetrad(i) = ( (zfluxds(i,1)-zfluxus(i,1)) - zfdss(i)) + ( (zfluxdl(i,1)-zfluxul(i,1)) - zfnsi(i) )
         zfns(i) =  zfnsi(i)+zfdss(i)

         ! iv represente le flux entrant au sommet de l'atmosphere
         ! if below ensures iv is zero when sun is set

         if (thold(i)) then
            ziv(i) = CONSOL2 * r0r * rmu0(i)
         else
            ziv(i) = 0.0
            zsalb6z(i)=0.0
         endif

         if (ziv(i) > 1.0) then
            zap(i) = zev(i) / ziv(i)
         else
            zap(i) = 0.
         endif

         p1(i) = ziv(i) - zev(i) - zei(i)
      enddo

      ! extraction pour diagnostics

      IF_SERIES: if (series_isstep()) then
         call series_xst(zti    , 'ti',   trnch)
         call series_xst(zt2    , 't2',   trnch)
         call series_xst(zctp   , 'bp',   trnch)
         call series_xst(zctt   , 'be',   trnch)
         call series_xst(ztopthw, 'w3',   trnch)
         call series_xst(ztopthi, 'w4',   trnch)
         call series_xst(ziv    , 'iv',   trnch)
         call series_xst(p1     , 'nr',   trnch)
         call series_xst(ztcc   , 'tcc',  trnch)
         call series_xst(znt    , 'nt', trnch)
         call series_xst(zecc   , 'ecc',  trnch)
         call series_xst(zeccl  , 'eccl', trnch)
         call series_xst(zeccm  , 'eccm', trnch)
         call series_xst(zecch  , 'ecch', trnch)
         call series_xst(zev    , 'ev',   trnch)
         call series_xst(zei    , 'ei',   trnch)
         call series_xst(zap    , 'ap',   trnch)
         call series_xst(zfdss  , 'fs',   trnch)
         call series_xst(zflusolis, 'fu', trnch)
         call series_xst(zfsd   , 'fsd',  trnch)
         call series_xst(zfsf   , 'fsf',  trnch)
         call series_xst(zfsv   , 'fsv',  trnch)
         call series_xst(zfsi   , 'fsi',  trnch)
         call series_xst(zparr  , 'parr', trnch)
         call series_xst(zclb   , 'clb',  trnch)
         call series_xst(zclt   , 'clt',  trnch)
         call series_xst(zcstt  , 'cst',  trnch)
         call series_xst(zcsb   , 'csb',  trnch)
         call series_xst(zcosas , 'co',   trnch)
         call series_xst(zcang  , 'cx',   trnch)
         call series_xst(zfdsi  , 'fi', trnch)
         call series_xst(zfnsi  , 'si', trnch)

         !PV for uv band fluxes
         ! zfatb, zfadb, zfafb, zfctb, zfcdb, zfcfb
         do iuv=1, RAD_NUVBRANDS
            write(niuv, '(i1)') iuv
            call series_xst(zfatb(:,iuv), 'fat'//trim(niuv), trnch)
            call series_xst(zfadb(:,iuv), 'fad'//trim(niuv), trnch)
            call series_xst(zfafb(:,iuv), 'faf'//trim(niuv), trnch)
            call series_xst(zfctb(:,iuv), 'fct'//trim(niuv), trnch)
            call series_xst(zfcdb(:,iuv), 'fcd'//trim(niuv), trnch)
            call series_xst(zfcfb(:,iuv), 'fcf'//trim(niuv), trnch)
         enddo
         call series_xst(zo3s, 'ozo', trnch)
      endif IF_SERIES

      ! tendances de la radiation
      ztrad = zti + zt2

 !  calculer avgcos(zcosas) et alvis_water(zcosas) pour le prochain pas de temps(kount+1) 
 !  alwater (alvis_water) sera donné au schéma water et ainsi correctement aggrege
      IF_NEXTSTEP: if (nextstepisrad) then
         step=kount+1
         
         julien = real(jdate_day_of_year(jdateo + step*int(tau) + MU_JDATE_HALFDAY))
         ncsec_now    = ncsec_deb + step*timestep
         day_reminder = mod(ncsec_now, csec_in_day)
         hz_8 = day_reminder / 360000.0d0
         hz  = hz_8

         call suncos2(rmu1, dummy1, dummy2, dummy3, dummy4, ni, &
                      zdlat, zdlon, hz, julien, .not.SLOPE_L)
         if (phy_error_L) return

         step=step+kntrad-1
         julien = real(jdate_day_of_year(jdateo + step*int(tau) + MU_JDATE_HALFDAY))
         ncsec_now    = ncsec_deb + step*timestep
         day_reminder = mod(ncsec_now, csec_in_day)
         hz_8 = day_reminder / 360000.0d0
         hz  = hz_8

         call suncos2(rmu2, dummy1, dummy2, dummy3, dummy4, ni, &
              zdlat, zdlon, hz, julien, .not.SLOPE_L)
         if (phy_error_L) return

         avgcos = (rmu1 + rmu2) * 0.5

         ! Local estimate of screen level wind
         IF_WS10m: if (kount == 0) then
            my_udiag = zudiag
            my_vdiag = zvdiag
         else
            i = sl_prelim(ztmoins(:, nkm1), zhumoins(:, nkm1), zumoins(:, nkm1), &
                 zvmoins(:, nkm1), zp0_plus, zzusl, spd_air=vmod2, dir_air=vdir, &
                 min_wind_speed=VAMIN)
            if (i /= SL_OK) then
               call msg_toall(MSG_ERROR, '(cccmarad) problem in sl_prelim()')
               phy_error_L = .true.
               return
            endif

            th_air = ztmoins(:, nkm1)*zsigt(:, nkm1)**(-cappa)
            i = sl_sfclayer(th_air, zhumoins(:, nkm1), vmod2, vdir, zzusl, zztsl, ztsrad, zqsurf, &
                 zz0_ag, zz0t_ag, zdlat, zfcor, hghtm_diag=10.,  &
                 u_diag=my_udiag, v_diag=my_vdiag)
            if (i /= SL_OK) then
               call msg_toall(MSG_ERROR, '(cccmarad) problem in sl_sfclayer()')
               phy_error_L = .true.
               return
            endif
         endif IF_WS10m

         !------------------------------------------------------------------------
         ! open water albedo adjusted for solar angle and white caps,
         ! fwcap is fraction of white caps, alwcap is albedo of white caps
         ! ws is the 10m wind speed, albrmu is albedo corrected for solar zenith angle
         ! white cap effect from : monahan et al., 1980, jpo, 10, 2094-2099
         ! solar angle dependence from  : taylor et al., 1996, qjrms, 122, 839-861
         !------------------------------------------------------------------------

         alwcap = 0.3
         zalwater  = 0.0
         do i = 1, ni
            if (avgcos(i) > seuil) then
               ws = (my_udiag(i)*my_udiag(i) + my_vdiag(i)*my_vdiag(i))**1.705
               fwcap      = amin1(3.84e-06 * ws, 1.0)
               albrmu     = 0.037 / (1.1 * (avgcos(i)**1.4) + 0.15)
               zalwater(i)  = (1.-fwcap) * albrmu + fwcap * alwcap
               zalwater(i)  = amax1(amin1(zalwater(i), 0.80), 0.03) ! this max comes from newrad!?!
            endif
         enddo
         zcosn=avgcos
      endif IF_NEXTSTEP

      call msg_toall(MSG_DEBUG, 'ccc2_cccmarad [BEGIN]')
      !----------------------------------------------------------------
      return
   end subroutine ccc2_cccmarad2

end module ccc2_cccmarad
