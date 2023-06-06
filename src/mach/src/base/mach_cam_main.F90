!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_cam_main.ftn90
! Creation       : S. Gong, W. Gong, V. Bouchet, S. Menard, S. Gravel, B. Pabla and
!                  P. Huang for GEM-MACH, June 2008
! Description    : Main entry of the CAM - CANADIAN AEROSOL MODULE VERSION 1.0
!
! Extra info     : - First version created by S. Gong Sept 05 1994 for CAM
!                  - Vectorized the whole program (S. Gong, Jan 19/1996)
!                  - Add sulphate aerosol (S. Gong, Jul 09, 1997)
!                  - Interfaced with adom gas drydeposition (S. Gong, Oct 01, 1997)
!                  - CAM version 1.0 with aerosol activation and cloud chemistry
!                    (S. Gong, Jan 01/1998)
!                  - Adapted for use in aurams. (S. Gong, Jun 20, 1999)
!                  - Modifications for aurams v0.25 (L. Zhang, V. Bouchet, Sept 2000)
!                  - Implement aerocld_new from v0.20+cldchem (W. Gong, Mar 2001)
!                  - Implement wet fluxes (modified scaveng.f) (W. Gong, May 2001)
!                  - Modification to implement het. chem. (V. Bouchet, Aug 2001)
!                  - Documentation/headers (S. Menard, July 2008)
!
! Arguments:  IN
!                  jlat    -> J Slice number
!                  throw   -> Temp
!                  rhrow   -> Relative humidity
!                  kount   -> Chemistry time step
!                  zmlwc   -> CWC content (bulk) (kg/kg)
!                  qr      -> rain water/snow
!                  pres    -> pressure (Pa)
!                  roarow  -> Air density (kg/m3)
!                  thlev   -> Layer thickness [m]
!                  rtso2   -> so2 oxidation rate
!                  fland   -> Landuse
!                  soa     -> Secondary organics aerosols
!                  iseasn  -> Assigned season descriptors
!                  ra      -> Aerodynamic resistance for species "species" (s/m)
!                  usi     -> friction velocity
!                  fctr    -> Cloud-to-rain conversion rate (s-1)
!                  frevp   -> evap. of strat. precip (consun)
!                  wetflx  -> Wet flux
!                  ccn     -> cloud droplet number density [1/m3]
!
!             OUT
!                  TRWTROW -> Aerosol liquid water content for each bin
!                  DRYFLX  -> Dry deposition flux
!
!             IN/OUT
!                  XROW    -> tracers conncentration
!                  TCLDCV  -> total cloud cover
!
!             LOCAL VARIABLE
!                  RTCOA   -> Coagulation rate
!                  RTBCLD  -> Rain Scavenging
!                  RTNUCL  -> Nucleation
!============================================================================
!
!!if_on
subroutine mach_cam_main(jlat, throw, rhrow, xrow, kount, zmlwc, qr, pres,  &
                         roarow, thlev, rtso2, fland, soa, iseasn, ra, usi, &
                         tcldcv, fctr, frevp, wetflx, dryflx, ccn, trwtrow, &
                         pni, pnk)
   use mach_cam_utils_mod,      only: isize, icom, nswdep, ntr
   use mach_drydep_mod,         only: lucprm
!!if_off
   use chm_utils_mod,           only: global_debug, chm_lun_out, chm_error_l, &
                                      dbg_id, dbg_jd
   use mach_cam_headers_mod,    only: mach_cam_aeroprop, mach_cam_sulfate, &
                                      mach_cam_condsoa,  mach_cam_aerocld, &
                                      mach_cam_coagd,    mach_cam_scaveng, &
                                      mach_cam_drydep_main
   use mach_aurams_headers_mod, only: mach_aurams_tsysms, mach_aurams_prntrow
   use mach_cam_utils_mod,      only: icob, iae1, iae2, chk_trc
   use chm_nml_mod,             only: chm_hetchem_s, chm_aqueous_s, chm_soa_s, &
                                      chm_pm_drydep_s, chm_diag_aero_opt_l,    &
                                      chm_pm_coag_step_intvl, chm_blcld_l,     &
                                      dbg_lev, dbg_step, dbg_itr
   implicit none
!!if_on
   integer(kind=4), intent   (in)                             :: jlat
   integer(kind=4), intent   (in)                             :: kount
   integer(kind=4), intent   (in)                             :: pni, pnk
   integer(kind=4), intent   (in), dimension(pni)             :: iseasn
   real(kind=4),    intent   (in), dimension(pni, pnk)        :: throw, rhrow, thlev
   real(kind=4),    intent   (in), dimension(pni, pnk)        :: soa, roarow, pres
   real(kind=4),    intent   (in), dimension(pni, pnk)        :: fctr, frevp, rtso2
   real(kind=4),    intent   (in), dimension(pni)             :: usi
   real(kind=4),    intent   (in), dimension(pni, pnk, 2)     :: qr
   real(kind=4),    intent   (in), dimension(pni, pnk)        :: zmlwc, tcldcv
   real(kind=4),    intent   (in), dimension(pni, lucprm)     :: ra, fland
   real(kind=4),    intent   (in), dimension(pni, pnk)        :: ccn
   real(kind=4),    intent  (out), dimension(pni, nswdep)     :: wetflx
   real(kind=4),    intent  (out), dimension(pni, icom)       :: dryflx
   real(kind=4),    intent  (out), dimension(pni, pnk, isize) :: trwtrow
   real(kind=4),    intent(inout), dimension(pni, pnk, ntr)   :: xrow
!!if_off
!
! local variables
!
   integer(kind=4)                               :: izon, l, i
   real(kind=4)                                  :: tmass
   real(kind=4), dimension(pni, pnk, nswdep)     :: flux
   real(kind=4), dimension(pni, pnk, isize)      :: pdiff, colef, rhop
   real(kind=4), dimension(pni, pnk, isize)      :: aeronum, pdepv, rtcond, rhsize
   real(kind=4), dimension(pni, pnk)             :: amu, amfp, rcrits, rtnucl
   real(kind=4), dimension(pni, pnk, ntr)        :: rtcoa, rtbcld
   logical(kind=4)                               :: local_dbg

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

   if (local_dbg) write (chm_lun_out, *) ' cam --> start'

!  section 0: initialization
   rtnucl = 0.0

   izon   = 3

   if (local_dbg .and. jlat == dbg_jd .and. kount == dbg_step) then
      call mach_aurams_prntrow(xrow(1, dbg_lev, dbg_itr) , 'xrow   ', 'in cam        ',  1, 1, &
                               1, pni, jlat)
      call mach_aurams_prntrow(tcldcv(1, dbg_lev)    , 'cloud  ', 'total         ',  1, 1, &
                               1, pni, jlat)
      call mach_aurams_prntrow(zmlwc(1, dbg_lev)     , 'clw    ', 'stratiform    ',  1, 1, &
                               1, pni, jlat)
      call mach_aurams_prntrow(qr(1, dbg_lev, 1:2)        , 'precip ', 'stratiform    ',  1, 1, &
                               1, pni, jlat)
   end if

!  section 2:  aerosol properties of clear sky
!     a - ambient radius and density
!     b - number density

   do l = 1, pnk
      do i = 1, pni
! air's dynamic viscosity
         amu(i, l)  = 145.8 * 1.0e-8 * throw(i, l) ** 1.5 / (throw(i, l) + 110.4)
!
! mean molecular free path. k.v. beard [1976], j atm. sci., 33
         amfp(i, l) = 6.54e-8 * (amu(i, l) / 1.818e-5) * (1.013e5 / pres(i, l)) *  &
                      (throw(i, l) / 293.15) ** (1.0 / 2.0)
      end do
   end do
!
   call mach_cam_aeroprop(rhsize, rhop, rhrow, throw, xrow, aeronum, &
                          pni, pnk, amu, amfp, roarow, pdiff, pdepv)

   if (jlat == dbg_jd .and. local_dbg) then
      write(chm_lun_out, *)'kt = ', kount
      call  mach_aurams_tsysms(xrow, tmass, 1, pni, izon, 1)
      write(chm_lun_out, *)' aeroprop -> aeronum  ', aeronum(dbg_id, dbg_lev, isize)
      write(chm_lun_out, *)' aeroprop -> xrow     ', xrow(dbg_id, dbg_lev, dbg_itr)
   end if

!  sections 3 and 4 have been swithed to allow for the proper addition of heterogeneous chemistry
!
!  section 4:  sulphur physics and chemistry [clear sky]
!              soa condensation
!              a - sulphur chemistry
!              b - nucleation and condensation

   call mach_cam_sulfate(aeronum, xrow, roarow, pres, throw, rhrow, rhsize, &
                         rtnucl, rtcond, rtso2, colef, 0, pni, pnk)

   if (jlat == dbg_jd .and. local_dbg) then
      write(chm_lun_out, *)'kt = ', kount
      call  mach_aurams_tsysms(xrow, tmass, 1, pni, izon, 3)
      write(chm_lun_out, *)' sulfate -> rtnucl  ', rtnucl(dbg_id, dbg_lev)
      write(chm_lun_out, *)' sulfate -> rtcon   ', rtcond(dbg_id, dbg_lev, isize)
      write(chm_lun_out, *)' sulfate -> rtso2   ', rtso2(dbg_id, dbg_lev)
      write(chm_lun_out, *)' sulfate -> xrow    ', xrow(dbg_id, dbg_lev, dbg_itr)
   end if

   if (chm_soa_s /= 'NIL') then
      call mach_cam_condsoa(aeronum, xrow, roarow, rtcond, colef, soa, &
                            pres, throw, pni, pnk)
   end if

!  section 3: aerosol-cloud interaction
!          a - aerosol activation
!          b - gas-aqueous phase transfer and in-cloud oxidation of siv to svi
!          c - cloud-to-rain conversion and removal
!  RTSO2 is not used in mach_cam_aerocld, WG mar 2001.

   if (chm_aqueous_s == 'GONG' .or. chm_hetchem_s /= 'NIL') then

      call mach_cam_aerocld(throw, xrow, rhsize, aeronum, thlev, roarow, &
                            pres, zmlwc, jlat, rcrits, tcldcv, flux,     &
                            wetflx, fctr, frevp, rhrow, ccn, kount, pni, pnk)
      if (chm_error_l) return
   end if

! Re-check the concentrations
   call chk_trc(xrow, iae1, iae2, jlat, 'aft_acld', kount, pni, pnk)
   if (chm_error_l .and. (.not.local_dbg)) return

   if (jlat == dbg_jd .and. local_dbg) then
      write(chm_lun_out, *)'kt = ', kount
      call  mach_aurams_tsysms(xrow, tmass, 1, pni, izon, 2)
      write(chm_lun_out, *)' aerocld -> cld fr  ', tcldcv(dbg_id, dbg_lev)
      write(chm_lun_out, *)' aerocld -> cld wt  ', zmlwc(dbg_id, dbg_lev)
      write(chm_lun_out, *)' aerocld -> precip  ', qr(dbg_id, dbg_lev, 1:2)
      write(chm_lun_out, *)' aerocld -> xrow    ', xrow(dbg_id, dbg_lev, dbg_itr)
   end if

! Section 5: coagulation
   if (chm_pm_coag_step_intvl > 0 .and. &
       mod(kount, chm_pm_coag_step_intvl) == 0) then

      call mach_cam_coagd(throw, roarow, rtcoa, rhsize, xrow, pdepv, pdiff, &
                          0, rhop, amu, pni, pnk)

      if (jlat == dbg_jd .and. local_dbg) then
         write(chm_lun_out, *)'kt = ', kount, '  icob = ', icob
         write(chm_lun_out, *)' coagd -> rtcoa  ', rtcoa(dbg_id, dbg_lev, dbg_itr)
         write(chm_lun_out, *)' coagd -> xrow   ', xrow(dbg_id, dbg_lev, dbg_itr)
         call  mach_aurams_tsysms(xrow, tmass, 1, pni, izon, 4)
      end if
   end if

!  section 6: below-cloud scavenging of aerosols
!          a - rain or snow scavenging

   if (chm_blcld_l) then
      call mach_cam_scaveng(throw, xrow, tcldcv, frevp, pdepv, qr, rtbcld, &
                            rhsize, rhop, pdiff, thlev, roarow, wetflx,    &
                            flux, pres, amu, amfp, pni, pnk)

      if (jlat == dbg_jd .and. local_dbg) then
         write(chm_lun_out, *)'kt = ', kount
         call mach_aurams_tsysms(xrow, tmass, 1, pni, izon, 5)
         write(chm_lun_out, *)' scaveng -> rtbcld  ', rtbcld(dbg_id, dbg_lev, dbg_itr)
         write(chm_lun_out, *)' scaveng -> xrow    ', xrow(dbg_id, dbg_lev, dbg_itr)
      end if
   end if

!  section 7: sedimentation and dry deposition
!          a - particle dry deposition and gravitational settling

   if (chm_pm_drydep_s /= 'NIL') then

!  drygas subroutine is removed from drydepo, thus all drygas
!  related varibles are removed from drydepo and this subroutine (Sep. 2000, Leiming Zhang)

      call mach_cam_drydep_main(iseasn, ra, usi, thlev, roarow, rhsize, fland, &
                                xrow, amu, pdepv, pdiff, dryflx, pni, pnk)
!
      if (jlat == dbg_jd .and. local_dbg) then
         write(chm_lun_out, *)'kt = ', kount
         call mach_aurams_tsysms(xrow, tmass, 1, pni, izon, 6)
         write(chm_lun_out, *)' drydepo -> xrow   ', xrow(dbg_id, dbg_lev, dbg_itr)
      end if
   end if
!
!  Second call to mach_cam_aeroprop so that aerosol water matches new aerosol
!  distribution on output:
   if (chm_diag_aero_opt_l) then
      call mach_cam_aeroprop(rhsize, rhop, rhrow, throw, xrow, aeronum, &
                             pni, pnk, trwtrow=trwtrow)
   else
      trwtrow = 0.0
   end if

!  section 8: diagnostic outputs

   if (local_dbg .and. jlat == dbg_jd .and. kount == dbg_step) then
      call mach_aurams_prntrow(rhsize(1, dbg_lev, 1),     'rhsize ', 'in cam        ',  1, 1, &
                                1, pni, jlat)
      call mach_aurams_prntrow(aeronum(1, dbg_lev, 1),    'aeronum', 'in cam        ',  1, 1, &
                                1, pni, jlat)
      call mach_aurams_prntrow(rhop(1, dbg_lev, 1),       'rhop   ', 'in cam        ',  1, 1, &
                                1, pni, jlat)
      call mach_aurams_prntrow(rtso2(1, dbg_lev),         'rtso2  ', 'in cam        ',  1, 1, &
                                1, pni, jlat)
      call mach_aurams_prntrow(rtcond(1, dbg_lev, isize), 'rtcond ', 'in cam        ',  1, 1, &
                                1, pni, jlat)
      call mach_aurams_prntrow(rtnucl(1, dbg_lev),        'rtnucl ', 'in cam        ',  1, 1, &
                                1, pni, jlat)
      call mach_aurams_prntrow(rcrits(1, dbg_lev),        'rcrit  ', 'in cam        ',  1, 1, &
                                1, pni, jlat)
   end if

   return
end subroutine mach_cam_main
