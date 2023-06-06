!begin trap head
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
! Fichier/File   : mach_cam_headers_mod.ftn90
! Creation       : S. Menard July 2008
! Description    : Modules defining explicit interfaces for mach_cam* subroutines
!
! Extra info     :
!
!============================================================================

module mach_cam_headers_mod
   interface
!end trap head

subroutine gocart_so2so4(chem_tr, metvar3d)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV3D
   real(kind=4), intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4), intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
end subroutine gocart_so2so4

subroutine mach_cam_aeroact(q_bin, rhsize, rcrit, aeronum, zmlwc, &
                            roarow, clsize, ncp, ibulk, pni, pnk)
   use mach_cam_utils_mod, only: isize
   integer(kind=4), intent (in) :: ibulk
   integer(kind=4), intent (in) :: pni, pnk
   real(kind=4),    intent(out) :: q_bin  (pni, pnk, isize)
   real(kind=4),    intent (in) :: rhsize (pni, pnk, isize)
   real(kind=4),    intent(out) :: rcrit  (pni, pnk)
   real(kind=4),    intent (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent (in) :: zmlwc  (pni, pnk)
   real(kind=4),    intent (in) :: roarow (pni, pnk)
   real(kind=4),    intent (in) :: ncp    (pni, pnk)
   real(kind=4),    intent(out) :: clsize (pni, pnk, isize)
end subroutine mach_cam_aeroact

subroutine mach_cam_aerocld(throw, xrow, rhsize, aeronum, thlev, roarow, &
                            pres, zmlwc, jlat, rcrits, tcldcv, flux,     &
                            wetflx, fctr, frevp, rhrow, ccn, kount, pni, pnk)
   use mach_cam_utils_mod,      only: isize, ntr, nswdep
   integer(kind=4), intent   (in) :: jlat
   integer(kind=4), intent   (in) :: kount
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: throw  (pni, pnk)
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhsize (pni, pnk, isize)
   real(kind=4),    intent(inout) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent   (in) :: thlev  (pni, pnk)
   real(kind=4),    intent   (in) :: roarow (pni, pnk)
   real(kind=4),    intent   (in) :: pres   (pni, pnk)
   real(kind=4),    intent   (in) :: zmlwc  (pni, pnk)
   real(kind=4),    intent  (out) :: rcrits (pni, pnk)
   real(kind=4),    intent   (in) :: tcldcv (pni, pnk)
   real(kind=4),    intent  (out) :: flux   (pni, pnk, nswdep)
   real(kind=4),    intent  (out) :: wetflx (pni, nswdep)
   real(kind=4),    intent   (in) :: fctr   (pni, pnk)
   real(kind=4),    intent   (in) :: frevp  (pni, pnk)
   real(kind=4),    intent   (in) :: rhrow  (pni, pnk)
   real(kind=4),    intent   (in) :: ccn    (pni, pnk)
end subroutine mach_cam_aerocld

subroutine mach_cam_aeroprop(rhsize, rhop, rhrow, throw, rgrid, aeronum,  &
                             pni, pnk, amu, amfp, roarow, pdiff, pdepv, trwtrow)
   use mach_cam_utils_mod, only: isize, ntr
   integer(kind=4), intent (in)           :: pni, pnk
   real(kind=4),    intent(out)           :: rhsize (pni, pnk, isize)
   real(kind=4),    intent(out)           :: rhop   (pni, pnk, isize)
   real(kind=4),    intent (in)           :: rhrow  (pni, pnk)
   real(kind=4),    intent (in)           :: throw  (pni, pnk)
   real(kind=4),    intent (in)           :: rgrid  (pni, pnk, ntr)
   real(kind=4),    intent(out)           :: aeronum(pni, pnk, isize)
   real(kind=4),    intent (in), optional :: amu    (pni, pnk)
   real(kind=4),    intent (in), optional :: amfp   (pni, pnk)
   real(kind=4),    intent (in), optional :: roarow (pni, pnk)
   real(kind=4),    intent(out), optional :: pdiff  (pni, pnk, isize)
   real(kind=4),    intent(out), optional :: pdepv  (pni, pnk, isize)
   real(kind=4),    intent(out), optional :: trwtrow(pni, pnk, isize)
end subroutine mach_cam_aeroprop

subroutine mach_cam_cas(tp, colef, rscavg, rhop, roarow, rhsize, pres, &
                        qr, qr_vel, pdiff, pdepv, amu, amfp, pni, pnk)
   use mach_cam_utils_mod, only: isize
   integer(kind=4), intent (in) :: pni, pnk
   real(kind=4),    intent (in) :: tp    (pni, pnk)
   real(kind=4),    intent(out) :: colef (pni, pnk, isize, 2)
   real(kind=4),    intent(out) :: rscavg(pni, pnk, 4, 2)
   real(kind=4),    intent (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent (in) :: roarow(pni, pnk)
   real(kind=4),    intent (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent (in) :: qr    (pni, pnk, 2)
   real(kind=4),    intent(out) :: qr_vel(pni, pnk, 2)
   real(kind=4),    intent (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent (in) :: pres  (pni, pnk)
   real(kind=4),    intent (in) :: amu   (pni, pnk)
   real(kind=4),    intent (in) :: amfp  (pni, pnk)
end subroutine mach_cam_cas

subroutine mach_cam_cas_split(tp, colef, rscavg, rhop, roarow, rhsize, pres, &
                              qr, qr_vel, pdiff, pdepv, amu, amfp, pni, pnk)
   use mach_cam_utils_mod,   only: isize
   integer(kind=4), intent (in) :: pni, pnk
   real(kind=4),    intent (in) :: tp    (pni, pnk)
   real(kind=4),    intent(out) :: colef (pni, pnk, isize, 2)
   real(kind=4),    intent(out) :: rscavg(pni, pnk, 4, 2)
   real(kind=4),    intent (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent (in) :: roarow(pni, pnk)
   real(kind=4),    intent (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent (in) :: qr    (pni, pnk, 2)
   real(kind=4),    intent(out) :: qr_vel(pni, pnk, 2)
   real(kind=4),    intent (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent (in) :: pres  (pni, pnk)
   real(kind=4),    intent (in) :: amu   (pni, pnk)
   real(kind=4),    intent (in) :: amfp  (pni, pnk)
end subroutine mach_cam_cas_split

subroutine mach_cam_coagd(throw, roarow, rtcoa, rhsize, xrow, pdepv, pdiff, &
                          mae, rhop, amu, pni, pnk)
   use mach_cam_utils_mod, only: isize, ntr
   integer(kind=4), intent   (in) :: mae
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: throw (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent  (out) :: rtcoa (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
end subroutine mach_cam_coagd

subroutine mach_cam_condsoa(aeronum, xrow, roarow, rtcond, pcond, soa, pres, &
                            tp, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: roarow (pni, pnk)
   real(kind=4),    intent  (out) :: rtcond (pni, pnk, isize)
   real(kind=4),    intent   (in) :: pcond  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: soa    (pni, pnk)
   real(kind=4),    intent   (in) :: pres   (pni, pnk)
   real(kind=4),    intent   (in) :: tp     (pni, pnk)
end subroutine mach_cam_condsoa

subroutine mach_cam_drydep1(xrow, dryflx, pdepv, thlev, rho, pni, pnk)
   use mach_cam_utils_mod, only: isize, ntr, icom
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent  (out) :: dryflx(pni, icom)
   real(kind=4),    intent   (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: thlev (pni, pnk)
   real(kind=4),    intent   (in) :: rho   (pni, pnk)
end subroutine mach_cam_drydep1

subroutine mach_cam_drydep2(xrow, dryflx, pdepv, thlev, rho, pni, pnk)
   use mach_cam_utils_mod, only: isize, ntr, icom

   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent(out)   :: dryflx(pni, icom)
   real(kind=4),    intent(in)    :: pdepv (pni, pnk, isize)
   real(kind=4),    intent(in)    :: thlev (pni, pnk)
   real(kind=4),    intent(in)    :: rho   (pni, pnk)
end subroutine mach_cam_drydep2

 subroutine mach_cam_drydep_main(iseasn, ra, usi, thlev, roarow, rhsize, fland, &
                                 xrow, amu, pdepv, pdiff, dryflx, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr, icom
   use mach_drydep_mod,      only: lucprm
   integer(kind=4), intent   (in) :: pni, pnk
   integer(kind=4), intent   (in) :: iseasn (pni)
   real(kind=4),    intent   (in) :: usi   (pni)
   real(kind=4),    intent   (in) :: ra    (pni, lucprm)
   real(kind=4),    intent   (in) :: thlev (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent   (in) :: fland (pni, lucprm)
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent(inout) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent  (out) :: dryflx(pni, icom)
end subroutine mach_cam_drydep_main

subroutine mach_cam_drypar(thlev, roarow, pdiff, rhsize, amu, fland, &
                           xrow, pdepv, iseasn, ra, usi, pni, pnk)
   use mach_cam_utils_mod, only: isize, ntr
   use mach_drydep_mod,    only: lucprm
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: thlev (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
   real(kind=4),    intent   (in) :: fland (pni, lucprm)
   real(kind=4),    intent   (in) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent(inout) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: ra    (pni, lucprm)
   real(kind=4),    intent   (in) :: usi   (pni)
   integer(kind=4), intent   (in) :: iseasn(pni)
end subroutine mach_cam_drypar

subroutine mach_cam_flux(busper, busvol, metvar2d, fland)
   use chm_metvar_mod,       only: SIZE_MV2D
   use chm_ptopo_grid_mod,   only: chm_ni
   use mach_drydep_mod,      only: lucprm
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: fland   (chm_ni, lucprm)
end subroutine mach_cam_flux

 subroutine mach_cam_intrsec_inner(xrow, nn, rtcond, aeronum, mae, pni, pnk)
   use mach_cam_utils_mod, only: nbnd, nsb
   integer(kind=4), intent   (in) :: nn
   integer(kind=4), intent   (in) :: mae
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, nsb)
   real(kind=4),    intent   (in) :: rtcond (pni, pnk, nbnd)
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, nbnd)
end subroutine mach_cam_intrsec_inner

 subroutine mach_cam_intrsec_outer(xrow, nn, rtcond, aeronum, mae, pres, &
                                   tp, roarow, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
   integer(kind=4), intent   (in) :: mae
   integer(kind=4), intent   (in) :: nn
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rtcond (pni, pnk, isize)
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent   (in) :: pres   (pni, pnk)
   real(kind=4),    intent   (in) :: tp     (pni, pnk)
   real(kind=4),    intent   (in) :: roarow (pni, pnk)
end subroutine mach_cam_intrsec_outer

subroutine mach_cam_intrsec1_inner(xrow, daqchm, aeronum, rmass, pni, pnk)
   use mach_cam_utils_mod, only: nbnd, nsb
   integer(kind=4),  intent   (in) :: pni, pnk
   real(kind=4),     intent(inout) :: xrow   (pni, pnk, nsb)
   real(kind=4),     intent   (in) :: daqchm (pni, pnk, nbnd)
   real(kind=4),     intent   (in) :: aeronum(pni, pnk, nbnd)
   real(kind=4),     intent   (in) :: rmass  (pni, pnk, nbnd)
end subroutine mach_cam_intrsec1_inner

 subroutine mach_cam_intrsec1_outer(xrow, rhopd, daqchm, aeronum, q_bin, &
                                    rcrit, iswitch, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
   integer(kind=4), intent   (in) :: iswitch
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhopd  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: daqchm (pni, pnk, isize)
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent   (in) :: q_bin  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rcrit  (pni, pnk)
end subroutine mach_cam_intrsec1_outer

subroutine mach_cam_main(jlat, throw, rhrow, xrow, kount, zmlwc, qr, pres,  &
                         roarow, thlev, rtso2, fland, soa, iseasn, ra, usi, &
                         tcldcv, fctr, frevp, wetflx, dryflx, ccn, trwtrow, &
                         pni, pnk)
   use mach_cam_utils_mod,      only: isize, icom, nswdep, ntr
   use mach_drydep_mod,         only: lucprm
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
end subroutine mach_cam_main

subroutine mach_cam_rain(tp, qr, rhsize, pdepv, pres, roarow, pdiff, &
                         amu, amfp, xrow, rtbcld, rhop, cc2d, qr_vel, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: tp    (pni, pnk)
   real(kind=4),    intent   (in) :: qr    (pni, pnk, 2)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent   (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: pres  (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent  (out) :: rtbcld(pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: cc2d  (pni, pnk)
   real(kind=4),    intent  (out) :: qr_vel(pni, pnk, 2)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
   real(kind=4),    intent   (in) :: amfp  (pni, pnk)
end subroutine mach_cam_rain

subroutine mach_cam_scaveng(throw, xrow, cf, evpfac, pdepv, qr, rtbcld, &
                            rhsize, rhop, pdiff, thlev, roarow, wetflx, &
                            flux, pres, amu, amfp, pni, pnk)
   use mach_cam_utils_mod,      only: isize, nswdep, ntr
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: throw (pni, pnk)
   real(kind=4),    intent(inout) :: xrow  (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: cf    (pni, pnk)
   real(kind=4),    intent   (in) :: evpfac(pni, pnk)
   real(kind=4),    intent   (in) :: pres  (pni, pnk)
   real(kind=4),    intent   (in) :: pdepv (pni, pnk, isize)
   real(kind=4),    intent   (in) :: qr    (pni, pnk, 2)
   real(kind=4),    intent  (out) :: rtbcld(pni, pnk, ntr)
   real(kind=4),    intent   (in) :: rhsize(pni, pnk, isize)
   real(kind=4),    intent   (in) :: rhop  (pni, pnk, isize)
   real(kind=4),    intent   (in) :: pdiff (pni, pnk, isize)
   real(kind=4),    intent   (in) :: thlev (pni, pnk)
   real(kind=4),    intent   (in) :: roarow(pni, pnk)
   real(kind=4),    intent   (in) :: amu   (pni, pnk)
   real(kind=4),    intent   (in) :: amfp  (pni, pnk)
   real(kind=4),    intent(inout) :: wetflx(pni, nswdep)
   real(kind=4),    intent(inout) :: flux  (pni, pnk, nswdep)
end subroutine mach_cam_scaveng

subroutine mach_cam_sfss (tdiag, surfwd, rsfrow, fland, pni)
   use mach_cam_utils_mod, only: isize
   use mach_drydep_mod,    only: lucprm
   integer(kind=4), intent   (in) :: pni
   real(kind=4),    intent   (in) :: tdiag (pni)
   real(kind=4),    intent   (in) :: surfwd(pni)
   real(kind=4),    intent  (out) :: rsfrow(pni, isize)
   real(kind=4),    intent   (in) :: fland (pni, lucprm)
end subroutine mach_cam_sfss

subroutine mach_cam_sulfate(aeronum, xrow, roarow, pres, tp, rh, rhsize, &
                            rtnucl, rtcond, rtso2, pcond, mae, pni, pnk)
   use mach_cam_utils_mod,   only: isize, ntr
   integer(kind=4), intent   (in) :: mae
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent(inout) :: xrow   (pni, pnk, ntr)
   real(kind=4),    intent   (in) :: roarow (pni, pnk)
   real(kind=4),    intent   (in) :: pres   (pni, pnk)
   real(kind=4),    intent   (in) :: tp     (pni, pnk)
   real(kind=4),    intent   (in) :: rh     (pni, pnk)
   real(kind=4),    intent   (in) :: rhsize (pni, pnk, isize)
   real(kind=4),    intent  (out) :: rtnucl (pni, pnk)
   real(kind=4),    intent  (out) :: rtcond (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rtso2  (pni, pnk)
   real(kind=4),    intent  (out) :: pcond  (pni, pnk, isize)
end subroutine mach_cam_sulfate

subroutine mach_mie_opt(busper, busvol, aero_tr, trwtrow, thlev, rhoa, pni, pnk,&
                        nmod, kmod)
   use chm_ptopo_grid_mod,   only: chm_nk
   use mach_cam_utils_mod,   only: isize, ntr
   integer(kind=4), intent   (in) :: pni, pnk
   integer(kind=4), intent   (in) :: nmod   (pni)
   integer(kind=4), intent   (in) :: kmod   (pni, chm_nk)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: aero_tr(pni, pnk, ntr)
   real(kind=4),    intent   (in) :: trwtrow(pni, pnk, isize)
   real(kind=4),    intent   (in) :: thlev  (pni, pnk)
   real(kind=4),    intent   (in) :: rhoa   (pni, pnk)
end subroutine mach_mie_opt

subroutine mach_pm_chem(busvol, busper, chem_tr, metvar2d, pmetvar3d, fland, &
                        oldso4, iseasn, istep, f_j, pni, pnk, pnka, nmod, &
                        kmod, tracers_can, kcan)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk, nkc
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use mach_drydep_mod,      only: lucprm
   integer(kind=4), intent   (in) :: istep
   integer(kind=4), intent   (in) :: f_j
   integer(kind=4), intent   (in) :: pnk, pnka
   integer(kind=4), intent   (in) :: pni
   integer(kind=4), intent   (in) :: nmod     (pni)
   integer(kind=4), intent   (in) :: kmod     (pni, chm_nk)
   integer(kind=4), intent   (in) :: iseasn   (chm_ni)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent(inout) :: chem_tr  (chm_ni, chm_nk + 1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: oldso4   (pni, pnka)
   real(kind=4),    intent   (in) :: metvar2d (chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: pmetvar3d(pni, pnka, SIZE_MV3D)
   real(kind=4),    intent   (in) :: fland    (chm_ni, lucprm)
   real(kind=4),    intent(inout), optional :: tracers_can(pni, nkc, nb_dyn_tracers)
   integer(kind=4), intent   (in), optional :: kcan    (pni, nkc)
end subroutine mach_pm_chem

end interface
end module mach_cam_headers_mod
