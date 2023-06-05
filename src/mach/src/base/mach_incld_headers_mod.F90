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
! Fichier/File   : mach_incld_headers_mod.ftn90
! Creation       : H. Landry, Juillet 2008
! Description    : Modules defining explicit interfaces for
!                  mach in cloud subroutines
!
! Extra info     :
!
!============================================================================

module mach_incld_headers_mod
   use mach_incld_upaqr_mod,  only: mach_incld_upaqr
   use mach_incld_steady_mod, only: mach_incld_steady
   interface
!end trap headces

subroutine mach_incld_concmp (t, C, r, b, nptsnz)
   integer(kind=4), intent   (in) :: nptsnz
   real(kind=4),    intent   (in) :: t(nptsnz,5)
   real(kind=4),    intent(inout) :: c(nptsnz,12)
   real(kind=4),    intent   (in) :: r(nptsnz,25)
   real(kind=4),    intent   (in) :: b(nptsnz,5)
end subroutine mach_incld_concmp

subroutine mach_incld_diffun(gaz_conc, pg, dg, cdg, aq, paq, daq, cdaq, q, &
                             ideriv, baq, tempk, tempi, raq, psacw)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
   integer(kind=4), intent   (in) :: ideriv
   real(kind=4),    intent   (in) :: q(2)
   real(kind=4),    intent   (in) :: baq(5, 2)
   real(kind=4),    intent(inout) :: raq(25, 2)
   real(kind=4),    intent   (in) :: gaz_conc(maxnsg)
   real(kind=4),    intent  (out) :: pg      (maxnsg)
   real(kind=4),    intent  (out) :: dg      (maxnsg)
   real(kind=4),    intent  (out) :: cdg     (maxnsg)
   real(kind=4),    intent   (in) :: tempk
   real(kind=4),    intent   (in) :: tempi
   real(kind=4),    intent   (in) :: psacw
   real(kind=4),    intent  (out) :: paq     (maxnsaq, 2)
   real(kind=4),    intent  (out) :: daq     (maxnsaq, 2)
   real(kind=4),    intent  (out) :: cdaq    (maxnsaq, 2)
   real(kind=4),    intent   (in) :: aq      (maxnsaq, 2)
end subroutine mach_incld_diffun

subroutine mach_incld_dochem(gvec, aqvec, qvec, tempkvec, tempivec, bvec, &
                             raqvec, psacwvec, tin, tout, idrf, nptsnz)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq

   integer(kind=4), intent   (in)  :: idrf
   integer(kind=4), intent   (in)  :: nptsnz
   real(kind=4),    intent   (in)  :: tin, tout
   real(kind=4),    intent   (in)  :: qvec     (nptsnz, 2)
   real(kind=4),    intent   (in)  :: bvec     (nptsnz, 5, 2)
   real(kind=4),    intent   (in)  :: tempkvec (nptsnz)
   real(kind=4),    intent   (in)  :: tempivec (nptsnz)
   real(kind=4),    intent   (in)  :: raqvec   (nptsnz, 25, 2)
   real(kind=4),    intent   (in)  :: psacwvec (nptsnz)
   real(kind=4),   intent (inout)  :: gvec     (nptsnz, maxnsg)
   real(kind=4),   intent (inout)  :: aqvec    (nptsnz, maxnsaq, 2)
end subroutine mach_incld_dochem

subroutine mach_incld_fff (FNEW, t, c, r, b, nptsnz, nleft)
   
   integer(kind=4), intent (in) :: nptsnz
   real(kind=4),    intent(out) :: fnew(nptsnz)
   real(kind=4),    intent (in) :: t(nptsnz,5)
   real(kind=4),    intent (in) :: c(nptsnz)
   real(kind=4),    intent (in) :: r(nptsnz,25)
   real(kind=4),    intent (in) :: b(nptsnz,5)
   integer(kind=4), intent (in) :: nleft
end subroutine mach_incld_fff

subroutine mach_incld_findh (T, C, R, B, nptsnz)
   integer(kind=4), intent   (in) :: nptsnz
   real(kind=4),    intent(inout) :: t(nptsnz,5)
   real(kind=4),    intent(inout) :: c(nptsnz,12)
   real(kind=4),    intent(inout) :: r(nptsnz,25)
   real(kind=4),    intent(inout) :: b(nptsnz,5)
end subroutine mach_incld_findh

subroutine mach_incld_funeq (y, c, b, nptsnz)
   integer(kind=4), intent (in)  :: nptsnz
   real(kind=4),    intent (in)  :: b(nptsnz, 5)
   real(kind=4),    intent (in)  :: c(nptsnz, 12)
   real(kind=4),    intent(out)  :: y(nptsnz, 5)
end subroutine mach_incld_funeq

subroutine mach_incld_intrqf (aq, ppaq, ddaq, ccdaq, q, f13, ideriv, temp)
   use mach_cam_utils_mod, only:maxnsaq
   integer(kind=4), intent (in) :: ideriv
   real(kind=4),    intent (in) :: temp
   real(kind=4),    intent (in) :: f13
   real(kind=4),    intent (in) :: q(2)
   real(kind=4),    intent (in) :: aq(maxnsaq, 2)
   real(kind=4),    intent(out) :: ppaq(maxnsaq, 2)
   real(kind=4),    intent(out) :: ddaq(maxnsaq, 2)
   real(kind=4),    intent(out) :: ccdaq(maxnsaq, 2)
end subroutine mach_incld_intrqf

subroutine mach_incld_intsca_il(n1, n2, n3, n4, ideriv, kiter, dt, &
                                gaz_conc, aq, q, b, tempk, tempi, raq, psacw)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
   integer(kind=4), intent   (in) :: n1
   integer(kind=4), intent   (in) :: n2
   integer(kind=4), intent   (in) :: n3
   integer(kind=4), intent   (in) :: n4
   integer(kind=4), intent   (in) :: ideriv
   integer(kind=4), intent  (out) :: kiter
   real(kind=4),    intent(inout) :: dt
   real(kind=4),    intent(inout) :: gaz_conc(maxnsg)
   real(kind=4),    intent(inout) :: aq(maxnsaq, 2)
   real(kind=4),    intent   (in) :: q(2)
   real(kind=4),    intent   (in) :: b(5, 2)
   real(kind=4),    intent   (in) :: tempk
   real(kind=4),    intent   (in) :: tempi
   real(kind=4),    intent(inout) :: raq(25, 2)
   real(kind=4),    intent   (in) :: psacw
end subroutine mach_incld_intsca_il

subroutine mach_incld_main(gaz_conc, aerocon, q_bin, tempk, psacw, rad1, rcrit,&
                           roarow, ibulk, flux, fctr, aeronum, pni, pnk)
   use mach_cam_utils_mod,     only: icom, maxnsg, maxns, isize, nswdep
   integer(kind=4), intent   (in) :: pni, pnk
   integer(kind=4), intent   (in) :: ibulk
   real(kind=4),    intent(inout) :: gaz_conc(pni, pnk, maxns)
   real(kind=4),    intent(inout) :: aerocon (pni, pnk, icom, isize)
   real(kind=4),    intent(inout) :: q_bin   (pni, pnk, isize, 2)
   real(kind=4),    intent   (in) :: tempk   (pni, pnk)
   real(kind=4),    intent   (in) :: psacw   (pni, pnk)
   real(kind=4),    intent   (in) :: rad1    (pni, pnk, isize)
   real(kind=4),    intent   (in) :: rcrit   (pni, pnk)
   real(kind=4),    intent   (in) :: roarow  (pni, pnk)
   real(kind=4),    intent  (out) :: flux    (pni, pnk, nswdep)
   real(kind=4),    intent   (in) :: fctr    (pni, pnk)
   real(kind=4),    intent(inout) :: aeronum (pni, pnk, isize)
end subroutine mach_incld_main

subroutine mach_incld_rates(gaz_conc, dg, pg, cdg, aq, daq, paq, cdaq, &
                            r, b, q, ideriv)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
   integer(kind=4), intent (in) :: ideriv
   real(kind=4),    intent (in) :: q       (2)
   real(kind=4),    intent (in) :: b       (5, 2)
   real(kind=4),    intent (in) :: r       (25, 2)
   real(kind=4),    intent(out) :: dg      (maxnsg)
   real(kind=4),    intent(out) :: pg      (maxnsg)
   real(kind=4),    intent(out) :: cdg     (maxnsg)
   real(kind=4),    intent (in) :: aq      (maxnsaq, 2)
   real(kind=4),    intent(out) :: daq     (maxnsaq, 2)
   real(kind=4),    intent(out) :: paq     (maxnsaq, 2)
   real(kind=4),    intent(out) :: cdaq    (maxnsaq, 2)
   real(kind=4),    intent (in) :: gaz_conc(maxnsg)
end subroutine mach_incld_rates

subroutine mach_incld_soleq(gaz_conc, aq, b, r, nptsnz, iaq, ncw)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
   integer(kind=4), intent   (in) :: nptsnz
   integer(kind=4), intent   (in) :: ncw
   integer(kind=4), intent   (in) :: iaq     (nptsnz)
   real(kind=4),    intent(inout) :: b       (nptsnz, 5, 2)
   real(kind=4),    intent(inout) :: gaz_conc(nptsnz, maxnsg)
   real(kind=4),    intent(inout) :: aq      (nptsnz, maxnsaq, 2)
   real(kind=4),    intent   (in) :: r       (nptsnz, 25, 2)
end subroutine mach_incld_soleq

real function mach_incld_fn_dtnew (delt, k)
   integer(kind=4), intent(in) :: k
   real(kind=4),    intent(in) :: delt
end function mach_incld_fn_dtnew

end interface
end module mach_incld_headers_mod
