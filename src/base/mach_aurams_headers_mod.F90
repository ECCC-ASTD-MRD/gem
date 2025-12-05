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
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_aurams_headers_mod.ftn90
! Creation       : H. Landry, Aout 2008
! Description    : Explicit interfaces for mach_aurams subroutines
!
!
!============================================================================

module mach_aurams_headers_mod
   interface
!end trap heade

subroutine mach_aurams_cldcv2d (cc2d, fnuage, pni, pnk)
   integer(kind=4), intent (in) :: pni, pnk
   real(kind=4),    intent (in) :: fnuage(pni, pnk)
   real(kind=4),    intent(out) :: cc2d  (pni, pnk)
end subroutine mach_aurams_cldcv2d

subroutine mach_aurams_prntrow (fld, name, loc, nlev, loff, il1, il2, jlat)
   use chm_ptopo_grid_mod, only: chm_ni
   integer(kind=4),   intent(in) :: nlev
   real(kind=4),      intent(in) :: fld(chm_ni,nlev)
   character(len=7),  intent(in) :: name
   character(len=14), intent(in) :: loc
   integer(kind=4),   intent(in) :: loff
   integer(kind=4),   intent(in) :: il1
   integer(kind=4),   intent(in) :: il2
   integer(kind=4),   intent(in) :: jlat
end subroutine mach_aurams_prntrow

subroutine mach_aurams_tsysms (xrow, TMASS, il1, il2, nt, loc)
   use chm_ptopo_grid_mod, only: chm_ni, pm_nk
   use mach_cam_utils_mod, only: ntr
   real(kind=4),    intent(out) :: tmass
   integer(kind=4), intent (in) :: il1
   integer(kind=4), intent (in) :: il2
   integer(kind=4), intent (in) :: nt
   integer(kind=4), intent (in) :: loc
   real(kind=4),    intent (in) :: xrow(chm_ni,pm_nk,ntr)
end subroutine mach_aurams_tsysms

end interface
end module mach_aurams_headers_mod
