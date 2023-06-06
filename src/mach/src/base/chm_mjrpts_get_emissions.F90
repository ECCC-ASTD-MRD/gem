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
! Fichier/File   : chm_mjrpts_get_emissions.ftn90
! Creation       : A. Kallaur, H. Landry, S. Menard - janvier/fevrier 2007
! Description    : Read in the Major point source time dependant emission data.
!                  If model time step = 0 (gem_tstep_num=0), also read in
!                  point source time independant info. (lats, lons, stk height, etc ...)
!
! Extra info     :
!
! Arguments:
!            IN
!              nb_sources --> Total number of major point emissions sources
!              datev      --> major point emission validity date
!
!==============================================================================
!
!!if_on
subroutine chm_mjrpts_get_emissions(file_unit, Fstack_emis, nb_sources, datev, &
                                    err, extra1, extra2, ex1_nom, ex2_nom)
   use chm_mjrpts_sortinfo_mod, only: nb_me_species
!!if_off
   use chm_utils_mod,           only: global_debug, chm_lun_out, NOMV_LEN
   use chm_mjrpts_sortinfo_mod, only: me_species_index
   use chm_species_info_mod,    only: species_master
   implicit  none
!!if_on
   integer(kind=4),            intent (in) :: file_unit, nb_sources, datev
   integer(kind=4),            intent(out) :: err
   real(kind=4),               intent(out) :: Fstack_emis(nb_sources, nb_me_species)
   real(kind=4),     optional, intent(out) :: extra1(nb_sources), extra2(nb_sources)
   character(len=4), optional, intent (in) :: ex1_nom, ex2_nom
!!if_off
!
! Local variables
!
   character(len=NOMV_LEN)   :: sp_name
   integer(kind=4)           :: i, j0, j1, j2
   logical(kind=4)           :: local_dbg
   integer(kind=4), external :: fstlir, fstopc

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

!   Turn off FST information messages
   if (.not. local_dbg) err = fstopc('MSGLVL', 'SYSTEM', 0)
!
   err = 0
   do i = 1, nb_me_species
      sp_name = trim(species_master(me_species_index(i)) % me_name) 
      err = min( fstlir(Fstack_emis(1,i), file_unit, j0, j1, j2, datev, ' ', &
                        -1, -1, -1, ' ', sp_name), err)
   end do

   if (present(extra1)) then
      err = min(fstlir(extra1, file_unit, j0, j1, j2, datev, ' ', -1, -1,  &
                       -1, ' ', ex1_nom), err)
   end if
   if (present(extra2)) then
      err = min(fstlir(extra2, file_unit, j0, j1, j2, datev, ' ', -1, -1,  &
                       -1, ' ', ex2_nom), err)
   end if

   if (local_dbg) write(chm_lun_out, *) 'Exit chm_mjrpts_get_emissions'

   return
end subroutine chm_mjrpts_get_emissions
