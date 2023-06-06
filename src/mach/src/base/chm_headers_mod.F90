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
! Fichier/File   : chm_headers_mod.ftn90
! Creation       : H. Landry, Mai 2008
! Description    : Modules defining explicit interfaces for chm s/r
!
! Extra info     :
!
!============================================================================
!
module chm_headers_mod

interface
!end trap head

subroutine chm_businit(F_ni, F_nk)
   integer(kind=4), intent(in) ::  F_ni, F_nk
end subroutine chm_businit

subroutine chm_exe(busdyn     , busper        , busvol     ,   &
                    slab_index , step)
   integer(kind=4), intent   (in) :: slab_index, step
   real(kind=4), dimension(:), pointer, contiguous :: busdyn, busper, busvol
end subroutine chm_exe

subroutine chm_fst_closefile(file_unit)
   integer(kind=4), intent(in) :: file_unit
end subroutine chm_fst_closefile

subroutine chm_fst_openfile(FILE_UNIT, file_name, file_type, options, istatus)
   integer(kind=4), intent(inout) :: file_unit
   integer(kind=4), intent  (out) :: istatus
   character(*),    intent   (in) :: file_name
   character(*),    intent   (in) :: file_type
   character(*),    intent   (in) :: options
end subroutine chm_fst_openfile

subroutine chm_getphybus_struct( )
end subroutine chm_getphybus_struct

subroutine chm_load_emissions2(F_basedir_S, gem_tstep_num, inputobj, nbvar_input)
   use inputio_mod,          only: INPUTIO_T
   character(len=*)               :: F_basedir_S  !- base path for input data file
   integer(kind=4), intent(in)    :: gem_tstep_num
   type(INPUTIO_T), intent(inout) :: inputobj
   integer(kind=4), intent(inout) :: nbvar_input
end subroutine chm_load_emissions2

subroutine chm_load_metvar(busdyn, busper, busvol, metvar2d, metvar3d)
  use chm_metvar_mod
  use chm_ptopo_grid_mod, only: chm_ni, chm_nk
 real(kind=4),    dimension(:), pointer, contiguous :: busdyn
 real(kind=4),    dimension(:), pointer, contiguous :: busper
 real(kind=4),    dimension(:), pointer, contiguous :: busvol
 real(kind=4),    intent(out) :: metvar2d(chm_ni, SIZE_MV2D)
 real(kind=4),    intent(out) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
end subroutine chm_load_metvar

subroutine chm_load_store_tracers(busdyn, chem_tr, flag)
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   integer(kind=4), intent   (in) :: flag
   real(kind=4),    dimension(:), pointer, contiguous :: busdyn
   real(kind=4),    intent(inout) :: chem_tr(chm_ni, chm_nk + 1, nb_dyn_tracers)
end subroutine chm_load_store_tracers

subroutine chm_mjrpts_get_emissions(file_unit, Fstack_emis, nb_sources, datev, &
                                    err, extra1, extra2, ex1_nom, ex2_nom)
   use chm_mjrpts_sortinfo_mod, only: nb_me_species
   integer(kind=4),            intent (in) :: file_unit, nb_sources, datev
   integer(kind=4),            intent(out) :: err
   real(kind=4),               intent(out) :: Fstack_emis(nb_sources, nb_me_species)
   real(kind=4),     optional, intent(out) :: extra1(nb_sources), extra2(nb_sources)
   character(len=4), optional, intent (in) :: ex1_nom, ex2_nom
end subroutine chm_mjrpts_get_emissions

subroutine chm_mjrpts_get_size(file_unit, numsrcs, istatus)
   integer(kind=4), intent (in) :: file_unit
   integer(kind=4), intent(out) :: numsrcs
   integer(kind=4), intent(out) :: istatus
end subroutine chm_mjrpts_get_size

subroutine chm_mjrpts_get_stkinfo(file_unit, numsrcs, istatus, &
                                  lat, lon, tem, hgt, dia, vel, typ)
   integer(kind=4),        intent (in) :: numsrcs
   integer(kind=4),        intent (in) :: file_unit
   integer(kind=4),        intent(out) :: istatus
   real(kind=4),           intent(out) :: lat(numsrcs)
   real(kind=4),           intent(out) :: lon(numsrcs)
   real(kind=4), optional, intent(out) :: tem(numsrcs)
   real(kind=4), optional, intent(out) :: hgt(numsrcs)
   real(kind=4), optional, intent(out) :: dia(numsrcs)
   real(kind=4), optional, intent(out) :: vel(numsrcs)
   real(kind=4), optional, intent(out) :: typ(numsrcs)
end subroutine chm_mjrpts_get_stkinfo

function chm_init(F_path_S) result(F_istat)
   character(len=*), intent(in) :: F_path_S !# data/tables dir
   integer(kind=4) ::  F_istat
end function chm_init

logical function chm_pkg_fields_init()
end function chm_pkg_fields_init

integer function chm_nml (F_namelist, lun_out)
   character(len=*), intent(in) :: F_namelist
   integer(kind=4),  intent(in) :: lun_out
end function chm_nml

integer function chm_set_dbg_point()
end function chm_set_dbg_point

integer function chm_cffeps_init(F_basedir_S, my_pe, numproc, rpi0_j0)
   character(len=*),                       intent(in) :: F_basedir_S
   integer(kind=4),                        intent(in) :: my_pe, numproc
   integer(kind=4), dimension(5, numproc), intent(in) :: rpi0_j0
end function chm_cffeps_init

end interface
end module chm_headers_mod
