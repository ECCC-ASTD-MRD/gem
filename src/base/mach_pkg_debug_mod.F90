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


module mach_pkg_debug_mod
   use chm_utils_mod,        only: NOMV_LEN, LONG_VARNAME, DICTSTRING_LEN, MAX_DEBUG_VAR, pre_increment
   use chm_species_idx_mod
   use chm_species_info_mod, only: species_master

   save

   contains
!============================================================================
! Name           : pkg_debug_idxinit
!
! Description    : Initialize the indices of debug fields.  Starts with
!                  the index passed as a dummy argument.
!
! Arguments:  IN/OUT
!                idx  -> index to start from
!============================================================================
      integer function pkg_debug_idxinit(idx, chm_debug_2d_i, chm_debug_3d_i)
         implicit none

         integer(kind=4), intent(inout) :: idx ! the index from where to start
         integer(kind=4), intent(in)    :: chm_debug_2d_i, chm_debug_3d_i

         integer(kind=4) start_index, nb_fields, i

         start_index = idx

         do i = 1, chm_debug_2d_i
            dbg_2D(i) = pre_increment(idx)
         end do

         do i = 1, chm_debug_3d_i
            dbg_3D(i) = pre_increment(idx)
         end do

         nb_fields = idx - start_index
         pkg_debug_idxinit = nb_fields

      end function pkg_debug_idxinit

!============================================================================
! Name           : pkg_debug_metainit
!
! Description    : Initialize the meta information for each species
!
! Arguments:  None
!============================================================================
      subroutine pkg_debug_metainit(chm_debug_2d_i, chm_debug_3d_i)
         implicit none

         integer(kind=4), intent(in) :: chm_debug_2d_i, chm_debug_3d_i

         integer(kind=4) i

         character (len = NOMV_LEN)       :: var_name
         character (len = LONG_VARNAME)   :: var_lname
         character (len = DICTSTRING_LEN) :: var_string
         character (len = 1) , dimension(MAX_DEBUG_VAR) :: suffix

         suffix = (/'1', '2' ,'3', '4', '5', '6', '7', '8', '9', &
                    'A', 'B' ,'C' ,'D', 'E', 'F', 'G', 'H', 'I'/)

         do i = 1, chm_debug_2d_i

            var_name   = "2DB" // suffix(i)
            var_lname  = "DEBUG_2D_" // suffix(i)
            var_string = 'VN=' // trim(var_lname) // ';ON=' // var_name // '; VD=Debug 2D variable ' // suffix(i) // '; VS=A; VB=v0'
            print *, trim(var_string)

            species_master(dbg_2D(i)) % out_name   = var_name
            species_master(dbg_2D(i)) % out_string = trim(var_string)

         end do

         do i = 1, chm_debug_3d_i

            var_name   = "3DB" // suffix(i)
            var_lname  = "DEBUG_3D_" // suffix(i)
            var_string = "VN=" // trim(var_lname) // "; ON=" // var_name // "; VD=Debug 3D variable " // suffix(i) // "; VS=T ;VB=v0"

            species_master(dbg_3D(i)) % out_name   = var_name
            species_master(dbg_3D(i)) % out_string = trim(var_string)

         end do

      end subroutine pkg_debug_metainit

end module mach_pkg_debug_mod

