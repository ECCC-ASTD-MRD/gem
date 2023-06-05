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
! Fichier/File   : chm_species_info_mod.ftn90
! Creation       : H. Landry, Decembre 2007
! Description    : Modules defining species meta-information
!
! Extra info     :
!
!============================================================================

module chm_species_info_mod
   use chm_utils_mod, only: NOMV_LEN, LONG_VARNAME, DICTSTRING_LEN

   public 
   
   character, parameter :: UNASSIGNED = "*"

   type :: species_info

!  Identification information
!  *_name   => Short name, often same as output name
!  *_string => String used by gesdict to allocate memory on the bus
!  *_offset => Starting index on the physics bus
!
!    Entry on the dynamic bus
      character (len=NOMV_LEN)       :: dyn_name   = unassigned
      character (len=DICTSTRING_LEN) :: dyn_string = unassigned
      integer(kind=4)                :: dyn_offset = -1

!    Entry on the permanent bus
      character (len=NOMV_LEN)       :: per_name   = unassigned
      character (len=DICTSTRING_LEN) :: per_string = unassigned
      integer(kind=4)                :: per_offset = -1

!    Entry on the volatile bus
      character (len=NOMV_LEN)       :: out_name   = unassigned
      character (len=DICTSTRING_LEN) :: out_string = unassigned
      integer(kind=4)                :: out_offset = -1

!    Entry on permanent bus the for the area emissions
      character (len=NOMV_LEN)       :: ae_name   = unassigned
      character (len=DICTSTRING_LEN) :: ae_string = unassigned
      integer(kind=4)                :: ae_offset = -1

!    Entry for the modulated biogenic emissions
      character (len=NOMV_LEN)       :: be_name   = unassigned
      character (len=DICTSTRING_LEN) :: be_string = unassigned
      integer(kind=4)                :: be_offset = -1

!    Entry on permanent bus the for (aerosol) fugitive area emissions
      character (len=NOMV_LEN)       :: fae_name   = unassigned
      character (len=DICTSTRING_LEN) :: fae_string = unassigned
      integer(kind=4)                :: fae_offset = -1

!    Entry on permanent bus the for the mobile area emissions
      character (len=NOMV_LEN)       :: mae_name   = unassigned
      character (len=DICTSTRING_LEN) :: mae_string = unassigned
      integer(kind=4)                :: mae_offset = -1

!    Entry on the permanent bus for the major point sources emissions
      character (len=NOMV_LEN)       :: me_name   = unassigned

!    Entry on the volatile bus for bidirectional flux emissions
      character (len=NOMV_LEN)       :: bd_name   = unassigned
      character (len=DICTSTRING_LEN) :: bd_string = unassigned
      integer(kind=4)                :: bd_offset = -1

!    Entry on permanent bus the for the ground emissions potential
      character (len=NOMV_LEN)       :: gep_name   = unassigned
      character (len=DICTSTRING_LEN) :: gep_string = unassigned
      integer(kind=4)                :: gep_offset = -1 
      
!    Entry on the volatile bus for the vertical diffusion velocities
      character (len=NOMV_LEN)       :: vd_name   = unassigned
      character (len=DICTSTRING_LEN) :: vd_string = unassigned
      integer(kind=4)                :: vd_offset = -1

!    Entries on the volatile bus for the chemical resistances
!    Aerodynamic resistance
      character (len=NOMV_LEN)       :: ra_name   = unassigned
      character (len=DICTSTRING_LEN) :: ra_string = unassigned
      integer(kind=4)                :: ra_offset = -1
!    Molecular diffusion resistance
      character (len=NOMV_LEN)       :: rb_name   = unassigned
      character (len=DICTSTRING_LEN) :: rb_string = unassigned
      integer(kind=4)                :: rb_offset = -1
!    Total surface resistance
      character (len=NOMV_LEN)       :: rc_name   = unassigned
      character (len=DICTSTRING_LEN) :: rc_string = unassigned
      integer(kind=4)                :: rc_offset = -1
!    Vertical diffusion velocity for ground surface pathway
      character (len=NOMV_LEN)       :: vdg_name   = unassigned
      character (len=DICTSTRING_LEN) :: vdg_string = unassigned
      integer(kind=4)                :: vdg_offset = -1
!    Diagnostic dry deposition
      character (len=NOMV_LEN)       :: dd_name   = unassigned
      character (len=DICTSTRING_LEN) :: dd_string = unassigned
      integer(kind=4)                :: dd_offset = -1
!    Diagnostic wet deposition
      character (len=NOMV_LEN)       :: wd_name   = unassigned
      character (len=DICTSTRING_LEN) :: wd_string = unassigned
      integer(kind=4)                :: wd_offset = -1
!      
#if defined(MACH_TENDENCIES)
! Diffusion tendencies
      character (len=NOMV_LEN)       :: td_name   = unassigned
      character (len=DICTSTRING_LEN) :: td_string = unassigned
      integer(kind=4)                :: td_offset = -1
! PM tendencies
      character (len=NOMV_LEN)       :: tp_name   = unassigned
      character (len=DICTSTRING_LEN) :: tp_string = unassigned
      integer(kind=4)                :: tp_offset = -1
! Gas chemistry tendencies
      character (len=NOMV_LEN)       :: tg_name   = unassigned
      character (len=DICTSTRING_LEN) :: tg_string = unassigned
      integer(kind=4)                :: tg_offset = -1
#endif

!  Molecular weight
      real  :: mol_wt = -999.0

   end type

   type(species_info), allocatable, target, save :: species_master(:)
   type(species_info), pointer            , save :: sm(:)

   integer(kind=4), save :: nb_species, nb_dyn_tracers

#if defined(MACH_TENDENCIES)
   type :: tendency_idx
      integer(kind=4) process ! chemical process id
      integer(kind=4) tend  ! idx of the tendency in species_master
      integer(kind=4) spc   ! idx of the original species in species_master
   end type
   
   integer(kind=4), save :: nb_tracers
#endif

   type :: ent_vars
      character(len=LONG_VARNAME)     :: ent_name   = unassigned
      character(len=DICTSTRING_LEN)   :: ent_string = unassigned
   end type
   type(ent_vars), dimension(3), save :: chem_ent_vars
! Note that the dimension of chem_ent_vars is currently hard-coded
   integer(kind=4) :: ent_vars_num = 0

   type :: perm_vars
      character(len=NOMV_LEN)         :: per_name   = unassigned
      character(len=DICTSTRING_LEN)   :: per_string = unassigned
      integer(kind=4)                 :: per_offset = -1      
   end type
   type(perm_vars), dimension(:), allocatable, save :: chem_per_vars

   contains
!============================================================================
! Name           : zero_fields
!
! Description    : Reset content of an array of type species_info
!
! Arguments:  OUT
!                 array -> the array of species_info structure to reset
!
!              IN
!                 array_size -> size of array
!
!============================================================================
   subroutine zero_fields(array, array_size)
      implicit none
      integer(kind=4)   , intent   (in) :: array_size
      type(species_info), intent(inout) :: array(array_size)
!
      array(1:array_size) % mol_wt   = -999.0
      array(1:array_size) % dyn_name = UNASSIGNED
      array(1:array_size) % per_name = UNASSIGNED
      array(1:array_size) % out_name = UNASSIGNED
      array(1:array_size) % ae_name  = UNASSIGNED
      array(1:array_size) % be_name  = UNASSIGNED
      array(1:array_size) % fae_name = UNASSIGNED
      array(1:array_size) % mae_name = UNASSIGNED
      array(1:array_size) % me_name  = UNASSIGNED
      array(1:array_size) % bd_name  = UNASSIGNED
      array(1:array_size) % gep_name  = UNASSIGNED
      array(1:array_size) % vd_name  = UNASSIGNED
      array(1:array_size) % vdg_name = UNASSIGNED
      array(1:array_size) % ra_name  = UNASSIGNED
      array(1:array_size) % rb_name  = UNASSIGNED
      array(1:array_size) % rc_name  = UNASSIGNED
      array(1:array_size) % dd_name  = UNASSIGNED
      array(1:array_size) % wd_name  = UNASSIGNED
#if defined(MACH_TENDENCIES)
      array(1:array_size) % td_name  = UNASSIGNED
      array(1:array_size) % tg_name  = UNASSIGNED
      array(1:array_size) % tp_name  = UNASSIGNED
#endif

   end subroutine zero_fields
!
!============================================================================
! Name           : print_species_info
!
! Description    : Print information of one species_info structure
!
! Arguments:   IN
!                 id -> the index of the structure in species_master
!
!============================================================================
!
   subroutine print_species_info(id, iunit)
      implicit none
      integer(kind=4), intent(in)   :: id, iunit

      write (iunit, *) "---------------------------------------------"
      write (iunit, *) "Species # ", id
      write (iunit, *) "Molecular weight: ", species_master(id) % mol_wt


      if (species_master(id) % dyn_name /= UNASSIGNED) then
         write (iunit, *) "dyn Output name  : ", species_master(id) % dyn_name
         write (iunit, *) "dyn String       : ", species_master(id) % dyn_string
         write (iunit, *) "dyn Offset       : ", species_master(id) % dyn_offset
      end if

      if (species_master(id) % per_name /= UNASSIGNED) then
         write (iunit, *) "per Output name  : ", species_master(id) % per_name
         write (iunit, *) "per String       : ", species_master(id) % per_string
         write (iunit, *) "per Offset       : ", species_master(id) % per_offset
      end if

      if (species_master(id) % out_name /= UNASSIGNED) then
         write (iunit, *) "out Output name  : ", species_master(id) % out_name
         write (iunit, *) "out String       : ", species_master(id) % out_string
         write (iunit, *) "out Offset       : ", species_master(id) % out_offset
      end if

      if (species_master(id) % ae_name /= UNASSIGNED) then
         write (iunit, *) "ae  Output name  : ", species_master(id) % ae_name
         write (iunit, *) "ae  String       : ", species_master(id) % ae_string
         write (iunit, *) "ae  Offset       : ", species_master(id) % ae_offset
      end if

      if (species_master(id) % fae_name /= UNASSIGNED) then
         write (iunit, *) "fae  Output name  : ", species_master(id) % fae_name
         write (iunit, *) "fae  String       : ", species_master(id) % fae_string
         write (iunit, *) "fae  Offset       : ", species_master(id) % fae_offset
      end if

      if (species_master(id) % mae_name /= UNASSIGNED) then
         write (iunit, *) "mae  Output name  : ", species_master(id) % mae_name
         write (iunit, *) "mae  String       : ", species_master(id) % mae_string
         write (iunit, *) "mae  Offset       : ", species_master(id) % mae_offset
      end if

      if (species_master(id) % be_name /= UNASSIGNED) then
         write (iunit, *) "be  Output name  : ", species_master(id) % be_name
         write (iunit, *) "be  String       : ", species_master(id) % be_string
         write (iunit, *) "be  Offset       : ", species_master(id) % be_offset
      end if

      if (species_master(id) % me_name /= UNASSIGNED) then
         write (iunit, *) "me  Output name  : ", species_master(id) % me_name
      end if

      if (species_master(id) % bd_name /= UNASSIGNED) then
         write (iunit, *) "bd  Output name  : ", species_master(id) % bd_name
         write (iunit, *) "bd  String       : ", species_master(id) % bd_string
         write (iunit, *) "bd  Offset       : ", species_master(id) % bd_offset
      end if

      if (species_master(id) % gep_name /= UNASSIGNED) then
         write (iunit, *) "gep  Output name  : ", species_master(id) % gep_name
         write (iunit, *) "gep  String       : ", species_master(id) % gep_string
         write (iunit, *) "gep  Offset       : ", species_master(id) % gep_offset
      end if
      
      if (species_master(id) % vd_name /= UNASSIGNED) then
         write (iunit, *) "vd  Output name  : ", species_master(id) % vd_name
         write (iunit, *) "vd  String       : ", species_master(id) % vd_string
         write (iunit, *) "vd  Offset       : ", species_master(id) % vd_offset
      end if

      if (species_master(id) % vdg_name /= UNASSIGNED) then
         write (iunit, *) "vdg  Output name  : ", species_master(id) % vdg_name
         write (iunit, *) "vdg  String       : ", species_master(id) % vdg_string
         write (iunit, *) "vdg  Offset       : ", species_master(id) % vdg_offset
      end if

      if (species_master(id) % ra_name /= UNASSIGNED) then
         write (iunit, *) "ra  Output name  : ", species_master(id) % ra_name
         write (iunit, *) "ra  String       : ", species_master(id) % ra_string
         write (iunit, *) "ra  Offset       : ", species_master(id) % ra_offset
      end if

      if (species_master(id) % rb_name /= UNASSIGNED) then
         write (iunit, *) "rb  Output name  : ", species_master(id) % rb_name
         write (iunit, *) "rb  String       : ", species_master(id) % rb_string
         write (iunit, *) "rb  Offset       : ", species_master(id) % rb_offset
      end if

      if (species_master(id) % rc_name /= UNASSIGNED) then
         write (iunit, *) "rc  Output name  : ", species_master(id) % rc_name
         write (iunit, *) "rc  String       : ", species_master(id) % rc_string
         write (iunit, *) "rc  Offset       : ", species_master(id) % rc_offset
      end if

      if (species_master(id) % dd_name /= UNASSIGNED) then
         write (iunit, *) "dd  Output name  : ", species_master(id) % dd_name
         write (iunit, *) "dd  String       : ", species_master(id) % dd_string
         write (iunit, *) "dd  Offset       : ", species_master(id) % dd_offset
      end if

      if (species_master(id) % wd_name /= UNASSIGNED) then
         write (iunit, *) "wd  Output name  : ", species_master(id) % wd_name
         write (iunit, *) "wd  String       : ", species_master(id) % wd_string
         write (iunit, *) "wd  Offset       : ", species_master(id) % wd_offset
      end if
         
#if defined(MACH_TENDENCIES)
      if (species_master(id) % td_name /= UNASSIGNED) then
         write (iunit, *) "td  Output name  : ", species_master(id) % td_name
         write (iunit, *) "td  String       : ", species_master(id) % td_string
         write (iunit, *) "td  Offset       : ", species_master(id) % td_offset
      end if
         
      if (species_master(id) % tp_name /= UNASSIGNED) then
         write (iunit, *) "tp  Output name  : ", species_master(id) % tp_name
         write (iunit, *) "tp  String       : ", species_master(id) % tp_string
         write (iunit, *) "tp  Offset       : ", species_master(id) % tp_offset
      end if
         
      if (species_master(id) % tg_name /= UNASSIGNED) then
         write (iunit, *) "tg  Output name  : ", species_master(id) % tg_name
         write (iunit, *) "tg  String       : ", species_master(id) % tg_string
         write (iunit, *) "tg  Offset       : ", species_master(id) % tg_offset
      end if
#endif

      write (iunit, *) "---------------------------------------------"
      write (iunit, *) ""

   end subroutine print_species_info

!============================================================================
! Name           : print_species_info
!
! Description    : Print the information of all the species in the
!                  species_master array
!
! Arguments:   None
!
!============================================================================

   subroutine print_all_species_info(iunit_opt)
      use chm_utils_mod
      implicit none
      
      integer(kind=4) :: iunit_opt
!      integer(kind=4), optional, intent(in) :: unit_opt

      integer(kind=4) i, iunit
!      if(present(iunit_opt)) then
         iunit = iunit_opt
!      else
!         iunit = chm_lun_out
!      end if

      write (iunit, *) "           PRINT               "   
      write (iunit, *) "                             ,,"
      write (iunit, *) "                         ';;   "
      write (iunit, *) "                          ''   "
      write (iunit, *) "            ____          ||   "
      write (iunit, *) "           ;    \         ||   "
      write (iunit, *) "            \,---'-,-,    ||   "
      write (iunit, *) "            /     (  o)   ||   "
      write (iunit, *) "          (o )__,--'-' \  ||   "
      write (iunit, *) ",,,,       ;'uuuuu''   ) ;;    "
      write (iunit, *) "\   \      \ )      ) /\//     "
      write (iunit, *) " '--'       \'nnnnn' /  \      "
      write (iunit, *) "   \\      //'------'    \     "
      write (iunit, *) "    \\    //  \           \    "
      write (iunit, *) "     \\  //    )           )   "
      write (iunit, *) "      \\//     |           |   "
      write (iunit, *) "       \\     /            |   "
      write (iunit, *) "       ALL THE SPECIES !!!     "   
      
      write (iunit, *) "There are ", nb_species, " global species (without the different bins)"
      do i = 1, nb_species
         call print_species_info(i, iunit)
      end do

   end subroutine print_all_species_info

end module chm_species_info_mod

