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
   use chm_utils_mod, only: NOMV_LEN, LONG_VARNAME, DICTSTRING_LEN, chm_msg_debug

   public

#include <rmn/msg.h>
   character, parameter :: UNASSIGNED = "*"

   type :: species_info

!  Identification information
!  *_name   => Short name, often same as output name
!  *_string => String used by phymem_add to allocate memory on the bus
!  *_pvarid => phyvar variable index pvmetas(idxv)
!
!    Entry on the dynamic bus
      character (len=NOMV_LEN)       :: dyn_name   = unassigned
      character (len=DICTSTRING_LEN) :: dyn_string = unassigned
      integer(kind=4)                :: dyn_pvarid = -1
!    Entry on the permanent bus
      character (len=NOMV_LEN)       :: per_name   = unassigned
      character (len=DICTSTRING_LEN) :: per_string = unassigned
      integer(kind=4)                :: per_pvarid = -1
!    Entry on the volatile bus
      character (len=NOMV_LEN)       :: out_name   = unassigned
      character (len=DICTSTRING_LEN) :: out_string = unassigned
      integer(kind=4)                :: out_pvarid = -1
!    Entry on permanent bus the for the area emissions
      character (len=NOMV_LEN)       :: ae_name   = unassigned
      character (len=DICTSTRING_LEN) :: ae_string = unassigned
      integer(kind=4)                :: ae_pvarid = -1
!    Entry for the modulated biogenic emissions
      character (len=NOMV_LEN)       :: be_name   = unassigned
      character (len=DICTSTRING_LEN) :: be_string = unassigned
      integer(kind=4)                :: be_pvarid = -1
!    Entry on permanent bus the for (aerosol) fugitive area emissions
      character (len=NOMV_LEN)       :: fae_name   = unassigned
      character (len=DICTSTRING_LEN) :: fae_string = unassigned
      integer(kind=4)                :: fae_pvarid = -1
!    Entry on permanent bus the for the mobile area emissions
      character (len=NOMV_LEN)       :: mae_name   = unassigned
      character (len=DICTSTRING_LEN) :: mae_string = unassigned
      integer(kind=4)                :: mae_pvarid = -1
!    Entry on the permanent bus for the major point sources emissions
      character (len=NOMV_LEN)       :: me_name   = unassigned
!    Entry on the volatile bus for bidirectional flux emissions
      character (len=NOMV_LEN)       :: bd_name   = unassigned
      character (len=DICTSTRING_LEN) :: bd_string = unassigned
      integer(kind=4)                :: bd_pvarid = -1
!    Entry on the volatile bus for bidirectional flux time scale
      character (len=NOMV_LEN)       :: bdt_name   = unassigned
      character (len=DICTSTRING_LEN) :: bdt_string = unassigned
      integer(kind=4)                :: bdt_pvarid = -1
!    Entry on permanent bus for the static values of the ground emissions potential
      character (len=NOMV_LEN)       :: gep_name   = unassigned
      character (len=DICTSTRING_LEN) :: gep_string = unassigned
      integer(kind=4)                :: gep_pvarid = -1
!    Entry on permanent bus for the dynamic values of the ground emissions potential
      character (len=NOMV_LEN)       :: epd_name   = unassigned
      character (len=DICTSTRING_LEN) :: epd_string = unassigned
      integer(kind=4)                :: epd_pvarid = -1
!    Entry on volatile bus for the diagnostic atmospheric deposition potential
      character (len=NOMV_LEN)       :: epa_name   = unassigned
      character (len=DICTSTRING_LEN) :: epa_string = unassigned
      integer(kind=4)                :: epa_pvarid = -1
!    Entry on permanent bus for soil pH
      character (len=NOMV_LEN)       :: sph_name   = unassigned
      character (len=DICTSTRING_LEN) :: sph_string = unassigned
      integer(kind=4)                :: sph_pvarid = -1
!    Entry on the volatile bus for the vertical diffusion velocities
      character (len=NOMV_LEN)       :: vd_name   = unassigned
      character (len=DICTSTRING_LEN) :: vd_string = unassigned
      integer(kind=4)                :: vd_pvarid = -1
!    Entries on the volatile bus for the chemical resistances
!    Aerodynamic resistance
      character (len=NOMV_LEN)       :: ra_name   = unassigned
      character (len=DICTSTRING_LEN) :: ra_string = unassigned
      integer(kind=4)                :: ra_pvarid = -1
!    Molecular diffusion resistance
      character (len=NOMV_LEN)       :: rb_name   = unassigned
      character (len=DICTSTRING_LEN) :: rb_string = unassigned
      integer(kind=4)                :: rb_pvarid = -1
!    Total surface resistance
      character (len=NOMV_LEN)       :: rc_name   = unassigned
      character (len=DICTSTRING_LEN) :: rc_string = unassigned
      integer(kind=4)                :: rc_pvarid = -1
!    Vertical diffusion velocity for ground surface pathway
      character (len=NOMV_LEN)       :: vdg_name   = unassigned
      character (len=DICTSTRING_LEN) :: vdg_string = unassigned
      integer(kind=4)                :: vdg_pvarid = -1
!    Diagnostic dry deposition
      character (len=NOMV_LEN)       :: dd_name   = unassigned
      character (len=DICTSTRING_LEN) :: dd_string = unassigned
      integer(kind=4)                :: dd_pvarid = -1
!    Diagnostic wet deposition
      character (len=NOMV_LEN)       :: wd_name   = unassigned
      character (len=DICTSTRING_LEN) :: wd_string = unassigned
      integer(kind=4)                :: wd_pvarid = -1
!
#if defined(MACH_TENDENCIES)
! Diffusion tendencies
      character (len=NOMV_LEN)       :: td_name   = unassigned
      character (len=DICTSTRING_LEN) :: td_string = unassigned
      integer(kind=4)                :: td_pvarid = -1
! PM tendencies
      character (len=NOMV_LEN)       :: tp_name   = unassigned
      character (len=DICTSTRING_LEN) :: tp_string = unassigned
      integer(kind=4)                :: tp_pvarid = -1
! Gas chemistry tendencies
      character (len=NOMV_LEN)       :: tg_name   = unassigned
      character (len=DICTSTRING_LEN) :: tg_string = unassigned
      integer(kind=4)                :: tg_pvarid = -1
#endif

!  Molecular weight
      real  :: mol_wt = -999.0
   end type

   type(species_info), allocatable, target, save :: species_master(:)
   type(species_info), pointer            , save :: sm(:)

   integer(kind=4), save :: nb_species, nb_dyn_tracers, npidx

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

!  type :: perm_vars
!     character(len=NOMV_LEN)         :: per_name   = unassigned
!     character(len=DICTSTRING_LEN)   :: per_string = unassigned
!  end type
!  type(perm_vars), dimension(:), allocatable, save :: chem_per_vars

   external :: msg_toall

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
      call msg_toall(chm_msg_debug, 'zero_fields [BEGIN]')
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
      array(1:array_size) % bdt_name  = UNASSIGNED
      array(1:array_size) % gep_name  = UNASSIGNED
      array(1:array_size) % epd_name  = UNASSIGNED
      array(1:array_size) % epa_name  = UNASSIGNED
      array(1:array_size) % sph_name  = UNASSIGNED
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

      call msg_toall(chm_msg_debug, 'zero_fields [END]')

   end subroutine zero_fields

!============================================================================
! Description    : Print information of meta data from phybus

   integer function print_phymeta(idx, iunit) result(istat)
      use phymem,             only: phymeta, phymem_getmeta, PHY_NAMELEN
      implicit none
      integer(kind=4), intent(in)   :: idx, iunit
      type(phymeta), pointer :: vmeta

      istat = phymem_getmeta(vmeta, idx)
      if (istat < 0) then
         call msg(MSG_ERROR,'(print_phymeta) Cannot find vmeta.')
         return
      endif
      write(iunit,25) 'vname => ', vmeta%vname(1:PHY_NAMELEN)  ! vmeta%vname=wsoil
      write(iunit,35) 'ni    => ', vmeta%ni ! folded ni dim (p_runlenght)
      write(iunit,35) 'nk    => ', vmeta%nk ! number of atmospheric levels.
      write(iunit,35) 'ibus  => ', vmeta%ibus ! index of bus containing the field
      write(iunit,35) 'fmul  => ', vmeta%fmul ! number of arbitrarily-defined levels
      write(iunit,35) 'mosaic=> ', vmeta%mosaic ! number of surface sub-types
      write(iunit,35) 'size  => ', vmeta%size ! ni * nk * fmul * (mosaic+1)
      write(iunit,45) 'nlc(3)=> ', vmeta%nlcl ! local tile dims in "not folded" space
      write(iunit,35) 'idxb  => ', vmeta%idxb ! var index in the specified bus, pbuses(ibus)%meta(idxb)
      write(iunit,35) 'idxv  => ', vmeta%idxv ! var index in the specified bus, pvmetas(idxv)%meta
      write(iunit,35) 'i0    => ', vmeta%i0   ! index of first element in the bus pointer, pbuses(ibus)%bptr(i0:in,:)
      write(iunit,35) 'in    => ', vmeta%in   ! in=i0+size-1; index of first element in the bus pointer, pbuses(ibus)%bptr(i0:in,:)
      write(iunit,35) 'init  => ', vmeta%init ! 1 = init/mandatory, 0 otherwise
      write(iunit,25) 'bus   => ', vmeta%bus(1:PHY_NAMELEN) ! name of the bus containing the field
      write(iunit,25) 'iname => ', vmeta%iname(1:PHY_NAMELEN) ! input name
      write(iunit,25) 'oname => ', vmeta%oname(1:PHY_NAMELEN) ! output name
      write(iunit,25) 'sname => ', vmeta%sname(1:PHY_NAMELEN) ! series name
      write(iunit,25) 'desc  => ', trim(vmeta%desc)
      return
 25 format(A13," ",A,T32)
 35 format(A13," ",I10)
 45 format(A13," ",3(I10,", "))
   end function
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
      use phymem,             only: phymeta
      implicit none
      integer(kind=4), intent(in)   :: id, iunit
      integer(kind=4)        :: istat

      write (iunit, *) "---------------------------------------------"
      write (iunit, *) "Species # ", id
      write (iunit, *) "Molecular weight: ", species_master(id) % mol_wt

      if (species_master(id)%dyn_name /= UNASSIGNED) then
         write (iunit, *) "dyn Output name  : ", species_master(id)%dyn_name
         write (iunit, 25)"dyn String: ", species_master(id)%dyn_string
         istat = print_phymeta(species_master(id)%dyn_pvarid, iunit)
      end if

      if (species_master(id)%per_name /= UNASSIGNED) then
         write (iunit, *) "per Output name  : ", species_master(id)%per_name
         write (iunit, 25)"per String: ", species_master(id)%per_string
         istat = print_phymeta(species_master(id)%per_pvarid, iunit)
      end if

      if (species_master(id)%out_name /= UNASSIGNED) then
         write (iunit, *) "out Output name  : ", species_master(id)%out_name
         write (iunit, 25)"out String: ", species_master(id)%out_string
         istat = print_phymeta(species_master(id)%out_pvarid, iunit)
      end if

      if (species_master(id)%ae_name /= UNASSIGNED) then
         write (iunit, *) "ae  Output name  : ", species_master(id)%ae_name
         write (iunit, 25)"ae  String: ", species_master(id)%ae_string
         istat = print_phymeta(species_master(id)%ae_pvarid, iunit)
      end if

      if (species_master(id)%fae_name /= UNASSIGNED) then
         write (iunit, *) "fae  Output name  : ", species_master(id)%fae_name
         write (iunit, 25)"fae  String: ", species_master(id)%fae_string
         istat = print_phymeta(species_master(id)%fae_pvarid, iunit)
      end if

      if (species_master(id)%mae_name /= UNASSIGNED) then
         write (iunit, *) "mae  Output name  : ", species_master(id)%mae_name
         write (iunit, 25)"mae  String: ", species_master(id)%mae_string
         istat = print_phymeta(species_master(id)%mae_pvarid, iunit)
      end if

      if (species_master(id)%be_name /= UNASSIGNED) then
         write (iunit, *) "be  Output name  : ", species_master(id)%be_name
         write (iunit, 25)"be  String: ", species_master(id)%be_string
         istat = print_phymeta(species_master(id)%be_pvarid, iunit)
      end if

      if (species_master(id)%me_name /= UNASSIGNED) then
         write (iunit, *) "me  Output name  : ", species_master(id)%me_name
      end if

      if (species_master(id)%bd_name /= UNASSIGNED) then
         write (iunit, *) "bd  Output name  : ", species_master(id)%bd_name
         write (iunit, 25)"bd  String: ", species_master(id)%bd_string
         istat = print_phymeta(species_master(id)%bd_pvarid, iunit)
      end if

      if (species_master(id)%bdt_name /= UNASSIGNED) then
         write (iunit, *) "bdt  Output name  : ", species_master(id)%bdt_name
         write (iunit, 25)"bdt  String: ", species_master(id)%bdt_string
         istat = print_phymeta(species_master(id)%bdt_pvarid, iunit)
      end if

      if (species_master(id)%gep_name /= UNASSIGNED) then
         write (iunit, *) "gep  Output name  : ", species_master(id)%gep_name
         write (iunit, 25)"gep  String: ", species_master(id)%gep_string
         istat = print_phymeta(species_master(id)%gep_pvarid, iunit)
      end if

      if (species_master(id)%epd_name /= UNASSIGNED) then
         write (iunit, *) "epd  Output name  : ", species_master(id)%epd_name
         write (iunit, 25)"epd  String: ", species_master(id)%epd_string
         istat = print_phymeta(species_master(id)%epd_pvarid, iunit)
      end if

      if (species_master(id)%epa_name /= UNASSIGNED) then
         write (iunit, *) "epa  Output name  : ", species_master(id)%epa_name
         write (iunit, 25)"epa  String: ", species_master(id)%epa_string
         istat = print_phymeta(species_master(id)%epa_pvarid, iunit)
      end if

      if (species_master(id)%sph_name /= UNASSIGNED) then
         write (iunit, *) "sph  Output name  : ", species_master(id)%sph_name
         write (iunit, 25)"sph  String: ", species_master(id)%sph_string
         istat = print_phymeta(species_master(id)%sph_pvarid, iunit)
      end if

      if (species_master(id)%vd_name /= UNASSIGNED) then
         write (iunit, *) "vd  Output name  : ", species_master(id)%vd_name
         write (iunit, 25)"vd  String: ", species_master(id)%vd_string
         istat = print_phymeta(species_master(id)%vd_pvarid, iunit)
      end if

      if (species_master(id)%vdg_name /= UNASSIGNED) then
         write (iunit, *) "vdg  Output name  : ", species_master(id)%vdg_name
         write (iunit, 25)"vdg  String: ", species_master(id)%vdg_string
         istat = print_phymeta(species_master(id)%vdg_pvarid, iunit)
      end if

      if (species_master(id)%ra_name /= UNASSIGNED) then
         write (iunit, *) "ra  Output name  : ", species_master(id)%ra_name
         write (iunit, 25)"ra  String: ", species_master(id)%ra_string
         istat = print_phymeta(species_master(id)%ra_pvarid, iunit)
      end if

      if (species_master(id)%rb_name /= UNASSIGNED) then
         write (iunit, *) "rb  Output name  : ", species_master(id)%rb_name
         write (iunit, 25)"rb  String: ", species_master(id)%rb_string
         istat = print_phymeta(species_master(id)%rb_pvarid, iunit)
      end if

      if (species_master(id)%rc_name /= UNASSIGNED) then
         write (iunit, *) "rc  Output name  : ", species_master(id)%rc_name
         write (iunit, 25)"rc  String: ", species_master(id)%rc_string
         istat = print_phymeta(species_master(id)%rc_pvarid, iunit)
      end if

      if (species_master(id)%dd_name /= UNASSIGNED) then
         write (iunit, *) "dd  Output name  : ", species_master(id)%dd_name
         write (iunit, 25)"dd  String: ", species_master(id)%dd_string
         istat = print_phymeta(species_master(id)%dd_pvarid, iunit)
      end if

      if (species_master(id)%wd_name /= UNASSIGNED) then
         write (iunit, *) "wd  Output name  : ", species_master(id)%wd_name
         write (iunit, 25)"wd  String: ", species_master(id)%wd_string
         istat = print_phymeta(species_master(id)%wd_pvarid, iunit)
      end if

#if defined(MACH_TENDENCIES)
      if (species_master(id)%td_name /= UNASSIGNED) then
         write (iunit, *) "td  Output name  : ", species_master(id)%td_name
         write (iunit, 25)"td  String: ", species_master(id)%td_string
         istat = print_phymeta(species_master(id)%td_pvarid, iunit)
      end if

      if (species_master(id)%tp_name /= UNASSIGNED) then
         write (iunit, *) "tp  Output name  : ", species_master(id)%tp_name
         write (iunit, 25)"tp  String: ", species_master(id)%tp_string
         istat = print_phymeta(species_master(id)%tp_pvarid, iunit)
      end if

      if (species_master(id)%tg_name /= UNASSIGNED) then
         write (iunit, *) "tg  Output name  : ", species_master(id)%tg_name
         write (iunit, 25)"tg  String: ", species_master(id)%tg_string
         istat = print_phymeta(species_master(id)%td_pvarid, iunit)
      end if
#endif
      write (iunit, *) "---------------------------------------------"
      write (iunit, *) ""

 25 format(A13," ",A100)
   end subroutine print_species_info

!============================================================================
! Name           : print_all_species_info
!
! Description    : Print the information of all the species in the
!                  species_master array
!
! Arguments      : unit
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

      call msg_toall(chm_msg_debug, 'print_all_species_info [BEGIN]')

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

      write (iunit, *) "There are ", nb_species, " global species (nb_species)"
      write (iunit, *) "There are ", npidx, " variables on phybus (npidx)"
      write (iunit, *) "There are ", nb_dyn_tracers, " dynamic tracers (nb_dyn_tracers)"
      do i = 1, nb_species
         call print_species_info(i, iunit)
      end do

      call msg_toall(chm_msg_debug, 'print_all_species_info [END]')

   end subroutine print_all_species_info

end module chm_species_info_mod

