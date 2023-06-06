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

module mach_pkg_tendencies_mod
#if defined(MACH_TENDENCIES)
   use chm_utils_mod,        only: chm_lun_out, global_debug, NOMV_LEN, LONG_VARNAME
   use chm_species_idx_mod
   use chm_species_info_mod, only: sm, UNASSIGNED, nb_dyn_tracers
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk

   save

   integer(kind=4), parameter :: TEND_DIFFUSION  = 1   ! define for diffusion
   integer(kind=4), parameter :: TEND_GASCHEM    = 2   ! define for gas chemistry
   integer(kind=4), parameter :: TEND_PMCHEM     = 3   ! define for PM chemistry

   contains

!============================================================================
!>
!! Name           : pkg_tendencies_metainit
!!
!! Description    : Initialize the meta information for each species
!!
!! Arguments:  None
!============================================================================
   subroutine pkg_tendencies_metainit()
      implicit none

      integer(kind=4) :: i

      character (len = NOMV_LEN)         :: var_name
      character (len = LONG_VARNAME)     :: lname
      logical(kind=4)                    :: local_dbg

      local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

      if (local_dbg) then
         write (chm_lun_out, *) "Entering pkg_tendencies_metainit"
      end if

!         Since we're dealing with tracers, we replace the first letter T by the
!         prefix that represents the operation, there's no garanty that the forged name
!         isn't use elsewhere so be advised.
!         U = Diffusion (d was already taken)
!         G = gas phase chemistry
!         L = pm chemistry

      do i = 1, nb_dyn_tracers
         var_name = sm(i) % dyn_name
!        Diffusion
         sm(i) % td_name= 'U' // var_name(2:4)
         lname = var_name // "_tend_diff"
         sm(i) % td_string = 'VN=' // trim(lname) // ';ON=' // sm(i) % td_name // '; VD=' // "Tendencies of diffusion for " // var_name // '; VS=T; VB=v0'
!        Gas
         sm(i) % tg_name = 'G' // var_name(2:4)
         lname = var_name // "_tend_gas"
         sm(i) % tg_string = 'VN=' // trim(lname) // ';ON=' // sm(i) % tg_name // '; VD=' // "Tendencies of gas phase chem for " // var_name // '; VS=T; VB=v0'
!        PM
         sm(i) % tp_name = 'L' // var_name(2:4)
         lname = var_name // "_tend_pm"
         sm(i) % tp_string = 'VN=' // trim(lname) // ';ON=' // sm(i) % tp_name // '; VD=' // "Tendencies of pm chem for " // var_name // '; VS=T; VB=v0'
      end do

      if (local_dbg) then
         write (chm_lun_out) "Leaving pkg_tendencies_metainit"
      end if

   end subroutine pkg_tendencies_metainit


   subroutine tendency_store(process_id, chem_tr, busvol)
!! Description    : Store the value of all species' concentrations in
!!                  order to diff it after the process' completion.
!!                  Basically: tend = conc
!
!! @param IN process_id defines which of the defined process is involved
!! @param IN chem_tr array which holds the current tracers concentrations
!! @param INOUT busvol volatile bus
!============================================================================
      integer(kind=4), intent(in)    :: process_id
      real(kind=4)   , intent(in)    :: chem_tr(chm_ni, chm_nk+1, nb_dyn_tracers)
      real(kind=4)   , dimension(:), pointer, contiguous :: busvol

      integer(kind=4) i, k, indx
      logical(kind=4) local_dbg

      local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

      if (local_dbg) then
         write (chm_lun_out, *) "Entering tendency_store"
      end if

      select case (process_id)
         case (TEND_DIFFUSION)
            do sp = 1, nb_dyn_tracers
               indx = 0
               do k = 1, chm_nk
                  do i = 1, chm_ni
                     busvol(sm(sp) % td_offset + indx) = chem_tr(i, k, sp)
                     indx = indx + 1
                  end do
               end do
            end do
         case (TEND_GASCHEM)
            do sp = 1, nb_dyn_tracers
               indx = 0
               do k = 1, chm_nk
                  do i = 1, chm_ni
                     busvol(sm(sp) % tg_offset + indx) = chem_tr(i, k, sp)
                     indx = indx + 1
                  end do
               end do
            end do
         case (TEND_PMCHEM)
            do sp = 1, nb_dyn_tracers
               indx = 0
               do k = 1, chm_nk
                  do i = 1, chm_ni
                     busvol(sm(sp) % tp_offset + indx) = chem_tr(i, k, sp)
                     indx = indx + 1
                  end do
               end do
            end do
         case default
            print *, "tendency_store: The case that should not be."

      end select

      if (local_dbg) then
         write (chm_lun_out, *) "Leaving tendency_store"
      end if

   end subroutine tendency_store


   subroutine tendency_delta(process_id, chem_tr, busvol)
!============================================================================
!! Description    : Store the value of all species' concentrations in
!!                  order to diff it after the process' completion.
!!                  Basically: tend = conc - tend
!
!! @param IN process_id defines which of the defined process is involved
!! @param IN chem_tr array which holds the current tracers concentrations
!! @param INOUT busvol volatile bus
!============================================================================
      integer(kind=4), intent(in)    :: process_id
      real(kind=4)   , intent(in)    :: chem_tr(chm_ni, chm_nk+1, nb_dyn_tracers)
      real(kind=4)   , dimension(:), pointer, contiguous :: busvol

      integer(kind=4) i, k, indx
      logical(kind=4) local_dbg

      local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))

      if (local_dbg) then
         write (chm_lun_out, *) "Entering tendency_delta"
      end if

      select case (process_id)
         case (TEND_DIFFUSION)
            do sp = 1, nb_dyn_tracers
               indx = 0
               do k = 1, chm_nk
                  do i = 1, chm_ni
                     busvol(sm(sp) % td_offset + indx) = chem_tr(i, k, sp) - &
                                              busvol(sm(sp) % td_offset + indx)
                     indx = indx + 1
                  end do
               end do
            end do
         case (TEND_GASCHEM)
            do sp = 1, nb_dyn_tracers
               indx = 0
               do k = 1, chm_nk
                  do i = 1, chm_ni
                     busvol(sm(sp) % tg_offset + indx) = chem_tr(i, k, sp) - &
                                              busvol(sm(sp) % tg_offset + indx)
                     indx = indx + 1
                  end do
               end do
            end do
         case (TEND_PMCHEM)
            do sp = 1, nb_dyn_tracers
               indx = 0
               do k = 1, chm_nk
                  do i = 1, chm_ni
                     busvol(sm(sp) % tp_offset + indx) = chem_tr(i, k, sp) - &
                                              busvol(sm(sp) % tp_offset + indx)
                     indx = indx + 1
                  end do
               end do
            end do
         case default
            print *, "tendency_delta: The case that should not be."
      end select


      if (local_dbg) then
         write (chm_lun_out, *) "Leaving tendency_delta"
      end if
   end subroutine tendency_delta


#else

! Here there must be a stub version of all subroutines declared above
!    subroutine pkg_tendencies_metainit()
!          implicit none
!          print *,"mach_pkg_tendencies_mod.pkg_tendencies_metainit stub"
!    end subroutine pkg_tendencies_metainit
!
!    subroutine tendency_store(process_id)
!       implicit none
!       print *,"mach_pkg_tendencies_mod.tendency_store stub"
!    end subroutine tendency_store
!
!    subroutine tendency_delta(process_id)
!       implicit none
!       print *,"mach_pkg_tendencies_mod.tendency_store stub"
!    end subroutine tendency_delta

#endif

end module mach_pkg_tendencies_mod
