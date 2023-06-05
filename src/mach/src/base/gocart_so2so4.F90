!---------------------------------- licence begin -------------------------------
! gem-mach - atmospheric chemistry library for the gem numerical atmospheric model
! copyright (c) 2007-2013 - air quality research division &
!                           national prediction operations division
!                           environnement canada
! this library is free software; you can redistribute it and/or
! modify it under the terms of the gnu lesser general public
! license as published by the free software foundation; either
! version 2.1 of the license, or (at your option) any later version.
!
! this library is distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose.  see the gnu
! lesser general public license for more details.
!
! you should have received a copy of the gnu lesser general public
! license along with this library; if not, write to the free software
! foundation, inc., 51 franklin street, fifth floor, boston, ma  02110-1301  usa
!---------------------------------- licence end ---------------------------------

!============================================================================!
!         environnement canada         |        environment canada           !
!                                      |                                     !
! - service meteorologique du canada   | - meteorological service of canada  !
! - direction generale des sciences    | - science and technology branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! projet/project : gem-mach
! fichier/file   : gocart_so2so4.ftn90
! creation       : from wrfchem subroutine so2so4 (GOCART option)
! description    : simple oxidation parameterization for so2 in absence of full
!                  aerosol module
!
!
!============================================================================
!
!!if_on
subroutine gocart_so2so4(chem_tr, metvar3d)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV3D
!!if_off
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_H2O2, sp_SO2, sp_SO4
   use chm_consphychm_mod,   only: rgasi, rgasd, delta
   use chm_metvar_mod,       only: MV3D_HUPLUS, MV3D_TPLUS, MV3D_FTOT
   implicit none
!
!  declaration of subroutine arguments
!
!!if_on
   real(kind=4), intent(inout) :: chem_tr (chm_ni, chm_nk+1, nb_dyn_tracers)
   real(kind=4), intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
!!if_off
!
!  declaration of local variables
!
   integer(kind=4) :: i, k
   real   (kind=4) :: tc2, tc3, h2o2, cldf, conv1
!
! CODE BEGINS
!
   do k = 1, chm_nk
      do i = 1, chm_ni
         cldf = metvar3d(i, k, MV3D_FTOT)
         conv1 = rgasi / (rgasd * (1.0 + delta * metvar3d(i, k, MV3D_HUPLUS)))
         tc2 = chem_tr(i, k, sp_SO2) * conv1 / sm(sp_SO2) % mol_wt
         if (cldf > 0.0 .and. tc2 > 0.0 .and. &
             metvar3d(i, k, MV3D_TPLUS) > 258.0) then
!            tc3 = chem_tr(i, k, sp_SO4) * conv1 / sm(sp_SO4 ) % mol_wt
            h2o2 = chem_tr(i, k, sp_H2O2) * conv1 / sm(sp_H2O2) % mol_wt
!
! ****************************************************************************
! *  Update SO2 concentration after cloud chemistry.                         *
! *  SO2 chemical loss rate  = SO4 production rate (MixingRatio/timestep).   *
! ****************************************************************************
           ! Cloud chemistry (above 258K): 

            if (tc2 > h2o2) then
               cldf = cldf * (h2o2 / tc2)
               h2o2 = h2o2 * (1.0 - cldf)
            else
               h2o2 = h2o2 * (1.0 - cldf * tc2 / h2o2)
            end if
!            tc3  = max(1.e-16,(tc3 + tc2*cldf))
            tc3  = 1.e-16
            tc2  = max(1.e-16, tc2 * (1.0 - cldf))
            h2o2 = max(1.e-16, h2o2)
            chem_tr(i, k, sp_SO2)  = tc2 / conv1 * sm(sp_SO2 ) % mol_wt
            chem_tr(i, k, sp_SO4)  = tc3 / conv1 * sm(sp_SO4 ) % mol_wt
            chem_tr(i, k, sp_H2O2) = h2o2 / conv1 * sm(sp_H2O2) % mol_wt
         end if
      end do
   end do

   return
end subroutine gocart_so2so4
