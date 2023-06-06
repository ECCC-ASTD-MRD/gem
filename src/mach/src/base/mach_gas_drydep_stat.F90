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
! Projet / Project : GEM-MACH
! Fichier / File   : mach_gas_drydep_stat.ftn90
! Creation         : Alexander Kallaur DEC 2007
! Description      : Provide statistics and values to all the parameters
!                    passed to mach_gas_drydep_solver.
!
!
! Arguments:  IN
!
!      metvar2D(:,MV2d_ILMO)  -> Monin-Obukhov length
!
!      metvar2D(:,MV2d_TSURF) -> Surface temperature (K)
!                NOTE: Two choices here:
!                      1) Take TSURF; averaged skin Temp from PHY (ISBA)
!                      2) Take TT(1) 1st level model temp
!                      In discussion with B. Bilodeau Nov 22 2007, he mentioned
!                      that there is a way to provide profile between TSURF & TT(1),
!                      if need be.
!
!      metvar2D(:,MV2d_UE) -> Surface friction velocity (m/s) taken from Phyvar(UE)
!
!      metvar2D(:,MV2d_QDIAG) -> Specific humidity
!
!      metvar2D(:,MV2d_FLUSOLIS) -> Downwards visible solar flux
!                Phyvar(FLUSOLIS)  (W/m2)
!
!      metvar2D(:,MV2d_RAINRATE) -> Precipitation rate (m/s)
!                Phyvar (U1 or RAINRATE from volatile bus)
!
!         lfu -> Land Form Use
!                CHEMVar(LAND_USE_15) from Chemical permanent bus
!                Derived from the CMC26 category data set into 15 category
!                set in subroutine "mach_landuse"
!
!          metvar2D(:,MV2d_PPLUS)  -> Surface presssure [Pa]
!                PhyVar(2p -> pplus from Dyn bus)
!
!          vd -> Deposition velocity for species "species"   (m/s)
!
! aero_resist -> Aerodynamic resistance for species "species" (s/m)
!
! diff_resist -> Molecular diffusion resistance for species "species" (s/m)
!
! surf_resist -> Total surface resistance for species "species" (s/m)
!
!============================================================================!
!
!!if_on
subroutine mach_gas_drydep_stat(vd, aero_resist, diff_resist, surf_resist, &
                                lfu, metvar2d)
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm, nb_gas_depo
   use chm_ptopo_grid_mod,   only: chm_ni
!!if_off
   use chm_metvar_mod,       only: MV2D_UE,    MV2D_ILMO,      MV2D_PPLUS,    &
                                   MV2D_TSURF, MV2D_FLUSOLIS, MV2D_RAINRATE, &
                                   MV2D_SNODP, MV2D_QDIAG
   use chm_utils_mod,        only: chm_lun_out, global_debug
   use chm_species_info_mod, only: species_master
   use mach_drydep_mod,      only: gas_depo
   implicit none
!!if_on
   real(kind=4),    intent(in) :: vd         (nb_gas_depo, chm_ni)
   real(kind=4),    intent(in) :: aero_resist(chm_ni, lucprm)
   real(kind=4),    intent(in) :: diff_resist(nb_gas_depo, chm_ni)
   real(kind=4),    intent(in) :: surf_resist(lucprm, nb_gas_depo, chm_ni)
   real(kind=4),    intent(in) :: lfu        (chm_ni, lucprm)
   real(kind=4),    intent(in) :: metvar2d   (chm_ni, SIZE_MV2D)
!!if_off
!
!  Local variables
!
   logical(kind=4)  :: local_dbg
   integer(kind=4)  :: i, sp, specie

   local_dbg = (.false. .or. global_debug)
   if( chm_lun_out < 0) return

   if (local_dbg) then
      write(chm_lun_out, *) 'begin data field analysis for mach_gas_drydep_solver (post mortem)'
   end if
   
   call calc_prn_stats_1d('metvar2d(:,            MV2D_ILMO) ' ,           metvar2d(:, MV2D_ILMO)          , chm_ni)
   call calc_prn_stats_1d('metvar2d(:,            MV2D_QDIAG)',            metvar2d(:, MV2D_QDIAG)         , chm_ni)
   call calc_prn_stats_1d('metvar2d(:,            MV2D_SNODP) ',           metvar2d(:, MV2D_SNODP)         , chm_ni)
   call calc_prn_stats_1d('metvar2d(:,            MV2D_PPLUS) ',           metvar2d(:, MV2D_PPLUS)         , chm_ni)
   call calc_prn_stats_1d('metvar2d(:,            MV2D_TSURF) ',           metvar2d(:, MV2D_TSURF)         , chm_ni)
   call calc_prn_stats_1d('metvar2d(:,            MV2D_RAINRATE) ',        metvar2d(:, MV2D_RAINRATE)      , chm_ni)
   call calc_prn_stats_1d('metvar2d(:,            MV2D_UE) ',              metvar2d(:, MV2D_UE)            , chm_ni)
   call calc_prn_stats_1d('metvar2d(:,            MV2D_FLUSOLIS) ',        metvar2d(:, MV2D_FLUSOLIS)      , chm_ni)
   do i = 1, lucprm
      call calc_prn_stats_1d('lfu        ', lfu(:,i),         chm_ni)
      call calc_prn_stats_1d('aero_resist', aero_resist(:,i), chm_ni)
   end do

   do sp = 1, nb_gas_depo
      specie = gas_depo(sp) % sp_id

      write (chm_lun_out,1) species_master(specie) % dyn_name, chm_ni
   1  format('For specie: ', a10, '  size: ', i5)

      call calc_prn_stats_1d('vd         ', vd(sp,:)         , chm_ni)
      call calc_prn_stats_1d('diff_resist', diff_resist(sp, :), chm_ni)

      do i = 1, lucprm
         call calc_prn_stats_1d('surf_resist', surf_resist(i, sp, :), chm_ni)
      end do
  
   end do

   if (local_dbg) then
      write(chm_lun_out, *) ' end data field analysis for mach_gas_drydep_solver (post mortem)'
      write(chm_lun_out, *) ' '
      write(chm_lun_out, *) ' '
   end if

   return

   contains

      subroutine calc_prn_stats_1d(arnam, f, fsiz)

         implicit none

         integer,              intent(in) :: fsiz
         real,dimension(fsiz), intent(in) :: f
         character(len=*),     intent(in) :: arnam

         write (chm_lun_out, 12) arnam, sum(f) / fsiz, minval(f), maxval(f), f(1), f(fsiz / 2), f(fsiz)
         12 format (a15, ' mean: ',e12.4, ' min: ', e12.4, ' max: ', e12.4, ' selected values (1,mid,end): ', 3e12.4)

      end subroutine calc_prn_stats_1d

end subroutine mach_gas_drydep_stat
