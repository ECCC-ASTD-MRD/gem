!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/r itf_phy_step - Apply the physical processes: CMC/RPN package

      subroutine itf_phy_step ( F_step_kount, F_lctl_step )
      use iso_c_binding
      use phy_itf, only: phy_input,phy_step,phy_snapshot
      use itf_phy_filter, only: ipf_smooth_fld
      use init_options
      use lun
      use tr3d
      use rstr
      use path
      use wb_itf_mod
      use ptopo
      use gem_timing
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_step_kount, F_lctl_step

!arguments
!  Name                    Description
!----------------------------------------------------------------
! F_step_kount             step count
! F_lctl_step              step number
!----------------------------------------------------------------

      include "rpn_comm.inc"

      integer,external :: itf_phy_prefold_opr

      integer err_geom, err_input, err_step, err_smooth, err
!
!     ---------------------------------------------------------------
!
      if ( Rstri_user_busper_L .and. (F_step_kount == 0) ) return

      call gemtime_start ( 40, 'PHYSTEP', 1 )

      if (Lun_out > 0) write (Lun_out,1001) F_lctl_step

      if (F_step_kount == 0) then
         call itf_phy_geom (err_geom)
         if (NTR_Tr3d_ntr > 0) then
            err = wb_put('itf_phy/READ_TRACERS', &
                                NTR_Tr3d_name_S(1:NTR_Tr3d_ntr))
         end if
      end if

      !call pw_glbstat('PW_BEF')

      call itf_phy_copy ()
      !if (F_step_kount == 0) call itf_phy_glbstat('befinp')

      call gemtime_start ( 45, 'PHY_input', 40 )
      err_input = phy_input ( itf_phy_prefold_opr, F_step_kount, &
            Path_phyincfg_S, Path_phy_S, 'GEOPHY/Gem_geophy.fst' )

      !if (F_step_kount == 0) call itf_phy_glbstat('Aftinp')
      call gem_error (err_input,'itf_phy_step','Problem with phy_input')
      call gemtime_stop  ( 45 )

      ! Smooth the thermodynamic state variables on request
      err_smooth = min(&
           ipf_smooth_fld('PW_TT:M','TTMS'), &
           ipf_smooth_fld('TR/HU:M','HUMS'), &
           ipf_smooth_fld('PW_TT:P','TTPS'), &
           ipf_smooth_fld('TR/HU:P','HUPS') &
           )
      call gem_error (err_smooth,'itf_phy_step','Problem with ipf_smooth_fld')

      call set_num_threads ( Ptopo_nthreads_phy, F_step_kount )
      !call itf_phy_glbstat('befphy')

      call gemtime_start ( 46, 'PHY_step', 40 )
      err_step = phy_step ( F_step_kount, F_lctl_step )
      call rpn_comm_barrier (RPN_COMM_ALLGRIDS, err)
      call gemtime_stop  ( 46 )

      call gem_error (err_step,'itf_phy_step','Problem with phy_step')
      !call itf_phy_glbstat('aftphy')

      call set_num_threads ( Ptopo_nthreads_dyn, F_step_kount )

      call gemtime_start ( 47, 'PHY_update', 40 )
      call itf_phy_update3 ( F_step_kount > 0 )
      call gemtime_stop  ( 47 )
      !call pw_glbstat('PW_AFT')

      call gemtime_start ( 48, 'PHY_output', 40 )
      call itf_phy_output2 ( F_lctl_step )
      call gemtime_stop  ( 48 )

      if ( Init_mode_L ) then
         if (F_step_kount ==   Init_halfspan) then
            err = phy_snapshot('W')
         else if (F_step_kount == 2*Init_halfspan) then
            err = phy_snapshot('R')
         end if
      end if

      call gemtime_stop ( 40 )

 1001 format(/,'PHYSICS : PERFORMING TIMESTEP #',I9, &
             /,'========================================')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine itf_phy_step
