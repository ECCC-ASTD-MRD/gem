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

!/@*
module phy_snapshot_mod
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use series_mod, only: series_pause, series_resume
   use phy_status, only: phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
   use cpl_itf, only: cpl_snapshot
   implicit none
   private
   !@Public Functions
   public :: phy_snapshot
   !@Public Parameters
   character(len=*), parameter, public :: PHY_SNAPSHOT_STORE  = 'W'
   character(len=*), parameter, public :: PHY_SNAPSHOT_RESUME = 'R'
   !*@/

contains

   !/@*
   function phy_snapshot(F_mode) result(F_istat)
      implicit none
      !@Arguments
      character(len=*), intent(in) :: F_mode  !Snapshot mode: PHY_SNAPSHOT_STORE/RESUME
      !@Rerturns
      integer :: F_istat  !Return status (RMN_OK or RMN_ERR)
      !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2015
      !@revision
      !*@/

#include <rmnlib_basics.hf>
#include <mu_gmm.hf>
#include <msg.h>

      type(gmm_metadata) :: meta_busper
      integer :: istat
      real, pointer :: BUSPER_3d_digf(:,:), BUSPER_3d(:,:)
      ! ------------------------------------------------------------------
      F_istat = RMN_ERR
      if (phy_init_ctrl == PHY_NONE) then
         F_istat = PHY_NONE
         return
      else if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_snapshot) Physics not properly initialized.')
         return
      endif
      
      nullify(BUSPER_3d_digf, BUSPER_3d)

      select case (F_mode)

      case (PHY_SNAPSHOT_STORE)

         istat= gmm_get('BUSPER_3d', BUSPER_3d, meta_busper)
         if (.not.associated(BUSPER_3d)) then
            call msg(MSG_ERROR,'(phy_snapshot) Cannot find BUSPER')
            return
         endif
         istat= gmm_create('BUSPER_3d_digf', BUSPER_3d_digf, meta_busper, GMM_FLAG_RSTR)
         istat= gmm_get('BUSPER_3d_digf', BUSPER_3d_digf)
         if (associated(BUSPER_3d_digf)) then
            call msg(MSG_INFO,'(phy_snapshot) Saving BUSPER')
            BUSPER_3d_digf = BUSPER_3d
         else
            call msg(MSG_ERROR,'(phy_snapshot) Cannot save BUSPER')
         endif

         istat = series_pause()

      case (PHY_SNAPSHOT_RESUME)

         istat= gmm_get('BUSPER_3d_digf', BUSPER_3d_digf)
         istat= gmm_get('BUSPER_3d', BUSPER_3d)
         if (associated(BUSPER_3d_digf) .and. associated(BUSPER_3d)) then
            call msg(MSG_INFO,'(phy_snapshot) Restoring BUSPER')
            BUSPER_3d = BUSPER_3d_digf
         else
            call msg(MSG_ERROR,'(phy_snapshot) Cannot restore BUSPER')
         endif

         istat = series_resume()

      case default
         call msg(MSG_ERROR,'(phy_snapshot) Unknown mode.')
         return
      end select

      ! coupling may have something to do for snapshot

      call cpl_snapshot(F_mode)

      F_istat = RMN_OK
      ! ------------------------------------------------------------------

      return
   end function phy_snapshot

end module phy_snapshot_mod
