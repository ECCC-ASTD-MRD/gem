
!/@*
module phy_snapshot_mod
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use series_mod, only: series_pause, series_resume
   use phy_status, only: phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE

#ifdef HAVE_NEMO
   use cpl_itf, only: cpl_snapshot
#endif

   use phymem, only: phymem_gmmname, PHY_PBUSIDX
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
      use rmn_gmm
      implicit none
      !@Arguments
      character(len=*), intent(in) :: F_mode  !Snapshot mode: PHY_SNAPSHOT_STORE/RESUME
      !@Rerturns
      integer :: F_istat  !Return status (RMN_OK or RMN_ERR)
      !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2015
      !@revision
      !*@/

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

      type(gmm_metadata) :: gmmmeta
      integer :: istat
      real, pointer :: busprt_digf(:,:), busptr(:,:)
      character(len=32) :: gmmname, gmmname_digf
      ! ------------------------------------------------------------------
      F_istat = RMN_ERR
      if (phy_init_ctrl == PHY_NONE) then
         F_istat = PHY_NONE
         return
      else if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_snapshot) Physics not properly initialized.')
         return
      endif
      
      gmmname = phymem_gmmname(PHY_PBUSIDX)
      gmmname_digf = trim(gmmname)//'_digf'
      nullify(busprt_digf, busptr)

      select case (F_mode)

      case (PHY_SNAPSHOT_STORE)

         istat= gmm_get(gmmname, busptr, gmmmeta)
         if (.not.associated(busptr)) then
            call msg(MSG_ERROR,'(phy_snapshot) Cannot find BUSPER')
            return
         endif
         istat= gmm_create(gmmname_digf, busprt_digf, gmmmeta, GMM_FLAG_RSTR)
         istat= gmm_get(gmmname_digf, busprt_digf)
         if (associated(busprt_digf)) then
            call msg(MSG_INFO,'(phy_snapshot) Saving BUSPER')
            busprt_digf = busptr
         else
            call msg(MSG_ERROR,'(phy_snapshot) Cannot save BUSPER')
         endif

         istat = series_pause()

      case (PHY_SNAPSHOT_RESUME)

         istat= gmm_get(gmmname_digf, busprt_digf)
         istat= gmm_get(gmmname, busptr)
         if (associated(busprt_digf) .and. associated(busptr)) then
            call msg(MSG_INFO,'(phy_snapshot) Restoring BUSPER')
            busptr = busprt_digf
         else
            call msg(MSG_ERROR,'(phy_snapshot) Cannot restore BUSPER')
         endif

         istat = series_resume()

      case default
         call msg(MSG_ERROR,'(phy_snapshot) Unknown mode.')
         return
      end select

#ifdef HAVE_NEMO
      ! coupling may have something to do for snapshot

      call cpl_snapshot(F_mode)
#endif

      F_istat = RMN_OK
      ! ------------------------------------------------------------------

      return
   end function phy_snapshot

end module phy_snapshot_mod
