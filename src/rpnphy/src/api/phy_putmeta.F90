
module phy_putmeta_mod
   use phy_status, only: phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
   use phyfold, only: phyfoldmeta1
   use phymem, only: phymeta, phymem_find, phymem_updatemeta, PHY_NPATH_DEFAULT, PHY_BPATH_DEFAULT
   private
   public :: phy_putmeta

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

  integer, parameter :: PATHLENGTH = 8

contains

   !/@*
   function phy_putmeta(F_meta, F_name, F_npath, F_bpath, F_quiet) result(F_istat)
      implicit none
      !@object Update physics var metadata
      !@arguments
      type(phymeta),    intent(in)           :: F_meta     !Physics field metadata for first/only matching var
      character(len=*), intent(in)           :: F_name     !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath    !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath    !Bus path to search ['PVD']
      logical,          intent(in), optional :: F_quiet    !Quiet mode [.false.]
     !@return
      integer :: F_istat                                   !Return status (RMN_OK or RMN_ERR)
      !@author 
      !@revision
      !*@/
      integer :: istat, idxv1(1)
      character(len=PATHLENGTH) :: npath, bpath
      logical :: quiet_L
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_putmeta) Physics not properly initialized.')
         return
      endif

      ! Set default values
      npath = PHY_NPATH_DEFAULT
      if (present(F_npath)) npath = F_npath
      if (len_trim(npath) == 0) npath = PHY_NPATH_DEFAULT
      bpath = PHY_BPATH_DEFAULT
      if (present(F_bpath)) bpath = F_bpath
      if (len_trim(bpath) == 0) bpath = PHY_BPATH_DEFAULT

      quiet_L = .false.
      if (present(F_quiet)) quiet_L = F_quiet

      ! Retrieve matching record information
      istat = phymem_find(idxv1, F_name, npath, bpath, &
           F_quiet=quiet_L, F_shortmatch=.false.)
      if (istat <= 0) then
         if (.not. quiet_L) &
              call msg(MSG_WARNING,'(phy_putmeta) Cannot retrieve metadata for '//trim(F_name))
         return
      endif
      F_istat = phymem_updatemeta(F_meta, idxv1(1))
      if (istat <= 0) then
         if (.not. quiet_L) &
              call msg(MSG_WARNING,'(phy_putmeta) problem updating metadata for '//trim(F_name))
         return
      endif
      ! ---------------------------------------------------------------------
      return
   end function phy_putmeta


end module phy_putmeta_mod
