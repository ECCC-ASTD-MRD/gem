
module phy_getmeta_mod
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use clib_itf_mod, only: clib_toupper
   use phy_status, only: phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
   use phy_typedef, only: NPATH_DEFAULT, BPATH_DEFAULT, PHY_MAXNAMELENGTH
   use phymem, only: phymeta, nphyvars, phymem_find, phymem_getmeta_copy
   private
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   public :: phy_getmeta, phy_getmeta_single, phy_getmeta_list

   integer, parameter, public :: PATHLENGTH = 8

   interface phy_getmeta
      module procedure phy_getmeta_single
      module procedure phy_getmeta_list
   end interface phy_getmeta

contains

   !/@*
   function phy_getmeta_single(F_meta, F_name, F_npath, F_bpath, F_quiet) result(F_istat)
      implicit none
      !@objective Retreive physic var metadata for first/only matching var
      !@arguments
      type(phymeta),    intent(out)          :: F_meta   !Physics field metadata for first/only matching var
      character(len=*), intent(in)           :: F_name   !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath  !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath  !Bus path to search ['PVD']
      logical,          intent(in), optional :: F_quiet  !Do not emit warning for unmatched entry [.false.]
      !@return
      integer :: F_istat !# RMN_OK if matching var, RMN_ERR on error or no matching var
      !@author Ron McTaggart-Cowan - Spring 2014
      !*@/
      character(len=PHY_MAXNAMELENGTH) :: npath, bpath
      logical :: quiet
      type(phymeta), target  :: meta_tbl(1)
      type(phymeta), pointer :: meta_ptr(:)

      !# Set default values
      npath = NPATH_DEFAULT
      bpath = BPATH_DEFAULT
      quiet = .false.
      if (present(F_npath)) then
         if (len_trim(F_npath) /= 0) npath = F_npath
      endif
      if (present(F_bpath)) then
         if (len_trim(F_bpath) /= 0) bpath = F_bpath
      endif
      if (present(F_quiet)) quiet = F_quiet

      !# Retrieve metadata into temporary space
      meta_ptr => meta_tbl
      F_istat = phy_getmeta_list(meta_ptr, F_name, npath, bpath, &
           F_maxmeta=1, F_quiet=quiet)
      if (F_istat == 0) F_istat = RMN_ERR
      if (RMN_IS_OK(F_istat)) F_meta = meta_tbl(1)

      return
   end function phy_getmeta_single


   !/@*
   function phy_getmeta_list(F_meta,F_name,F_npath,F_bpath,F_maxmeta,F_quiet,F_shortmatch) result(F_istat)
      implicit none
      !@objective Retreive list of physic var metadata for all matching var
      !@arguments
      type(phymeta),    pointer              :: F_meta(:)    !Physics field metadata for all matching vars
      character(len=*), intent(in)           :: F_name       !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath      !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath      !Bus path to search ['PVD']
      integer,          intent(in), optional :: F_maxmeta    !Maximum number of vars to retrieve
      logical,          intent(in), optional :: F_quiet      !Do not emit warning for unmatched entry [.false.]
      logical,          intent(in), optional :: F_shortmatch  !if true, Match F_name against only the first len_trim(F_name) of input, variable or output name
      !@return
      integer :: F_istat !# number of matching vars (>=0), RMN_ERR on error
      !@author Ron McTaggart-Cowan - Spring 2014
      !*@/
      character(len=PHY_MAXNAMELENGTH) :: npath, bpath
      integer :: istat, maxmeta, i, nmatch
      integer, target :: idxv_target(nphyvars)
      logical :: quiet, shortmatch, to_alloc
      integer, pointer :: idxv(:)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_getmeta) Physics not properly initialized.')
         return
      endif

      !# Set default values
      npath = NPATH_DEFAULT
      bpath = BPATH_DEFAULT
      quiet = .false.
      maxmeta = size(idxv_target)
      shortmatch = .false.
      if (present(F_npath)) then
         if (len_trim(F_npath) /= 0) then
            npath = F_npath
            istat = clib_toupper(npath)
         endif
      endif
      if (present(F_bpath)) then
         if (len_trim(F_bpath) /= 0) then
            bpath = F_bpath
            istat = clib_toupper(bpath)
         endif
      endif
      if (present(F_quiet)) quiet = F_quiet
      if (present(F_maxmeta)) then
         if (F_maxmeta > 0) maxmeta = min(maxmeta,F_maxmeta)
      endif
      if (present(F_shortmatch)) shortmatch = F_shortmatch

      !# Retrieve metadata into temporary space
      idxv => idxv_target(1:maxmeta)
      nmatch = phymem_find(idxv, F_name, npath, bpath, quiet, shortmatch)

      !# Extract public part of the metadata
      if (nmatch > 0) then
         to_alloc = .true.
         if (associated(F_meta)) then
            if (size(F_meta) >= nmatch) to_alloc = .false.
         endif
         if (to_alloc) then
            nullify(F_meta)  !#NOTE: This may lead to memory leak but we cannot deallocate since it may not be allocated mem and changing behavious would be an API change
            allocate(F_meta(nmatch))
         endif
         do i = 1, nmatch
            istat = phymem_getmeta_copy(F_meta(i), idxv(i))
         enddo
      endif

      F_istat = nmatch
      return
   end function phy_getmeta_list

end module phy_getmeta_mod
