!-------------------------------------- LICENCE BEGIN ------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------

module phy_put_mod
   use phy_status, only: phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
   use phyfold, only: phyfoldmeta1
   use phymem, only: phymeta, phyvar, phymem_find, PHY_NPATH_DEFAULT, PHY_BPATH_DEFAULT
   private
   public :: phy_put

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

  interface phy_put
     module procedure phy_put_2d
     module procedure phy_put_3d
     module procedure phy_put_4d
  end interface phy_put

  integer, parameter :: PATHLENGTH = 8

contains

   !/@*
   function phy_put_2d(F_fld, F_name, F_npath, F_bpath, F_start, &
        F_end, F_quiet) result(F_istat)
      implicit none
      !@object
      !	Transfer data to the physic space
      !@arguments
      real,             pointer              :: F_fld(:,:) !Field for the physics
      character(len=*), intent(in)           :: F_name     !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath    !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath    !Bus path to search ['PVD']
      integer,          intent(in), optional :: F_start(3) !Start index (lbound) for each dimension [(/1,1,1/)]
      integer,          intent(in), optional :: F_end(3)   !End index (ubound) for each dimension [(/ni,nj,1/)]
      logical,          intent(in), optional :: F_quiet    !Quiet mode [.false.]
     !@return
      integer :: F_istat                                   !Return status (RMN_OK or RMN_ERR)
      !@author Ron McTaggart-Cowan - Spring 2014
      !@revision
      !*@/
      integer :: istat
      integer, dimension(3) :: istart, iend
      character(len=PATHLENGTH) :: npath, bpath
      character(len=128) :: msg_S
      logical :: quiet_L
      type(phymeta) :: meta
      type(phyvar) :: myphyvar(1)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_put) Physics not properly initialized.')
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
      istat = phymem_find(myphyvar, F_name, npath, bpath, &
           F_quiet=quiet_L, F_shortmatch=.false.)
      if (istat <= 0) then
         if (.not. quiet_L) &
              call msg(MSG_WARNING,'(phy_put) Cannot retrieve metadata for '//trim(F_name))
         return
      endif
      meta = myphyvar(1)%meta

      ! Set automatic grid dimensions
      istart = 1
      if (present(F_start)) then
         where (F_start > 0)
            istart = F_start
         end where
      endif
      iend = (/meta%nlcl(1:2),1/)
      if (present(F_end)) then
         where (F_end > 0)
            iend = F_end
         end where
      endif

      ! Check bound and Set up space for return array
      if (iend(3) /= istart(3)) then
         call msg(MSG_WARNING,'(phy_put_2d) Can only put 1 level for '//trim(F_name))
         return
      endif

      if (associated(F_fld)) then
         if (.not.(size(F_fld,dim=1) == iend(1)-istart(1)+1 .and. &
              size(F_fld,dim=2) == iend(2)-istart(2)+1)) then
            write(msg_S,"(' :: target: ',3i4,5x,'phy: ',3i4)") shape(F_fld), meta%nlcl
            call msg(MSG_WARNING,'(phy_put) Invalid input pointer shape for '//trim(F_name)//trim(msg_S))
            return
         endif
      else
         call msg(MSG_WARNING,'(phy_put) Invalid input pointer for '//trim(F_name))
         return
      endif

      ! Fold information onto the bus
      !#TODO: check 2d to 3d through the phyfoldmeta itf
      F_istat = phyfoldmeta1(F_fld, istart, iend, meta)
     if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phy_put) Cannot fold '//trim(meta%vname))
      endif
      ! ---------------------------------------------------------------------
      return
   end function phy_put_2d


   !/@*
   function phy_put_3d(F_fld, F_name, F_npath, F_bpath, F_start, &
        F_end, F_quiet) result(F_istat)
      implicit none
      !@object
      !	Transfer data to the physic space
      !@arguments
      real,             pointer              :: F_fld(:,:,:) !Field for the physics
      character(len=*), intent(in)           :: F_name     !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath    !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath    !Bus path to search ['PVD']
      integer,          intent(in), optional :: F_start(3) !Start index (lbound) for each dimension [(/1,1,1/)]
      integer,          intent(in), optional :: F_end(3)   !End index (ubound) for each dimension [automatic]
      logical,          intent(in), optional :: F_quiet    !Quiet mode [.false.]
      !@return
      integer :: F_istat                                   !Return status (RMN_OK or RMN_ERR)
      !@author Ron McTaggart-Cowan - Spring 2014
      !@revision
      !*@/
      integer :: istat, lijk(3)
      integer, dimension(3) :: istart, iend
      character(len=PATHLENGTH) :: npath, bpath
      character(len=128) :: msg_S
      logical :: quiet_L
      type(phymeta) :: meta
      type(phyvar) :: myphyvar(1)
      real, pointer :: fld1(:,:,:)

      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_put) Physics not properly initialized.')
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
      istat = phymem_find(myphyvar, F_name, npath, bpath, &
           F_quiet=quiet_L, F_shortmatch=.false.)
      if (istat <= 0) then
         if (.not. quiet_L) &
              call msg(MSG_WARNING,'(phy_put) Cannot retrieve metadata for '//trim(F_name))
         return
      endif
      meta = myphyvar(1)%meta

      ! Set automatic grid dimensions
      istart = 1
      if (present(F_start)) then
         where (F_start > 0)
            istart = F_start
         end where
      endif

      iend = meta%nlcl
      if (present(F_end)) then
         where (F_end > 0)
            iend = F_end
         end where
      endif

      ! Check bound and Set up space for return array
      if (associated(F_fld)) then
         if (.not.(size(F_fld,dim=1) == iend(1)-istart(1)+1 .and. &
              size(F_fld,dim=2) == iend(2)-istart(2)+1 .and. &
              size(F_fld,dim=3) >= iend(3)-istart(3)+1)) then
            write(msg_S,"(' :: target: ',3i4,5x,'phy: ',3i4)") shape(F_fld), meta%nlcl
            call msg(MSG_WARNING,'(phy_put) Invalid input pointer shape for '//trim(F_name)//trim(msg_S))
            return
         endif
      else
         call msg(MSG_WARNING,'(phy_put) Invalid input pointer for '//trim(F_name))
         return
      endif

      ! Fold information onto the bus
      lijk = lbound(F_fld)
      fld1(istart(1):,istart(2):,istart(3):) => F_fld(lijk(1):,lijk(2):,lijk(3):)
      F_istat = phyfoldmeta1(fld1, istart, iend, meta)
     if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phy_put) Cannot fold '//trim(meta%vname))
      endif
      ! ---------------------------------------------------------------------
      return
   end function phy_put_3d


   !/@*
   function phy_put_4d(F_fld, F_name, F_npath, F_bpath, F_start, &
        F_end, F_quiet) result(F_istat)
      implicit none
      !@object
      !	Transfer data to the physic space
      !@arguments
      real,             pointer              :: F_fld(:,:,:,:) !Field for the physics
      character(len=*), intent(in)           :: F_name     !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath    !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath    !Bus path to search ['PVD']
      integer,          intent(in), optional :: F_start(4) !Start index (lbound) for each dimension [(/1,1,1/)]
      integer,          intent(in), optional :: F_end(4)   !End index (ubound) for each dimension [automatic]
      logical,          intent(in), optional :: F_quiet    !Quiet mode [.false.]
      !@return
      integer :: F_istat                                   !Return status (RMN_OK or RMN_ERR)
      !@author Ron McTaggart-Cowan - Spring 2014
      !@revision
      !*@/
      integer :: istat
      integer, dimension(4) :: istart, iend, lijk, iend0
      character(len=PATHLENGTH) :: npath, bpath
      character(len=128) :: msg_S
      logical :: quiet_L
      type(phymeta) :: meta
      type(phyvar) :: myphyvar(1)
      real, pointer :: fld1(:,:,:,:)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_put) Physics not properly initialized.')
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
      istat = phymem_find(myphyvar, F_name, npath, bpath, &
           F_quiet=quiet_L, F_shortmatch=.false.)
      if (istat <= 0) then
         if (.not. quiet_L) &
              call msg(MSG_WARNING,'(phy_put) Cannot retrieve metadata for '//trim(F_name))
         return
      endif
      meta = myphyvar(1)%meta

      ! Set automatic grid dimensions
      istart = 1
      if (present(F_start)) then
         where (F_start > 0)
            istart = F_start
         end where
      endif

      iend0(1:2) = meta%nlcl(1:2)
      iend0(3)   = meta%nk
      iend0(4)   = meta%fmul * (meta%mosaic + 1)
      iend(1:4)  = iend0(1:4)
      if (present(F_end)) then
         where (F_end > 0)
            iend = F_end
         end where
      endif

      ! Check bound and Set up space for return array
      if (associated(F_fld)) then
         if (.not.(size(F_fld,dim=1) == iend(1)-istart(1)+1 .and. &
              size(F_fld,dim=2) == iend(2)-istart(2)+1 .and. &
              size(F_fld,dim=3) >= iend(3)-istart(3)+1 .and. &
              size(F_fld,dim=4) >= iend(4)-istart(4)+1)) then
            write(msg_S,"(' :: target: ',3i4,5x,'phy: ',3i4)") shape(F_fld), iend0
            call msg(MSG_WARNING,'(phy_put) Invalid input pointer shape for '//trim(F_name)//trim(msg_S))
            return
         endif
      else
         call msg(MSG_WARNING,'(phy_put) Invalid input pointer for '//trim(F_name))
         return
      endif

      ! Fold information onto the bus
      lijk = lbound(F_fld)
      fld1(istart(1):,istart(2):,istart(3):,istart(4):) => &
           F_fld(lijk(1):,lijk(2):,lijk(3):,lijk(4):)
      F_istat = phyfoldmeta1(fld1, istart, iend, meta)
     if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phy_put) Cannot fold '//trim(meta%vname))
      endif
      ! ---------------------------------------------------------------------
      return
   end function phy_put_4d

end module phy_put_mod
