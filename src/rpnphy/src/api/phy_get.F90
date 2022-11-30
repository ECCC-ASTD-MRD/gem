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

module phy_get_mod
  use phy_status, only: phy_init_ctrl, PHY_CTRL_INI_OK, PHY_NONE
  use phyunfold, only: phyunfoldmeta1
  use phymem, only: phymeta, phyvar, phymem_find, PHY_NPATH_DEFAULT, PHY_BPATH_DEFAULT
  private
  public :: phy_get

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

  interface phy_get
     module procedure phy_get_2d
     module procedure phy_get_3d
     module procedure phy_get_4d
  end interface phy_get

  integer, parameter :: PATHLENGTH = 8

contains

   !/@*
   function phy_get_2d(F_fld, F_name, F_npath, F_bpath, F_start, &
        F_end, F_meta, F_realloc, F_quiet) result(F_istat)
      implicit none
      !@object Get a data copy from then physic space
      !@arguments
      real,              pointer             :: F_fld(:,:)   !Copy of physics field
      character(len=*), intent(in)           :: F_name       !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath      !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath      !Bus path to search ['PVD']
      integer,          intent(in), optional :: F_start(3)   !Start index (lbound) for each dimension [(/1,1,1/)]
      integer,          intent(in), optional :: F_end(3)     !End index (ubound) for each dimension [(/ni,nj,1/)]
      logical,          intent(in), optional :: F_realloc    !Allow reallocation [.false.]
      logical,          intent(in), optional :: F_quiet      !Quiet mode [.false.]
      type(phymeta),    intent(out),optional :: F_meta       !Physics field metadata for returned field
      !@return
      integer :: F_istat                                     !Return status (RMN_OK or RMN_ERR)
      !@author Ron McTaggart-Cowan and Stephane Chamberland - Winter 2015
      !*@/
      integer :: istat
      integer, dimension(3) :: istart, iend
      character(len=PATHLENGTH) :: npath, bpath
      character(len=128) :: msg_S
      logical :: to_alloc, allow_realloc, quiet_L
      type(phymeta) :: meta
      type(phyvar) :: myphyvar(1)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR

      if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_get) Physics not properly initialized.')
         return
      endif

      ! Set default values
      npath = PHY_NPATH_DEFAULT
      if (present(F_npath)) npath = F_npath
      if (len_trim(npath) == 0) npath = PHY_NPATH_DEFAULT

      bpath = PHY_BPATH_DEFAULT
      if (present(F_bpath)) bpath = F_bpath
      if (len_trim(bpath) == 0) bpath = PHY_BPATH_DEFAULT

      allow_realloc = .false. ; quiet_L = .false.
      if (present(F_realloc)) allow_realloc = F_realloc
      if (present(F_quiet)) quiet_L = F_quiet

      ! Retrieve matching record information
      istat = phymem_find(myphyvar, F_name, npath, bpath, &
           F_quiet=quiet_L, F_shortmatch=.false.)
      if (istat <= 0) then
         if (.not. quiet_L) &
              call msg(MSG_WARNING,'(phy_get) Cannot retrieve metadata for '//trim(F_name))
         return
      endif
      meta = myphyvar(1)%meta
      if (present(F_meta)) F_meta = meta

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
         call msg(MSG_WARNING,'(phy_get_2d) Can only get 1 level for '//trim(F_name))
         return
      endif

      to_alloc = .true.
      if (associated(F_fld)) then
         if (.not.(size(F_fld,dim=1) == iend(1)-istart(1)+1 .and. &
              size(F_fld,dim=2) == iend(2)-istart(2)+1)) then
            if (allow_realloc) then
               call msg(MSG_INFOPLUS,'(phy_get) reallocating output array')
               deallocate(F_fld)
            else
               write(msg_S,"(' :: target: ',3i4,5x,'phy: ',3i4)") shape(F_fld), meta%nlcl
               call msg(MSG_WARNING,'(phy_get) Invalid input pointer shape for '//trim(F_name)//trim(msg_S))
               return
            endif
         else
            to_alloc = .false.
         endif
      endif
      if (to_alloc) allocate(F_fld(istart(1):iend(1),istart(2):iend(2)))

      ! Unfold physics field and copy into output array
      !#TODO: check 2d to 3d through the phyunfoldmeta itf
      F_istat = phyunfoldmeta1(F_fld, istart, iend, meta)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phy_get) Cannot unfold '//trim(meta%vname))
      endif
      ! ---------------------------------------------------------------------
      return
   end function phy_get_2d


   !/@*
   function phy_get_3d(F_fld, F_name, F_npath, F_bpath, F_start, &
        F_end, F_meta, F_realloc, F_quiet) result(F_istat)
      implicit none
      !@object Get a data copy from then physic space
      !@arguments
      real,              pointer             :: F_fld(:,:,:) !Copy of physics field
      character(len=*), intent(in)           :: F_name       !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath      !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath      !Bus path to search ['PVD']
      integer,          intent(in), optional :: F_start(3)   !Start index (lbound) for each dimension [(/1,1,1/)]
      integer,          intent(in), optional :: F_end(3)     !End index (ubound) for each dimension [automatic]
      logical,          intent(in), optional :: F_realloc    !Allow reallocation [.false.]
      logical,          intent(in), optional :: F_quiet      !Quiet mode [.false.]
      type(phymeta),    intent(out),optional :: F_meta       !Physics field metadata for returned field
      !@return
      integer :: F_istat                                     !Return status (RMN_OK or RMN_ERR)
      !@author Ron McTaggart-Cowan - Spring 2014
      !*@/
      integer :: istat
      integer, dimension(3) :: istart, iend
      character(len=PATHLENGTH) :: npath, bpath
      character(len=128) :: msg_S
      logical :: to_alloc, allow_realloc, quiet_L
      type(phymeta) :: meta
      type(phyvar) :: myphyvar(1)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR

      if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_get) Physics not properly initialized.')
         return
      endif

      ! Set default values
      npath = PHY_NPATH_DEFAULT
      if (present(F_npath)) npath = F_npath
      if (len_trim(npath) == 0) npath = PHY_NPATH_DEFAULT

      bpath = PHY_BPATH_DEFAULT
      if (present(F_bpath)) bpath = F_bpath
      if (len_trim(bpath) == 0) bpath = PHY_BPATH_DEFAULT

      allow_realloc = .false. ; quiet_L = .false.
      if (present(F_realloc)) allow_realloc = F_realloc
      if (present(F_quiet)) quiet_L = F_quiet

      ! Retrieve matching record information
      istat = phymem_find(myphyvar, F_name, npath, bpath, &
           F_quiet=quiet_L, F_shortmatch=.false.)
      if (istat <= 0) then
         if (.not. quiet_L) &
              call msg(MSG_WARNING,'(phy_get) Cannot retrieve metadata for '//trim(F_name))
         return
      endif
      meta = myphyvar(1)%meta
      if (present(F_meta)) F_meta = meta

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
      to_alloc = .true.
      if (associated(F_fld)) then
         if (.not.(size(F_fld,dim=1) == iend(1)-istart(1)+1 .and. &
              size(F_fld,dim=2) == iend(2)-istart(2)+1 .and. &
              size(F_fld,dim=3) >= iend(3)-istart(3)+1)) then
            if (allow_realloc) then
               call msg(MSG_INFOPLUS,'(phy_get) reallocating output array')
               deallocate(F_fld)
            else
               write(msg_S,"(' :: target: ',3i4,5x,'phy: ',3i4)") shape(F_fld), meta%nlcl
               call msg(MSG_WARNING,'(phy_get) Invalid input pointer shape for '//trim(F_name)//trim(msg_S))
               return
            endif
         else
            to_alloc = .false.
         endif
      endif
      if (to_alloc) allocate(F_fld(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)))

      ! Unfold physics field into output array
      F_istat = phyunfoldmeta1(F_fld, istart, iend, meta)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phy_get) Cannot unfold '//trim(meta%vname))
      endif
      ! ---------------------------------------------------------------------
      return
   end function phy_get_3d


   !/@*
   function phy_get_4d(F_fld, F_name, F_npath, F_bpath, F_start, &
        F_end, F_meta, F_realloc, F_quiet) result(F_istat)
      implicit none
      !@object Get a data copy from then physic space
      !@arguments
      real,              pointer             :: F_fld(:,:,:,:) !Copy of physics field
      character(len=*), intent(in)           :: F_name       !Name of field to retrieve (input, variable or output name)
      character(len=*), intent(in), optional :: F_npath      !Name path to search ['VOI']
      character(len=*), intent(in), optional :: F_bpath      !Bus path to search ['PVD']
      integer,          intent(in), optional :: F_start(4)   !Start index (lbound) for each dimension [(/1,1,1,1/)]
      integer,          intent(in), optional :: F_end(4)     !End index (ubound) for each dimension [automatic]
      logical,          intent(in), optional :: F_realloc    !Allow reallocation [.false.]
      logical,          intent(in), optional :: F_quiet      !Quiet mode [.false.]
      type(phymeta),    intent(out),optional :: F_meta       !Physics field metadata for returned field
      !@return
      integer :: F_istat                                     !Return status (RMN_OK or RMN_ERR)
      !@author Ron McTaggart-Cowan - Spring 2014
      !*@/
      integer :: istat
      integer, dimension(4) :: istart, iend, iend0
      character(len=PATHLENGTH) :: npath, bpath
      character(len=128) :: msg_S
      logical :: to_alloc, allow_realloc, quiet_L
      type(phymeta) :: meta
      type(phyvar) :: myphyvar(1)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR

      if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
         call msg(MSG_ERROR,'(phy_get) Physics not properly initialized.')
         return
      endif

      ! Set default values
      npath = PHY_NPATH_DEFAULT
      if (present(F_npath)) npath = F_npath
      if (len_trim(npath) == 0) npath = PHY_NPATH_DEFAULT

      bpath = PHY_BPATH_DEFAULT
      if (present(F_bpath)) bpath = F_bpath
      if (len_trim(bpath) == 0) bpath = PHY_BPATH_DEFAULT

      allow_realloc = .false. ; quiet_L = .false.
      if (present(F_realloc)) allow_realloc = F_realloc
      if (present(F_quiet)) quiet_L = F_quiet

      ! Retrieve matching record information
      istat = phymem_find(myphyvar, F_name, npath, bpath, &
           F_quiet=quiet_L, F_shortmatch=.false.)
      if (istat <= 0) then
         if (.not. quiet_L) &
              call msg(MSG_WARNING,'(phy_get) Cannot retrieve metadata for '//trim(F_name))
         return
      endif
      meta = myphyvar(1)%meta
      if (present(F_meta)) F_meta = meta

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
      to_alloc = .true.
      if (associated(F_fld)) then
         if (.not.(size(F_fld,dim=1) == iend(1)-istart(1)+1 .and. &
              size(F_fld,dim=2) == iend(2)-istart(2)+1 .and. &
              size(F_fld,dim=3) >= iend(3)-istart(3)+1 .and. &
              size(F_fld,dim=4) >= iend(4)-istart(4)+1)) then
            if (allow_realloc) then
               call msg(MSG_INFOPLUS,'(phy_get) reallocating output array')
               deallocate(F_fld)
            else
               write(msg_S,"(' :: target: ',4i4,5x,'phy: ',4i4)") shape(F_fld),iend0
               call msg(MSG_WARNING,'(phy_get) Invalid input pointer shape for '//trim(F_name)//trim(msg_S))
               return
            endif
         else
            to_alloc = .false.
         endif
      endif
      if (to_alloc) allocate(F_fld(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3),istart(4):iend(4)))

      ! Unfold physics field into output array
      F_istat = phyunfoldmeta1(F_fld, istart, iend, meta)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(phy_get) Cannot unfold '//trim(meta%vname))
      endif
      ! ---------------------------------------------------------------------
      return
   end function phy_get_4d

end module phy_get_mod
