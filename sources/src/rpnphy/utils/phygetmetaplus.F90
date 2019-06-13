!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

module phygetmetaplus_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use clib_itf_mod, only: clib_toupper
   use phybus, only: entbus, perbus, dynbus, volbus
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj
   use phy_typedef
   implicit none
   private
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
#include <mu_gmm.hf>

   include "buses.cdk"

   public :: phygetmetaplus, phygetmetaplus_single, phygetmetaplus_list, phymetaplus

   integer, parameter, public :: PATHLENGTH = 8

   type phymetaplus
      type(phymeta) :: meta
      real, pointer :: vptr(:,:)
      integer :: index
   end type phymetaplus

   interface phygetmetaplus
      module procedure phygetmetaplus_single
      module procedure phygetmetaplus_list
   end interface phygetmetaplus

contains

   !/@*
   function phygetmetaplus_single(F_meta, F_name, F_npath, F_bpath, &
        F_quiet, F_shortmatch) result(F_istat)
      implicit none
      !@objective Retreive physic var metadata for first/only matching var
      !@arguments
      type(phymetaplus), intent(out) :: F_meta    !Physics field metadata
      character(len=*),  intent(in)  :: F_name    !Name of field to retrieve (input, variable or output name)
      character(len=*),  intent(in)  :: F_npath   !Name path to search ['VOI']
      character(len=*),  intent(in)  :: F_bpath   !Bus path to search ['PVD']
      logical,           intent(in)  :: F_quiet   !Do not emit warning for unmatched entry [.false.]
      logical,           intent(in)  :: F_shortmatch
      !@return
      integer :: F_istat !# RMN_OK if matching var, RMN_ERR on error or no matching var
      !@author Ron McTaggart-Cowan - Spring 2014
      !*@/
      type(phymetaplus), dimension(1), target  :: meta_tbl
      type(phymetaplus), dimension(:), pointer :: meta_ptr

      meta_ptr => meta_tbl
      F_istat = phygetmetaplus_list(meta_ptr, F_name, F_npath, F_bpath, &
           1, F_quiet, F_shortmatch)
      if (F_istat == 0) F_istat = RMN_ERR
      if (RMN_IS_OK(F_istat)) F_meta = meta_tbl(1)

      return
   end function phygetmetaplus_single


   !/@*
   function phygetmetaplus_list(F_meta, F_name, F_npath, F_bpath, &
        F_maxmeta, F_quiet, F_shortmatch) result(F_istat)
      implicit none
      !@objective Retreive list of physic var metadata for all matching var
      !@arguments
      type(phymetaplus), pointer    :: F_meta(:)  !Physics field metadata for all matching vars
      character(len=*),  intent(in) :: F_name     !Name of field to retrieve (input, variable or output name)
      character(len=*),  intent(in) :: F_npath    !Name path to search ['VOI']
      character(len=*),  intent(in) :: F_bpath    !Bus path to search ['PVD']
      integer,           intent(in) :: F_maxmeta  !Maximum number of vars to retrieve
      logical,           intent(in) :: F_quiet    !Do not emit warning for unmatched entry [.false.]
      logical,           intent(in) :: F_shortmatch  !if true, Match F_name against only the first len_trim(F_name) of input, variable or output name
      !@return
      integer :: F_istat !# number of matching vars (>=0), RMN_ERR on error
      !@author Ron McTaggart-Cowan - Spring 2014
      !*@/
      integer, external :: phygetvarlist2, phygetindx

      integer :: np, bp, v, istat, index, cnt, maxmeta, nmatch, ngetmax, nktot, niktot
      integer :: param(BUSPAR_MAXPAR)
      real    :: vmin, vmax
      character(len=PHY_MAXNAMELENGTH) :: name, npath, bpath, bus
      character(len=PHY_MAXNAMELENGTH) :: iname, vname, oname, iname_v, oname_v
      character(len=PHY_MAXNAMELENGTH) :: vlist(MAXBUS)
      logical :: full, to_alloc
      type(phymetaplus) :: meta_tmp(MAXBUS)
      real, pointer :: busptr(:,:)
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR

      name  = F_name  ; istat = clib_toupper(name)
      npath = F_npath ; istat = clib_toupper(npath)
      bpath = F_bpath ; istat = clib_toupper(bpath)
      maxmeta = min(size(meta_tmp),F_maxmeta)
      if (associated(F_meta)) &
           maxmeta = min(size(F_meta),maxmeta)

      !# Loop to search for matching record
      full = .false.
      cnt  = 0
      DO_BUS: do bp=1,len_trim(bpath)
         bus = bpath(bp:bp)
         DO_NAME: do np=1,len_trim(npath)
            iname = ' '
            vname = ' '
            oname = ' '
            select case (npath(np:np))
            case ('V')
               vname = F_name
            case ('I')
               iname = F_name
            case ('O')
               oname = F_name
            case DEFAULT
               call msg(MSG_WARNING,'(phy_getmeta) Ignoring unknown variable path entry '//npath(np:np))
               cycle
            end select
            ngetmax = max(1,min(size(vlist),maxmeta-cnt))
            nmatch = phygetvarlist2(vlist, vname, oname, iname, bus, -1, ngetmax, F_shortmatch)
            DO_FOUND: do v=1,nmatch
               oname_v = ' '
               iname_v = ' '
               istat = phygetindx(vlist(v), oname_v, iname_v, bus, index, param, size(param))
               if (istat < 0 .or. index < 1) cycle

               ! Record is found - fill metadata structure
               nullify(busptr)
               IF_NOTMAX: if (cnt < maxmeta) then
                  select case(bus)
                  case ('D')
                     busptr => dynbus
                  case ('P')
                     busptr => perbus
                  case ('V')
                     busptr => volbus
                  case ('E')
                     busptr => entbus
                  case DEFAULT
                     call msg(MSG_WARNING,'(phygetmetaplus) Unknown bus: '//trim(bus))
                     cycle
                  end select
                  cnt = cnt+1
                  vmin = transfer(param(BUSPAR_VMIN), vmin)
                  vmax = transfer(param(BUSPAR_VMAX), vmax)
                  nktot = param(BUSPAR_NK) * param(BUSPAR_FMUL) &
                       * (param(BUSPAR_MOSAIC)+1)
                  niktot = phy_lcl_ni * nktot
                  meta_tmp(cnt)%vptr  => busptr(index:index-1+niktot,:)
                  meta_tmp(cnt)%index =  index
                  meta_tmp(cnt)%meta%iname  = iname_v
                  meta_tmp(cnt)%meta%oname  = oname_v
                  meta_tmp(cnt)%meta%vname  = vlist(v)
                  meta_tmp(cnt)%meta%bus    = bus
                  meta_tmp(cnt)%meta%n(1:3) = (/phy_lcl_ni, phy_lcl_nj, nktot/)
                  meta_tmp(cnt)%meta%init   = param(BUSPAR_INIT)
                  meta_tmp(cnt)%meta%stag   = param(BUSPAR_STAG)
                  meta_tmp(cnt)%meta%esp    = param(BUSPAR_ESP)
                  meta_tmp(cnt)%meta%fmul   = param(BUSPAR_FMUL)
                  meta_tmp(cnt)%meta%nk     = param(BUSPAR_NK)
                  meta_tmp(cnt)%meta%mosaic = param(BUSPAR_MOSAIC)
                  meta_tmp(cnt)%meta%monot  = param(BUSPAR_MONOT)
                  meta_tmp(cnt)%meta%massc  = param(BUSPAR_MASSC)
                  meta_tmp(cnt)%meta%vmin   = vmin
                  meta_tmp(cnt)%meta%vmax   = vmax
                  meta_tmp(cnt)%meta%wload  = (param(BUSPAR_WLOAD) > 0)
                  meta_tmp(cnt)%meta%hzd    = (param(BUSPAR_HZD)   > 0)
               else
                  if (maxmeta > 1 .and. .not.F_quiet) &
                       call msg(MSG_WARNING,'(phy_getmeta) Metadata buffer overflow')
                  full = .true.
                  exit
               endif IF_NOTMAX
            enddo DO_FOUND
            if (full) exit
         enddo DO_NAME
         if (full) exit
      enddo DO_BUS

      !# Check that a match was found and allocate space for return
      if (cnt == 0) then
         if (.not.F_quiet) call msg(MSG_WARNING,'(phy_getmeta) No matching entry found in physics for '// &
              trim(F_name)//' on '//trim(bpath))
      else
         to_alloc = .true.
         if (associated(F_meta)) then
            if (size(F_meta) < cnt) then
               deallocate(F_meta, stat=istat)
            else
               to_alloc = .false.
            endif
         endif
         if (to_alloc) allocate(F_meta(cnt))
         F_meta(1:cnt) = meta_tmp(1:cnt)
      endif

      F_istat = cnt
      return
   end function phygetmetaplus_list

end module phygetmetaplus_mod
