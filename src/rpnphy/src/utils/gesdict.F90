!-------------------------------------- LICENCE BEGIN --------------------------
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
!-------------------------------------- LICENCE END ----------------------------

module gesdictmod
   use clib_itf_mod, only: clib_toupper
   use str_mod, only: str_normalize
   use phymem, only: pbuslist, PHY_BUSID, PHY_NBUSES, PHY_NAMELEN, PHY_MAXFLAGS, PHY_MAXVARS
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj
   implicit none
   private
   public :: gesdictadd, gesdictcheck

#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   
contains

   !/@*
   function gesdictadd(varname, outname, inname, sername, vardesc, &
        shape, dynini, stagg, fmul, nmosaic, vmin, vmax, &
        wload, hzd, monot, massc, flags, &
        busidx, ni ,nk, memgap) result(F_istat)
      implicit none
      integer :: F_istat
      character(len=*), intent(in) :: varname, outname, inname, sername, vardesc, shape, flags(:)
      integer,          intent(in) :: dynini, stagg, fmul, nmosaic, wload, hzd, monot, massc, busidx, ni ,nk, memgap
      real :: vmin, vmax
      !*@/

      character(len=64) :: prefix_S, basename_S, time_S, ext_S
      integer :: idxb, esp, fmosaik, nk1, wload1, istat
      !-------------------------------------------------------------------
      F_istat = RMN_ERR

      if (nmosaic > 0) then
         call msg(MSG_ERROR, '(gesdict) mosaic not supported; VN='//trim(varname))
         return
      endif

      fmosaik = nmosaic + 1

      istat = gesdictcheck(varname, outname, inname, sername, vardesc)
      !#Note: possibility to double specify the same var was removed.
      if (istat /= RMN_OK) return

      call gmmx_name_parts(varname, prefix_S, basename_S, time_S, ext_S)
      wload1 = wload
      if (time_S /= ':P') wload1 = 0  !# For backward compat
      nk1 = nk
      if (shape == "A") nk1 = 1

      if (pbuslist(busidx)%nvars >= PHY_MAXVARS) then
         call msg(MSG_ERROR, &
              '(gesdict) Too many vars - increase PHY_MAXVAR in phymem')
         return
      endif
      
      pbuslist(busidx)%nvars = pbuslist(busidx)%nvars + 1
      idxb = pbuslist(busidx)%nvars

      pbuslist(busidx)%meta(idxb)%ni = ni
      pbuslist(busidx)%meta(idxb)%nk = nk1
      pbuslist(busidx)%meta(idxb)%fmul = fmul
      pbuslist(busidx)%meta(idxb)%mosaic = nmosaic
      pbuslist(busidx)%meta(idxb)%size = ni*nk1*fmul*fmosaik

      pbuslist(busidx)%meta(idxb)%ibus = busidx
      
      pbuslist(busidx)%meta(idxb)%idxb = idxb
      pbuslist(busidx)%meta(idxb)%idxv = -1
      pbuslist(busidx)%meta(idxb)%i0 = pbuslist(busidx)%nsize + 1 + memgap
      pbuslist(busidx)%meta(idxb)%in = pbuslist(busidx)%meta(idxb)%i0 + pbuslist(busidx)%meta(idxb)%size - 1
      pbuslist(busidx)%meta(idxb)%nlcl = (/phy_lcl_ni, phy_lcl_nj, nk1*fmul*fmosaik/)
      
      pbuslist(busidx)%meta(idxb)%init = dynini
      pbuslist(busidx)%meta(idxb)%stag = stagg
      pbuslist(busidx)%meta(idxb)%wload = (wload1 > 0)
      pbuslist(busidx)%meta(idxb)%hzd = (hzd > 0)
      pbuslist(busidx)%meta(idxb)%monot = monot
      pbuslist(busidx)%meta(idxb)%massc = massc
      pbuslist(busidx)%meta(idxb)%vmin = vmin
      pbuslist(busidx)%meta(idxb)%vmax = vmax

      nullify(pbuslist(busidx)%meta(idxb)%bptr)

      pbuslist(busidx)%meta(idxb)%bus = PHY_BUSID(busidx)
      pbuslist(busidx)%meta(idxb)%iname = inname
      pbuslist(busidx)%meta(idxb)%oname = outname
      pbuslist(busidx)%meta(idxb)%vname = varname
      pbuslist(busidx)%meta(idxb)%sname = sername
      
      call str_normalize(pbuslist(busidx)%meta(idxb)%iname)
      call str_normalize(pbuslist(busidx)%meta(idxb)%oname)
      call str_normalize(pbuslist(busidx)%meta(idxb)%vname)
      call str_normalize(pbuslist(busidx)%meta(idxb)%sname)

      pbuslist(busidx)%meta(idxb)%flags(:) = flags(:)
      pbuslist(busidx)%meta(idxb)%desc = vardesc
      
      pbuslist(busidx)%nsize = pbuslist(busidx)%meta(idxb)%i0 + pbuslist(busidx)%meta(idxb)%size - 1
      
      F_istat = pbuslist(busidx)%meta(idxb)%i0
      !-------------------------------------------------------------------
      return
   end function gesdictadd

   
   !/@*
   function gesdictcheck(vname, iname, oname, sname, desc) result(F_istat)
      implicit none
      character(len=*), intent(in) :: vname, iname, oname, sname, desc
      integer :: F_istat
      !*@/
      integer :: ibus, idxb
      !-------------------------------------------------------------------
      F_istat = RMN_ERR
      
      DO_BUS: do ibus = 1, PHY_NBUSES
         DO_VAR: do idxb = 1, pbuslist(ibus)%nvars
            if (vname == pbuslist(ibus)%meta(idxb)%vname) then
               call msg(MSG_ERROR, '(gesdict) Duplicate entry (VN) for: '//trim(vname))
               return
            endif
            if (oname == pbuslist(ibus)%meta(idxb)%oname) then
               call msg(MSG_ERROR, '(gesdict) Duplicate entry (ON) for: '//trim(oname))
               return
            endif
            if (iname == pbuslist(ibus)%meta(idxb)%iname) then
               call msg(MSG_ERROR, '(gesdict) Duplicate entry (IN) for: '//trim(iname))
               return
            endif
            if (sname == pbuslist(ibus)%meta(idxb)%sname) then
               call msg(MSG_ERROR, '(gesdict) Duplicate entry (SN) for: '//trim(sname))
               return
            endif
            if (desc == pbuslist(ibus)%meta(idxb)%desc) then
               call msg(MSG_WARNING, '(gesdict) Duplicate entry (VD) for: '//trim(desc))
            endif
         enddo DO_VAR
      enddo DO_BUS
      
      F_istat = RMN_OK
      !-------------------------------------------------------------------
      return
   end function gesdictcheck

end module gesdictmod


!/@*
subroutine gesdict(ni, nk, lindex, lachaine)
   use phy_options, only: debug_mem_L
   use phymem, only: phymem_busidx, phymem_init, phymem_isalloc, PHY_NAMELEN, PHY_MAXFLAGS
   use gesdictmod, only: gesdictadd, gesdictcheck
   use splitst, only: splitst4
   implicit none
!!!#include <arch_specific.hf>
   !@Arguments
   ! ni         horizontal dimension
   ! nk         vertical dimension
   ! lindex     starting index on the bus
   ! lachaine   string identifying the variable attributes
   character(len=*), intent(in) :: lachaine
   integer, intent(in)  :: ni, nk
   integer, intent(out) :: lindex

   !@Object
   !    Manages the dictionary describing the 4 main buses of the unified
   !    CMC-RPN physics package interface (BUSENT, BUSDYN, BUSPER and BUSVOL).
   !    Each variable has a formal name <bus>nm(*) and a formal
   !    description <bus>dc(*) along with 4 attributes <bus>par(*,4).
   !    The first and second attributes are respectively the starting
   !    index on <bus> and the length of the variable. The third
   !    attribute is the multiplicity factor. The fourth attribute is
   !    the a flag to identify variables that are defined on staggered levels.
   !    The recognized token in "lachaine" are:
   !         VN=  ;       ===> formal name
   !         ON=  ;       ===> output name (4 letters only)
   !         IN=  ;       ===> input  name (4 letters only)
   !         SN=  ;       ===> series name (4 letters only)
   !         VD=  ;       ===> formal description
   !         VS=  ;       ===> variable shape (accepted shapes are M, T, E and
   !                           A with +, - or * followed by an integer)
   !         VB=  ;       ===> bus identification (D, P and V)
   !         MIN= ;       ===> minimum value of the field
   !         MAX= ;       ===> maximum value of the field
   !       WLOAD= ;       ===> water load flag (default=0)
   !         HZD= ;       ===> Horizontal diffusion (default=0)
   !       MASSC= ;       ===> mass conserv (default=0)
   !       MONOT= ;       ===> monotone interpolation (default=1)
   !       FLAGS= ;       ===> list of keywords '+' separated
   !*@/
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   integer, parameter :: MAXVALPERCHAR = 34

   logical, save :: init_L = .false.

   integer :: memgap, busidx
   character(len=1) ::   bus
   character(len=4) ::   outname,inname,sername
   character(len=3) ::   shape
   character(len=7) ::   struc
   character(len=PHY_NAMELEN) ::  varname, flags(PHY_MAXFLAGS)
   character(len=256) :: vdescrp
   character(len=256) :: vardesc
   character(len=256) :: string
   integer :: nmosaic, fmul, dynini, stagg
   integer :: i, wload, hzd, monot, massc, istat
   real :: vmin, vmax
   !-------------------------------------------------------------------
   if (.not.init_L) then
      init_L = .true.
      istat = phymem_init()
   endif

   if (phymem_isalloc()) then
      call physeterror('gesdict', 'Call to gesdict after buses have been allocated')
      return
   endif

   call low2up(lachaine,string)
   istat = splitst4(varname, outname, inname, sername, &
        vdescrp, struc, shape, nmosaic, fmul, bus, dynini,&
        stagg, vmin, vmax, wload, hzd, monot, massc, flags, string)
   if (.not.RMN_IS_OK(istat)) then
      call physeterror('gesdict', 'Invalid gesdict string: '//trim(string))
      return
   endif

   !# Mosaic field have an extra row for the averaged layer (1st row)
   vardesc = trim(vdescrp)//';VS='//struc

   i = MAXVALPERCHAR ** (len(outname) - max(len_trim(outname), len_trim(inname)))
   if (any(shape == (/'M', 'T', 'E'/)) .and. fmul > i) then
      call msg(MSG_WARNING,'Varname '//trim(varname)//' is declared as a 3D var with multiple categories, fmul>1, I/O names must be 2/3 char max.')
   endif

   lindex  = RMN_ERR
   busidx = phymem_busidx(bus)
   if (busidx <= 0) then
      call physeterror('gesdict', 'Unknown bus: '//trim(string))
      return
   endif

   memgap = 0
   if (debug_mem_L) memgap = max(ni/5,1)+2*(busidx-1)
   lindex = gesdictadd( &
        varname, outname, inname, sername, &
        vardesc, shape, dynini, stagg, fmul, nmosaic, vmin, vmax, &
        wload, hzd, monot, massc, flags, &
        busidx, ni, nk, memgap)
   if (.not.RMN_IS_OK(lindex)) then
      call physeterror('gesdict', 'Problem getting lindex')
      return
   endif
   !-------------------------------------------------------------------
   return
end subroutine gesdict


