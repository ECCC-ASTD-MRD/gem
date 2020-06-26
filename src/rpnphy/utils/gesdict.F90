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
   implicit none
   private
   public :: gesdictlock, gesdictislock, gesdictadd, gesdictcheck

#include <rmnlib_basics.hf>
#include <msg.h>
   include "buses.cdk"

   logical, save :: buslck = .false.
   
contains

   !/@*
   subroutine gesdictlock()
      implicit none
      !*@/
      buslck = .true.
      return
   end subroutine gesdictlock

   !/@*
   function gesdictislock() result(F_istat)
      implicit none
      logical :: F_istat
      !*@/
      F_istat = buslck
      return
   end function gesdictislock

   !/@*
   function gesdictadd(varname, outname, inname, sername, vardesc, &
        shape, dynini, stagg, fmul, nmosaic, ivmin, ivmax, &
        wload, hzd, monot, massc, &
        bustop, busnm, busdc, buspar, busspc, ni ,nk, memgap) result(F_istat)
      implicit none
      integer :: F_istat
      character(len=*), intent(in) :: varname, outname, inname, sername, vardesc, shape
      integer,          intent(in) :: dynini, stagg, fmul, nmosaic, ivmin, ivmax, wload, hzd, monot, massc, ni ,nk, memgap
      character(len=*), intent(inout) :: busnm(:,:), busdc(:)
      integer,          intent(inout) :: bustop, buspar(:,:), busspc
      !*@/

      character(len=64) :: prefix_S, basename_S, time_S, ext_S
      integer :: i, lindex, esp, fmosaik
      ! real :: vmin, vmin0
      !-------------------------------------------------------------------
      F_istat = RMN_ERR

      !TODO: should we discontinue support for mosaic since it's only half supported?
!!$      if (nmosaic > 0) then
!!$         call msg(MSG_ERROR, '(gesdict) mosaic not supported; VN='//trim(varname))
!!$         return
!!$      endif

      fmosaik = nmosaic + 1
      do i=1,bustop

         if (inname == busnm(i, BUSNM_IN)) then
            call msg(MSG_ERROR, '(gesdict) input name conflict, IN='//trim(inname)//' (VN='//trim(varname)//') is already used for VN='//trim(busnm(i, BUSNM_VN)))
            return
         endif

         !# verifier si la meme description existe deja
         if (vardesc == busdc(i)) then
            if (varname /= busnm(i, BUSNM_VN)) then
               call msg(MSG_ERROR, '(gesdict) name conflict, same description for VN='//trim(varname)//' and VN='//trim(busnm(i, BUSNM_VN)))
               return
            endif
         endif

         if (varname == busnm(i, BUSNM_VN)) then
            if (vardesc /= busdc(i)) then
               call msg(MSG_ERROR, '(gesdict) name conflict, same description for VN='//trim(varname)//' and VN='//trim(busnm(i, BUSNM_VN)))
               return
            endif
            esp = ni*nk
            if (shape == "A") esp = ni
            if (buspar(i, BUSNM_ON) /= (esp*fmul*fmosaik)) then
               call msg(MSG_ERROR, '(gesdict) shape confilct, VN='//trim(varname)//' and VN='//trim(busnm(i, BUSNM_VN)))
               return
            endif
            lindex = buspar(i, BUSPAR_I0)
            return
         endif

      enddo

      bustop = bustop + 1
      esp = ni*nk
      buspar(bustop, BUSPAR_NK) = nk
      if (shape == "A") then
         esp = ni
         buspar(bustop, BUSPAR_NK) = 1
      endif
      buspar(bustop, BUSPAR_ESP) = esp
      esp = esp*fmul*fmosaik
      busnm(bustop, BUSNM_VN) = varname
      busnm(bustop, BUSNM_ON) = outname
      busnm(bustop, BUSNM_IN) = inname
      busnm(bustop, BUSNM_SN) = sername
      busdc(bustop) = vardesc
      buspar(bustop, BUSPAR_I0)   = busspc + 1 + memgap
      buspar(bustop, BUSPAR_NIK)  = esp
      buspar(bustop, BUSPAR_INIT) = dynini
      buspar(bustop, BUSPAR_STAG) = stagg
      buspar(bustop, BUSPAR_FMUL) = fmul
      buspar(bustop, BUSPAR_MOSAIC) = nmosaic
      buspar(bustop, BUSPAR_WLOAD) = wload
      buspar(bustop, BUSPAR_HZD)   = hzd
      buspar(bustop, BUSPAR_MONOT) = monot
      buspar(bustop, BUSPAR_MASSC) = massc
      buspar(bustop, BUSPAR_VMIN)  = ivmin
      buspar(bustop, BUSPAR_VMAX)  = ivmax
      busspc = buspar(bustop, BUSPAR_I0) + esp - 1
      lindex = buspar(bustop, BUSPAR_I0)
      call gmmx_name_parts(varname, prefix_S, basename_S, time_S, ext_S)
      if (time_S /= ':P') buspar(bustop, BUSPAR_WLOAD) = 0  !# For backward compat, TODO:should be removed and handled by the dyn
      F_istat = lindex
      !-------------------------------------------------------------------
      return
   end function gesdictadd


   !/@*
   function gesdictcheck(inbus, varname, outname, &
        bustop, busnm, busname) result(F_istat)
      implicit none
      integer :: F_istat
      character(len=*), intent(in) :: inbus, varname, outname, busname
      character(len=*), intent(inout) :: busnm(:,:)
      integer,          intent(inout) :: bustop
      !*@/
      integer :: i
      !-------------------------------------------------------------------
      F_istat = RMN_ERR

      do i=1,bustop

         !# verifier que le nom de la variable est unique
         if (inbus /= busname) then
            if (varname == busnm(i, BUSNM_VN)) then
               call msg(MSG_ERROR, '(gesdict) varname conflict: VN='//trim(varname)//' VB='//trim(inbus)//'; already declared in VB='//trim(busname))
               return
            endif
         endif

         !# verifier que le nom de 4 lettres est unique
         if (outname == busnm(i, BUSNM_ON) .and. varname /= busnm(i, BUSNM_VN)) then
            call msg(MSG_ERROR, '(gesdict) output name conflict: ON='//trim(outname)//' VN= '//trim(varname)//'; already used by VN='//trim(busnm(i, BUSNM_VN)))
            return
         endif

         !# verifier qu'une variable ne porte qu'un seul nom de 4 lettres
         if (varname == busnm(i, BUSNM_VN) .and. outname /= busnm(i, BUSNM_ON)) then
            call msg(MSG_ERROR, '(gesdict) name conflict: VN='//trim(varname)//' has 2 output names ON= '//trim(outname)//'; ON='//trim(busnm(i, BUSNM_ON)))
            return
         endif

         !# verifier qu'un "varname" ne soit pas identique a un "outname"
         if (varname /= busnm(i, BUSNM_VN) .and. varname == busnm(i, BUSNM_ON)) then
            call msg(MSG_ERROR, '(gesdict) name conflict: VN='//trim(varname)//' already used as ON= '//trim(busnm(i, BUSNM_ON))//' for VN='//trim(busnm(i, BUSNM_VN)))
            return
         endif
         !# verifier qu'un "outname" ne soit pas identique a un "varname"
         if (outname == busnm(i, BUSNM_VN) .and. varname /= busnm(i, BUSNM_VN)) then
            call msg(MSG_ERROR, '(gesdict) name conflict: ON='//trim(outname)//' already used as VN= '//trim(busnm(i, BUSNM_VN))//' with ON='//trim(busnm(i, BUSNM_ON)))
            return
         endif
      end do
      F_istat = RMN_OK
      !-------------------------------------------------------------------
      return
   end function gesdictcheck

end module gesdictmod


!/@*
subroutine gesdict(ni, nk, lindex, lachaine)
   use phy_options, only: debug_mem_L
   use gesdictmod, only: gesdictislock, gesdictadd, gesdictcheck
   use splitst, only: splitst4
   implicit none
!!!#include <arch_specific.hf>
   !@Arguments
   ! n          horizontal dimension
   ! nk         vertical dimension
   ! lindex     starting index on the bus
   ! lachaine   string identifying the variable attributes
   character(len=*), intent(in) :: lachaine
   integer, intent(in)  :: ni, nk
   integer, intent(out) :: lindex

   !@Author M. Desgagne (Oct 1995)
   !@Revision
   ! 001      B. Bilodeau (Jan 1996) - Check name conflicts for
   !                                   a given description
   ! 002      B. Bilodeau (Sep 1996) - Add 2-letter names
   ! 003      B. Bilodeau (Aug 1998) - Add staggered levels
   ! 004      B. Bilodeau (Dec 1998) - Add entry bus
   ! 005      B. Bilodeau (Feb 1999) - Add fmul to entpar, dynpar,
   !                                   perpar and volpar
   ! 006      G. Bergeron (Oct 1999) - Test if top < maxbus
   ! 007      B. Bilodeau (Mar 2000) - Test conflicting output names
   !                                   for a given variable
   ! 008      B. Bilodeau (Feb 2004) - 4-letter output names and
   !                                  16-letter names
   ! 009      B. Bilodeau (Mar 2005) - Test conflicting variable names
   !                                   and output names
   ! 010      B. Bilodeau (Jun 2005) - Forbid SLB*n and SLS*n for n > 1
   !                                   Add mosaic capability for CLASS
   ! 011      V. Lee      (Mar 2011) - nmosaic=real number of mosaic tiles
   !                                 - entpar(*,8),perpar(*,8),..=0 if no mosaic
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
   !*@/
#include <rmnlib_basics.hf>
#include <msg.h>
   include "buses.cdk"

   integer, parameter :: MAXVALPERCHAR = 34

   logical, save :: init_L = .false.

   integer :: memgap
   character(len=1) ::   bus
   character(len=4) ::   outname,inname,sername
   character(len=3) ::   shape
   character(len=7) ::   struc
   character(len=16) ::  varname
   character(len=48) ::  vdescrp
   character(len=60) ::  vardesc
   character(len=256) :: string
   integer :: nmosaic, fmul, dynini, stagg
   integer :: i, ivmin, ivmax, wload, hzd, monot, massc, istat
   real :: vmin, vmax
   !-------------------------------------------------------------------
   if (.not.init_L) then
      init_L = .true.
      enttop = 0; dyntop = 0; pertop = 0; voltop = 0
      entspc = 0; dynspc = 0; perspc = 0; volspc = 0
   endif

   if (gesdictislock()) then
      call physeterror('gesdict', 'Call to gesdict after buses have been allocated')
      return
   endif

   call low2up(lachaine,string)
   istat = splitst4(varname, outname, inname, sername, &
        vdescrp, struc, shape, nmosaic, fmul, bus, dynini,&
        stagg, vmin, vmax, wload, hzd, monot, massc, string)
   if (.not.RMN_IS_OK(istat)) then
      call physeterror('gesdict', 'Invalid gesdict string: '//trim(string))
      return
   endif

   !# Mosaic field have an extra row for the averaged layer (1st row)
   vardesc = vdescrp//';VS='//struc

   i = MAXVALPERCHAR ** (len(outname) - max(len_trim(outname), len_trim(inname)))
   if (any(shape == (/'M', 'T', 'E'/)) .and. fmul > i) then
      call msg(MSG_WARNING,'Varname '//trim(varname)//' is declared as a 3D var with multiple categories, fmul>1, I/O names must be 2/3 char max.')
   endif

   ivmin = transfer(vmin,ivmin)
   ivmax = transfer(vmax,ivmax)

   memgap = 0
   lindex  = RMN_ERR
   select case(bus)
   case("E")
      if (debug_mem_L) memgap = max(ni/5,1)+8
      !#TODO: should all var have a different gap? possible using enttop
      lindex = gesdictadd( &
        varname, outname, inname, sername, &
        vardesc, shape, dynini, stagg, fmul, nmosaic, ivmin, ivmax, &
        wload, hzd, monot, massc, &
        enttop, entnm, entdc, entpar, entspc, ni, nk, memgap)
   case("D")
      if (debug_mem_L) memgap = max(ni/5,1)+4
      lindex = gesdictadd( &
        varname, outname, inname, sername, &
        vardesc, shape, dynini, stagg, fmul, nmosaic, ivmin, ivmax, &
        wload, hzd, monot, massc, &
        dyntop, dynnm, dyndc, dynpar, dynspc, ni, nk, memgap*2)
   case("P")
      if (debug_mem_L) memgap = max(ni/5,1)+2
      lindex = gesdictadd( &
        varname, outname, inname, sername, &
        vardesc, shape, dynini, stagg, fmul, nmosaic, ivmin, ivmax, &
        wload, hzd, monot, massc, &
        pertop, pernm, perdc, perpar, perspc, ni, nk, memgap*3)
   case("V")
      if (debug_mem_L) memgap = max(ni/5,1)
      lindex = gesdictadd( &
        varname, outname, inname, sername, &
        vardesc, shape, dynini, stagg, fmul, nmosaic, ivmin, ivmax, &
        wload, hzd, monot, massc, &
        voltop, volnm, voldc, volpar, volspc, ni, nk, memgap*4)
   case default
      call physeterror('gesdict', 'Unknown bus: '//trim(string))
      return
   end select
   if (.not.RMN_IS_OK(lindex)) then
      call physeterror('gesdict', 'Problem getting lindex')
      return
   endif

   istat = gesdictcheck(bus, varname, outname, &
        enttop, entnm, 'E')
   istat = min(gesdictcheck(bus, varname, outname, &
        dyntop, dynnm, 'D'), istat)
   istat = min(gesdictcheck(bus, varname, outname, &
        pertop, pernm, 'P'), istat)
   istat = min(gesdictcheck(bus, varname, outname, &
        voltop, volnm, 'V'), istat)
   if (.not.RMN_IS_OK(istat)) then
      call physeterror('gesdict', 'Problem while checking values')
      return
   endif
   !-------------------------------------------------------------------
   return
end subroutine gesdict


