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

subroutine inichamp4(kount, trnch, ni, nk)
   use sfc_options
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: kount, trnch , ni, nk

   !@Author B. Bilodeau (July 1997)
   !@Revisions
   ! 001 K. Winger      (Feb 2017) - Add initialization for lake fraction and rivers (m.a.)
   ! 002 M. Mackay      (Oct 2018/Sep 2022) - Code added for CSLM
   !@Object initialize arrays.
   !@Arguments
   !          - Input -
   ! kount    timestep number
   ! trnch    row number
   ! ni       horizontal dimension
   ! nk       vertical dimension

#include <rmn/msg.h>
   include "sfcinput.cdk"

#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define PTR1D(NAME2) busptr(vd%NAME2%i)%ptr(1,trnch)

   integer :: i
   real, dimension(ni) :: land
   real, pointer, dimension(:)   :: zglacier, zglsea, zh, zlst, &
        zmg, zsnoden, ztglacier, ztmice, ztsrad, ztwater
   real, pointer, dimension(:,:) :: zhst, ztmoins
   !***********************************************************************
   call msg_toall(MSG_DEBUG, 'inichamp [BEGIN]')
   if (timings_L) call timing_start_omp(401, 'inichamp', 46)

   nullify(zglacier, zglsea, zh, zmg, zsnoden, ztglacier, &
        ztmice, ztsrad, ztwater,zhst, ztmoins)

   PRE_INIT: if (kount == 0) then

      MKPTR1D(zglacier, glacier)
      MKPTR1D(zglsea, glsea)
      MKPTR1D(zh, h)
      MKPTR1D(zlst, lst)
      MKPTR1D(zmg, mg)
      MKPTR1D(ztglacier, tglacier)
      MKPTR1D(ztmice, tmice)
      MKPTR1D(ztsrad, tsrad)
      MKPTR1D(ztwater, twater)

      MKPTR2D(zhst, hst)
      MKPTR2D(ztmoins, tmoins)
      
   endif PRE_INIT

   ! Initialization of surface fields
   call inisurf4(kount, ni, nk, trnch)

   ! for slope only
!VDIR NODEP
   if (radslope) then
      if (any('fsa' == phyinread_list_s(1:phyinread_n)) .or. &
           any('sla' == phyinread_list_s(1:phyinread_n))) then
         call radcons2(ni, trnch)
      endif
   endif

   POST_INIT: if (kount == 0) then

!vdir nodep
      do i=1,ni
         ! hauteur de la couche limite
         zhst(i,indx_soil   ) = 300.
         zhst(i,indx_glacier) = 300.
         zhst(i,indx_water  ) = 300.
         zhst(i,indx_ice    ) = 300.
         if (schmurb.ne.'NIL') then
            zhst(i,indx_urb ) = 300.
         endif
         if (schmlake.ne.'NIL') then
            zhst(i,indx_lake ) = 300.
         endif
         if (schmriver.ne.'NIL') then
            zhst(i,indx_river ) = 300.
         endif
         zh   (i) = 300.
         ! Initial radiative surface temperature estimate for the radiation scheme
         if (any('tsoil' == phyinread_list_s(1:phyinread_n))) then

           if (schmlake.ne.'NIL') then
            ztsrad(i) = zmg(i)*ztsrad(i) + (1.-zmg(i))*zlst(i)
           else
            ztsrad(i) = zmg(i)*ztsrad(i) + (1.-zmg(i))*ztwater(i)
           endif
            ztsrad(i) = (1.-zglsea(i)*(1.-zmg(i)))*ztsrad(i) + &
                 (zglsea(i)*(1.-zmg(i)))*ztmice(i)
            ztsrad(i) = (1.-zglacier(i)*zmg(i))*ztsrad(i) + &
                 (zglacier(i)*zmg(i))*ztglacier(i)
         endif
      end do

   endif POST_INIT

   NEW_TOPO: if (any('dhdxdy' == phyinread_list_s(1:phyinread_n)) .or. &
        any('dhdx' == phyinread_list_s(1:phyinread_n)) .or. &
        any('dhdy' == phyinread_list_s(1:phyinread_n)) .or. &
        any('mg' == phyinread_list_s(1:phyinread_n)) .or. &
        any('lhtg' == phyinread_list_s(1:phyinread_n))) then

      do i=1,ni
         land(i)  = - abs( nint( zmg(i) ) )
      enddo

      call equivmount1(land, PTR1D(dhdx), &
           PTR1D(dhdy), PTR1D(dhdxdy), &
           ni, 1, ni, PTR1D(slope), PTR1D(xcent), PTR1D(mtdir))

   endif NEW_TOPO

   if (timings_L) call timing_stop_omp(401)
   call msg_toall(MSG_DEBUG, 'inichamp [END]')

   return
end subroutine inichamp4
