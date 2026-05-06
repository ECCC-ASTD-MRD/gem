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

      subroutine cplocn_step (F_stepcount, F_stepdriver, tag)
      use cplocn_mod, only : cplao_xchg_mode
      implicit none

#include <rmn/msg.h>

      integer unout
      logical print_L

      integer, intent(in) :: F_stepcount, F_stepdriver
      character(len=5), intent(in) :: tag

!     locals
      integer nsend_offset

      integer, external :: msg_getUnit
!     ________________________________________________________________
!

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)

      select case( tag )
      case( 'first' )  ! step=0 special operations for cplao_xchg_mode=0,1
         if (cplao_xchg_mode < 2) then
            call cplocn_step_frst_call( F_stepcount, F_stepdriver )
         endif
      case( 'ssend' ) ! send at the start of the phy step, only used by cplao_xchg_mode=0
         if (cplao_xchg_mode < 1) then
            nsend_offset = -1
            call cplocn_step_send_call( F_stepcount, F_stepdriver, nsend_offset )
         endif
      case( 'esend' ) ! send at the end of the phy step
         if (cplao_xchg_mode > 0) then
            nsend_offset = 0
            call cplocn_step_send_call( F_stepcount, F_stepdriver, nsend_offset )
         endif
      case( 'srecv' ) ! receive outside step=0
         call cplocn_step_recv_call( F_stepcount, F_stepdriver ) ! receive only
      case default
         if (print_L) write(unout,*) 'cplocn_step : operation not recognized: ',tag
         goto 998
      end select

      return
 998  call handle_error(-1,'cplocn_step','Problems')
      end subroutine cplocn_step


      subroutine cplocn_step_frst_call (F_stepcount, F_stepdriver)
      use rmn_gmm
      use cpl_mod
      use cplocn_mod
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_stepcount, F_stepdriver

!authors    Francois Roy - Spring 2015
! 
!revision
! v4_80 - Roy F.          - initial version

!Purpose
! Performs coupled stepping with ocean model

#include <rmn/msg.h>

      integer unout
      logical print_L

      character(len=16) :: datev,date_in,date_ou
      logical flag
      integer err,send,recv

      real(kind=8) ::  dayfrac,one,sid,rsid,secs_elapsed
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)

      integer, external :: msg_getUnit
!     ________________________________________________________________
!

      if ( F_stepdriver > 0 ) return ! only used at phy step=0

      err = gmm_get(gmmk_ocn_busin_s , ocn_busin)

      secs_elapsed=dble(F_stepcount)*cpl_drv_delt

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)

! Performs stepping

      !*         First exchange for ocean to advance in time
      !*              0
      !* (atm)   ----------->
      !*                    |
      !*                    |
      !*                    |
      !*                    |
      !*                    |
      !*                    v
      !* (oce)              ----------------------------->

      ! Import date

! Pack information to send to ocean model

      call cplocn_busou ()

      dayfrac = secs_elapsed*rsid
      call incdatsd  (datev,cplocn_runstrt_S,dayfrac)
      call cplocn_send(datev, send, F_stepdriver)

      ! reset to zero
      cplocn_it = 0
      ocn_busou(:,:,:) = 0.

      if (send < 0) then
        if (print_L) write(unout,*) 'cplocn_step=0 : error sending'
        goto 998
      endif

      !*         Second exchange with time offset for full ice-ocean initialization
      !*         (ocean fields used for the entire first atm sequence)
      !*           0        1
      !* (atm)   ----->----------->
      !*              ^
      !*                 -
      !*                      -
      !*                          -
      !*                              -
      !*                                   -
      !*
      !* (oce)     ----------------------------->
      !*                         1


      ! Import date at end of ocean time step

      call cplocn_recv(datev, recv, F_stepdriver)

      if (recv < 0) then
        if (print_L) write(unout,*) 'cplocn_step=0 : error receiving'
        goto 998
      endif

      goto 999

 998  call handle_error(err,'cplocn_step_frst_call','Problems')
 999  continue

!
!     ________________________________________________________________
!
      return
      end subroutine cplocn_step_frst_call

      subroutine cplocn_step_recv_call (F_stepcount, F_stepdriver)
      use rmn_gmm
      use cpl_mod
      use cplocn_mod
      implicit none

      integer, intent(in) :: F_stepcount, F_stepdriver

      integer unout
      logical print_L

      character(len=16) :: datev,date_in,date_ou
      logical flag
      integer err,send,recv
      integer noffset

      real(kind=8) ::  dayfrac,one,sid,rsid,secs_elapsed
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)

      integer, external :: msg_getUnit
!     ________________________________________________________________
!

! Performs stepping

      if (F_stepdriver <= 1 .or. mod(F_stepdriver-1,cplocn_rap_dt) /= 0) return

      err = gmm_get(gmmk_ocn_busin_s , ocn_busin)

      noffset = -1
      secs_elapsed=dble(F_stepcount+noffset)*cpl_drv_delt

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)


      !*       Two-way main exchange at the end of atm time step
      !*       NEMO sbc module is called at the beginning of full exchange time step
      !*
      !*                       n
      !*  (atm)             ------->
      !*                           ^
      !*                           |
      !*                           |
      !*                           |
      !*                           |
      !*                           |
      !*  (oce)                   ------------------->
      !*                                   n

      ! Import date (end of atm time)

      dayfrac = secs_elapsed*rsid
      call incdatsd  (datev,cplocn_runstrt_S,dayfrac)

      call cplocn_recv(datev, recv, F_stepdriver)
      if (recv < 0) then
        if (print_L) write(unout,*) 'cplocn_step=',F_stepdriver,' : error receiving'
        goto 998
      endif

      goto 999

 998  call handle_error(err,'cplocn_step_recv_call','Problems')
 999  continue

!
!     ________________________________________________________________
!
      return
      end subroutine cplocn_step_recv_call

      subroutine cplocn_step_send_call (F_stepcount, F_stepdriver, noffset)
      use rmn_gmm
      use cpl_mod
      use cplocn_mod
      implicit none

      integer, intent(in) :: F_stepcount, F_stepdriver, noffset

#include <rmn/msg.h>

      integer unout
      logical print_L

      character(len=16) :: datev,date_in,date_ou
      logical flag
      integer err,i,j,send,recv

      real(kind=8) ::  dayfrac,one,sid,rsid,secs_elapsed
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)

      integer, external :: msg_getUnit
 !     ________________________________________________________________
 !

      ! Performs stepping

      if (F_stepdriver <  1 .and. cplao_xchg_mode == 1 ) return
      if (F_stepdriver <= 1 .and. cplao_xchg_mode == 0 ) return

      err = gmm_get(gmmk_ocn_busin_s , ocn_busin)

      secs_elapsed=dble(F_stepcount+noffset)*cpl_drv_delt

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)


      call cplocn_busou ()

      !*       Two-way main exchange at the end of atm time step
      !*       NEMO sbc module is called at the beginning of full exchange time step
      !*
      !*                       n
      !*  (atm)             ------->
      !*                          |
      !*                          |
      !*                          |
      !*                          |
      !*                          |
      !*                          v
      !*  (oce)                   ------------------->
      !*                                   n

      ! Import date (end of atm time)

      dayfrac = secs_elapsed*rsid
      call incdatsd  (datev,cplocn_runstrt_S,dayfrac)

      if ( cplocn_it == 0 .or. F_stepdriver == 0 ) then ! we are ready to do the sending

         call cplocn_send(datev, send, F_stepdriver)
         if (send < 0) then
            if (print_L) write(unout,*) 'cplocn_step=',F_stepdriver,' : error sending'
            goto 998
         endif

         ! reset array for next averaging
         ocn_busou(:,:,:) = 0
         cplocn_it = 0 ! only needed for cplao_xchg_mode=2

      endif ! cplocn_it == 0

      return

 998  call handle_error(err,'cplocn_step_scnd_step','Problems')

      end subroutine cplocn_step_send_call
