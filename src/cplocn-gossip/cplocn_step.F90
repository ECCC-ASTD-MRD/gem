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

      subroutine cplocn_step (F_stepcount, F_stepdriver)
      use rmn_gmm
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

      include "cpl.cdk"
      include "cplocn.cdk"

      integer unout
      logical print_L

      character*16 datev,date_in,date_ou
      logical flag
      integer err,i,j,send,recv,n,ni,nj,ivar
      integer nstp_st
      real strt_dt

      real*8  dayfrac,one,sid,rsid,secs_elapsed
      parameter(one=1.0d0, sid=86400.0d0, rsid=one/sid)

      integer, external :: msg_getUnit
!     ________________________________________________________________
!

      cplocn_mystep=F_stepdriver

      err = gmm_get(gmmk_ocn_busin_s , ocn_busin)

      secs_elapsed=dble(F_stepcount)*cpl_drv_delt

      ni = cpl_drv_lni
      nj = cpl_drv_lnj

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)

! Pack information to send to ocean model

      call cplocn_busou ()

! Performs stepping

      nstp_st = cplocn_oc_dt/cpl_drv_delt
      strt_dt = cplocn_oc_dt
      if (cplocn_oc_dt.lt.cpl_drv_delt) then
        nstp_st = 1
        strt_dt = cpl_drv_delt
      endif

      if (F_stepdriver.eq.0) then

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

! ddeacu: modification for IAU
          datev   = cplocn_runstrt_S
         dayfrac = secs_elapsed*rsid
         call incdatsd  (datev,cplocn_runstrt_S,dayfrac)
         date_in = datev

         ! Export date
         date_ou = datev

         call cplocn_xchng (date_ou, date_in, send, recv, F_stepdriver, err)

         if (err < 0) then
           if (print_L) &
            write(unout,*) 'cplocn_step 1: error received from cplocn_xchng'
           goto 998
         endif

         if (recv /= 0) then
           if (print_L) write(unout,*) &
            'cplocn_step 1: After cplocn_xchng recv =',recv, &
            ' -> Dummy communication step. No update of ocn_busin'
         endif
         if (recv == 0) then
           if (print_L) write(unout,*) &
            'cplocn_step 1: After cplocn_xchng recv =',recv, &
            ' -> Coupling communication should be dummy: will abort'
           err = -1
           goto 998
         endif

         !*         Second exchange with time offset for full ice-ocean initialization
         !*         (ocean fields used for the entire first atm sequence based on strt_dt)
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

! ddeacu: modification for IAU
          dayfrac = strt_dt*rsid
         dayfrac = (secs_elapsed + strt_dt)*rsid
         call incdatsd  (datev,cplocn_runstrt_S,dayfrac)
         date_in = datev

         ! Export date
         date_ou(1:15) = '00000000.000000'

         call cplocn_xchng (date_ou, date_in, send, recv, F_stepdriver, err)

         if (err < 0) then
           if (print_L) write(unout,*) &
            'cplocn_step 2: error received from cplocn_xchng'
           goto 998
         endif

         if (recv /= 0) then
           if (print_L) write(unout,*) &
            'cplocn_step 2: After cplocn_xchng recv =',recv, &
                      ' -> Coupling communication should not be dummy: will abort'
           err = -1
           goto 998
         endif
         if (recv == 0) then
           if (print_L) write(unout,*) &
            'cplocn_step 2: After cplocn_xchng recv =',recv, &
            ' -> Useful communication step. Using ocn_busin'
         endif

      elseif (F_stepdriver.gt.nstp_st) then

         !*       Two-way main exchange at the end of atm time step
         !*       NEMO sbc module is called at the beginning of full exchange time step
         !*       (receive may be dummy)
         !*
         !*                       n
         !*  (atm)             ------->
         !*                          |^
         !*                          ||
         !*                          ||
         !*                          ||
         !*                          ||
         !*                          v|
         !*  (oce)                   ------------------->
         !*                                   n

         ! Import date (end of atm time)

! ddeacu: modification for IAU
!         dayfrac = dble(F_stepdriver)*cpl_drv_delt*rsid
         dayfrac = secs_elapsed*rsid
         call incdatsd  (datev,cplocn_runstrt_S,dayfrac)
         date_in = datev

         ! Export date

         date_ou = date_in

         call cplocn_xchng (date_ou, date_in, send, recv, F_stepdriver, err)

         if (err < 0) then
           if (print_L) write(unout,*) &
            'cplocn_step: error received from cplocn_xchng'
           goto 998
         endif

         if (recv /= 0) then
           if (print_L) write(unout,*) &
            'cplocn_step: After cplocn_xchng recv =',recv, &
            ' -> Dummy communication step. Re-using ocn_busin'
         endif
         if (recv == 0) then
           if (print_L) write(unout,*) &
            'cplocn_step: After cplocn_xchng recv =',recv, &
            ' -> Useful communication step. Using ocn_busin'
         endif

      endif

      do ivar = 1, cplocn_n_fldin
        if ( cplocn_cvin_S(ivar) == 'MCP' ) then
          do j=1,nj
          do i=1,ni
            ! Avoid negative values of tau later - Fred Dupont
            ocn_busin(i,j,ivar) = min ( max ( ocn_busin(i,j,ivar), 0.), 1.)
            if ( ocn_mwgt(i,j) .eq. 0 ) ocn_busin(i,j,ivar) = 0.
               ! only coupling zone but keep fractional part of ocn_busin
               ! if so
            if ( cplocn_off_L ) then
              if (nint(ocn_busin(i,j,ivar)).ne.0) then
                if (print_L) write(unout,*) &
                 'cplocn_step - MASK Inconsistent - cplocn_off_L - ABORT -'
                err = -1
                goto 998
              endif
            endif
          enddo
          enddo
          goto 997
        endif
      enddo
997   continue

      goto 999

 998  call handle_error(err,'cplocn_step','Problems')
 999  continue

!
!     ________________________________________________________________
!
      return
      end subroutine cplocn_step
