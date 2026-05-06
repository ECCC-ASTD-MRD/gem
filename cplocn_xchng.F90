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

SUBROUTINE cplocn_send(F_date, F_send, F_stepdriver)
    use app
    use iris_mod
    use rmn_gmm
    use cpl_timers_mod
    use iso_c_binding
    use cpl_mod
    use cplocn_mod
    implicit none
#include <arch_specific.hf>
!
!authors    Francois Roy - Fall 2014
! 
!revision
! v4_7  - Roy F.  - initial version

! Purpose: Performs exchange with ocean model
!          through gossip interface and apply 
!          distributed interpolation weights 
!          (including angles and vector rotation)

#include <rmn/msg.h>

      ! arguments
      character(len=16), intent(in) :: F_date
      integer, intent(in ) :: F_stepdriver
      integer, intent(out) :: F_send

      ! locals
      integer :: stamp

      integer unout
      logical print_L
      integer(kind(IRIS_PROCESS_NONE)) :: iris_process

      integer, external :: msg_getUnit

      integer ier
      real(C_FLOAT) :: level

!     ________________________________________________________________
!
      if(cplocn_myproc == 0) then
          call cplocn_iris_share_timer%start()
      endif


      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)
      iris_process = IRIS_PROCESS_NONE
      if (cplocn_iweight_L) then
         iris_process = IRIS_PROCESS_IWEIGHT
      endif
      level = 0.0

! Send step
      F_send = -1
      call datp2f (stamp,F_date)
      ier = iris%model_share(iris_grid, cplocn_iris_provides, level, int(stamp, C_LONG), C_LOC(ocn_busou))
      if(ier /= 0) then
          write(app_msg,'("Error in iris%model_share: ",i3)') ier
          call app_log(APP_ERROR, app_msg)
          call gem_error(ier)
      endif

      F_send = 0

      if(cplocn_myproc == 0) then
          call cplocn_iris_share_timer%stop()
          write (*,*) "TIMING, cplocn_iris_share, ", cplocn_iris_share_timer%get_latest_time_ms(), " (ms), step=", F_stepdriver, ",", F_date(1:15)
      endif

END SUBROUTINE cplocn_send

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE cplocn_recv(F_date, F_recv, F_stepdriver)
    use app
    use iris_mod
    use rmn_gmm
    use cpl_timers_mod
    use cpl_mod
    use cplocn_mod
    use iso_c_binding
    implicit none

#include <rmn/msg.h>

      ! arguments
      character(len=16), intent(in) :: F_date
      integer, intent(in ) :: F_stepdriver
      integer, intent(out) :: F_recv

      ! locals
      integer :: stamp
      integer unout
      logical print_L
      integer(kind(IRIS_PROCESS_NONE)) :: iris_process
      real(C_FLOAT) :: level

      integer, external :: msg_getUnit

      integer ivar,ni,nj,i,j
      integer ier

!     ________________________________________________________________
!
      if(cplocn_myproc == 0) then
          call cplocn_iris_need_timer%start()
      endif

      ier = gmm_get(gmmk_ocn_busin_s , ocn_busin)

      ni=cpl_drv_lni
      nj=cpl_drv_lnj

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)
      iris_process = IRIS_PROCESS_NONE
      if (cplocn_iweight_L) then
         iris_process = IRIS_PROCESS_IWEIGHT
      endif
      level = 0.0

! Receive step
      F_recv = -1
      call datp2f (stamp,F_date)
      ier = iris%model_need(iris_grid, cplocn_iris_consumes, level, int(stamp, C_LONG), C_LOC(ocn_busin), iris_process)
      if(ier /= 0) then
          write(app_msg,'("Error in iris%model_need: ",i3)') ier
          call app_log(APP_ERROR, app_msg)
          call gem_error(ier)
      endif

      F_recv = 0

      ! do some clean up of the ocn_busin
      do ivar = 1, cplocn_n_fldin
        do j=1,nj
           do i=1,ni
              if(isnan(ocn_busin(i,j,ivar))) then
                 ocn_busin(i,j,ivar) = 0.
              endif
           enddo
        enddo
        if ( cplocn_cvin_S(ivar) == 'MCP' ) then
          do j=1,nj
            do i=1,ni
              ! Avoid negative values of tau later - Fred Dupont
              ! coupling mask should be always bounded between 0 and 1, can be fractional
              ! used in cplocn_update.F90 when updating the water/ice/surface conditions
              ocn_busin(i,j,ivar) = min ( max ( ocn_busin(i,j,ivar), 0.), 1.)
            enddo
          enddo
        endif
      enddo

      if(cplocn_myproc == 0) then
          call cplocn_iris_need_timer%stop()
          write (*,*) "TIMING, cplocn_iris_need, ", cplocn_iris_need_timer%get_latest_time_ms(), " (ms), ", F_stepdriver, ",", F_date(1:15)
      endif

END SUBROUTINE cplocn_recv
