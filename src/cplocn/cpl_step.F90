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

module cpl_step_mod

   private
   public :: cpl_step

contains

      subroutine cpl_step (F_stepcount, F_stepdriver)
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

      integer unout
      logical print_L

      logical first
      integer err

      data first/.true./
      save first

      integer, external :: msg_getUnit
!     ________________________________________________________________
!

      if ( (.not. cpl_ocn_L .and. .not. cpl_wav_L) &
           .or. F_stepdriver < 0 .or. cpl_dgflt_H ) return

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)

      if ( cpl_ocn_L ) &
        call cplocn_step (F_stepcount, F_stepdriver)

!      if ( cpl_wav_L ) &
!        call cplwav_step (F_stepcount, F_stepdriver)

      goto 999

 998  call handle_error(err,'cpl_step','Problems')
 999  continue

      first = .false.
!
!     ________________________________________________________________
!
      return
      end subroutine cpl_step

end module cpl_step_mod


