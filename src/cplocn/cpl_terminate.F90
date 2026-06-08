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

module cpl_terminate_mod

   private
   public :: cpl_terminate

contains

      subroutine cpl_terminate (F_pcomm_L)
      use iris_mod
      use cpl_mod
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_pcomm_L

!authors    Francois Roy -- spring 2015
! 
!revision
! v4_80 - Roy, F.  - initial version

!Purpose
! Terminate coupler

!     ---------------------------------------------------------------
!
      if (cpl_ocn_L) then
          call iris%model_finalize()
      endif

!      if (cpl_wav_L) then
!         call cpl_aw_terminate (F_pcomm_L)
!      endif
!
!     ---------------------------------------------------------------
!
      return
      end subroutine cpl_terminate

end module cpl_terminate_mod
