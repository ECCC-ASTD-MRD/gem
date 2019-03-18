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

!**s/r tracers_step

      subroutine tracers_step (F_before_psadj_L)
      use adv_options
      use lun
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_before_psadj_L

! local variables
      logical :: verbose_L
!     _________________________________________________________________
!

      verbose_L = Adv_verbose/=0

      if (Lun_debug_L) write (Lun_out,1000) F_before_psadj_L

      if (verbose_L .and. .not.F_before_psadj_L) &
         call stat_mass_tracers (1,"BEFORE ADVECTION")

      if (.not. F_before_psadj_L) Adv_pos_reset(1:2) = 1

      call adv_tracers (F_before_psadj_L)

      if (verbose_L .and. .not.F_before_psadj_L) &
         call stat_mass_tracers (0,"AFTER ADVECTION")

1000  format(3X,'ADVECT TRACERS: (S/R TRACERS_STEP) H2O only=',L2)
!     _________________________________________________________________
!
      return
      end
