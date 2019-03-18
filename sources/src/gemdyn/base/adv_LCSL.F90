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

!**s/r adv_LCSL: Various options of Local Conservative Semi-Lagrangian (LCSL) Advection

      subroutine adv_LCSL (F_cub, F_ni, F_nj, F_nk)

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: F_ni,F_nj,F_nk                                          !Dimension of position fields
      real,   dimension(F_ni,F_nj,F_nk),intent(out):: F_cub                          !High-order SL solution

      !object
      !============================================================================
      !     Various options of Local Conservative Semi-Lagrangian (LCSL) Advections
      !============================================================================
      F_cub = 0. ! shut up compiler error with -warn all
      call handle_error(-1,'ADV_LCSL','F_conserv_local NOT AVAILABLE')

      return

      end subroutine adv_LCSL
