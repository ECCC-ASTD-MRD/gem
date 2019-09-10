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
!**s/r bubble_cfg - Configure bubble case (Height or pressure coord)
      integer function bubble_cfg()
      use dynkernel_options
      use bubble_options
      use HORgrid_options
      use VERgrid_options

      implicit none

#include <arch_specific.hf>

      integer,external :: bubble_H_cfg, bubble_P_cfg
!     ---------------------------------------------------------------
!
      bubble_cfg = -1
      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' .or. &
          trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H') then
          bubble_cfg = bubble_H_cfg()
      else
          bubble_cfg = bubble_P_cfg()
      end if

      if (bubble_ictr < 0) bubble_ictr = int(float(Grd_ni-1)*0.5)+1
      if (bubble_kctr < 0) bubble_kctr = G_nk - bubble_rad - 1

      bubble_cfg = 1

      return
      end function bubble_cfg
