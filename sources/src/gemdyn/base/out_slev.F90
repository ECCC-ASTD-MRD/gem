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

!**s/r out_slev
!
      subroutine out_slev2 (F_level,F_levmax,F_nk,F_indo,F_nko,F_nearsfc_L)
      implicit none
#include <arch_specific.hf>

      integer F_nk,F_indo(*),F_nko, F_levmax
      real F_level(*)
      logical F_nearsfc_L

!AUTHOR   Michel Desgagne     July 2004 (MC2)
!
!REVISION
! v3_20 - Lee V.            - Adapted for GEMDM
! v4_40 - Lee V.            - change in argument call for this routine
!                             in order to select the "bot" levels
!ARGUMENTS
!    NAMES       I/O  TYPE  DESCRIPTION

      integer k,kk
!
!----------------------------------------------------------------------
!
!     initialize output level array to 0
!     output in model levels for dynamic variables
      kk=0
      do k=1,F_levmax
         if (F_level(k) == 0.0) then
!     ground zero (surface or lowest possible level in model)
             kk=kk+1
             F_indo(kk)=F_nk
         else if (F_level(k) < 0.0) then
             kk=kk+1
!     ground zero with surfaces above it
             F_indo(kk)=F_nk + nint(F_level(k))
         else
!     model level as defined by user. Discard levels that do not exist
             if (F_level(k) <= F_nk)then
                 kk=kk+1
                 F_indo(kk)=nint(F_level(k))
             end if
         end if
      end do
      F_nko=kk

      ! Diagnostic level additions
      F_nearsfc_L = .false.
      if ( (F_nko == F_nk) .or. &
              any(F_indo(1:F_nko)==F_nk)) F_nearsfc_L = .true.

      if (F_nearsfc_L .and. all(F_level(1:F_levmax)==0.0)) F_nko = 0
!
!----------------------------------------------------------------------
!
      return
      end

