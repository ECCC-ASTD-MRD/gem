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

!**s/r theo_bndry

      subroutine theo_bndry ()
      use mtn_options
      use theo_options
      implicit none
#include <arch_specific.hf>
!
!----------------------------------------------------------------------
!
      if ( mtn_zblen_thk > 0. ) call height_sponge()
      
      call slabsym()

      if( Theo_case_S == 'BUBBLE' ) call mirror()
!
!----------------------------------------------------------------------
!
      return
      end subroutine theo_bndry
