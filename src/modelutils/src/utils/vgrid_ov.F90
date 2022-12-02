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

!/@*
module vgrid_ov
   use iso_c_binding
   use vGrid_Descriptors
   implicit none
   !@author Stephane Chamberland, 2018-10
   !@Objective Some vgrid overrides
   private
   public :: vgrid_nullify, vgrid_associated
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

contains

   !/@*
   subroutine vgrid_nullify(F_vgd)
      implicit none
      !@arguments
      type(vgrid_descriptor),intent(out) :: F_vgd
      !*@/
      !---------------------------------------------------------------
      F_vgd%cptr = C_NULL_PTR
      !---------------------------------------------------------------
      return
   end subroutine vgrid_nullify


   !/@*
   function vgrid_associated(F_vgd) result(F_istat_L)
      implicit none
      !@objective Store a new vgrid
      !@arguments
      type(vgrid_descriptor),intent(in) :: F_vgd
      !@returns
      logical :: F_istat_L
      !*@/
      !---------------------------------------------------------------
      F_istat_L = c_associated(F_vgd%cptr)
      !---------------------------------------------------------------
      return
   end function vgrid_associated

end module vgrid_ov
