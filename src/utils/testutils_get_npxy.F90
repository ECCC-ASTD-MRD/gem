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
subroutine testutils_get_npxy(F_npex,F_npey)
   use testutils
   implicit none
   !@objective init fn to be passed to rpn_comm_init_multigrid
   !@author  Stephane Chamberland, 2010-05
   !@arguments
   integer,intent(out) :: F_npex,F_npey
   !*@/
#include "rmnlib_basics.hf"
   integer :: err
   !---------------------------------------------------------------------
   err = testutils_getenv_int('MPI_NPEX',F_npex)
   if (.not.RMN_IS_OK(err)) F_npex = 1
   err = testutils_getenv_int('MPI_NPEY',F_npey)
   if (.not.RMN_IS_OK(err)) F_npey = F_npex
   !---------------------------------------------------------------------
   return
end subroutine testutils_get_npxy


!/@*
subroutine testutils_get_doms(F_ndomains,F_offset,F_istat)
   use testutils
   implicit none
   !@objective init fn to be passed to rpn_comm_mydomain
   !@author  Stephane Chamberland, 2012-01
   !@arguments
   integer, intent(out) :: F_ndomains,F_offset,F_istat
   !*@/
#include "rmnlib_basics.hf"
   integer :: err
   ! ---------------------------------------------------------------------
   err = testutils_getenv_int('MPI_NDOM',F_ndomains)
   if (.not.RMN_IS_OK(err)) F_ndomains = 1
   err = testutils_getenv_int('MPI_IDOM',F_offset)
   if (.not.RMN_IS_OK(err)) F_offset = 0
   F_istat = 0
   ! ---------------------------------------------------------------------
   return
end subroutine testutils_get_doms
