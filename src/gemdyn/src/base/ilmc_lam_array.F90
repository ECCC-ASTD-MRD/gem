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

module ilmc_lam_array
   implicit none
   public
   save
   type a_ilmc
      ! Comment the following line to let the compiler reorganize structure and avoid misalignments
      !sequence
      integer, pointer, dimension(:):: i_rd ! I indice of neighbor in a given sweep
      integer, pointer, dimension(:):: j_rd ! J indice of neighbor in a given sweep
      integer, pointer, dimension(:):: k_rd ! J indice of neighbor in a given sweep
      integer                       :: cell ! # neighbors in a given sweep
   end type a_ilmc

   type (a_ilmc), dimension(:,:,:,:), pointer :: sweep_rd
end module ilmc_lam_array
