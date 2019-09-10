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

module outgrid
   implicit none
   public
   save
!
!revision
! v2_30 - V. Lee                 - moved definition of OutGrid_MAXGRID1 to dimout.cdk
! v2_30                          - Added Grid_phi_ig2
! v3_30 - R. McTaggart-Cowan     - Added user defined tag extension
!______________________________________________________________________
!                                                                      |
!  VARIABLES FOR DEFINITION OF THE OUTPUT GRIDS (set_grid)             |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! OutGrid_MAXGRID1   | maximum number of grids that can be defined     |
! OutGrid_sets       | total number of sets of defined output grids    |
!  The following variables carry values for each defined output grid   |
! OutGrid_id         | OutGrid_id(i) are the id of each defined grid   |
! OutGrid_x0         | x origin of output grid                         |
! OutGrid_y0         | y origin of output grid                         |
! OutGrid_x1         | x at outermost corner of output grid            |
! OutGrid_y1         | y at outermost corner of output grid            |
! OutGrid_stride     | every ith point to be outputted                 |
!----------------------------------------------------------------------
!

   integer, parameter :: OutGrid_MAXGRID1 = 4

   logical :: OutGrid_reduc(OutGrid_MAXGRID1)
   integer :: OutGrid_x0 (OutGrid_MAXGRID1), OutGrid_x1    (OutGrid_MAXGRID1)
   integer :: OutGrid_y0 (OutGrid_MAXGRID1), OutGrid_y1    (OutGrid_MAXGRID1)
   integer :: OutGrid_id (OutGrid_MAXGRID1), OutGrid_stride(OutGrid_MAXGRID1)
   integer :: OutGrid_sets

end module outgrid
