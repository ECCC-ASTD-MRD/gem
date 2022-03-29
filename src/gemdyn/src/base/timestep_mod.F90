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

module timestep
   use dimout, only: MAXSET
   implicit none
   public
   save
!
!______________________________________________________________________
!                                                                      |
!  VARIABLES FOR DEFINITION OF THE FREQUENCY OF OUTPUTS (set_step)     |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! MAXSTEP            | maximum number of output timesteps per set      |
! Timestep_sets      | The total number of sets of defined output      |
!                    |                                      timesteps  |
! Timestep_id        | Timestep_id(i) are the ID of each defined set   |
! Timestep_tbl       | A table containing the defined output timesteps |
!                    | for each individual set.                        |
! Timestep_max       | Timestep_max(i) contains the number of timesteps|
!                    | defined in each set                             |
! Timestep_init_L    | if Timestep_init_L(i) is TRUE,this set of output|
!                    | timesteps only apply to initialization period   |
!----------------------------------------------------------------------
!
   integer, parameter :: MAXSTEP = 50000
   integer :: Timestep_tbl(MAXSTEP,MAXSET)
   integer :: Timestep_max(MAXSET),Timestep_id(MAXSET)
   integer :: Timestep_sets
   logical :: Timestep_init_L(MAXSET)
end module timestep
