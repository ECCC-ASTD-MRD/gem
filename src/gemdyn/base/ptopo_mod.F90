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

module ptopo
   implicit none
   public
   save
!
!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH LOGICAL PARALLEL PROCESSOR TOPOLOGY       |
!                                                                      |
!     along Y                                                          |
!        .           .                    .                            |
!        .           .                    .                            |
!   +-----------+-----------+     +---------------+                    |
!   | (0,myrow) | (1,myrow) |.....| (mycol,myrow) |.....               |
!   +-----------+-----------+     +---------------+                    |
!        .           .                    .                            |
!        .           .                    .                            |
!   +-----------+-----------+     +---------------+                    |
!   |   (0,2)   |   (1,2)   |.....|   (mycol,2)   |.....               |
!   +-----------+-----------+     +---------------+                    |
!   |   (0,1)   |   (1,1)   |.....|   (mycol,1)   |.....               |
!   +-----------+-----------+     +---------------+                    |
!   |   (0,0)   |   (1,0)   |.....|   (mycol,0)   |..... along X       |
!   +-----------+-----------+     +---------------+                    |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Ptopo_myproc       | local processor number (zero based numbering)   |
!                    | 0,1,2 ... (Ptopo_npex * Ptopo_npey -1)          |
! Ptopo_myrow        | local row    number in processor topology       |
! Ptopo_mycol        | local column number in processor topology       |
! Ptopo_numproc      | total number of processors used                 |
! Ptopo_npex         | number of processors along X                    |
! Ptopo_npey         | number of processors along Y                    |
! Ptopo_npeOpenMP    | number of processors requested for OpenMp       |
! Ptopo_nthreads_dyn | number of threads for the dynamics              |
! Ptopo_nthreads_phy | number of threads for the physics               |
! Ptopo_gindx        | contains global indices that represents:        |
!                    | (1,*)-the minimum I indices on each local PE    |
!                    | (2,*)-the maximum I indices on each local PE    |
!                    | (3,*)-the minimum J indices on each local PE    |
!                    | (4,*)-the maximum J indices on each local PE    |
!                    | (5,*)-the minimum K indices on each local PE    |
!                    | (6,*)-the maximum K indices on each local PE    |
! Ptopo_ncolors      | Number of colors (or sub grids)                 |
! Ptopo_couleur      | 0 for Yin, 1 for Yan (Yang) domain              |
! Ptopo_tag          | tag number within the inter-communicator        |
! Ptopo_intracomm    | intra communicator number for Yin, Yan          |
! Ptopo_intercomm    | inter-communicator number between Yin and Yan   |
!----------------------------------------------------------------------
!
   integer Ptopo_myproc    , Ptopo_myrow       , Ptopo_mycol
   integer Ptopo_numproc   , Ptopo_npex        , Ptopo_npey
   integer Ptopo_npeOpenMP , Ptopo_nthreads_dyn, Ptopo_nthreads_phy
   integer Ptopo_nodes     , Ptopo_ncolors     , Ptopo_couleur
   integer Ptopo_intracomm , Ptopo_intercomm   , Ptopo_tag
   integer Ptopo_world_myproc, Ptopo_world_numproc
   logical Ptopo_last_domain_L

   integer, dimension(:,:), allocatable :: Ptopo_gindx

end module ptopo
