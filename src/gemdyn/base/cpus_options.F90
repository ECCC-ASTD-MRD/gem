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

module cpus_options
   implicit none
   public
   save

   !# number of processors along X
   integer :: Cpus_npex = 1
   namelist /cpus/ Cpus_npex

   !# number of processors along Y
   integer :: Cpus_npey = 1
   namelist /cpus/ Cpus_npey

   !# number of threads for the dynamics
   integer :: Cpus_nthreads_dyn = 0
   namelist /cpus/ Cpus_nthreads_dyn

   !# number of threads for the physics
   integer :: Cpus_nthreads_phy = 0
   namelist /cpus/ Cpus_nthreads_phy

end module cpus_options

