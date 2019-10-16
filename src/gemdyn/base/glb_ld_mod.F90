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

module glb_ld
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   ! Variables for global(G) and local(L) reference

   !> .true. if periodic on X
   logical :: G_periodx
   !> .true. if periodic on Y
   logical :: G_periody
   !> Global # of grid points along X  (scalar grid)
   integer :: G_ni
   !> Global # of grid points along Y  (scalar grid)
   integer :: G_nj
   !> Total number of computational vertical levels
   integer :: G_nk
   !> Global # of grid points along X (U grid)
   integer :: G_niu
   !> Global # of grid points along Y (V grid)
   integer :: G_njv
   !> Largest l_ni of all subdomain
   integer :: G_lnimax
   !> Largest l_nj of all subdomain
   integer :: G_lnjmax
   !> Longitudes of the scalar grid in radians
   real(kind=REAL64), dimension(:), allocatable :: G_xg_8
   !> Latitudes  of the scalar grid in radians
   real(kind=REAL64), dimension(:), allocatable :: G_yg_8

   !> .true. if subdomain owns global north boundary
   logical :: l_north
   !> .true. if subdomain owns global south boundary
   logical :: l_south
   !> .true. if subdomain owns global east  boundary
   logical :: l_east
   !> .true. if subdomain owns global west  boundary
   logical :: l_west
   !> NOT USED
   logical :: l_mesg_proc

   !> Local # of grid points on X (scalar grid)
   integer :: l_ni
   !> Local # of grid points on Y (scalar grid)
   integer :: l_nj
   !> Total number of computational vertical levels
   integer :: l_nk
   !> Local # of grid points on X (U grid)
   integer :: l_niu
   !> Local # of grid points on Y (V grid)
   integer :: l_njv
   !> left global index of local subdomain
   integer :: l_i0
   !> Bottom global index of local subdomain
   integer :: l_j0

   !> Mminimum value for first  index of main 3D var
   integer :: l_minx
   !> Maximum value for first  index of main 3D var
   integer :: l_maxx
   !> Minimum value for second index of main 3D var
   integer :: l_miny
   !> Maximum value for second index of main 3D var
   integer :: l_maxy

   !> NOT USED
   integer :: l_dimmsg
   !> NOT USED
   integer :: l_dim2d
   !> NOT USED
   integer :: l_dim3d

   !> # of points on global north boundary for pilot
   integer :: pil_n
   !> # of points on global south boundary for pilot
   integer :: pil_s
   !> # of points on global west  boundary for pilot
   integer :: pil_w
   !> # of points on global east  boundary for pilot
   integer :: pil_e

   !> 1 if touching north boundary
   integer :: north
   !> 1 if touching south boundary
   integer :: south
   !> 1 if touching east  boundary
   integer :: east
   !> 1 if touching west  boundary
   integer :: west

end module glb_ld
