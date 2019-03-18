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
   implicit none
   public
   save

!______________________________________________________________________
!                                                                      |
!  VARIABLES FOR GLOBAL(G) and LOCAL(L) reference                      |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! G_periodx          | .true. if periodic on X                         |
! G_periody          | .true. if periodic on Y                         |
! G_ni               | global # of grid points along X  (scalar grid)  |
! G_nj               | global # of grid points along Y  (scalar grid)  |
! G_nk               | total number of computational vertical levels   |
! G_niu              | global # of grid points along X  (U      grid)  |
! G_njv              | global # of grid points along Y  (V      grid)  |
! G_lnimax           | largest l_ni of all subdomain                   |
! G_lnjmax           | largest l_nj of all subdomain                   |
! G_halox            | number of points for the halo on X              |
! G_haloy            | number of points for the halo on y              |
! G_xg_8             | longitudes of the scalar grid in radians        |
! G_yg_8             | latitudes  of the scalar grid in radians        |
! l_north            | .true. if subdomain owns global north boundary  |
! l_south            | .true. if subdomain owns global south boundary  |
! l_east             | .true. if subdomain owns global east  boundary  |
! l_west             | .true. if subdomain owns global west  boundary  |
! l_mesg_proc        | NOT USED                                        |
! l_ni               | local # of grid points on X (scalar grid)       |
! l_nj               | local # of grid points on Y (scalar grid)       |
! l_nk               | total number of computational vertical levels   |
! l_niu              | local # of grid points on X (U      grid)       |
! l_njv              | local # of grid points on Y (V      grid)       |
! l_i0               | left   global index of local subdomain          |
! l_j0               | bottom global index of local subdomain          |
! l_minx             | minimum value for first  index of main 3D var.  |
! l_maxx             | maximum value for first  index of main 3D var.  |
! l_miny             | minimum value for second index of main 3D var.  |
! l_maxy             | maximum value for second index of main 3D var.  |
! l_dimmsg           | NOT USED                                        |
! l_dim2d            | NOT USED                                        |
! l_dim3d            | NOT USED                                        |
! pil_n              | # of points on global north boundary for pilot  |
! pil_s              | # of points on global south boundary for pilot  |
! pil_w              | # of points on global west  boundary for pilot  |
! pil_e              | # of points on global east  boundary for pilot  |
! north              | =1 if touching north boundary                   |
! south              | =1 if touching south boundary                   |
! east               | =1 if touching east  boundary                   |
! west               | =1 if touching west  boundary                   |
!----------------------------------------------------------------------

   logical G_periodx, G_periody
   integer G_ni, G_nj, G_nk, G_niu, G_njv, G_lnimax, G_lnjmax
!      integer G_halox, G_haloy
   real*8, dimension(:), allocatable :: G_xg_8, G_yg_8

   logical l_north, l_south, l_east, l_west, l_mesg_proc
   integer l_ni, l_nj, l_nk, l_niu, l_njv, l_i0, l_j0
   integer l_minx, l_maxx, l_miny, l_maxy
   integer l_dimmsg, l_dim2d, l_dim3d
   integer pil_n,pil_s,pil_w,pil_e,north,south,east,west

end module glb_ld
