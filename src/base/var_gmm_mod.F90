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

module var_gmm
   use rmn_gmm
   implicit none
   public
   save

   type(gmm_metadata) :: meta2d      ! (l_minx:l_maxx, l_miny:l_maxy          )
   type(gmm_metadata) :: meta3d_nk   ! (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk  )
   type(gmm_metadata) :: meta3d_nk1  ! (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk+1)
   type(gmm_metadata) :: meta3d_0nk  ! (l_minx:l_maxx, l_miny:l_maxy, 0:l_nk  )
   type(gmm_metadata) :: meta3d_0nk1 ! (l_minx:l_maxx, l_miny:l_maxy, 0:l_nk+1)

end module var_gmm
