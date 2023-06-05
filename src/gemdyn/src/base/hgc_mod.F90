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

module hgc
   implicit none
   public
   save
!______________________________________________________________________
!                                                                      |
!  VARIABLES FOR HORIZONTAL DESCRIPTORS OF OUTPUT GRIDS (set_hgc)      |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Hgc_gxtyp_s        | grid type for positional records (>> and ^^)    |
!----------------------------------------------------------------------|
! The following four grid descriptors are used when                    |
! generating >> and ^^ records. They convey grid rotation information. |
! They are the same for phi, U and V grids.                            |
!----------------------------------------------------------------------|
! Hgc_ig1ro          | first grid descriptor grid                      |
! Hgc_ig2ro          | second grid descriptor grid                     |
! Hgc_ig3ro          | third grid descriptor grid                      |
! Hgc_ig4ro          | fourth grid descriptor grid                     |
!----------------------------------------------------------------------
!
!
   integer :: Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro
   integer :: dstf_gid, dstu_gid, dstv_gid, dstp_gid
   integer :: dstf_gis, dstu_gis, dstv_gis
   character(len=4) :: Hgc_gxtyp_s

end module hgc
