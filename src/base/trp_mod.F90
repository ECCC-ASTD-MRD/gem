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

module trp
   implicit none
   public
   save
!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH TRANSPOSE (set_transpose)                 |
!   (see subroutine set_transpose for documentation)                   |
!______________________________________________________________________|
!
   integer Trp_12dmin, Trp_12dmax, Trp_12dn, Trp_12dn0
   integer Trp_p12dmin,Trp_p12dmax,Trp_p12dn,Trp_p12dn0
   integer Trp_12smin, Trp_12smax, Trp_12sn, Trp_12sn0
   integer Trp_12emin, Trp_12emax, Trp_12en, Trp_12en0
   integer Trp_22min , Trp_22max , Trp_22n , Trp_22n0
   integer Trp_22emin, Trp_22emax, Trp_22en, Trp_22en0

end module trp
