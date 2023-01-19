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

      subroutine itf_phy_diag ()
      use phy_itf, only: phy_get
      use HORgrid_options
      use outp
      use wb_itf_mod
      implicit none

   !@objective
   ! To compute diagnostic level values using the surface layer module.
   !@arguments
   !@author  Ron McTaggart-Cowan - Winter 2015
   !*@/

   ! Local variable declarations
      logical, save :: init_L= .false., dodiag_L= .true.
      integer :: istat,i,j
      real :: zu,zt
      real, dimension(:,:  ), pointer :: ptr2d
!
!----------------------------------------------------------------------
!     
! Retrieve physics information only if physics is active
      if (.not. init_L) then
         istat = wb_get('phy/zu', zu)
         istat = wb_get('phy/zt', zt)
         dodiag_L = (zu >= 0. .and. zt >= 0.)
         init_L = .true.
      end if
! Retrieve physics internal diagnostics
      if (dodiag_L) then
         ptr2d => tdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
         istat = phy_get(ptr2d,'tdiag',F_npath='V')
         ptr2d => udiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
         istat = phy_get(ptr2d,'udiag',F_npath='V')
         ptr2d => vdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
         istat = phy_get(ptr2d,'vdiag',F_npath='V')
         ptr2d => qdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
         istat = phy_get(ptr2d,'qdiag',F_npath='V')
      endif
!
!----------------------------------------------------------------------
!
      return
      end subroutine itf_phy_diag
