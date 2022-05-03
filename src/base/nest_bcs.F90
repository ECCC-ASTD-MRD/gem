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

!**s/r nest_bcs

      subroutine nest_bcs ( F_dt_8, F_rhsu, F_rhsv, &
                            Minx, Maxx, Miny, Maxy, Nk )
      use lam_options
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: &
                                               F_rhsu, F_rhsv
      real(kind=REAL64), intent(IN) :: F_dt_8
!
!----------------------------------------------------------------------
!
      if (.not. Lam_ctebcs_L) call nest_intt

      call nest_bcs_t0 ( F_dt_8, F_rhsu, F_rhsv, Minx,Maxx,Miny,Maxy,Nk )
!
!----------------------------------------------------------------------
!
      return
      end
