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

module hzd_ctrl
  use hzd_exp
  use gem_options
  use glb_ld
  implicit none
#include <arch_specific.hf>
  private
  public :: hzd_ctrl4

  interface hzd_ctrl4
     module procedure hzd_ctrl_scalar
     module procedure hzd_ctrl_vector
  end interface

contains

      subroutine hzd_ctrl_scalar ( F_f2hzd, F_type_S, Minx,Maxx,Miny,Maxy,Nk )
      use glb_ld
      implicit none

      character(len=*), intent(in) :: F_type_S
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent (inout) :: F_f2hzd
!
!-------------------------------------------------------------------
!
         call hzd_exp_visco (F_f2hzd, F_type_S, l_minx,l_maxx,l_miny,l_maxy, Nk)
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_ctrl_scalar

      subroutine hzd_ctrl_vector ( F_u, F_v, Minx,Maxx,Miny,Maxy,Nk )
      use glb_ld
      implicit none

      integer Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: F_u, F_v
!
!-------------------------------------------------------------------
!
         call hzd_exp_visco(F_u, 'U', l_minx,l_maxx,l_miny,l_maxy, Nk)
         call hzd_exp_visco(F_v, 'V', l_minx,l_maxx,l_miny,l_maxy, Nk)
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_ctrl_vector

end module hzd_ctrl
