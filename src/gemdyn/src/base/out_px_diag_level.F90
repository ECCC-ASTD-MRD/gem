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

!**s/p out_px_diag_level - interpolated pressure at given height
!

subroutine out_px_diag_level (F_f,F_height,minx,maxx,miny,maxy,nk, &
     f_px_m,F_zm)
  use glb_ld
  implicit none
#include <arch_specific.hf>
  !
  integer, intent(in) :: minx,maxx,miny,maxy,nk
  real, dimension(minx:maxx,miny:maxy), intent(out) :: F_f
  real, intent(in) :: F_height
  real, dimension(minx:maxx,miny:maxy,nk+1), intent(in) :: F_px_m,F_zm
  !
  !author
  !     Andre Plante mars 2021.
  !
  !object
  !	Interpolate the pressure to given height e.g. diag level height.
  !     F_height and F_zm,F_zt must have the same units.
  !
  
  !arguments
  !  Name               Description
  !---------------------------------------------------
  ! F_f                 Pressure a height F_height
  ! F_height            Height or geopotential (AGL) at which pression must be interpolated
  ! F_px_m              Pressure on mometum levels
  ! F_gzm               Height or geopotential on mometum levels

  integer :: i,j,k
  real :: ff, wp, wm, den
  !__________________________________________________________________
  !
  !
  !     F_zm(k)    ------------------
  ! 
  !                       ...
  !
  !     F_zm(nk-1) -----------------
  !                        
  !     F_zm(nk)   -----------------
  !                        
  !     F_zm(nk+1) ======================= surface
  !
  !____________________________________________________________________

  ! Note: level k just above diag level may not be to same for all i,j.
  ! This is why there is a loop on nested in the i j loops.
  do j=1,l_nj
     do i=1,l_ni
        ff = F_height + F_zm(i,j,nk+1)
        do k=nk,1,-1
           if(F_zm(i,j,k) > ff )then
              den = F_zm(i,j,k) - F_zm(i,j,k+1)
              wm = (ff - F_zm(i,j,k+1)) / den
              wp = (F_zm(i,j,k) - ff) / den
              F_f(i,j) = wm*F_px_m(i,j,k) + wp*f_px_m(i,j,k+1)
              exit
           endif
        enddo
     end do
  end do
  !     __________________________________________________________________
  !
  return
end subroutine out_px_diag_level
