!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_marine_halo.ftn90
! Creation       : Diane Pendlebury
! Description    : Parameterization of effect of halogen chemistry on ozone over oceans
! Extra info     :  
!
! Arguments:
!            IN
!               LU14          -> "sea" vegetation type 
!               Pressure      -> 3d pressure field
!
!           IN/OUT
!               Ozone         -> 3d ozone field
!               
!
!================================================================================================
!
!!if_on
subroutine mach_marine_halo(ozone, pressure, lfu, imod, ni, nk)
   use chm_ptopo_grid_mod,   only: chm_ni
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_utils_mod,        only: chm_timestep
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: ni, nk
   integer(kind=4), intent   (in) :: imod(ni)
   real(kind=4),    intent(inout) :: ozone(ni, nk)
   real(kind=4),    intent   (in) :: lfu(chm_ni, lucprm)
   real(kind=4),    intent   (in) :: pressure(ni, nk)
!!if_off
!
!  Local variables
!
   integer(kind=4)  :: i, k, ii
   real   (kind=8)  :: k_ozone

!------------------ parameterized marine halogen chemistry---------------------!!!!!!!!
!  Note that anything in the stratosphere should revert right back to linoz

   do i = 1, ni
      ii = imod(i)
      if  (lfu(ii,14) >  0.01)  then
         do k = 1, nk
            !loss rate of ozone (s^(-1))
            k_ozone = 1d-40*exp(7.74d-4*dble(pressure(i,k))) + 4.0582d-9*exp(5.7451d-5*dble(pressure(i,k)))
            !loss of ozone for this timestep due to marine halogen chemistry
            ozone(i,k) = ozone(i,k) * (1.0 - chm_timestep * real(k_ozone) * lfu(ii,14))
         enddo
      endif
   enddo

end subroutine mach_marine_halo


