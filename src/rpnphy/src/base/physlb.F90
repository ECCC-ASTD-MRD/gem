!---------------------------------- LICENCE BEGIN -----------------------------
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
!---------------------------------- LICENCE END -------------------------------

module physlb
   use phy_status, only: phy_error_L
   use phy_options
   use phyexe, only: phyexe1
   private
   public :: physlb1

!!!#include <arch_specific.hf>

contains

   !/@*
   subroutine physlb1(busdyn3D ,busper3D, busvol3D, &
        dsiz, fsiz, vsiz, kount, ni, nj, nk, pslic)
      implicit none
      !@Author L. Spacek (May 2010)
      !@Object The main physics subroutine
      !@Arguments
      
      integer, intent(in) :: dsiz, fsiz, vsiz, kount, ni, nj, nk
      integer, intent(inout) :: pslic
      real, dimension(:,:), pointer, contiguous :: busdyn3D, busper3D, busvol3D
      
      !          - Input/Output -
      ! busdyn3D  - dynamics input field
      ! busper3D  - historic variables for the physics
      !
      !          - Output -
      ! busvol3D  - physics tendencies and other output fields from the physics
      !          - Input -
      ! dsiz     dimension of d
      ! fsiz     dimension of f
      ! vsiz     dimension of v
      ! kount    timestep number
      ! ni       horizontal running length
      ! nj       number of slices
      ! nk       vertical dimension
      !*@/

      integer :: jdo
      real, dimension(:), pointer, contiguous :: dbus, fbus, vbus
      !     ---------------------------------------------------------------

100   continue

!$omp critical
      pslic = pslic+1
      jdo   = pslic
!$omp end critical

      if (jdo > nj) return

      dbus(1:dsiz) => busdyn3D(:,jdo)
      fbus(1:fsiz) => busper3D(:,jdo)
      vbus(1:vsiz) => busvol3D(:,jdo)
      
      if (test_phy) then
         call physeterror('physlb1', 'testphy_phyexe needs to be updated')
!!$       call testphy_phyexe(dbus, fbus, vbus, &
!!$         jdo, kount, ni, nk)
      else
         call phyexe1(dbus, fbus, vbus, &
              jdo, kount, ni, nk)
      endif
      if (phy_error_L) return

      goto 100
      !     ---------------------------------------------------------------
   end subroutine physlb1

end module physlb
