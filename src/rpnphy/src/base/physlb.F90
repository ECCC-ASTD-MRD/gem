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
   use phymem, only: phyvar, phymem_get_slabvars
   use phy_status, only: phy_error_L
   use phy_options
   use phyexe, only: phyexe1
   use testphy_phyexe, only: testphy_phyexe1
   private
   public :: physlb1

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

contains

   !/@*
   subroutine physlb1(kount, ni, nj, nk, pslic)
      implicit none
      !@Author L. Spacek (May 2010)
      !@Object The main physics subroutine
      !@Arguments
      
      integer, intent(in) :: kount, ni, nj, nk
      integer, intent(inout) :: pslic
      
      !          - Input -
      ! kount    timestep number
      ! ni       horizontal running length
      ! nj       number of slices
      ! nk       vertical dimension
      !*@/

      integer :: jdo, istat
      type(phyvar), pointer, contiguous :: pvars(:)
      !     ---------------------------------------------------------------

100   continue

!$omp critical
      pslic = pslic+1
      jdo   = pslic
!$omp end critical

      if (jdo > nj) return

      nullify(pvars)
      istat = phymem_get_slabvars(pvars, jdo)
      if (.not.(RMN_IS_OK(istat) .and. associated(pvars))) then
         call physeterror('physlb1', 'Problem getting slab vars pointers')
         return
      endif
           
      if (test_phy) then
         call testphy_phyexe1(pvars, kount, ni, nk, jdo)
      else
         call phyexe1(pvars, kount, ni, nk, jdo)
      endif
      deallocate(pvars)
      if (phy_error_L) return

      goto 100
      !     ---------------------------------------------------------------
   end subroutine physlb1

end module physlb
