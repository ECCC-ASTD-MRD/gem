!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module radslop
   implicit none
   private
   public :: radslop3

contains

   !/@*
   subroutine radslop3(fbus, vbus, ni, hz, julien, trnch)
      use debug_mod, only: init2nan
      use phy_options
      use phybus
      use series_mod, only: series_xst
      implicit none
!!!#include <arch_specific.hf>

      !@Object add the effects of the sloping terrain to the radiation computation.
      !@Arguments
      !          - Input/Output -
      ! f        field of permanent physics variables
      !          - Input
      ! v        field of volatile physics variables
      ! ni       horizontal dimension
      ! hz       Greenwich hour (0 to 24)
      ! julien   Julien days
      ! trnch    number of the slice

      integer, intent(in) :: ni, trnch
      real, pointer, contiguous :: fbus(:), vbus(:)
      real, intent(in) :: julien
      real, intent(in) :: hz

      !@Authors J. Mailhot, A. Erfani, J. Toviessi
      !*@/
#include "phymkptr.hf"

      real, pointer, dimension(:), contiguous :: zc1slop, zc2slop, zc3slop, zc4slop, zc5slop, zdlat, zdlon, zfluslop, zfsd0, zfsf0, zvv1, zap
      real, dimension(ni) :: bcos, bsin, stan, ssin, scos
      real :: dire, difu, albe
      integer :: i
      !----------------------------------------------------------------

      MKPTR1D(zc1slop, c1slop, fbus)
      MKPTR1D(zc2slop, c2slop, fbus)
      MKPTR1D(zc3slop, c3slop, fbus)
      MKPTR1D(zc4slop, c4slop, fbus)
      MKPTR1D(zc5slop, c5slop, fbus)
      MKPTR1D(zdlat, dlat, fbus)
      MKPTR1D(zdlon, dlon, fbus)
      MKPTR1D(zfluslop, fluslop, fbus)
      MKPTR1D(zfsd0, fsd0, fbus)
      MKPTR1D(zfsf0, fsf0, fbus)
      MKPTR1D(zvv1, vv1, fbus)
      MKPTR1D(zap, ap, vbus)

      if (.not.radslope) then
         zfluslop(1:ni) = 0.0
         return
      endif

      call init2nan(bcos, bsin, stan, ssin, scos)

      call suncos2(scos, ssin, stan, bsin, bcos, ni, zdlat, zdlon, hz, &
           julien, radslope)

      do i = 1, ni

         dire = zfsd0(i) * ( &
              zc1slop(i) + &
              stan(i) * (zc2slop(i)*bcos(i) + zc3slop(i)*bsin(i)) &
              ) * zvv1(i)
         dire = max(0.0, dire)

         difu = (zfsf0(i) * zc4slop(i)) * zvv1(i)

         albe = ((zfsf0(i) + zfsd0(i)) * zc5slop(i) * zap(i)) * zvv1(i)

         zfluslop(i) = dire + difu + albe

      enddo

      call series_xst(zfluslop, 'fusl', trnch)
      !----------------------------------------------------------------
      return
   end subroutine radslop3

end module radslop
