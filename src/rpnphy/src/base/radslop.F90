
module radslop
   implicit none
   private
   public :: radslop3

contains

   !/@*
   subroutine radslop3(pvars, hz, julien, ni, trnch)
      use debug_mod, only: init2nan
      use phy_options
      use phybusidx
      use phymem, only: phyvar
      use series_mod, only: series_xst
      implicit none
!!!#include <arch_specific.hf>

      !@Object add the effects of the sloping terrain to the radiation computation.
      !@Arguments
      !          - Input/Output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - Input -
      ! ni       horizontal dimension
      ! hz       Greenwich hour (0 to 24)
      ! julien   Julien days
      ! trnch    number of the slice

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni, trnch
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

      MKPTR1D(zc1slop, c1slop, pvars)
      MKPTR1D(zc2slop, c2slop, pvars)
      MKPTR1D(zc3slop, c3slop, pvars)
      MKPTR1D(zc4slop, c4slop, pvars)
      MKPTR1D(zc5slop, c5slop, pvars)
      MKPTR1D(zdlat, dlat, pvars)
      MKPTR1D(zdlon, dlon, pvars)
      MKPTR1D(zfluslop, fluslop, pvars)
      MKPTR1D(zfsd0, fsd0, pvars)
      MKPTR1D(zfsf0, fsf0, pvars)
      MKPTR1D(zvv1, vv1, pvars)
      MKPTR1D(zap, ap, pvars)

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
