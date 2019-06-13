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

module diagno_cw_rad
   implicit none
   private
   public :: diagno_cw_rad1

contains

   !/@*
   subroutine diagno_cw_rad1(f, fsiz, d,dsiz, v, vsiz, &
        liqwcin, icewcin, liqwp, icewp, cloud, &
        trnch, ni, nk)
      use series_mod, only: series_xst
      use phybus
      implicit none
!!!#include <arch_specific.hf>

      integer, intent(in) :: fsiz, dsiz, vsiz, ni, nk, trnch
      real, intent(inout), target ::  f(fsiz), d(dsiz), v(vsiz)
      real, intent(inout) :: liqwcin(ni,nk), icewcin(ni,nk)
      real, intent(inout) :: liqwp(ni,nk-1), icewp(ni,nk-1) !#TODO: check this
      real, intent(inout) :: cloud(ni,nk)

      !@Author L. Spacek (Apr 2005)
      !@Object Calculate diagnostic for the radiation package
      !@Arguments
      !     - input -
      !     dsiz     Dimension of d
      !     fsiz     Dimension of f
      !     vsiz     Dimension of v
      !     liqwcin  in-cloud liquid water content
      !     icewcin  in-cloud ice    water content
      !     liqwp    in-cloud liquid water path
      !     icewp    in-cloud ice    water path
      !     cloud    cloudiness passed to radiation
      !     kount    index of timestep
      !     trnch    number of the slice
      !     ni       horizontal Dimension
      !     nk       number of layers
      !     - output -
      !     tlwp     total integrated liquid water path
      !     tiwp     total integrated ice    water path
      !     tlwpin   total integrated in-cloud liquid water path
      !     tiwpin   total integrated in-cloud ice    water path
      !     lwcrad   liquid water content passed to radiation
      !     iwcrad   ice    water content passed to radiation
      !     cldrad  cloudiness passed to radiation
      !*@/
#include "phymkptr.hf"

      integer :: i, k

      real, pointer, dimension(:)   :: ztlwp, ztiwp, ztlwpin, ztiwpin
      real, pointer, dimension(:,:) :: zlwcrad, ziwcrad, zcldrad
      !----------------------------------------------------------------

      MKPTR1D(ztlwp, tlwp, f)
      MKPTR1D(ztiwp, tiwp, f)
      MKPTR1D(ztlwpin, tlwpin, f)
      MKPTR1D(ztiwpin, tiwpin, f)
      MKPTR2D(zlwcrad, lwcrad, v)
      MKPTR2D(ziwcrad, iwcrad, v)
      MKPTR2D(zcldrad, cldrad, v)

      do i=1,ni
         ztlwp(i)      = 0.0
         ztiwp(i)      = 0.0
         ztlwpin(i)    = 0.0
         ztiwpin(i)    = 0.0
         zlwcrad(i,nk) = 0.0
         ziwcrad(i,nk) = 0.0
         zcldrad(i,nk) = 0.0
      enddo

      do k=1,nk-1
         do i=1,ni
            ztlwp(i)   = ztlwp(i)   + liqwp(i,k)*cloud(i,k)
            ztiwp(i)   = ztiwp(i)   + icewp(i,k)*cloud(i,k)
            ztlwpin(i) = ztlwpin(i) + liqwp(i,k)
            ztiwpin(i) = ztiwpin(i) + icewp(i,k)
         enddo
      enddo

      !# conversion d'unites : tlwp et tiwp en kg/m2
      do i=1,ni
         ztlwp(i)   = ztlwp(i) * 0.001
         ztiwp(i)   = ztiwp(i) * 0.001
         ztlwpin(i) = ztlwpin(i) * 0.001
         ztiwpin(i) = ztiwpin(i) * 0.001
      enddo

      do k=1,nk-1
         do i=1,ni
            zlwcrad(i,k) = liqwcin(i,k)*cloud(i,k)
            ziwcrad(i,k) = icewcin(i,k)*cloud(i,k)
            zcldrad(i,k) = cloud(i,k)
         enddo
      enddo

      !# extraction pour diagnostics
      call series_xst(ztlwp, 'icr', trnch)
      call series_xst(ztiwp, 'iir', trnch)
      call series_xst(ztlwpin, 'w1', trnch)
      call series_xst(ztiwpin, 'w2', trnch)
      call series_xst(ziwcrad, 'iwcr', trnch)
      call series_xst(zlwcrad, 'lwcr', trnch)
      call series_xst(zcldrad, 'cldr', trnch)

      !----------------------------------------------------------------
      return
   end subroutine diagno_cw_rad1

end module diagno_cw_rad
