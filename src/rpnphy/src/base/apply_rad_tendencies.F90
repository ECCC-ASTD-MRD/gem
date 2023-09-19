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

module apply_rad_tendencies
   implicit none
   private
   public :: apply_rad_tendencies1

contains

   !/@*
   subroutine apply_rad_tendencies1(pvars, ni, nk, nkm1)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use phybudget, only: pb_compute, pb_conserve, pb_residual
      use phy_status, only: PHY_OK
      use phy_options
      use phybusidx
      use phymem, only: phyvar
      use tendency, only: apply_tendencies
      use ens_perturb, only: ens_nc2d, ens_spp_get
      implicit none
!!!#include <arch_specific.hf>
      !@Object Apply radiative tendencies
      !@Arguments
      !          - input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - Input -
      ! ni       horizontal running length
      ! nk       vertical dimension
      ! nkm1     vertical scope of the operator

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni, nk, nkm1

      !@Author Ron McTaggart-Cowan, Summer 2017
      !*@/
#include <rmn/msg.h>
#include <rmnlib_basics.hf>
#include "phymkptr.hf"

      ! Local variables
      integer :: k
      real, dimension(ni) :: tradmult
      real, dimension(ni,nkm1) :: qrad, mtrad
      real, dimension(:), pointer, contiguous :: zconerad, zconqrad, zps, ztdmask, znetrad
      real, dimension(:,:), pointer, contiguous :: ztplus, zqplus, zqcplus, zgztherm, zsigt, ztrad, zmrk2
      real(REAL64), dimension(ni) :: l_en0, l_pw0
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'apply_rad_tendencies [BEGIN]')
      if (timings_L) call timing_start_omp(422, 'apply_rad_td', 46)
      
      ! Bus pointer association
      MKPTR1D(zconerad, conerad, pvars)
      MKPTR1D(zconqrad, conqrad, pvars)
      MKPTR1D(znetrad, netrad, pvars)
      MKPTR1D(zps, pmoins, pvars)
      MKPTR1D(ztdmask, tdmask, pvars)

      MKPTR2Dm1(zgztherm, gztherm, pvars)
      MKPTR2Dm1(zqcplus, qcplus, pvars)
      MKPTR2Dm1(zqplus, huplus, pvars)
      MKPTR2Dm1(zsigt, sigt, pvars)
      MKPTR2Dm1(ztplus, tplus, pvars)
      MKPTR2Dm1(ztrad, trad, pvars)
      
      MKPTR2D(zmrk2, mrk2, pvars)

      ! Early exit if radiation is not used
      if (radia == 'NIL') then
         if (associated(zconerad)) zconerad = 0.
         if (associated(zconqrad)) zconqrad = 0.
         if (timings_L) call timing_stop_omp(422)
         return
      endif

      call init2nan(tradmult)
      call init2nan(qrad, mtrad)
      call init2nan(l_en0, l_pw0)

      ! Pre-radiation state for budget
      if (pb_compute(zconerad, zconqrad, l_en0, l_pw0, &
           pvars, nkm1) /= PHY_OK) then
         call physeterror('apply_rad_tendencies', &
              'Problem computing preliminary budget for '//trim(radia))
         return
      endif

      ! Adjust radiation tendencies to impose conservation
      qrad = 0.
      if (pb_conserve(rad_conserve, ztrad, qrad, pvars, &
           F_rad=znetrad) /= PHY_OK) then
         call physeterror('apply_rad_tendencies', &
              'Cannot correct conservation for '//trim(radia))
         return
      endif

      ! Apply tendency perturbation on request
      tradmult = ens_spp_get('trad_mult', zmrk2, default=1.)
      do k=1,nkm1
         mtrad(:,k) = tradmult(:) * ztrad(:,k)
      enddo
      
      ! Apply radiative tendencies
      call apply_tendencies(ztplus, mtrad, ztdmask, ni, nk, nkm1)

      ! Post-radiation budget analysis
      if (pb_residual(zconerad, zconqrad, l_en0, l_pw0, pvars, &
           delt, nkm1, F_rad=znetrad) /= PHY_OK) then
         call physeterror('apply_rad_tendencies', &
              'Problem computing final budget for '//trim(radia))
         return
      endif

      if (timings_L) call timing_stop_omp(422)
      call msg_toall(MSG_DEBUG, 'apply_rad_tendencies [END]')
      !----------------------------------------------------------------
      return
   end subroutine apply_rad_tendencies1

end module apply_rad_tendencies
