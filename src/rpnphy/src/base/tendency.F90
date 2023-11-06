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

module tendency
   implicit none
   private
   public :: tendency5, apply_tendencies

#include "phymkptr.hf"

   interface apply_tendencies
      module procedure apply_tendencies_ptr
   end interface apply_tendencies

contains

   !/@*
   subroutine tendency5(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, &
        pvars, rcdt1, kount, ni, nk)
      use phy_options
      use phybusidx
      use phymem, only: phyvar
      use ens_perturb, only: ptp_L
      implicit none
!!!#include <arch_specific.hf>
      !@Author  L. Spacek (Oct 2011)
      !@Object Calculates tendencies in physics

      !@Arguments

      integer, intent(in) :: kount, ni, nk
      real, dimension(ni,nk), intent(in) :: uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0
      type(phyvar), pointer, contiguous :: pvars(:)
      real, intent(in)  :: rcdt1

      ! ni       horizontal running length
      ! nk       vertical dimension
      ! rcdt1    1/cdt1
      ! uplus0   initial value of pvars(uplus)
      ! vplus0   initial value of pvars(vplus)
      ! tplus0   initial value of pvars(tplus)
      ! huplus0  initial value of pvars(huplus)
      ! qcplus0  initial value of pvars(qcplus)
      ! pvars    list of all phy vars (meta + slab data)
      !*@/
#include <rmn/msg.h>

      integer :: i, k
      real, pointer, dimension(:,:), contiguous :: zhuphytd, zhuplus, zqcphytd, zqcplus, zqdifv, ztdifv, ztphytd, ztplus, zuphytd, zudifv, zuplus, zvphytd, zvdifv, zvplus, zwphytd, zwplus
      !-------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'tendency [BEGIN]')
      if (timings_L) call timing_start_omp(450, 'tendency', 46)

      MKPTR2D(zhuphytd, huphytd, pvars)
      MKPTR2D(zhuplus , huplus, pvars)
      MKPTR2D(zqcphytd, qcphytd, pvars)
      MKPTR2D(zqcplus , qcplus, pvars)
      MKPTR2D(zqdifv  , qdifv, pvars)
      MKPTR2D(ztdifv  , tdifv, pvars)
      MKPTR2D(ztphytd , tphytd, pvars)
      MKPTR2D(ztplus  , tplus, pvars)
      MKPTR2D(zuphytd , uphytd, pvars)
      MKPTR2D(zudifv  , udifv, pvars)
      MKPTR2D(zuplus  , uplus, pvars)
      MKPTR2D(zvphytd , vphytd, pvars)
      MKPTR2D(zvdifv  , vdifv, pvars)
      MKPTR2D(zvplus  , vplus, pvars)
      MKPTR2D(zwphytd , wphytd, pvars)
      MKPTR2D(zwplus  , wplus, pvars)

      do k=1,nk
         do i=1,ni
            zuphytd (i,k) = (zuplus (i,k) - uplus0 (i,k)) * rcdt1
            zvphytd (i,k) = (zvplus (i,k) - vplus0 (i,k)) * rcdt1
            ztphytd (i,k) = (ztplus (i,k) - tplus0 (i,k)) * rcdt1
            zhuphytd(i,k) = (zhuplus(i,k) - huplus0(i,k)) * rcdt1
         enddo
      enddo

      if (associated(zqcphytd)) then
         do k=1,nk
            do i=1,ni
               zqcphytd(i,k) = (zqcplus(i,k) - qcplus0(i,k)) * rcdt1
            enddo
         enddo
      endif

      if (diffuw) then
         do k=1,nk
            do i=1,ni
               zwphytd(i,k)  = (zwplus(i,k) - wplus0(i,k)) * rcdt1
            enddo
         enddo
      endif

      do i=1,ni
         zuphytd (i,nk) = zudifv(i,nk)
         zvphytd (i,nk) = zvdifv(i,nk)
         ztphytd (i,nk) = ztdifv(i,nk)
         zhuphytd(i,nk) = zqdifv(i,nk)
      end do

      if (ptp_L .and. kount >= 1) then
         do k=1,nk
            do i=1,ni
               zuplus (i,k) = uplus0 (i,k)
               zvplus (i,k) = vplus0 (i,k)
               ztplus (i,k) = tplus0 (i,k)
               !       zhuplus(i,k) = huplus0(i,k)
               !       zqcplus(i,k) = qcplus0(i,k)
            enddo
         enddo
      endif

      if (timings_L) call timing_stop_omp(450)
      call msg_toall(MSG_DEBUG, 'tendency [END]')
      !-------------------------------------------------------------
      return
   end subroutine tendency5

   
   !/@*
   subroutine apply_tendencies_ptr(zivar, ziten, ztdmask, ni, nk, F_nkscope, F_minval, F_maxval)
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      integer, intent(in) :: ni, nk
      real, dimension(ni, nk), intent(in) :: ziten
      real, dimension(ni), intent(in)  :: ztdmask
      real, dimension(ni, nk), intent(inout) :: zivar
      integer, intent(in), optional :: F_nkscope
      real, intent(in), optional :: F_minval, F_maxval
      !*@/
      integer :: k, nkscope
     !----------------------------------------------------------------
      nkscope = nk
      if (present(F_nkscope)) nkscope = F_nkscope
      if (present(F_minval)) then
         if (present(F_maxval)) then
            do k=1,nkscope
               zivar(:,k) = max(F_minval, max(F_minval, &
                    zivar(:,k) + ztdmask(:)*ziten(:,k)*delt))
            end do
         else
            do k=1,nkscope
               zivar(:,k) = max(F_minval, &
                    zivar(:,k) + ztdmask(:)*ziten(:,k)*delt)
            end do
         endif
      else if (present(F_maxval)) then
         do k=1,nkscope
            zivar(:,k) = min(F_maxval, &
                 zivar(:,k) + ztdmask(:)*ziten(:,k)*delt)
         enddo
      else
         do k=1,nkscope
            zivar(:,k) = zivar(:,k) + ztdmask(:)*ziten(:,k)*delt
         enddo
      endif
      !----------------------------------------------------------------
      return
   end subroutine apply_tendencies_ptr


end module tendency
