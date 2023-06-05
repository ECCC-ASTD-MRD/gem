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
      module procedure apply_tendencies_bus
      module procedure apply_tendencies_bus2
      module procedure apply_tendencies_ptr
   end interface apply_tendencies

contains

   !/@*
   subroutine tendency5(uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0, &
        vbus, dbus, rcdt1, kount, ni, nk)
      use phy_options
      use phybus
      use ens_perturb, only: ptp_L
      implicit none
!!!#include <arch_specific.hf>
      !@Author  L. Spacek (Oct 2011)
      !@Object Calculates tendencies in physics

      !@Arguments

      integer, intent(in) :: kount, ni, nk
      real, dimension(ni,nk), intent(in) :: uplus0, vplus0, wplus0, tplus0, huplus0, qcplus0
      real, dimension(:), pointer, contiguous :: dbus, vbus
      real, intent(in)  :: rcdt1

      ! ni       horizontal running length
      ! nk       vertical dimension
      ! rcdt1    1/cdt1
      ! uplus0   initial value of dbus(uplus)
      ! vplus0   initial value of dbus(vplus)
      ! tplus0   initial value of dbus(tplus)
      ! huplus0  initial value of dbus(huplus)
      ! qcplus0  initial value of dbus(qcplus)
      ! dbus     dynamics input field
      ! vbus     physics tendencies and other output fields from the physics
      !*@/
#include <rmn/msg.h>

      integer :: i, k
      real, pointer, dimension(:,:), contiguous :: zhuphytd, zhuplus, zqcphytd, zqcplus, zqdifv, ztdifv, ztphytd, ztplus, zuphytd, zudifv, zuplus, zvphytd, zvdifv, zvplus, zwphytd, zwplus
      !-------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'tendency [BEGIN]')
      if (timings_L) call timing_start_omp(450, 'tendency', 46)

      MKPTR2D(zhuphytd, huphytd, vbus)
      MKPTR2D(zhuplus , huplus, dbus)
      MKPTR2D(zqcphytd, qcphytd, vbus)
      MKPTR2D(zqcplus , qcplus, dbus)
      MKPTR2D(zqdifv  , qdifv, vbus)
      MKPTR2D(ztdifv  , tdifv, vbus)
      MKPTR2D(ztphytd , tphytd, vbus)
      MKPTR2D(ztplus  , tplus, dbus)
      MKPTR2D(zuphytd , uphytd, vbus)
      MKPTR2D(zudifv  , udifv, vbus)
      MKPTR2D(zuplus  , uplus, dbus)
      MKPTR2D(zvphytd , vphytd, vbus)
      MKPTR2D(zvdifv  , vdifv, vbus)
      MKPTR2D(zvplus  , vplus, dbus)
      MKPTR2D(zwphytd , wphytd, vbus)
      MKPTR2D(zwplus  , wplus, dbus)

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
   subroutine apply_tendencies_bus(dbus,dsiz,vbus,vsiz,fbus,fsiz,ivar,iten,ni,nk,nkscope)
      use phy_options
      use phybus
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      !          - input -
      ! ivar     variable  index
      ! iten     tendency  index
      ! ni       horizonal index
      ! nk       vertical  index
      ! nkscope  vertical  operator scope
      ! vbus     volatile bus
      ! fbus     permanent bus
      !
      !          - input/output -
      ! dbus     dynamics bus

      integer, intent(in)    :: dsiz,vsiz,fsiz,ivar,iten,ni,nk,nkscope
      real, target, intent(in)    :: vbus(vsiz), fbus(fsiz)
      real, target, intent(inout) :: dbus(dsiz)

      !@Author L. Spacek (Oct 2011)
      !*@/

      integer :: k

      real, pointer, dimension(:), contiguous   :: ztdmask
      real, pointer, dimension(:,:), contiguous :: ziten, zivar
      !----------------------------------------------------------------

      MKPTR1D(ztdmask, tdmask, fbus)
      MKPTR2D(ziten, iten, vbus)
      MKPTR2D(zivar, ivar, dbus)

      if (.not.(associated(ztdmask) .and. associated(ziten) .and. associated(zivar))) then
         call physeterror('apply_tendencies_bus', 'Problem getting pointers')
         return
      endif

      do k=1,nkscope
         zivar(:,k) = zivar(:,k) + ztdmask(:)*ziten(:,k)*delt
     enddo

      !----------------------------------------------------------------
      return
   end subroutine apply_tendencies_bus


   !/@*
   subroutine apply_tendencies_bus2(dbus, vbus, fbus, ivar, iten, ni, nk, nkscope)
      use phy_options
      use phybus
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      !          - input -
      ! ivar     variable  index
      ! iten     tendency  index
      ! ni       horizonal index
      ! nk       vertical  index
      ! nkscope  vertical  operator scope
      ! vbus     volatile bus
      ! fbus     permanent bus
      !
      !          - input/output -
      ! dbus     dynamics bus

      real, target, intent(inout) :: dbus(:)
      real, target, intent(in)    :: vbus(:),fbus(:)
      integer, intent(in) :: ivar,iten,ni,nk
      integer, intent(in), optional :: nkscope

      !@Author L. Spacek (Oct 2011)
      !*@/

      integer :: k, nkscope1

      real, pointer, dimension(:), contiguous   :: ztdmask
      real, pointer, dimension(:,:), contiguous :: ziten, zivar
      !----------------------------------------------------------------

      MKPTR1D(ztdmask, tdmask, fbus)
      MKPTR2D(ziten, iten, vbus)
      MKPTR2D(zivar, ivar, dbus)

      nkscope1 = nk
      if (present(nkscope)) nkscope1 = nkscope

      if (.not.(associated(ztdmask) .and. associated(ziten) .and. associated(zivar))) then
         call physeterror('apply_tendencies_bus2', 'Problem getting pointers')
         return
      endif

      do k=1,nkscope1
         zivar(:,k) = zivar(:,k) + ztdmask(:)*ziten(:,k)*delt
     enddo

      !----------------------------------------------------------------
      return
   end subroutine apply_tendencies_bus2

   
   !/@*
   subroutine apply_tendencies_ptr(zivar, ziten, ztdmask, ni, nk, nkscope, F_minval, F_maxval)
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      integer, intent(in) :: ni, nk, nkscope
      real, dimension(ni)  :: ztdmask
      real, dimension(ni, nk) :: ziten, zivar
      real, intent(in), optional :: F_minval, F_maxval
      !*@/
      integer :: k
     !----------------------------------------------------------------
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
