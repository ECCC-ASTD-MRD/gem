
module tendency
   implicit none
   private
   public :: tendency5, apply_tendencies, apply_tendencies_mm
   public :: apply_tendencies1, apply_tendencies2, apply_tendencies3, apply_tendencies4
   
#include "phymkptr.hf"

   interface apply_tendencies
      module procedure apply_tendencies1
      module procedure apply_tendencies2
      module procedure apply_tendencies3
      module procedure apply_tendencies4
   end interface apply_tendencies

contains

   !/@*
   subroutine tendency5(pvars, dt, kount, ni, nk)
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
      type(phyvar), pointer, contiguous :: pvars(:)
      real, intent(in)  :: dt

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
      real :: rcdt1
      real, pointer, dimension(:,:), contiguous :: zhuphytd, zhuplus, &
           zqcphytd, zqcplus, zqdifv, ztdifv, ztphytd, ztplus, zuphytd, &
           zudifv, zuplus, zvphytd, zvdifv, zvplus, zwphytd, zwplus, &
           zuplus0, zvplus0, zwplus0, ztplus0, zhuplus0, zqcplus0
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

      MKPTR2D(zqcplus0, qcplus0, pvars)
      MKPTR2D(zhuplus0, huplus0, pvars)
      MKPTR2D(zuplus0,  uplus0, pvars)
      MKPTR2D(zvplus0,  vplus0, pvars)
      MKPTR2D(ztplus0,  tplus0, pvars)
      MKPTR2D(zwplus0,  wplus0, pvars)

      rcdt1 = 1. / dt
      
      if (associated(zuphytd) .and. associated(zuplus0) .and. (ISREQSTEP('UPHY') .or. ISREQOUT('UPHM'))) then
         zuphytd = (zuplus - zuplus0) * rcdt1
         zuphytd(:,nk) = zudifv(:,nk)
      endif
      if (associated(zvphytd) .and. associated(zvplus0) .and. (ISREQSTEP('VPHY') .or. ISREQOUT('VPHM'))) then
         zvphytd = (zvplus - zvplus0) * rcdt1
         zvphytd(:,nk) = zvdifv(:,nk)
      endif
      if (associated(ztphytd) .and. associated(ztplus0) .and. (ISREQSTEP('TPHY') .or. ISREQOUT('TPHM'))) then
         ztphytd = (ztplus - ztplus0) * rcdt1
         ztphytd (:,nk) = ztdifv(:,nk)
      endif
      if (associated(zhuphytd) .and. associated(zhuplus0) .and. (ISREQSTEP('QPHY') .or. ISREQOUT('QPHM'))) then
         zhuphytd = (zhuplus - zhuplus0) * rcdt1
         zhuphytd(:,nk) = zqdifv(:,nk)
      endif
      if (associated(zqcphytd) .and. associated(zqcplus0)) &
           zqcphytd = (zqcplus - zqcplus0) * rcdt1
      if (diffuw .and. associated(zwphytd) .and. associated(zwplus0) .and. ISREQSTEP('WPHY')) &
           zwphytd  = (zwplus - zwplus0) * rcdt1

      if (ptp_L .and. kount >= 1) then
         if (associated(zuplus0)) zuplus = zuplus0
         if (associated(zvplus0)) zvplus = zvplus0
         if (associated(ztplus0)) ztplus = ztplus0
         ! if (associated(zhuplus0)) zhuplus = zhuplus0
         ! if (associated(zqcplus0)) zcqplus = zcqplus0
      endif

      if (timings_L) call timing_stop_omp(450)
      call msg_toall(MSG_DEBUG, 'tendency [END]')
      !-------------------------------------------------------------
      return
   end subroutine tendency5

   
   !/@*
   subroutine apply_tendencies_mm(zivar, ziten, ztdmaskxdt, ni, nk, F_nkscope, F_minval, F_maxval)
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      integer, intent(in) :: ni, nk
      real, dimension(ni, nk), intent(in) :: ziten
      real, dimension(ni), intent(in)  :: ztdmaskxdt
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
                    zivar(:,k) + ztdmaskxdt(:)*ziten(:,k)))
            end do
         else
            do k=1,nkscope
               zivar(:,k) = max(F_minval, &
                    zivar(:,k) + ztdmaskxdt(:)*ziten(:,k))
            end do
         endif
      else if (present(F_maxval)) then
         do k=1,nkscope
            zivar(:,k) = min(F_maxval, &
                 zivar(:,k) + ztdmaskxdt(:)*ziten(:,k))
         enddo
      else
         do k=1,nkscope
            zivar(:,k) = zivar(:,k) + ztdmaskxdt(:)*ziten(:,k)
         enddo
      endif
      !----------------------------------------------------------------
      return
   end subroutine apply_tendencies_mm

   
   !/@*
   subroutine apply_tendencies1(zivar1, ziten1, ztdmaskxdt, ni, nk, F_nkscope)
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      integer, intent(in) :: ni, nk
      real, dimension(ni, nk), intent(in) :: ziten1
      real, dimension(ni), intent(in)  :: ztdmaskxdt
      real, dimension(ni, nk), intent(inout) :: zivar1
      integer, intent(in), optional :: F_nkscope
      !*@/
      integer :: i, k, nkscope
      !----------------------------------------------------------------
      nkscope = nk
      if (present(F_nkscope)) nkscope = F_nkscope
      do k=1,nkscope
         do i=1,ni
            zivar1(i,k) = zivar1(i,k) + ztdmaskxdt(i)*ziten1(i,k)
         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine apply_tendencies1
   
   !/@*
   subroutine apply_tendencies2(zivar1, zivar2, ziten1, ziten2, ztdmaskxdt, ni, nk, F_nkscope)
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      integer, intent(in) :: ni, nk
      real, dimension(ni, nk), intent(in) :: ziten1, ziten2
      real, dimension(ni), intent(in)  :: ztdmaskxdt
      real, dimension(ni, nk), intent(inout) :: zivar1, zivar2
      integer, intent(in), optional :: F_nkscope
      !*@/
      integer :: i, k, nkscope
      !----------------------------------------------------------------
      nkscope = nk
      if (present(F_nkscope)) nkscope = F_nkscope
      do k=1,nkscope
         do i=1,ni
            zivar1(i,k) = zivar1(i,k) + ztdmaskxdt(i)*ziten1(i,k)
            zivar2(i,k) = zivar2(i,k) + ztdmaskxdt(i)*ziten2(i,k)
         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine apply_tendencies2

   !/@*
   subroutine apply_tendencies3(zivar1, zivar2, zivar3, ziten1, ziten2, ziten3, ztdmaskxdt, ni, nk, F_nkscope)
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      integer, intent(in) :: ni, nk
      real, dimension(ni, nk), intent(in) :: ziten1, ziten2, ziten3
      real, dimension(ni), intent(in)  :: ztdmaskxdt
      real, dimension(ni, nk), intent(inout) :: zivar1, zivar2, zivar3
      integer, intent(in), optional :: F_nkscope
      !*@/
      integer :: i, k, nkscope
      !----------------------------------------------------------------
      nkscope = nk
      if (present(F_nkscope)) nkscope = F_nkscope
      do k=1,nkscope
         do i=1,ni
            zivar1(i,k) = zivar1(i,k) + ztdmaskxdt(i)*ziten1(i,k)
            zivar2(i,k) = zivar2(i,k) + ztdmaskxdt(i)*ziten2(i,k)
            zivar3(i,k) = zivar3(i,k) + ztdmaskxdt(i)*ziten3(i,k)
         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine apply_tendencies3
   
   !/@*
   subroutine apply_tendencies4(zivar1, zivar2, zivar3, zivar4, ziten1, ziten2, ziten3, ziten4, ztdmaskxdt, ni, nk, F_nkscope)
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      !@Object Linear combination of two arrays
      !@Arguments
      integer, intent(in) :: ni, nk
      real, dimension(ni, nk), intent(in) :: ziten1, ziten2, ziten3, ziten4
      real, dimension(ni), intent(in)  :: ztdmaskxdt
      real, dimension(ni, nk), intent(inout) :: zivar1, zivar2, zivar3, zivar4
      integer, intent(in), optional :: F_nkscope
      !*@/
      integer :: i, k, nkscope
      !----------------------------------------------------------------
      nkscope = nk
      if (present(F_nkscope)) nkscope = F_nkscope
      do k=1,nkscope
         do i=1,ni
            zivar1(i,k) = zivar1(i,k) + ztdmaskxdt(i)*ziten1(i,k)
            zivar2(i,k) = zivar2(i,k) + ztdmaskxdt(i)*ziten2(i,k)
            zivar3(i,k) = zivar3(i,k) + ztdmaskxdt(i)*ziten3(i,k)
            zivar4(i,k) = zivar4(i,k) + ztdmaskxdt(i)*ziten4(i,k)
         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine apply_tendencies4

end module tendency
