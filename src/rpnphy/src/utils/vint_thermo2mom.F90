
module vintphy
   implicit none
   private
   public :: vint_thermo2mom, vint_mom2thermo
   public :: vint_thermo2mom1, vint_mom2thermo1
   public :: vint_thermo2mom2, vint_mom2thermo2
   public :: vint_thermo2mom3, vint_mom2thermo3

   interface vint_thermo2mom
      module procedure vint_thermo2mom1
      module procedure vint_thermo2mom2
      module procedure vint_thermo2mom3
   end interface vint_thermo2mom
   
   interface vint_mom2thermo
      module procedure vint_mom2thermo1
      module procedure vint_mom2thermo2
      module procedure vint_mom2thermo3
   end interface vint_mom2thermo

contains

   !/@*
   subroutine vint_thermo2mom1(fldout,fldin,vcoef,ni,nk)
      use phygridmap, only: phydim_nk
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      integer, intent(in)  :: ni, nk         !# dimensions
      real,    intent(in)  :: vcoef(ni,phydim_nk,2) !# interp coef
      real,    intent(in)  :: fldin(ni,nk)   !# var on thermo levels
      real,    intent(out) :: fldout(ni,nk)  !# var on momentum levels
      !@Author L. Spacek (Dec 2007)
      !@Object
      ! Vertical interpolation of T/Q to momentum or thermo levels
      !*@/
      integer :: k, km1
      do k=nk,1,-1
         km1 = max(1,k-1)
         fldout(:,k) = fldin(:,k) - vcoef(:,k,1) * (fldin(:,k) - fldin(:,km1))
      enddo
      return
   end subroutine vint_thermo2mom1

   
   !/@*
   subroutine vint_thermo2mom2(fldout1,fldout2,fldin1,fldin2,vcoef,ni,nk)
      use phygridmap, only: phydim_nk
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      integer, intent(in)  :: ni, nk          !# dimensions
      real,    intent(in)  :: vcoef(ni,phydim_nk,2) !# interp coef
      real,    intent(in)  :: fldin1(ni,nk)   !# var on thermo levels
      real,    intent(in)  :: fldin2(ni,nk)   !# var on thermo levels
      real,    intent(out) :: fldout1(ni,nk)  !# var on momentum levels
      real,    intent(out) :: fldout2(ni,nk)  !# var on momentum levels
      !@Author L. Spacek (Dec 2007)
      !@Object
      ! Vertical interpolation of T/Q to momentum or thermo levels
      !*@/
      integer :: i, k, km1
      do k=nk,1,-1
         km1 = max(1,k-1)
         do i=1,ni
            fldout1(i,k) = fldin1(i,k) - vcoef(i,k,1) * (fldin1(i,k) - fldin1(i,km1))
            fldout2(i,k) = fldin2(i,k) - vcoef(i,k,1) * (fldin2(i,k) - fldin2(i,km1))
         enddo
      enddo
      return
   end subroutine vint_thermo2mom2


   !/@*
   subroutine vint_thermo2mom3(fldout1,fldout2,fldout3,fldin1,fldin2,fldin3,vcoef,ni,nk)
      use phygridmap, only: phydim_nk
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      integer, intent(in)  :: ni, nk          !# dimensions
      real,    intent(in)  :: vcoef(ni,phydim_nk,2) !# interp coef
      real,    intent(in)  :: fldin1(ni,nk)   !# var on thermo levels
      real,    intent(in)  :: fldin2(ni,nk)   !# var on thermo levels
      real,    intent(in)  :: fldin3(ni,nk)   !# var on thermo levels
      real,    intent(out) :: fldout1(ni,nk)  !# var on momentum levels
      real,    intent(out) :: fldout2(ni,nk)  !# var on momentum levels
      real,    intent(out) :: fldout3(ni,nk)  !# var on momentum levels
      !@Author L. Spacek (Dec 2007)
      !@Object
      ! Vertical interpolation of T/Q to momentum or thermo levels
      !*@/
      integer :: i, k, km1
      do k=nk,1,-1
         km1 = max(1,k-1)
         do i=1,ni
            fldout1(i,k) = fldin1(i,k) - vcoef(i,k,1) * (fldin1(i,k) - fldin1(i,km1))
            fldout2(i,k) = fldin2(i,k) - vcoef(i,k,1) * (fldin2(i,k) - fldin2(i,km1))
            fldout3(i,k) = fldin3(i,k) - vcoef(i,k,1) * (fldin3(i,k) - fldin3(i,km1))
         enddo
      enddo
      return
   end subroutine vint_thermo2mom3


   !/@*
   subroutine vint_mom2thermo1(fldout,fldin,vcoef,ni,nk)
      use phygridmap, only: phydim_nk
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      integer, intent(in)  :: ni, nk         !#  dimensions
      real,    intent(in)  :: vcoef(ni,phydim_nk,2) !# interp coef
      real,    intent(in)  :: fldin(ni,nk)   !# var on momentum levels
      real,    intent(out) :: fldout(ni,nk)  !# var on thermo levels
      !@Author L. Spacek (Dec 2007)
      !@Object
      ! Vertical interpolation of T/Q to momentum or thermo levels
      !*@/
      integer :: k, kp1
      do k=1,nk
         kp1 = min(nk,k+1)
         fldout(:,k) = fldin(:,k) + vcoef(:,k,2) * (fldin(:,kp1) - fldin(:,k))
      enddo
      return
   end subroutine vint_mom2thermo1


   !/@*
   subroutine vint_mom2thermo2(fldout1,fldout2,fldin1,fldin2,vcoef,ni,nk)
      use phygridmap, only: phydim_nk
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      integer, intent(in)  :: ni, nk          !#  dimensions
      real,    intent(in)  :: vcoef(ni,phydim_nk,2) !# interp coef
      real,    intent(in)  :: fldin1(ni,nk)   !# var on momentum levels
      real,    intent(in)  :: fldin2(ni,nk)   !# var on momentum levels
      real,    intent(out) :: fldout1(ni,nk)  !# var on thermo levels
      real,    intent(out) :: fldout2(ni,nk)  !# var on thermo levels
      !@Author L. Spacek (Dec 2007)
      !@Object
      ! Vertical interpolation of T/Q to momentum or thermo levels
      !*@/
      integer :: i, k, kp1
      do k=1,nk
         kp1 = min(nk,k+1)
         do i=1,ni
            fldout1(:,k) = fldin1(:,k) + vcoef(:,k,2) * (fldin1(:,kp1) - fldin1(:,k))
            fldout2(:,k) = fldin2(:,k) + vcoef(:,k,2) * (fldin2(:,kp1) - fldin2(:,k))
         enddo
      enddo
      return
   end subroutine vint_mom2thermo2

   
   !/@*
   subroutine vint_mom2thermo3(fldout1,fldout2,fldout3,fldin1,fldin2,fldin3,vcoef,ni,nk)
      use phygridmap, only: phydim_nk
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      integer, intent(in)  :: ni, nk          !#  dimensions
      real,    intent(in)  :: vcoef(ni,phydim_nk,2) !# interp coef
      real,    intent(in)  :: fldin1(ni,nk)   !# var on momentum levels
      real,    intent(in)  :: fldin2(ni,nk)   !# var on momentum levels
      real,    intent(in)  :: fldin3(ni,nk)   !# var on momentum levels
      real,    intent(out) :: fldout1(ni,nk)  !# var on thermo levels
      real,    intent(out) :: fldout2(ni,nk)  !# var on thermo levels
      real,    intent(out) :: fldout3(ni,nk)  !# var on thermo levels
      !@Author L. Spacek (Dec 2007)
      !@Object
      ! Vertical interpolation of T/Q to momentum or thermo levels
      !*@/
      integer :: i, k, kp1
      do k=1,nk
         kp1 = min(nk,k+1)
         do i=1,ni
            fldout1(:,k) = fldin1(:,k) + vcoef(:,k,2) * (fldin1(:,kp1) - fldin1(:,k))
            fldout2(:,k) = fldin2(:,k) + vcoef(:,k,2) * (fldin2(:,kp1) - fldin2(:,k))
            fldout3(:,k) = fldin3(:,k) + vcoef(:,k,2) * (fldin3(:,kp1) - fldin3(:,k))
         enddo
      enddo
      return
   end subroutine vint_mom2thermo3

end module vintphy
