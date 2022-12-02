!-------------------------------------- LICENCE BEGIN --------------------------
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
!-------------------------------------- LICENCE END ----------------------------

module turbulence
   implicit none
   private
   public :: turbulence2

contains

   !/@*
   subroutine turbulence2(dbus, fbus, vbus, ficebl, seloc, cdt1, &
        kount, trnch, ni, nk)
      use debug_mod, only: init2nan
      use tdpack_const, only: RGASD, CPD
      use boundary_layer, only: boundary_layer4
      use integrals, only: int_profile, INT_OK
      use form_drag, only: form_drag1
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use tendency, only: apply_tendencies
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Object
      !@Arguments
      !          - Input -
      ! dbus     dynamics input field
      !          - Input/Output -
      ! fbus     historic variables for the physics
      !          - Output -
      ! vbus     physics tendencies and other output fields from the physics
      !          - Input -
      ! dt       timestep (sec.)
      ! trnch    slice number
      ! kount    timestep number
      ! icpu     cpu number executing slice "trnch"
      ! n        horizontal running length
      ! nk       vertical dimension

      integer, intent(in) :: trnch, kount, ni, nk
      real, dimension(ni,nk), intent(inout) :: ficebl
      real, dimension(ni,nk), intent(in) :: seloc
      real, intent(in) :: cdt1
      real, pointer, contiguous :: dbus(:), fbus(:), vbus(:)

      !@Author L. Spacek (Nov 2011)
      !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"
      include "surface.cdk"

      ! Internal variables
      integer :: istat, istat1, j, k
      real dsig
      real, dimension(ni,nk) :: dkem, dket, rho
      real, dimension(:), pointer, contiguous :: ztsrad, zpplus, zutautofd, zvtautofd, zsigs, zgzmom1, zgzmom2, ztdmask
      real, dimension(:,:), pointer, contiguous :: zutofd, zvtofd, zgzmom, zumoins, zvmoins, zz0, zuplus, zvplus, zsigt, ztplus, zttofd, zsigm
      real, pointer, dimension(:,:,:), contiguous :: zvcoef
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'turbulence [BEGIN]')
      if (timings_L) call timing_start_omp(430, 'turbulence', 46)

      ! Reshape bus entries
      MKPTR1D(zpplus, pplus, vbus)
      MKPTR1D(zsigs, sigs, fbus)
      MKPTR1D(ztdmask, tdmask, fbus)
      MKPTR1D(ztsrad, tsrad, fbus)
      MKPTR1D(zutautofd, utautofd, vbus)
      MKPTR1D(zvtautofd, vtautofd, vbus)
      MKPTR2D(zsigt, sigt, dbus)
      MKPTR2D(zsigm, sigm, dbus)

      MKPTR2D(zgzmom, gzmom, vbus)
      MKPTR2D(zumoins, umoins, dbus)
      MKPTR2D(zuplus, uplus, dbus)
      MKPTR2D(zutofd, utofd, vbus)
      MKPTR2D(zvmoins, vmoins, dbus)
      MKPTR2D(zvplus, vplus, dbus)
      MKPTR2D(zvtofd, vtofd, vbus)
      MKPTR2D(ztplus, tplus, dbus)
      MKPTR2D(zttofd, ttofd, vbus)
      MKPTR2DN(zz0, z0, ni, nagrege, fbus)

      MKPTR3D(zvcoef, vcoef, 2, vbus)
      
      call init2nan(dkem,dket,rho)

      ! Turbulent orographic form drag using the "distributed drag" scheme
      if (tofd == 'BELJAARS04') then

         ! Compute and apply tendencies
         istat = form_drag1(zutofd, zvtofd, zumoins, zvmoins, zgzmom, zsigs, zz0, cdt1, ni, nk)
         if (.not.RMN_IS_OK(istat)) then
            call physeterror('turbulence', 'Problem in form_drag')
            return
         endif

         call apply_tendencies(zuplus, zutofd, ztdmask, ni, nk, nk)
         call apply_tendencies(zvplus, zvtofd, ztdmask, ni, nk, nk)
         
         ! dissipative heating calculations
         
         dkem = zuplus(:,1:nk) * zutofd(:,1:nk) + &
                zvplus(:,1:nk) * zvtofd(:,1:nk)
         call vint_mom2thermo(dket, dkem, zvcoef, ni, nk)
         do j=1,ni
            dsig = (zsigt(j,nk) - 1.)/(zsigt(j,nk-1) - 1.)
            dket(j,nk) = dsig * dket(j,nk-1)
         end do
         zttofd = - (1./CPD) * dket         
         call apply_tendencies(ztplus, zttofd, ztdmask, ni, nk, nk)
         
         ! Diagnose surface stress
         zgzmom1 => zgzmom(:,1)
         zgzmom2 => zgzmom(:,nk-1)
         do k=1,nk
            rho(:,k) = zpplus(:) * zsigm(:,k) / (RGASD * ztplus(:,k))
         enddo
         istat  = int_profile(zutautofd, -rho*zutofd, zgzmom, zgzmom2, zgzmom1)
         istat1 = int_profile(zvtautofd, -rho*zvtofd, zgzmom, zgzmom2, zgzmom1)
         if (istat /= INT_OK .or. istat1 /= INT_OK) then
            call physeterror('turbulence', 'Problem in int_profile')
            return
         endif

      endif

      ! Turbulent vertical diffusion using the PBL scheme
      if (.not.(fluvert == 'NIL' .or. fluvert == 'SURFACE')) then
         call boundary_layer4(dbus, fbus, vbus, ficebl, seloc, cdt1, kount, trnch, ni, nk)
         if (phy_error_L) return
      endif

      if (timings_L) call timing_stop_omp(430)
      call msg_toall(MSG_DEBUG, 'turbulence [END]')
      !----------------------------------------------------------------
      return
   end subroutine turbulence2

end module turbulence
