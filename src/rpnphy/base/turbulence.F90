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
   subroutine turbulence2(d, f, v, dsiz, fsiz, vsiz, ficebl, seloc, cdt1, &
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
      ! e        entry    input field
      ! d        dynamics input field
      !          - Input/Output -
      ! f        historic variables for the physics
      !          - Output -
      ! v        physics tendencies and other output fields from the physics
      !          - Input -
      ! esiz     dimension of e
      ! dsiz     dimension of d
      ! fsiz     dimension of f
      ! vsiz     dimension of v
      ! dt       timestep (sec.)
      ! trnch    slice number
      ! kount    timestep number
      ! icpu     cpu number executing slice "trnch"
      ! n        horizontal running length
      ! nk       vertical dimension

      integer, intent(in) :: dsiz, fsiz, vsiz, trnch, kount, ni, nk
      real, dimension(ni,nk), intent(inout) :: ficebl
      real, dimension(ni,nk), intent(in) :: seloc
      real, intent(in) :: cdt1
      real, target, intent(inout) :: d(dsiz), f(fsiz), v(vsiz)

      !@Author L. Spacek (Nov 2011)
      !*@/
#include <msg.h>
#include "phymkptr.hf"
      include "surface.cdk"

      ! Internal variables
      integer :: istat, istat1, j
      real dsig
      real, dimension(ni) :: rhos
      real, dimension(ni,nk) :: dkem, dket
      real, dimension(:), pointer :: ztsrad, zpmoins, zutautofd, zvtautofd, zsigs, zgzmom1, zgzmom2, ztdmask
      real, dimension(:,:), pointer :: zutofd, zvtofd, zgzmom, zumoins, zvmoins, zz0, zuplus, zvplus, zsigt, ztplus, zttofd
      real, pointer, dimension(:,:,:) :: zvcoef
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'turbulence [BEGIN]')
      if (timings_L) call timing_start_omp(430, 'turbulence', 46)

      ! Reshape bus entries
      MKPTR1D(zpmoins, pmoins, f)
      MKPTR1D(zsigs, sigs, f)
      MKPTR1D(ztdmask, tdmask, f)
      MKPTR1D(ztsrad, tsrad, f)
      MKPTR1D(zutautofd, utautofd, v)
      MKPTR1D(zvtautofd, vtautofd, v)
      MKPTR2D(zsigt, sigt, d)

      MKPTR2D(zgzmom, gzmom, v)
      MKPTR2D(zumoins, umoins, d)
      MKPTR2D(zuplus, uplus, d)
      MKPTR2D(zutofd, utofd, v)
      MKPTR2D(zvmoins, vmoins, d)
      MKPTR2D(zvplus, vplus, d)
      MKPTR2D(zvtofd, vtofd, v)
      MKPTR2D(ztplus, tplus, d)
      MKPTR2D(zttofd, ttofd, v)
      MKPTR2DN(zz0, z0, ni, nagrege, f)

      MKPTR3D(zvcoef, vcoef, 2, v)
      
      call init2nan(rhos)
      call init2nan(dkem,dket)

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
         istat  = int_profile(zutautofd, zutofd, zgzmom, zgzmom2, zgzmom1)
         istat1 = int_profile(zvtautofd, zvtofd, zgzmom, zgzmom2, zgzmom1)
         if (istat /= INT_OK .or. istat1 /= INT_OK) then
            call physeterror('turbulence', 'Problem in int_profile')
            return
         endif
         rhos = zpmoins / (RGASD*ztsrad)
         zutautofd = zutautofd * rhos
         zvtautofd = zvtautofd * rhos

      endif

      ! Turbulent vertical diffusion using the PBL scheme
      if (any(fluvert == (/&
           'MOISTKE', &
           'CLEF   ', &
           'SIMPLE '/))) then
         call boundary_layer4(d, f, v, dsiz, fsiz, vsiz, ficebl, seloc, cdt1, kount, trnch, ni, nk)
         if (phy_error_L) return
      endif

      if (timings_L) call timing_stop_omp(430)
      call msg_toall(MSG_DEBUG, 'turbulence [END]')
      !----------------------------------------------------------------
      return
   end subroutine turbulence2

end module turbulence
