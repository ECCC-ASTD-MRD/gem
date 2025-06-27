
module turbulence
   implicit none
   private
   public :: turbulence2

contains

   !/@*
   subroutine turbulence2(pvars, cdt1, kount, ni, nk, trnch)
      use debug_mod, only: init2nan
      use tdpack_const, only: RGASD, CPD
      use boundary_layer, only: boundary_layer4
      use integrals, only: int_profile, INT_OK
      use form_drag, only: form_drag1
      use phy_options
      use phy_status, only: phy_error_L
      use phybusidx
      use phymem, only: phyvar
      use tendency, only: apply_tendencies
      use vintphy, only: vint_mom2thermo
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Object
      !@Arguments
      !          - Input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - Input -
      ! dt       timestep (sec.)
      ! trnch    slice number
      ! kount    timestep number
      ! n        horizontal running length
      ! nk       vertical dimension

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: trnch, kount, ni, nk
      real, intent(in) :: cdt1

      !@Author L. Spacek (Nov 2011)
      !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"
      include "surface.cdk"

      ! Internal variables
      integer :: istat, istat1, j, k, i
      real dsig
      real, dimension(ni,nk) :: dkem, dket, rho
      real, dimension(:), pointer, contiguous :: ztsrad, zpplus, zutautofd, zvtautofd, zsigs, zgzmom1, zgzmom2, ztdmaskxdt
      real, dimension(:,:), pointer, contiguous :: zutofd, zvtofd, zgzmom, zumoins, zvmoins, zz0, zuplus, zvplus, zsigt, ztplus, zttofd, zsigm
      real, pointer, dimension(:,:,:), contiguous :: zvcoef
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'turbulence [BEGIN]')
      if (timings_L) call timing_start_omp(430, 'turbulence', 46)

      ! Reshape bus entries
      MKPTR1D(zpplus, pplus, pvars)
      MKPTR1D(zsigs, sigs, pvars)
      MKPTR1D(ztdmaskxdt, tdmaskxdt, pvars)
      MKPTR1D(ztsrad, tsrad, pvars)
      MKPTR1D(zutautofd, utautofd, pvars)
      MKPTR1D(zvtautofd, vtautofd, pvars)
      MKPTR2D(zsigt, sigt, pvars)
      MKPTR2D(zsigm, sigm, pvars)

      MKPTR2D(zgzmom, gzmom, pvars)
      MKPTR2D(zumoins, umoins, pvars)
      MKPTR2D(zuplus, uplus, pvars)
      MKPTR2D(zutofd, utofd, pvars)
      MKPTR2D(zvmoins, vmoins, pvars)
      MKPTR2D(zvplus, vplus, pvars)
      MKPTR2D(zvtofd, vtofd, pvars)
      MKPTR2D(ztplus, tplus, pvars)
      MKPTR2D(zttofd, ttofd, pvars)
      MKPTR2D(zz0, z0, pvars)

      MKPTR3D(zvcoef, vcoef, pvars)
      
      call init2nan(dkem,dket,rho)

      ! Turbulent orographic form drag using the "distributed drag" scheme
      if (tofd == 'BELJAARS04') then

         ! Compute and apply tendencies
         istat = form_drag1(zutofd, zvtofd, zumoins, zvmoins, zgzmom, zsigs, zz0, cdt1, ni, nk)
         if (.not.RMN_IS_OK(istat)) then
            call physeterror('turbulence', 'Problem in form_drag')
            return
         endif

         call apply_tendencies(zuplus, zvplus, &
              &                zutofd, zvtofd, ztdmaskxdt, ni, nk, nk)
         ! dissipative heating calculations
         
         dkem = zuplus(:,1:nk) * zutofd(:,1:nk) + &
                zvplus(:,1:nk) * zvtofd(:,1:nk)
         call vint_mom2thermo(dket, dkem, zvcoef, ni, nk)
         do j=1,ni
            dsig = (zsigt(j,nk) - 1.)/(zsigt(j,nk-1) - 1.)
            dket(j,nk) = dsig * dket(j,nk-1)
         end do
         zttofd = - (1./CPD) * dket         
         call apply_tendencies(ztplus, zttofd, ztdmaskxdt, ni, nk, nk)
         
         ! Diagnose surface stress
         if (ISREQSTEPL((/"UTFD", "VTFD"/))) then
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

      endif

      ! Turbulent vertical diffusion using the PBL scheme
      if (.not.(fluvert == 'NIL' .or. fluvert == 'SURFACE')) then
         call boundary_layer4(pvars, cdt1, kount, ni, nk, trnch)
         if (phy_error_L) return
      endif

      if (timings_L) call timing_stop_omp(430)
      call msg_toall(MSG_DEBUG, 'turbulence [END]')
      !----------------------------------------------------------------
      return
   end subroutine turbulence2

end module turbulence
