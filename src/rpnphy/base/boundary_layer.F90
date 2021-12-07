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

module boundary_layer
   implicit none
   private
   public :: boundary_layer4

contains

   !/@*
   subroutine boundary_layer4(dbus, fbus, vbus, &
        ficebl, seloc, cdt1, kount, trnch, ni, nk)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use difver, only: difver8
      use energy_budget, only: eb_en,eb_pw,eb_residual_en,eb_residual_pw,eb_conserve_en,eb_conserve_pw,EB_OK
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use tendency, only: apply_tendencies
      use turbul, only: turbul2
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Arguments
      !          - Input -
      ! Dbus     dynamics input field
      !
      !          - Input/Output -
      ! Fbus     historic variables for the physics
      !
      !          - Output -
      ! Vbus     physics tendencies and other output fields from the physics
      !
      !          - Input -
      ! CDT1     timestep (sec.)
      ! TRNCH    slice number
      ! KOUNT    timestep number
      ! NI       horizontal running length
      ! NK       vertical dimension
      integer, intent(in) :: trnch, kount, ni, nk
      real, pointer, contiguous :: dbus(:), fbus(:), vbus(:)
      real,    intent(in) :: seloc(ni,nk), cdt1
      real,    intent(inout) :: ficebl(ni,nk)
      !@Author M. Desgagne summer 2011
      !*@/

      ! External definitions
#include "phymkptr.hf"
      include "surface.cdk"

      ! Local variables
      integer :: i, istat, istat2, nkm1
      real    :: wk1(ni,nk), wk2(ni,nk), rcdt1
      real, dimension(ni,nk), target :: zero
      real(REAL64), dimension(ni) :: l_en0, l_en, l_pw0, l_pw, l_enr, l_pwr
      ! Pointers to busdyn
      real, pointer, dimension(:,:), contiguous :: zqplus, ztplus, zuplus, zvplus, zumoins, zvmoins, zsigt, zqcplus, zwplus
      ! Pointers to busper
      real, pointer, dimension(:), contiguous   :: zqdiag, ztdiag, zudiag, zvdiag, zz0, zps, ztdmask
      real, pointer, dimension(:,:), contiguous :: zqtbl
      ! Pointers to busvol
      real, pointer, dimension(:), contiguous   :: zconepbl, zconqpbl, zflw
      real, pointer, dimension(:) :: zfc  !#TODO: should be contiguous
      real, pointer, dimension(:,:), contiguous :: zqdifv, ztdifv, zudifv, zvdifv, zkm, zkt, zgzmom, zgztherm, zldifv, &
           zwdifv, zqcdifv

      ! External symbols
      integer, external :: pbl_simple

      call init2nan(wk1,wk2)
      call init2nan(l_en0,l_en,l_pw0,l_pw,l_enr,l_pwr)

      ! Pointers to busdyn
      MKPTR2D(zqplus, huplus, dbus)
      MKPTR2D(zsigt, sigt, dbus)
      MKPTR2D(ztplus, tplus, dbus)
      MKPTR2D(zumoins, umoins, dbus)
      MKPTR2D(zuplus, uplus, dbus)
      MKPTR2D(zvmoins, vmoins, dbus)
      MKPTR2D(zvplus, vplus, dbus)
      MKPTR2D(zqcplus, qcplus, dbus)
      MKPTR2D(zwplus, wplus, dbus)
      ! Pointers to busper
      MKPTR1D(zps, pmoins, fbus)
      MKPTR1D(zqdiag, qdiag, fbus)
      MKPTR2D(zqtbl, qtbl, fbus)
      MKPTR1D(ztdiag, tdiag, fbus)
      MKPTR1D(zudiag, udiag, fbus)
      MKPTR1D(zvdiag, vdiag, fbus)
      MKPTR1D(zz0, z0, fbus)
      MKPTR1D(ztdmask, tdmask, fbus)
     ! Pointers to busvol
      MKPTR1D(zconepbl, conepbl, vbus)
      MKPTR1D(zconqpbl, conqpbl, vbus)
      MKPTR1DK(zfc, fc, indx_agrege, vbus)
      MKPTR1D(zflw, flw, vbus)
      MKPTR2D(zkm, km, vbus)
      MKPTR2D(zkt, kt, vbus)
      MKPTR2D(zgzmom, gzmom, vbus)
      MKPTR2D(zgztherm, gztherm, vbus)
      MKPTR2D(zqcdifv, qcdifv, vbus)
      MKPTR2D(zqdifv, qdifv, vbus)
      MKPTR2D(ztdifv, tdifv, vbus)
      MKPTR2D(zudifv, udifv, vbus)
      MKPTR2D(zvdifv, vdifv, vbus)
      MKPTR2D(zldifv, ldifv, vbus)
      MKPTR2D(zwdifv, wdifv, vbus)

      if (ldifv <= 0) zldifv => zero

      ! Initialization
      zero = 0.
      rcdt1 = 1./cdt1
      nkm1 = nk-1

      ! Turbulence closures used to compute diffusion coefficients
      if (any(fluvert == (/ &
           'MOISTKE', &
           'CLEF   '/))) then
         call turbul2(dbus, fbus, vbus, seloc, kount, trnch, ni, nk, nkm1)
         if (phy_error_L) return
      elseif (fluvert == 'SIMPLE') then
         istat = pbl_simple(zkm,zkt,zumoins,zvmoins,zgzmom,zgztherm,zz0,ni,nk)
         if (.not.RMN_IS_OK(istat)) then
            call physeterror('boundary_layer', 'Problem in phy_simple')
            return
         endif
      endif

      ! Pre-scheme state for energy budget
      if (associated(zconepbl)) then
         istat  = eb_en(l_en0,ztplus,zqplus,zqtbl,zsigt,zps,nkm1,F_inttype='linear')
         istat2 = eb_pw(l_pw0,zqplus,zqtbl,zsigt,zps,nkm1,F_inttype='linear')
         if (istat /= EB_OK .or.istat2  /= EB_OK) then
            call physeterror('boundary_layer', 'Problem computing preliminary energy budget for '//trim(fluvert))
            return
         endif
      endif

      ! Compute diagnostic-level tendencies as the change in diagnostic values (there
      ! appears to be no good reason for why the PBL scheme should be responsible
      ! for this)
      do i=1,ni
         zqdifv(i,nk) = (zqdiag(i) - zqplus(i,nk))*rcdt1
         ztdifv(i,nk) = (ztdiag(i) - ztplus(i,nk))*rcdt1
         zudifv(i,nk) = (zudiag(i) - zuplus(i,nk))*rcdt1
         zvdifv(i,nk) = (zvdiag(i) - zvplus(i,nk))*rcdt1
         zqcdifv(i,nk) = 0.
         zqplus(i,nk) =  zqdiag(i)
         ztplus(i,nk) =  ztdiag(i)
         zuplus(i,nk) =  zudiag(i)
         zvplus(i,nk) =  zvdiag(i)
      end do

      ! Apply diffusion operator to compute PBL tendencies
      call difver8(dbus, fbus, vbus, seloc, &
           cdt1, kount, trnch, ni, nk, nkm1)
      if (phy_error_L) return

      ! Impose conservation by adjustment of moisture and temperature tendencies
      TENDENCY_ADJUSTMENT: if (pbl_conserve == 'TEND') then

         ! Apply humidity tendency correction for total water conservation
         istat = eb_conserve_pw(zqdifv,zqdifv,ztplus,zqplus,zsigt,zps,nkm1, &
              F_dqc=zldifv,F_inttype='linear',F_lhf=zflw)
         ! Apply temperature tendency correction for liquid water static energy conservation
         istat2 = eb_conserve_en(ztdifv,ztdifv,zqdifv,ztplus,zqplus,zqtbl,zsigt,zps,nkm1, &
              F_dqc=zldifv,F_inttype='linear',F_shf=zfc,F_lhf=zflw)
         if (istat /= EB_OK .or.istat2  /= EB_OK) then
            call physeterror('boundary_layer', 'Problem correcting for liquid water static energy conservation in boundary_layer')
            return
         endif

      endif TENDENCY_ADJUSTMENT

      ! Apply computed tendencies to state variables
      call apply_tendencies(zqcplus, zqcdifv, ztdmask, ni, nk, nkm1)
      call apply_tendencies(zqplus, zqdifv, ztdmask, ni, nk, nkm1)
      call apply_tendencies(ztplus, ztdifv, ztdmask, ni, nk, nkm1)
      call apply_tendencies(zuplus, zudifv, ztdmask, ni, nk, nkm1)
      call apply_tendencies(zvplus, zvdifv, ztdmask, ni, nk, nkm1)
      if (diffuw) call apply_tendencies(zwplus, zwdifv, ztdmask, ni, nk, nkm1)

      ! Post-scheme energy budget analysis
      if (associated(zconepbl)) then
         ! Compute post-scheme state
         istat  = eb_en(l_en,ztplus,zqplus,zqtbl,zsigt,zps,nkm1,F_inttype='linear')
         istat2 = eb_pw(l_pw,zqplus,zqtbl,zsigt,zps,nkm1,F_inttype='linear')
         if (istat == EB_OK .and. istat2 == EB_OK) then
            ! Compute residuals
            istat  = eb_residual_en(l_enr,l_en0,l_en,ztplus,zqplus,zqcplus,delt,nkm1,F_shf=zfc,F_lhf=zflw)
            istat2 = eb_residual_pw(l_pwr,l_pw0,l_pw,ztplus,delt,nkm1,F_lhf=zflw)
         endif
         if (istat /= EB_OK .or.istat2  /= EB_OK) then
            call physeterror('boundary_layer', 'Problem computing final energy budget for '//trim(fluvert))
            return
         endif
         zconepbl(:) = real(l_enr)
         zconqpbl(:) = real(l_pwr)
      endif

      ! Compute ice fraction for later use (in cloud water section)
      if (fluvert == 'MOISTKE') then
         call ficemxp2(ficebl, wk1, wk2, ztplus, ni, nkm1)
         if (phy_error_L) return
      endif

      ! End of subprogram
      return
   end subroutine boundary_layer4

end module boundary_layer
