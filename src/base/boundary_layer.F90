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
   subroutine boundary_layer4(pvars, &
        ficebl, seloc, cdt1, kount, ni, nk, trnch)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use difver, only: difver8
      use phybudget, only: pb_compute, pb_conserve, pb_residual
      use phy_options
      use phy_status, only: phy_error_L, PHY_OK
      use phybusidx
      use phymem, only: phyvar
      use tendency, only: apply_tendencies
      use turbul, only: turbul2
      use pbl_ysu, only: pbl_ysu1
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Arguments
      !          - Input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - Input -
      ! CDT1     timestep (sec.)
      ! TRNCH    slice number
      ! KOUNT    timestep number
      ! NI       horizontal running length
      ! NK       vertical dimension
      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: trnch, kount, ni, nk
      real,    intent(in) :: seloc(ni,nk), cdt1
      real,    intent(inout) :: ficebl(ni,nk)
      !@Author M. Desgagne summer 2011
      !*@/

      ! External definitions
#include "phymkptr.hf"
      include "surface.cdk"

      ! Local variables
      integer :: i, istat, nkm1
      real    :: wk1(ni,nk), wk2(ni,nk), rcdt1
      real, dimension(ni,nk), target :: zero
      real(REAL64), dimension(ni) :: l_en0, l_pw0
      real, pointer, dimension(:,:), contiguous :: zqplus, ztplus, zuplus, zvplus, zumoins, zvmoins, zsigt, zqcplus, zwplus
      real, pointer, dimension(:), contiguous   :: zqdiag, ztdiag, zudiag, zvdiag, zz0, zps, ztdmask
      real, pointer, dimension(:,:), contiguous :: zqtbl
      real, pointer, dimension(:), contiguous   :: zconepbl, zconqpbl, zflw
      real, pointer, dimension(:) :: zfc  !#TODO: should be contiguous
      real, pointer, dimension(:,:), contiguous :: zqdifv, ztdifv, zudifv, zvdifv, zkm, zkt, zgzmom, zgztherm, zldifv, &
           zwdifv, zqcdifv, ztmoins, zqmoins, ztve

      ! External symbols
      integer, external :: pbl_simple

      call init2nan(wk1,wk2)
      call init2nan(l_en0,l_pw0)

      MKPTR2D(zqmoins, humoins, pvars)
      MKPTR2D(zqplus, huplus, pvars)
      MKPTR2D(zsigt, sigt, pvars)
      MKPTR2D(ztmoins, tmoins, pvars)
      MKPTR2D(ztplus, tplus, pvars)
      MKPTR2D(zumoins, umoins, pvars)
      MKPTR2D(zuplus, uplus, pvars)
      MKPTR2D(zvmoins, vmoins, pvars)
      MKPTR2D(zvplus, vplus, pvars)
      MKPTR2D(zqcplus, qcplus, pvars)
      MKPTR2D(zwplus, wplus, pvars)
      MKPTR1D(zps, pmoins, pvars)
      MKPTR1D(zqdiag, qdiag, pvars)
      MKPTR2D(zqtbl, qtbl, pvars)
      MKPTR1D(ztdiag, tdiag, pvars)
      MKPTR1D(zudiag, udiag, pvars)
      MKPTR1D(zvdiag, vdiag, pvars)
      MKPTR1D(zz0, z0, pvars)
      MKPTR1D(ztdmask, tdmask, pvars)
      MKPTR1D(zconepbl, conepbl, pvars)
      MKPTR1D(zconqpbl, conqpbl, pvars)
      MKPTR1DK(zfc, fc, indx_agrege, pvars)
      MKPTR1D(zflw, flw, pvars)
      MKPTR2D(zkm, km, pvars)
      MKPTR2D(zkt, kt, pvars)
      MKPTR2D(zgzmom, gzmom, pvars)
      MKPTR2D(zgztherm, gztherm, pvars)
      MKPTR2D(zqcdifv, qcdifv, pvars)
      MKPTR2D(zqdifv, qdifv, pvars)
      MKPTR2D(ztve, tve, pvars)      
      MKPTR2D(ztdifv, tdifv, pvars)
      MKPTR2D(zudifv, udifv, pvars)
      MKPTR2D(zvdifv, vdifv, pvars)
      MKPTR2D(zldifv, ldifv, pvars)
      MKPTR2D(zwdifv, wdifv, pvars)

      if (ldifv <= 0) zldifv => zero

      ! Initialization
      zero = 0.
      rcdt1 = 1./cdt1
      nkm1 = nk-1
      
      ! Compute thermodynamic quantities on energy levels
      call mfotvt(ztve, ztmoins, zqmoins, ni, nkm1, ni)

      ! Turbulence closures used to compute diffusion coefficients
      if (any(fluvert == (/ &
           'MOISTKE', &
           'CLEF   '/))) then
         call turbul2(pvars, seloc, kount, ni, nk, nkm1, trnch)
         if (phy_error_L) return
      elseif (fluvert == 'SIMPLE') then
         istat = pbl_simple(zkm,zkt,zumoins,zvmoins,zgzmom,zgztherm,zz0,ni,nk)
         if (.not.RMN_IS_OK(istat)) then
            call physeterror('boundary_layer', 'Problem in phy_simple')
            return
         endif
      elseif (fluvert == 'YSU') then
         istat = pbl_ysu1(pvars, cdt1, ni, nk, nkm1, trnch)
         if (phy_error_L) return
      endif

      ! Pre-PBL state for budget
      if (pb_compute(zconepbl, zconqpbl, l_en0, l_pw0, &
           pvars, nkm1, F_inttype='linear') /= PHY_OK) then
         call physeterror('boundary_layer', &
              'Problem computing preliminary budget for '//trim(fluvert))
         return
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
      call difver8(pvars, seloc, cdt1, kount, ni, nk, nkm1, trnch)
      if (phy_error_L) return
         
      ! Adjust PBL tendencies to impose conservation
      if (pb_conserve(pbl_conserve, ztdifv, zqdifv, pvars, &
           F_dqc=zldifv, F_shf=zfc, F_wvf=zflw, F_inttype='linear') /= PHY_OK) then
         call physeterror('boundary_layer', &
              'Cannot correct conservation for '//trim(fluvert))
         return
      endif

      ! Apply computed tendencies to state variables
      call apply_tendencies(zqplus, zqdifv, ztdmask, ni, nk, nkm1)
      call apply_tendencies(ztplus, ztdifv, ztdmask, ni, nk, nkm1)
      call apply_tendencies(zuplus, zudifv, ztdmask, ni, nk, nkm1)
      call apply_tendencies(zvplus, zvdifv, ztdmask, ni, nk, nkm1)
      if (associated(zqcplus)) &
           call apply_tendencies(zqcplus, zqcdifv, ztdmask, ni, nk, nkm1)
      if (diffuw) call apply_tendencies(zwplus, zwdifv, ztdmask, ni, nk, nkm1)
      
      ! Post-PBL budget analysis
      if (pb_residual(zconepbl, zconqpbl, l_en0, l_pw0, pvars, &
           delt, nkm1, F_shf=zfc, F_wvf=zflw, F_inttype='linear') /= PHY_OK) then
         call physeterror('boundary_layer', &
              'Problem computing final budget for '//trim(fluvert))
         return
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
