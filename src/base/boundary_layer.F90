
module boundary_layer
   implicit none
   private
   public :: boundary_layer4

contains

   !/@*
   subroutine boundary_layer4(pvars, cdt1, kount, ni, nk, trnch)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use phybudget, only: pb_compute, pb_conserve, pb_residual, INT_TYPE_LINEAR
      use phy_options
      use phy_status, only: phy_error_L, PHY_OK
      use phybusidx
      use phymem, only: phyvar
      use tendency, only: apply_tendencies
      use pbl_maintke, only: maintke
      use pbl_ysu, only: pbl_ysu1
      use pbl_sim, only: simplepbl
      use pbl_ri_diffuse, only: diffuseall
      use pbl_diffuse, only: difver
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
      real,    intent(in) :: cdt1
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
      real, pointer, dimension(:,:), contiguous :: zqplus, ztplus, zuplus, zvplus, zumoins, zvmoins, zqcplus, zwplus
      real, pointer, dimension(:), contiguous   :: zqdiag, ztdiag, zudiag, zvdiag, zz0, zps, ztdmaskxdt
      real, pointer, dimension(:), contiguous   :: zconepbl, zconqpbl, zflw
      real, pointer, dimension(:) :: zfc  !#TODO: should be contiguous
      real, pointer, dimension(:,:), contiguous :: zqdifv, ztdifv, zudifv, zvdifv, zkm, zkt, zgzmom, zgztherm, &
           zwdifv, zqcdifv, ztmoins, zqmoins, ztve, zvcoef
      logical :: dqc_applied

      ! External symbols
      integer, external :: pbl_simple

      call init2nan(wk1,wk2)
      call init2nan(l_en0,l_pw0)

      MKPTR2D(zqmoins, humoins, pvars)
      MKPTR2D(zqplus, huplus, pvars)
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
      MKPTR1D(ztdiag, tdiag, pvars)
      MKPTR1D(zudiag, udiag, pvars)
      MKPTR1D(zvdiag, vdiag, pvars)
      MKPTR1D(zz0, z0, pvars)
      MKPTR1D(ztdmaskxdt, tdmaskxdt, pvars)
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
      MKPTR2D(zvcoef, vcoef, pvars)
      MKPTR2D(zwdifv, wdifv, pvars)
      
      ! Initialization
      zero = 0.
      rcdt1 = 1./cdt1
      nkm1 = nk-1
      
      ! Compute thermodynamic quantities on energy levels
      call mfotvt(ztve, ztmoins, zqmoins, ni, nkm1, ni)

      ! Turbulence closures used to compute diffusion coefficients
      if (any(fluvert == (/ &
           'MOISTKE', &
           'RPNINT '/))) then
         call maintke(pvars, kount, ni, nk, nkm1, trnch)
         if (phy_error_L) return
      elseif (fluvert == 'SIMPLE') then
         istat = simplepbl(zkm,zkt,zumoins,zvmoins,zgzmom,zgztherm,zz0,ni,nk)
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
           pvars, nkm1, F_inttype=INT_TYPE_LINEAR) /= PHY_OK) then
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
      if (fluvert == 'RPNINT') then
         call diffuseall(pvars, dqc_applied, cdt1, ni, nk, nkm1)
      else
         call difver(pvars, cdt1, kount, ni, nk, nkm1, trnch)
         dqc_applied = .false.
      endif
      if (phy_error_L) return
         
      ! Adjust PBL tendencies to impose conservation
      if (pb_conserve(pbl_conserve, ztdifv, zqdifv, pvars, &
           F_dqc=zqcdifv, F_shf=zfc, F_wvf=zflw, F_inttype=INT_TYPE_LINEAR) /= PHY_OK) then
         call physeterror('boundary_layer', &
              'Cannot correct conservation for '//trim(fluvert))
         return
      endif

      ! Apply computed tendencies to state variables
      call apply_tendencies(zqplus, ztplus, zuplus, zvplus, &
           &                zqdifv, ztdifv, zudifv, zvdifv, &
           &                ztdmaskxdt, ni, nk, nkm1)
      if (associated(zqcplus) .and. .not.dqc_applied) &
           call apply_tendencies(zqcplus, zqcdifv, ztdmaskxdt, ni, nk, nkm1)
      if (diffuw) call apply_tendencies(zwplus, zwdifv, ztdmaskxdt, ni, nk, nkm1)
      
      ! Post-PBL budget analysis
      if (pb_residual(zconepbl, zconqpbl, l_en0, l_pw0, pvars, &
           delt, nkm1, F_shf=zfc, F_wvf=zflw, F_inttype=INT_TYPE_LINEAR) /= PHY_OK) then
         call physeterror('boundary_layer', &
              'Problem computing final budget for '//trim(fluvert))
         return
      endif

      ! End of subprogram
      return
   end subroutine boundary_layer4

end module boundary_layer
