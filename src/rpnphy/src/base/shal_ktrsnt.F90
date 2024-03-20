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

module shal_ktrsnt
   implicit none
   private
   public :: sc_ktrsnt

contains

   !/@*
  subroutine sc_ktrsnt(pvars, dt, ni, nkm1)
    use, intrinsic :: iso_fortran_env, only: REAL64
    use debug_mod, only: init2nan
    use phy_status, only: phy_error_L, PHY_OK
    use phybusidx
    use phymem, only: phyvar
    use tdpack_const, only: RAUW
    use phy_options, only: conv_shal
    use cnv_options, only: shal_conserve
    use phybudget, only: pb_compute, pb_conserve, pb_residual
    use shallconv, only: shallconv5
    use tendency, only: apply_tendencies
    implicit none
#include <rmnlib_basics.hf>
    
    !@Objective Call to Kuo-transient shallow convection
    
    !@Arguments
    type(phyvar), pointer, contiguous :: pvars(:)  !all phy vars (meta + slab data)
    real, intent(in) :: dt                           !Time step (s)
    integer, intent(in) :: ni                        !Row length
    integer, intent(in) :: nkm1                      !Levels (not including diagnostic)
    
    !@Author Ron McTaggart-Cowan, Feb 2022
    
    !@Revisions
    !*@/
#include "phymkptr.hf"

    ! Local variables
    real, pointer, dimension(:), contiguous :: psm, ztdmaskxdt, &
         zkshal, zconesc, zconqsc, ztlcs, ztscs
    real, pointer, dimension(:,:), contiguous :: qqm, sigma, ttp, &
         zgztherm, zhushal, ztshal, zfsc, zqlsc, zqssc, zqcz, &
         zqdifv, qqp
    real(kind=REAL64), dimension(ni) :: l_en0, l_pw0
    
    ! Return if not running Kuo Transient scheme
    if (conv_shal /= 'KTRSNT') return

    ! Basic setup
    call init2nan(l_en0, l_pw0)

    ! Extract bus information
    MKPTR1D(psm, pmoins, pvars)
    MKPTR1D(ztdmaskxdt, tdmaskxdt, pvars)
    MKPTR1D(zkshal, kshal, pvars)
    MKPTR1D(zconesc, conesc, pvars)
    MKPTR1D(zconqsc, conqsc, pvars)
    MKPTR1D(ztlcs, tlcs, pvars)
    MKPTR1D(ztscs, tscs, pvars)
    MKPTR2Dm1(qqm, humoins, pvars)
    MKPTR2Dm1(qqp, huplus, pvars)
    MKPTR2Dm1(sigma, sigw, pvars)
    MKPTR2Dm1(ttp, tplus, pvars)
    MKPTR2Dm1(zgztherm, gztherm, pvars)
    MKPTR2Dm1(zhushal, hushal, pvars)
    MKPTR2Dm1(ztshal, tshal, pvars)
    MKPTR2Dm1(zfsc, fsc, pvars)
    MKPTR2Dm1(zqlsc, qlsc, pvars)
    MKPTR2Dm1(zqssc, qssc, pvars)
    MKPTR2Dm1(zqcz, qcz, pvars)
    MKPTR2Dm1(zqdifv, qdifv, pvars)

    ! Pre-shallow state for energy budget
    if (pb_compute(zconesc, zconqsc, l_en0, l_pw0, &
         pvars, nkm1) /= PHY_OK) then
       call physeterror('shal_ktrsnt::sc_ktrsnt', &
            'Problem computing preliminary budget for '//trim(conv_shal))
       return
    endif

    ! Shallow convective scheme (Kuo transient types)
    call shallconv5(ttp, qqm, sigma, zgztherm, psm, &
         ztshal, zhushal, ztlcs, ztscs, zqlsc, zqssc, zkshal, &
         zfsc, zqcz, zqdifv, dt, ni, nkm1, nkm1)
    if (phy_error_L) return

    ! Adjust tendencies to impose conservation
    if (pb_conserve(shal_conserve, ztshal, zhushal, pvars, &
         F_rain=ztlcs*RAUW, F_snow=ztscs*RAUW) /= PHY_OK) then
       call physeterror('shal_ktrsnt::sc_ktrsnt', &
            'Cannot correct conservation for '//trim(conv_shal))
       return
    endif

    ! Apply shallow convective tendencies to state variables
    call apply_tendencies(ttp,    qqp, &
         &                ztshal, zhushal, ztdmaskxdt, ni, nkm1, nkm1)

    ! Post-scheme energy budget analysis
    if (pb_residual(zconesc, zconqsc, l_en0, l_pw0, pvars, &
         dt, nkm1, F_rain=ztlcs*RAUW, F_snow=ztscs*RAUW) /= PHY_OK) then
       call physeterror('shal_ktrsnt::sc_ktrsnt', &
            'Problem computing final budget for '//trim(conv_shal))
       return
    endif

  end subroutine sc_ktrsnt
  
end module shal_ktrsnt
