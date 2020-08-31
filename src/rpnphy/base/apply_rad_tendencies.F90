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

module apply_rad_tendencies
   implicit none
   private
   public :: apply_rad_tendencies1

contains

   !/@*
   subroutine apply_rad_tendencies1(d, dsiz, v, vsiz, f, fsiz, ni, nk, nkm1)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use energy_budget, only: eb_en, eb_pw, eb_residual_en, eb_residual_pw, eb_conserve_en, EB_OK
      use phy_options
      use phybus
      use tendency, only: apply_tendencies
      use ens_perturb, only: ens_nc2d, ens_spp_get
      implicit none
!!!#include <arch_specific.hf>
      !@Object Apply radiative tendencies
      !@Arguments
      !          - Input -
      ! dsiz     dimension of dbus
      ! vsiz     dimension of vbus
      ! fsiz     dimension of fbus
      ! ni       horizontal running length
      ! nk       vertical dimension
      ! nkm1     vertical scope of the operator

      integer, intent(in) :: fsiz, vsiz, dsiz, ni, nk, nkm1
      real, intent(inout), target :: d(dsiz), v(vsiz), f(fsiz)

      !@Author Ron McTaggart-Cowan, Summer 2017
      !*@/
#include <msg.h>
#include <rmnlib_basics.hf>
#include "phymkptr.hf"

      ! Local variables
      integer ::  istat, istat2, k
      real, dimension(ni) :: tradmult
      real, dimension(ni,nkm1) :: qrad, mtrad
      real, dimension(:), pointer :: zconerad,zconqrad, zps, zfusi, zfdsi, zfdss, ziv, &
           zev, zei, ztdmask, znetrad
      real, dimension(:,:), pointer :: ztplus, zqplus, zqcplus, zgztherm, zsigt, ztrad, &
           zmrk2
      real(REAL64), dimension(ni) :: l_en0, l_pw0, l_en, l_pw, l_enr, l_pwr
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'apply_rad_tendencies [BEGIN]')
      if (timings_L) call timing_start_omp(422, 'apply_rad_td', 46)

      ! Bus pointer association
      MKPTR1D(zconerad, conerad, v)
      MKPTR1D(zconqrad, conqrad, v)
      MKPTR1D(zei, ei, f)
      MKPTR1D(zev, ev, f)
      MKPTR1D(zfdsi, fdsi, f)
      MKPTR1D(zfdss, fdss, f)
      MKPTR1D(zfusi, fusi, f)
      MKPTR1D(ziv, iv, v)
      MKPTR1D(znetrad, netrad, v)
      MKPTR1D(zps, pmoins, f)
      MKPTR1D(ztdmask, tdmask, f)

      MKPTR2Dm1(zgztherm, gztherm, v)
      MKPTR2Dm1(zqcplus, qcplus, d)
      MKPTR2Dm1(zqplus, huplus, d)
      MKPTR2Dm1(zsigt, sigt, d)
      MKPTR2Dm1(ztplus, tplus, d)
      MKPTR2Dm1(ztrad, trad, v)
      
      MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, f)
      
      ! Early exit if radiation is not used
      if (radia == 'NIL') then
         if (associated(zconerad)) zconerad = 0.
         if (associated(zconqrad)) zconqrad = 0.
         if (timings_L) call timing_stop_omp(422)
         return
      endif

      call init2nan(qrad)
      call init2nan(l_en0, l_pw0, l_en, l_pw, l_enr, l_pwr)

      ! Pre-scheme state for energy budget
      if (associated(zconerad)) then
         istat  = eb_en(l_en0, ztplus, zqplus, zqcplus, zsigt, zps, nkm1)
         istat2 = eb_pw(l_pw0, zqplus, zqcplus, zsigt, zps, nkm1)
         if (istat /= EB_OK .or. istat2 /= EB_OK) then
            call physeterror('apply_rad_td', 'Problem computing preliminary energy budget for '//trim(radia))
            return
         endif
      endif

      ! Compute source of column integral of liquid water static energy:
      !  (iv-ev): TOA SW net
      !  fdss: SFC SW absorption (equals net by balance)
      !  (0.-ei): TOA LW net
      !  (fdsi-fusi): SW LW net
      znetrad = ((ziv-zev) - zfdss) + ((0.-zei) - (zfdsi-zfusi))

      ! Impose conservation by adjustment of moisture and temperature tendencies
      TENDENCY_ADJUSTMENT: if (rad_conserve == 'TEND') then

         ! Apply temperature tendency correction for liquid water static energy conservation
         qrad = 0.
         istat = eb_conserve_en(ztrad, ztrad, qrad, ztplus, zqplus, zqcplus, zsigt, zps, nkm1, F_rad=znetrad)
         if (istat /= EB_OK) then
            call physeterror('apply_rad_td', 'Problem correcting for liquid water static energy conservation')
            return
         endif

      endif TENDENCY_ADJUSTMENT

      ! Apply tendency perturbation on request
      tradmult = ens_spp_get('trad_mult', zmrk2, default=1.)
      do k=1,nkm1
         mtrad(:,k) = tradmult(:) * ztrad(:,k)
      enddo
      
      ! Apply radiative tendencies
      call apply_tendencies(ztplus, mtrad, ztdmask, ni, nk, nkm1)

      ! Post-scheme energy budget analysis
      if (associated(zconerad)) then
         ! Compute post-scheme state
         istat  = eb_en(l_en, ztplus, zqplus, zqcplus, zsigt, zps, nkm1)
         istat2 = eb_pw(l_pw, zqplus, zqcplus, zsigt, zps, nkm1)
         if (istat == EB_OK .and. istat2 == EB_OK) then
            ! Compute residuals
            istat  = eb_residual_en(l_enr, l_en0, l_en, ztplus, zqplus, zqcplus, delt, nkm1, F_rad=znetrad)
            istat2 = eb_residual_pw(l_pwr, l_pw0, l_pw, ztplus, delt, nkm1)
         endif
         if (istat /= EB_OK .or. istat2 /= EB_OK) then
            call physeterror('apply_rad_td', 'Problem computing final energy budget for '//trim(radia))
            return
         endif
         zconerad(:) = real(l_enr)
         zconqrad(:) = real(l_pwr)
      endif

      if (timings_L) call timing_stop_omp(422)
      call msg_toall(MSG_DEBUG, 'apply_rad_tendencies [END]')
      !----------------------------------------------------------------
      return
   end subroutine apply_rad_tendencies1

end module apply_rad_tendencies
