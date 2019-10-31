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

module prep_cw
   use debug_mod, only: init2nan
   use phy_options
   use phybus
   implicit none
   private
   public :: prep_cw3

#include <msg.h>
#include "phymkptr.hf"

contains

  !/@*
   subroutine prep_cw3(f, fsiz, d, dsiz, v, vsiz, ficebl, ni, nk)
      implicit none
!!!#include <arch_specific.hf>

      !@Object Save water contents and cloudiness in the permanent bus
      !@Arguments
      !          - Input -
      ! dsiz     dimension of d
      ! fsiz     dimension of f
      ! vsiz     dimension of v
      ! ficebl   fraction of ice
      ! ni       horizontal dimension
      ! nk       vertical dimension
      !          - Input/Output -
      ! d        dynamic             bus
      ! f        permanent variables bus
      ! v        volatile (output)   bus
      !*@/
      integer, intent(in) :: fsiz, dsiz, vsiz, ni, nk
      real, intent(inout), target :: f(fsiz), d(dsiz), v(vsiz)
      real, intent(in) :: ficebl(ni,nk)
      !*@/
      integer :: nkm1
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'prep_cw [BEGIN]')
      if (timings_L) call timing_start_omp(445, 'prep_cw', 46)
      if (stcond(1:3) == 'MP_') then
         nkm1 = nk-1
         call prep_cw_MP(f, fsiz, v, vsiz, ficebl, ni, nk, nkm1)
      else
         call prep_cw_noMP(f, fsiz, d, dsiz, v, vsiz, ni, nk)
      endif
      if (timings_L) call timing_stop_omp(445)
      call msg_toall(MSG_DEBUG, 'prep_cw [END]')
      !----------------------------------------------------------------
      return
   end subroutine prep_cw3


   !/@*
   subroutine prep_cw_noMP(f, fsiz, d, dsiz, v, vsiz, ni, nk)
      implicit none
!!!#include <arch_specific.hf>

      !@Object Save water contents and cloudiness in the permanent bus
      !@Arguments
      !          - Input -
      ! dsiz     dimension of d
      ! fsiz     dimension of f
      ! vsiz     dimension of v
      ! ni       horizontal dimension
      ! nk       vertical dimension
      !          - Input/Output -
      ! d        dynamic             bus
      ! f        permanent variables bus
      ! v        volatile (output)   bus
      !*@/

      integer, intent(in) :: fsiz, dsiz, vsiz, ni, nk
      real, intent(inout), target :: f(fsiz), d(dsiz), v(vsiz)

      !@Author L. Spacek (Oct 2004)
      !*@/

      integer :: i, k
      real :: cfblxp(ni,nk)
      real, target :: zero(ni,nk)
      real, pointer, dimension(:,:) :: zfbl, zfdc, zfsc, zftot, zfxp, ziwc, zlwc, &
           zqcplus, zqiplus, zqldi, zqlsc, zqlmi, zqsmi, zfmc, &
           zqsdi, zqssc, zqtbl, zsnow, zqi_cat1, zqi_cat2, zqi_cat3, zqi_cat4
      !----------------------------------------------------------------
      MKPTR2D(zfbl, fbl, f)
      MKPTR2D(zfdc, fdc, f)
      MKPTR2D(zfmc, fmc, f)
      MKPTR2D(zfsc, fsc, f)
      MKPTR2D(zftot, ftot, f)
      MKPTR2D(zfxp, fxp, f)
      MKPTR2D(ziwc, iwc, f)
      MKPTR2D(zlwc, lwc, f)
      MKPTR2D(zqcplus, qcplus, d)
      MKPTR2D(zqi_cat1, qti1plus, d)
      MKPTR2D(zqi_cat2, qti2plus, d)
      MKPTR2D(zqi_cat3, qti3plus, d)
      MKPTR2D(zqi_cat4, qti4plus, d)
      MKPTR2D(zqiplus, qiplus, d)
      MKPTR2D(zqldi, qldi, f)
      MKPTR2D(zqlmi, qlmi, v)
      MKPTR2D(zqlsc, qlsc, v)
      MKPTR2D(zqsdi, qsdi, f)
      MKPTR2D(zqsmi, qsmi, v)
      MKPTR2D(zqssc, qssc, v)
      MKPTR2D(zqtbl, qtbl, f)
      MKPTR2D(zsnow, qnplus, d)

      zero(:,:) =  0.0

      if (any(convec == (/ &
           'NIL   ', &
           'KUOSTD', &
           'OLDKUO'  &
           /))) then
         zqldi => zero(1:ni,1:nk)
         zqsdi => zero(1:ni,1:nk)
      endif

      if (conv_shal == 'NIL') then
         zqlsc => zero(1:ni,1:nk)
         zqssc => zero(1:ni,1:nk)
         zfsc => zero(1:ni,1:nk)
      endif

      if (conv_mid == 'NIL') then
         zqlmi => zero(1:ni,1:nk)
         zqsmi => zero(1:ni,1:nk)
         zfmc => zero(1:ni,1:nk)
      endif

      call init2nan(cfblxp)

      ! ------------------------------------------
      ! Cloud water
      ! ------------------------------------------

      if (stcond /= 'NIL') then

         ! qcplus may contain total water content from consun/newsund and detrained explicit clouds from kfc/bech
         ! if (stcond = 'MP') qcplus is liquid clouds from expicit scheme + detrained explicit liquid clouds from kfc/bech

         do k=1,nk
            do i=1,ni
               zlwc(i,k) = zqcplus(i,k)
               cfblxp (i,k) = zfxp(i,k)
            enddo
         enddo
      else
         cfblxp = 0.
      endif

      ! If we have the CONSUN scheme, the cloud water from MoisTKE has
      ! priority over the cloud water from the grid-scale scheme.

      if (stcond == 'CONSUN' .and. fluvert == 'MOISTKE') then
         do k=1,nk
            do i=1,ni
               if (zqtbl(i,k) > zqcplus(i,k))then
                  zlwc(i,k) = zqtbl(i,k)
                  cfblxp(i,k) = zfbl(i,k)
               endif
            enddo
         enddo
      endif

      ! Add the cloud water (liquid and solid) coming from shallow and deep
      ! cumulus clouds (only for the Kuo Transient and Kain-Fritsch schemes).
      ! Note that no conditions are used for these calculations ...
      ! qldi, qsdi, and qlsc, qssc are zero if these schemes are not used.
      ! Also note that qldi, qsdi, qlsc and qssc are NOT IN-CLOUD values
      ! (multiplication done in kfcp4 and ktrsnt)
      !
      ! For Sundqvist schemes all
      ! the cloud water is put in LWC (and will be partitioned later in prep_cw_rad)

      zlwc = zlwc + zqldi + zqsdi + zqlsc + zqssc + zqlmi + zqsmi

      ! Combine explicit and Implicit clouds using the random overlap
      ! approximation:
      !     CFBLXP is either fxp or fbl depending on choice of merging between explicit and moistke
      !     FDC is for deep convection clouds
      !          (always defined as necessary for condensation too)
      !     FSC is for the shallow convection clouds
      !     FMC is for mid-level convective clouds

      do k=1,nk
         do i=1,ni
            zftot(i,k) = min(1., max(0., &
                 1. - (1.-cfblxp(i,k))*(1.-zfdc(i,k))*(1.-zfsc(i,k))*(1.-zfmc(i,k)) &
                 ))
         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine prep_cw_noMP


   !/@*
   subroutine prep_cw_MP(f, fsiz,  v, vsiz, ficebl, ni, nk, nkm1)
      implicit none
!!!#include <arch_specific.hf>

      !@Object
      ! When MP schemes are used; merge water contents from "implicit cloud
      ! sources"  for radiation
      !
      !@Arguments
      !
      !      - Input -
      ! dsiz     dimension of d
      ! fsiz     dimension of f
      ! vsiz     dimension of v
      ! ficebl   fraction of ice
      ! ni       horizontal dimension
      ! nk       vertical dimension
      ! nkm1     vertical scope of the operator
      !
      !      - Input/Output -
      ! d        dynamic             bus
      ! f        permanent variables bus
      ! v        volatile (output)   bus

      integer, intent(in) :: fsiz, vsiz, ni, nk, nkm1
      real, intent(inout), target :: f(fsiz), v(vsiz)
      real, intent(in) :: ficebl(ni,nk)

      !@Author
      ! D. Paquin-Ricard (June 2017)
      ! P. Vaillancourt (July  2016)
      !*@/

      integer :: i, k
      real, target :: zero(ni,nk)
      real, pointer, dimension(:,:) :: zfbl, zfdc, zfsc, zftot, zfxp, zfmp,  &
           ziwcimp, zlwcimp, zqldi, zqlsc, zqsdi, zqssc, zqtbl, zfmc, &
           zqlmi, zqsmi
      !----------------------------------------------------------------
      MKPTR2D(zfbl, fbl, f)
      MKPTR2D(zfdc, fdc, f)
      MKPTR2D(zfmc, fmc, f)
      MKPTR2D(zfmp, fmp, f)
      MKPTR2D(zfsc, fsc, f)
      MKPTR2D(zftot, ftot, f)
      MKPTR2D(zfxp, fxp, f)
      MKPTR2D(zlwcimp, lwcimp, f)
      MKPTR2D(ziwcimp, iwcimp, f)
      MKPTR2D(zqldi, qldi, f)
      MKPTR2D(zqlmi, qlmi, v)
      MKPTR2D(zqlsc, qlsc, v)
      MKPTR2D(zqsdi, qsdi, f)
      MKPTR2D(zqsmi, qsmi, v)
      MKPTR2D(zqssc, qssc, v)
      MKPTR2D(zqtbl, qtbl, f)

      zero(1:ni,1:nkm1) =  0.0

      if (any(convec == (/ &
           'NIL   ', &
           'KUOSTD', &
           'OLDKUO'  &
           /))) then
         zqldi => zero(1:ni,1:nkm1)
         zqsdi => zero(1:ni,1:nkm1)
         if (convec == 'NIL') zfdc => zero(1:ni,1:nkm1)
      endif

      if (conv_shal == 'NIL') then
         zqlsc => zero(1:ni,1:nkm1)
         zqssc => zero(1:ni,1:nkm1)
         zfsc => zero(1:ni,1:nkm1)
      endif

      if (conv_mid == 'NIL') then
         zqlmi => zero(1:ni,1:nkm1)
         zqsmi => zero(1:ni,1:nkm1)
         zfmc => zero(1:ni,1:nkm1)
      endif

      if (fluvert /= 'MOISTKE') then
         zfbl => zero(1:ni,1:nkm1)
         zqtbl => zero(1:ni,1:nkm1)
      endif

      ! PV-avril2016 : simplified version for MP schemes only
      !          - agregates implicit clouds only;
      !          - assumes that moistke is an implicit source of clouds amongst others (eliminate choice between mtke and exp clouds)
      !          - could choose a maximum overlap of clouds instead of random?
      !          - condensates are agregated assuming max overlap while fractions assume random ???
      !
      ! Add the cloud water (liquid and solid) coming from PBL, shallow  and deep cumulus clouds
      ! note that all condensates must be GRID-SCALE values (not in-cloud)

      if (fluvert == 'MOISTKE') then
         do k=1,nkm1
            do i=1,ni
               zlwcimp(i,k) = zqtbl(i,k) * (1.0 - ficebl(i,k))
               ziwcimp(i,k) = zqtbl(i,k) * ficebl(i,k)
            enddo
         enddo
      endif

      do k=1,nkm1
         do i=1,ni
            zlwcimp(i,k) = zlwcimp(i,k) + zqldi(i,k) + zqlsc(i,k) + zqlmi(i,k)
            ziwcimp(i,k) = ziwcimp(i,k) + zqsdi(i,k) + zqssc(i,k) + zqsmi(i,k)
         enddo
      enddo


      ! Agregate implicit cloud fractions
      !     FBL is for PBL clouds from moistke
      !     FDC is for deep convection clouds
      !     FSC is for the shallow convection clouds
      !     FMC is for mid-level convective clouds
      !     FMP is the sum over the implicit, FXP: the explicit, FTOT: total
      do k=1,nkm1
         do i=1,ni
            zfmp(i,k) = min(1., max(0., &     ! random overlap
                 1. - (1.-zfbl(i,k))*(1.-zfdc(i,k))*(1.-zfsc(i,k))*(1.-zfmc(i,k)) &
                 ))
            zftot(i,k) = min(1., max(0., &    ! maximum overlap
                 zfmp(i,k)+zfxp(i,k) &
                 ))

         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine prep_cw_MP

end module prep_cw
