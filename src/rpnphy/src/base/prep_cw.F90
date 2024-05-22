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
   use phybusidx
   use phymem, only: phyvar
   use pbl_utils, only: ficemxp
   implicit none
   private
   public :: prep_cw3

#include <rmn/msg.h>
#include "phymkptr.hf"

contains

  !/@*
   subroutine prep_cw3(pvars, ni, nk)
      implicit none
!!!#include <arch_specific.hf>

      !@Object Save water contents and cloudiness in the permanent bus
      !@Arguments
      !          - Input -
      ! ni       horizontal dimension
      ! nk       vertical dimension
      !          - Input/Output -
      ! pvars    list of all phy vars (meta + slab data)
      !*@/
      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni, nk
      !*@/
      integer :: nkm1
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'prep_cw [BEGIN]')
      if (timings_L) call timing_start_omp(445, 'prep_cw', 46)
      if (stcond(1:3) == 'MP_') then
         nkm1 = nk-1
         call prep_cw_MP(pvars, ni, nk, nkm1)
      else
         call prep_cw_noMP(pvars, ni, nk)
      endif
      if (timings_L) call timing_stop_omp(445)
      call msg_toall(MSG_DEBUG, 'prep_cw [END]')
      !----------------------------------------------------------------
      return
   end subroutine prep_cw3


   !/@*
   subroutine prep_cw_noMP(pvars, ni, nk)
      implicit none
!!!#include <arch_specific.hf>

      !@Object
      ! When  non MP schemes are used; merge cloud fractions and water contents from both implicit and  explicit cloud
      ! sources"  for radiation
      !@Arguments
      !          - Input -
      ! ni       horizontal dimension
      ! nk       vertical dimension
      !          - Input/Output -
      ! pvars    list of all phy vars (meta + slab data)
      !*@/

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni, nk

      !@Author L. Spacek (Oct 2004)
      !*@/

      integer :: i, k
      real :: cfblxp(ni,nk)
      real, target :: zero(ni,nk)
      real, pointer, dimension(:,:), contiguous :: zfbl, zfdc, zfsc, zftot, zfxp, zlwc, &
           zqcplus, zqldi, zqlsc, zqlmi, zqsmi, zfmc, &
           zqsdi, zqssc, zqtbl
      !----------------------------------------------------------------
      MKPTR2D(zfbl, fbl, pvars)
      MKPTR2D(zfdc, fdc, pvars)
      MKPTR2D(zfmc, fmc, pvars)
      MKPTR2D(zfsc, fsc, pvars)
      MKPTR2D(zftot, ftot, pvars)
      MKPTR2D(zfxp, fxp, pvars)
      MKPTR2D(zlwc, lwc, pvars)
      MKPTR2D(zqcplus, qcplus, pvars)
      MKPTR2D(zqldi, qldi, pvars)
      MKPTR2D(zqlmi, qlmi, pvars)
      MKPTR2D(zqlsc, qlsc, pvars)
      MKPTR2D(zqsdi, qsdi, pvars)
      MKPTR2D(zqsmi, qsmi, pvars)
      MKPTR2D(zqssc, qssc, pvars)
      
      if (advecqtbl) then
         MKPTR2D(zqtbl, qtblplus, pvars)
      else
         MKPTR2D(zqtbl, qtbl, pvars)
      endif

      zero(:,:) =  0.0

      if (convec == 'NIL') then
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

         ! qcplus may contain total water content from consun and detrained explicit clouds from kfc/bech
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

      
      ! Recette 1 - La traditionnelle
      ! the cloud water from MoisTKE has
      ! priority over the cloud water from the grid-scale scheme.
      ! water contents are combined assuming no overlap
      ! cloud fractions are combined assuming random overlap
      
      if (fluvert == 'MOISTKE' .or. fluvert == 'RPNINT') then
         if (associated(zqcplus)) then
            where (zqtbl(:,:) > zqcplus(:,:))
               zlwc(:,:) = zqtbl(:,:)
               cfblxp(:,:) = zfbl(:,:)
            endwhere
         else
            zlwc(:,:) = zqtbl(:,:)
            cfblxp(:,:) = zfbl(:,:)
         endif
      endif
      zlwc = zlwc + zqldi + zqsdi + zqlsc + zqssc + zqlmi + zqsmi
      do k=1,nk
         do i=1,ni
            zftot(i,k) = min(1., max(0., &
                 1. - (1.-cfblxp(i,k))*(1.-zfdc(i,k))*(1.-zfsc(i,k))*(1.-zfmc(i,k)) &
                 ))
         enddo
      enddo

!!$      ! Recette 2 - Min overlap for all sources
!!$      if(.false.) then
!!$        zlwc =  zqcplus + zqtbl + zqldi + zqsdi + zqlsc + zqssc + zqlmi + zqsmi
!!$        do k=1,nk
!!$           do i=1,ni
!!$              zftot(i,k) = zfxp(i,k)+zfbl(i,k)+zfdc(i,k)+zfsc(i,k)+zfmc(i,k)  
!!$              zftot(i,k) = min(1., max(0., zftot(i,k) )) 
!!$           enddo
!!$        enddo
!!$      endif

!!$      ! Recette 3 - Min overlap for all implicit sources(including pbl) ; random overlap to combine explicit and implicit
!!$      if(.false.) then
!!$         zlwc =  zqtbl + zqldi + zqsdi + zqlsc + zqssc + zqlmi + zqsmi
!!$         do k=1,nk
!!$            do i=1,ni
!!$               cfblxp(i,k) = min(1., max(0., zfbl(i,k)+zfdc(i,k)+zfsc(i,k)+zfmc(i,k) )) 
!!$               zftot(i,k) = min(1., max(0., 1. - (1.-cfblxp(i,k))*(1.-zfxp(i,k)) ))  
!!$               zlwc(i,k) = zlwc(i,k) + zqcplus(i,k)*(1.-cfblxp(i,k))
!!$            enddo
!!$         enddo
!!$      endif

      !----------------------------------------------------------------
      return
   end subroutine prep_cw_noMP


   !/@*
   subroutine prep_cw_MP(pvars, ni, nk, nkm1)
      implicit none
!!!#include <arch_specific.hf>

      !@Object
      ! When MP schemes are used; merge cloud fractions and water contents from "implicit cloud
      ! sources"  for radiation
      !
      !@Arguments
      !
      !      - Input -
      ! ni       horizontal dimension
      ! nk       vertical dimension
      ! nkm1     vertical scope of the operator
      !
      !      - Input/Output -
      ! pvars list of all phy vars (meta + slab data)

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni, nk, nkm1

      !@Author
      ! D. Paquin-Ricard (June 2017)
      ! P. Vaillancourt (July  2016)
      !*@/

      integer :: i, k
      real, dimension(ni,nk) :: unused, ficebl
      real, target :: zero(ni,nk)
      real, pointer, dimension(:,:), contiguous :: zfbl, zfdc, zfsc, zftot, zfxp, zfmp,  &
           ziwcimp, zlwcimp, zqldi, zqlsc, zqsdi, zqssc, zqtbl, zfmc, &
           zqlmi, zqsmi, ztplus
      !----------------------------------------------------------------
      MKPTR2D(zfbl, fbl, pvars)
      MKPTR2D(zfdc, fdc, pvars)
      MKPTR2D(zfmc, fmc, pvars)
      MKPTR2D(zfmp, fmp, pvars)
      MKPTR2D(zfsc, fsc, pvars)
      MKPTR2D(zftot, ftot, pvars)
      MKPTR2D(zfxp, fxp, pvars)
      MKPTR2D(zlwcimp, lwcimp, pvars)
      MKPTR2D(ziwcimp, iwcimp, pvars)
      MKPTR2D(zqldi, qldi, pvars)
      MKPTR2D(zqlmi, qlmi, pvars)
      MKPTR2D(zqlsc, qlsc, pvars)
      MKPTR2D(zqsdi, qsdi, pvars)
      MKPTR2D(zqsmi, qsmi, pvars)
      MKPTR2D(zqssc, qssc, pvars)
      MKPTR2D(ztplus, tplus, pvars)

      if (advecqtbl) then
         MKPTR2D(zqtbl, qtblplus, pvars)
      else
         MKPTR2D(zqtbl, qtbl, pvars)
      endif
      
      zero(1:ni,1:nkm1) =  0.0

      if (convec == 'NIL') then
         zqldi => zero(1:ni,1:nkm1)
         zqsdi => zero(1:ni,1:nkm1)
         zfdc => zero(1:ni,1:nkm1)
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

      if (fluvert /= 'MOISTKE' .or. fluvert == 'RPNINT') then
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

      if (fluvert == 'MOISTKE' .or. fluvert == 'RPNINT') then
         call ficemxp(ficebl, unused, unused, ztplus, ni, nkm1)
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
            zfmp(i,k) = min(1., max(0., &     ! random overlap  == OPS
                 1. - (1.-zfbl(i,k))*(1.-zfdc(i,k))*(1.-zfsc(i,k))*(1.-zfmc(i,k)) &
                 ))
!            zfmp(i,k) = min(1., max(0., zfbl(i,k)+zfdc(i,k)+zfsc(i,k)+zfmc(i,k) )) ! no overlap
            zftot(i,k) = min(1., max(0., &    ! no overlap - should use same recipe as in  cldoppro_mp; this variable is diagnostic
                 zfmp(i,k)+zfxp(i,k) &
                 ))

         enddo
      enddo
      !----------------------------------------------------------------
      return
   end subroutine prep_cw_MP

end module prep_cw
