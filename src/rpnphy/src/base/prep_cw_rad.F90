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

module prep_cw_rad
   implicit none
   private
   public :: prep_cw_rad3

contains

   !/@*
   subroutine prep_cw_rad3(fbus, dbus, &
        tm, qm, ps, sigma, cloud, &
        liqwcin, icewcin, liqwpin, icewpin, &
        trav2d,  &
        kount, trnch, ni, nk, nkm1)
      use, intrinsic :: iso_fortran_env, only: INT64
      use debug_mod, only: init2nan
      use tdpack_const, only: GRAV, TCDK, RGASD
      use phy_options
      use phybus
      use series_mod, only: series_xst
      implicit none
      !@Author L. Spacek (Oct 2004)
      !@Object  Prepare liquid/ice water contents and cloudiness for the radiation
      !@Arguments
      !     - input -
      !     tm       temperature
      !     qm       specific humidity
      !     ps       surface pressure
      !     sigma    sigma levels
      !     kount    index of timestep
      !     trnch    number of the slice
      !     task     task number
      !     ni       horizontal dimension
      !     nk       number of layers, including the diag level
      !     nkm1     number of layers, excluding the diag level
      !     - output -
      !     liqwcin  in-cloud liquid water content
      !     icewcin  in-cloud ice    water content
      !     liqwpin  in-cloud liquid water path (g/m^2)
      !     icewpin  in-cloud ice    water path (g/m^2)
      !     cloud    cloudiness passed to radiation

      integer, intent(in) :: kount, trnch, ni, nk, nkm1
      real, pointer, contiguous :: dbus(:), fbus(:)
      real, intent(inout) :: tm(ni,nk), qm(ni,nk), ps(ni), sigma(ni,nk)
      real, intent(inout) :: liqwcin(ni,nk), icewcin(ni,nk)
      real, intent(inout) :: liqwpin(ni,nk), icewpin(ni,nk)
      real, intent(inout) :: cloud(ni,nk), trav2d(ni,nk)
      !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include "phymkptr.hf"

      include "phyinput.inc"
      include "surface.cdk"
      include "nocld.cdk"

      real, dimension(ni,nkm1) ::  frac, lwcth, tcel, vtcel, vliqwcin, twcdens
      real rhoa

      integer :: i, k
      real    :: dp, lwcm1, iwcm1, zz, rec_grav, press
      logical :: nostrlwc, hascond_L

      real, pointer, contiguous :: znt(:)
      real, pointer, dimension(:,:), contiguous :: zftot, ziwc, zlwc, zqcplus, &
           zqiplus, zsnow, zqi_cat1, zqi_cat2, zqi_cat3, zqi_cat4
      !----------------------------------------------------------------

      MKPTR1D(znt, nt, fbus)

      MKPTR2D(zftot, ftot, fbus)
      MKPTR2D(ziwc, iwc, fbus)
      MKPTR2D(zlwc, lwc, fbus)
      MKPTR2D(zqcplus, qcplus, dbus)
      MKPTR2D(zqi_cat1, qti1plus, dbus)
      MKPTR2D(zqi_cat2, qti2plus, dbus)
      MKPTR2D(zqi_cat3, qti3plus, dbus)
      MKPTR2D(zqi_cat4, qti4plus, dbus)
      MKPTR2D(zqiplus, qiplus, dbus)
      MKPTR2D(zsnow, qnplus, dbus)

      call init2nan(frac, lwcth, tcel, vtcel, vliqwcin)

      rec_grav = 1./GRAV
      nostrlwc = (climat.or.stratos)
      hascond_L = (stcond /= 'NIL')

      if (kount == 0 .and. inilwc .and. hascond_L) then
         !     initialiser le champ d'eau nuageuse ainsi que la
         !     fraction nuageuse pour l'appel a la radiation a kount=0
         !     seulement.
         !     ces valeurs seront remplacees par celles calculees dans
         !     les modules de condensation.

         call cldwin(zftot, zlwc, tm, qm, ps, &
              trav2d, sigma, ni, nkm1, satuco)
      endif

      !     Diagnostic initialization of QC has been suppressed (because
      !     QC has been read as an input).  ftot still comes from cldwin
      !     if inilwc==.true., but lwc is replaced by the input value here.

      if ((any(dyninread_list_s == 'qc') .or. &
           any(phyinread_list_s(1:phyinread_n) == 'tr/qc:p')) .and. &
           .not. any(phyinread_list_s(1:phyinread_n) == 'lwc') ) then
         call cldwin(zftot, zlwc, tm, qm, ps, &
              trav2d, sigma, ni, nkm1, satuco)
         do k = 1,nkm1
            do i = 1,ni
               zlwc(i,k) = zqcplus(i,k)
            enddo
         enddo
      endif

      !     For maximum of lwc (when using newrad) or Liquid water content when
      !     istcond=1 Always execute this part of the code to allow calculation
      !     of NT in calcNT

      call liqwc(lwcth, sigma, tm, ps, ni, nkm1, ni, satuco)

      !     extracted from newrad3

      do k = 1,nkm1
         do i = 1,ni
            cloud(i,k) = zftot(i,k)
         enddo
      enddo

      !...  "no stratospheric lwc" mode when CLIMAT or STRATOS = true
      !...  no clouds above TOPC (see nocld.cdk)

      if (nostrlwc) then
         do k = 1,nkm1
            do i = 1,ni
               press = sigma(i,k)*ps(i)
               if (topc > press) then
                  cloud(i,k) = 0.0
                  zlwc (i,k) = 0.0
                  ziwc (i,k) = 0.0
               endif
            enddo
         enddo
      endif

      !     ************************************************************
      !     one branch for radia /= CCCMARAD and a simplified branch for radia=CCCMARAD
      !     -----------------------------------------------------------
      !
      !      If(radia /= 'CCCMARAD') Then
      !
      !    Always execute this part of the code to allow calculation of NT in calcNT


      DO_K: do k = 1,nkm1
         do i = 1,ni
            liqwcin(i,k) = max(zlwc(i,k),0.)
            icewcin(i,k) = 0.0

            if (hascond_L) then
               if ((liqwcin(i,k)+icewcin(i,k)) > 1.e-6) then
                  cloud(i,k) = max(cloud(i,k) ,0.01)
               else
                  cloud(i,k) = 0.0
               endif
            endif

            if (cloud(i,k) < 0.01) then
               liqwcin(i,k) = 0.
               icewcin(i,k) = 0.
            endif

            !     Min,Max of cloud

            cloud(i,k) = min(cloud(i,k),1.)
            cloud(i,k) = max(cloud(i,k),0.)

            if (hascond_L) then

               !     Normalize water contents to get in-cloud values

               zz = max(cloud(i,k),0.05)
               lwcm1 = liqwcin(i,k)/zz
               iwcm1 = icewcin(i,k)/zz

               !     Consider diabatic lifting limit when Sundquist scheme only

               liqwcin(i,k) = min(lwcm1,lwcth(i,k))
               icewcin(i,k) = min(iwcm1,lwcth(i,k))

            endif

            !       calculation of argument for call vsexp
            tcel(i,k) = tm(i,k)-TCDK
            vtcel(i,k) = -.003102*tcel(i,k)*tcel(i,k)

         end do
      end do DO_K

      !     liquid/solid water partition when not provided by
      !     microphysics scheme
      !     as in cldoptx4 of phy4.2 - after Rockel et al, Beitr. Atmos. Phys, 1991,
      !     p.10 (depends on T only [frac = .0059+.9941*Exp(-.003102 * tcel*tcel)]

      call gem_vsexp(frac, vtcel, nkm1*ni)
      do k = 1,nkm1
         do i = 1,ni
            if (tcel(i,k) >= 0.) then
               frac(i,k) = 1.0
            else
               frac(i,k) = .0059+.9941*frac(i,k)
            endif
            if (frac(i,k) < 0.01) frac(i,k) = 0.

            icewcin(i,k) = (1.-frac(i,k))*liqwcin(i,k)
            liqwcin(i,k) = frac(i,k)*liqwcin(i,k)
         enddo
      enddo

      !     calculate in-cloud liquid and ice water paths in each layer
      !     note: the calculation of the thickness of the layers done here is not
      !     coherent with what is done elsewhere for newrad (radir and sun) or
      !     cccmarad this code was extracted from cldoptx4 for phy4.4 dp(nk) is wrong

      do i = 1,ni
         dp = 0.5*(sigma(i,1)+sigma(i,2))
         dp = max(dp*ps(i),0.)
         icewpin(i,1) = icewcin(i,1)*dp*rec_grav*1000.
         liqwpin(i,1) = liqwcin(i,1)*dp*rec_grav*1000.

         dp = 0.5*(1.-sigma(i,nkm1))
         dp = max(dp*ps(i),0.)
         icewpin(i,nkm1) = icewcin(i,nkm1)*dp*rec_grav*1000.
         liqwpin(i,nkm1) = liqwcin(i,nkm1)*dp*rec_grav*1000.
      end do

      do k = 2,nkm1-1
         do i = 1,ni
            dp = 0.5*(sigma(i,k+1)-sigma(i,k-1))
            dp = max(dp*ps(i),0.)
            icewpin(i,k) = icewcin(i,k)*dp*rec_grav*1000.
            liqwpin(i,k) = liqwcin(i,k)*dp*rec_grav*1000.
         end do
      end do


      !... cccmarad simplified branch
      !      Else

      IF_CCCMARAD: if (radia(1:8) == 'CCCMARAD') then

         !    Begin - Calculation of NT - reproduction of NT obtained with newrad
         !    code (see cldoptx4)

         call calcnt(liqwpin, icewpin, cloud, znt, ni, nkm1, nk)
         call series_xst(znt, 'nt', trnch)

         !     End - Calculation of NT - reproduction of NT obtained with newrad code
         !     (see cldoptx4)
         !     impose coherent thresholds to cloud fraction and content

         do k = 1,nkm1
            do i = 1,ni
               cloud  (i,k) = min(cloud(i,k),1.)
               cloud  (i,k) = max(cloud(i,k),0.)
               liqwcin(i,k) = max(zlwc (i,k),0.)
               icewcin(i,k) = max(ziwc (i,k),0.)

               !            If ((liqwcin(i,k)+icewcin(i,k)) <= 1.e-6) Then
               !               cloud(i,k) = 0.0
               !            Endif

               if ((liqwcin(i,k)+icewcin(i,k)) > 1.e-6) then
                  cloud(i,k) = max(cloud(i,k) ,0.01)
               else
                  cloud(i,k) = 0.0
               endif


               if (cloud(i,k) < 0.01) then
                  liqwcin(i,k) = 0.
                  icewcin(i,k) = 0.
                  cloud(i,k) = 0.0
               endif

               !     Normalize water contents to get in-cloud values

               if (hascond_L) then
                  !               zz = Max(cloud(i,k),0.01)
                  zz = max(cloud(i,k),0.05)
                  liqwcin(i,k) = liqwcin(i,k)/zz
                  icewcin(i,k) = icewcin(i,k)/zz
               endif
            end do
         enddo


         !    calculate liquid/solid water partition when not provided by
         !    microphysics scheme
         !    ioptpart=1 : as for newrad - after Rockel et al, Beitr. Atmos. Phys, 1991,
         !    p.10 (depends on T only) [frac = .0059+.9941*Exp(-.003102 * tcel*tcel)]
         !    ioptpart=2 : after Boudala et al. (2004), QJRMS, 130, pp. 2919-2931.
         !    (depends on T and twc) [frac=twc^(0.141)*exp(0.037*(tcel))]
         !    PV; feb 2016: simplification of code, eliminate ioptpart=1; if you want it DIY

         if (.not.prep_cw_rad_fix_l) then
            do k = 1,nkm1
               do i = 1,ni
                  tcel(i,k) = tm(i,k)-TCDK
                  vtcel(i,k) = .037*tcel(i,k)
               enddo
            enddo
            call gem_vspown1(vliqwcin, liqwcin, 0.141, nkm1*ni) !ancien-bug-defaut
         else
            do k = 1,nkm1
               do i = 1,ni
                  tcel(i,k) = tm(i,k)-TCDK
                  vtcel(i,k) = .037*tcel(i,k)
                  press = sigma(i,k)*ps(i)
                  rhoa = press/(tm(i,k)*RGASD)
                  twcdens(i,k) = rhoa*liqwcin(i,k)*1000.
               enddo
            enddo
            call gem_vspown1(vliqwcin, twcdens, 0.141, nkm1*ni) !nouveau-debug
         endif

         call gem_vsexp(frac, vtcel, nkm1*ni)
         do k = 1,nkm1
            do i = 1,ni
               frac(i,k) = vliqwcin(i,k)*frac(i,k)
            enddo
         enddo
         do k = 1,nkm1
            do I = 1,ni
               if (tcel(i,k) >= 0.) then
                  frac(i,k) = 1.0
               elseif (tcel(i,k) < -38.) then
                  frac(i,k) = 0.0
               endif
               if (frac(i,k) < 0.01) frac(i,k) = 0.

               icewcin(i,k) = (1.-frac(i,k))*liqwcin(i,k)
               liqwcin(i,k) = frac(i,k)*liqwcin(i,k)
            enddo
         enddo

         !    calculate in-cloud liquid and ice water paths in each layer
         !    note: the calculation of the thickness of the layers done here is
         !    not coherent with what is done elsewhere for newrad (radir and sun)
         !    or cccmarad this code was extracted from cldoptx4 for phy4.4
         !    dp(nk) is wrong

         do i = 1,ni
            dp = 0.5*(sigma(i,1)+sigma(i,2))
            dp = dp*ps(i)
            icewpin(i,1) = icewcin(i,1)*dp*rec_grav*1000.
            liqwpin(i,1) = liqwcin(i,1)*dp*rec_grav*1000.
            dp = 0.5*(1.-sigma(i,nkm1))
            dp = dp*ps(i)
            icewpin(i,nkm1) = icewcin(i,nkm1)*dp*rec_grav*1000.
            liqwpin(i,nkm1) = liqwcin(i,nkm1)*dp*rec_grav*1000.
         end do

         do k = 2,nkm1-1
            do i = 1,ni
               dp = 0.5*(sigma(i,k+1)-sigma(i,k-1))
               dp = dp*ps(i)
               icewpin(i,k) = icewcin(i,k)*dp*rec_grav*1000.
               liqwpin(i,k) = liqwcin(i,k)*dp*rec_grav*1000.
            end do
         end do

         !     to impose coherence between cloud fraction and thresholds on liqwpin
         !     and icewpin in calculations of cloud optical properties (cldoppro)
         !     in cccmarad

         do k = 1,nkm1
            do i = 1,ni
               if (liqwpin(i,k) <= 0.001 .and. icewpin(i,k) <= 0.001) &
                    cloud(i,k) = 0.0
            end do
         end do

      endif IF_CCCMARAD

      !     to simulate a clear sky radiative transfer, de-comment following lines
      !        do k = 1,nkm1
      !        do i = 1,ni
      !              liqwcin(i,k)   = 0.0
      !              icewcin(i,k)   = 0.0
      !              liqwpin(i,k)   = 0.0
      !              icewpin(i,k)   = 0.0
      !              cloud(i,k)     = 0.0
      !        end do
      !        end do

      return
   end subroutine prep_cw_rad3

end module prep_cw_rad
