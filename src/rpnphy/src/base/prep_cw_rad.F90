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
        kount, ni, nk, nkm1)
      use, intrinsic :: iso_fortran_env, only: INT64
      use debug_mod, only: init2nan
      use tdpack_const, only: GRAV, TCDK, RGASD, TRPL
      use phy_options
      use phybus
      implicit none

      !@Object  Prepare liquid/ice water contents and cloudiness for the radiation
      !@Arguments
      !     - input -
      !     tm       temperature
      !     qm       specific humidity
      !     ps       surface pressure
      !     sigma    sigma levels
      !     kount    index of timestep
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

      integer, intent(in) :: kount, ni, nk, nkm1
      real, pointer, contiguous :: dbus(:), fbus(:)
      real, intent(inout) :: tm(ni,nk), qm(ni,nk), ps(ni), sigma(ni,nk)
      real, intent(inout) :: liqwcin(ni,nk), icewcin(ni,nk)
      real, intent(inout) :: liqwpin(ni,nk), icewpin(ni,nk)
      real, intent(inout) :: cloud(ni,nk)
      !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include "phymkptr.hf"

      include "phyinput.inc"
      include "surface.cdk"
      include "nocld.cdk"

      real, dimension(ni,nkm1) ::  dp
      real, dimension(ni,nkm1) ::  tcel, frac
      real rhoa

      integer :: i, k
      real    :: rec_grav, press, vliqwcin, twcdens
      logical :: nostrlwc, hascond_L

      real, pointer, dimension(:,:), contiguous :: zftot, ziwc, zlwc, zqcplus

      !----------------------------------------------------------------

      MKPTR2D(zftot, ftot, fbus)
      MKPTR2D(ziwc, iwc, fbus)
      MKPTR2D(zlwc, lwc, fbus)
      MKPTR2D(zqcplus, qcplus, dbus)

      call init2nan(dp, tcel, frac)

      rec_grav = 1./GRAV
      nostrlwc = (climat.or.stratos)
      hascond_L = (stcond /= 'NIL')

      if (kount == 0 .and. inilwc .and. hascond_L) then
         !     initialiser le champ d'eau nuageuse ainsi que la
         !     fraction nuageuse pour l'appel a la radiation a kount=0
         !     seulement.
         !     ces valeurs seront remplacees par celles calculees dans
         !     les modules de condensation.

         call cldwin1(zftot, zlwc, tm, qm, ps, &
              sigma, ni, nkm1, satuco)
      endif

      !     Diagnostic initialization of QC has been suppressed (because
      !     QC has been read as an input).  ftot still comes from cldwin
      !     if inilwc==.true., but lwc is replaced by the input value here.

      if ((any(dyninread_list_s == 'qc') .or. &
           any(phyinread_list_s(1:phyinread_n) == 'tr/qc:p')) .and. &
           .not. any(phyinread_list_s(1:phyinread_n) == 'lwc') ) then
         call cldwin1(zftot, zlwc, tm, qm, ps, &
              sigma, ni, nkm1, satuco)
         if (associated(zqcplus)) then
            do k = 1,nkm1
               do i = 1,ni
                  zlwc(i,k) = zqcplus(i,k)
               enddo
            enddo
         else
            zlwc(:, 1:nkm1) = 0.
         endif
      endif

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

      do i = 1,ni
         dp(i,1) = 0.5*(sigma(i,1)+sigma(i,2))
         dp(i,1) = max(dp(i,1)*ps(i),0.)
         dp(i,nkm1) = 1. - 0.5*(sigma(i,nkm1) + sigma(i,nkm1-1))
         dp(i,nkm1) = max(dp(i,nkm1)*ps(i),0.)
      end do

      do k = 2,nkm1-1
         do i = 1,ni
            dp(i,k) = 0.5*(sigma(i,k+1)-sigma(i,k-1))
            dp(i,k) = max(dp(i,k)*ps(i),0.)
         end do
      end do

      !     impose same threshold on cloud fraction  has used in RT
      !     Normalize water contents to get in-cloud values

      do k = 1,nkm1
         do i = 1,ni
            cloud  (i,k) = min(cloud(i,k),1.)
            cloud  (i,k) = max(cloud(i,k),0.)
            liqwcin(i,k) = max(zlwc (i,k),0.)
            icewcin(i,k) = max(ziwc (i,k),0.)
            if (cloud(i,k) <= cldfth) then
               liqwcin(i,k) = 0.
               icewcin(i,k) = 0.
               cloud(i,k) = 0.0
            else
               liqwcin(i,k) = liqwcin(i,k)/cloud(i,k)
               icewcin(i,k) = icewcin(i,k)/cloud(i,k)
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

      do k = 1,nkm1
         do i = 1,ni
            tcel(i,k) = tm(i,k)-TCDK
         enddo
      enddo
      select case (rad_part_nomp)
      case ('BOUOPS')
         do k = 1,nkm1
            do i = 1,ni
               vliqwcin = liqwcin(i,k)**0.141  !# ancien-bug-default
               frac(i,k) = vliqwcin*exp(.037*tcel(i,k))
               if (tcel(i,k) >= 0.) then
                  frac(i,k) = 1.0
               elseif (tcel(i,k) < -38.) then
                  frac(i,k) = 0.0
               endif
               if (frac(i,k) < 0.01) frac(i,k) = 0.
            enddo
         enddo
      case ('BOUDALA')
         do k = 1,nkm1
            do i = 1,ni
               press = sigma(i,k)*ps(i)
               rhoa = press/(tm(i,k)*RGASD)
               twcdens = rhoa*liqwcin(i,k)*1000.
               vliqwcin = twcdens**0.141  !# nouveau-debug
               frac(i,k) = vliqwcin*exp(.037*tcel(i,k))
               if (tcel(i,k) >= 0.) then
                  frac(i,k) = 1.0
               elseif (tcel(i,k) < -38.) then
                  frac(i,k) = 0.0
               endif
               if (frac(i,k) < 0.01) frac(i,k) = 0.
            enddo
         enddo
      case ('ECMWF')
         do k = 1,nkm1
            do i = 1,ni
               if (tm(i,k) >= TRPL) then
                  frac(i,k) = 1.0
               elseif (tm(i,k) <= 250.16) then
                  frac(i,k) = 0.0
               else
                  frac(i,k) = ((tm(i,k)-250.16)/(TRPL-250.16))**2
               endif
            enddo
         enddo
      case ('ROCKEL')
         do k = 1,nkm1
            do i = 1,ni
               if (tcel(i,k) >= 0.0) then
                  frac(i,k) = 1.0
               else
                  frac(i,k) = .0059+.9941*Exp(-.003102 * tcel(i,k)**2)
               endif
               if (frac(i,k) < 0.01) frac(i,k) = 0.
            enddo
         enddo
      end select

      do k = 1,nkm1
         do i = 1,ni
            icewcin(i,k) = (1.-frac(i,k))*liqwcin(i,k)
            liqwcin(i,k) = frac(i,k)*liqwcin(i,k)
         enddo
      enddo


      !    calculate in-cloud liquid and ice water paths in each layer

      do k = 1,nkm1
         do i = 1,ni
            icewpin(i,k) = icewcin(i,k)*dp(i,k)*rec_grav*1000.
            liqwpin(i,k) = liqwcin(i,k)*dp(i,k)*rec_grav*1000.
         end do
      end do

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
