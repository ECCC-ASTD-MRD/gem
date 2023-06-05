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
module vintage_nt
   implicit none
   private
   public :: vintage_nt1

contains
   
   !/@*
   subroutine vintage_nt1(fbus,  &
        tm, ps, sigma,  &
        cloud,    &
        trnch, ni, nk, nkm1)
      use, intrinsic :: iso_fortran_env, only: INT64
      use debug_mod, only: init2nan
      use tdpack_const, only: GRAV, TCDK, RGASD
      use phy_options
      use phybus
      use series_mod, only: series_xst
      implicit none

      !@Object  Calculate 2d variable  cloud cover NT - reproduce old newrad style effective cloud cover
      !@Arguments
      !     - input -
      !     tm       temperature
      !     ps       surface pressure
      !     sigma    sigma levels
      !     cloud    cloud fraction
      !     trnch    number of the slice
      !     ni       horizontal dimension
      !     nk       number of layers, including the diag level
      !     nkm1     number of layers, excluding the diag level
      !     - output -
      !      znt - vintage effective total cloud cover
      !
      !      zlwc - total water content (consun) arrives from zlwc=qcplus (in prep_cw)


      integer, intent(in) :: trnch, ni, nk, nkm1
      real, pointer, contiguous :: fbus(:)
      real, intent(in) :: tm(ni,nk), ps(ni), sigma(ni,nk)
      real, intent(inout) :: cloud(ni,nkm1)
      !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include "phymkptr.hf"

      include "phyinput.inc"

      real, dimension(ni,nkm1) :: lwcth, vliqwcin
      real, dimension(ni,nk) :: liqwpin(ni,nk), icewpin(ni,nk)
      real, dimension(ni,nk) :: liqwcin(ni,nk), icewcin(ni,nk)

      integer :: i, k, kind
      real    :: lwcm1, zz, rec_grav, frac, tcel
      logical :: hascond_L
      real, dimension(ni,nkm1  ) :: emw,emi,eneb,dp
      real, dimension(ni) :: trmin, tmem
      real, dimension(ni,nk) :: ff

      real :: emiss,xnu
      ! diffusivity factor of Elsasser
      real, parameter :: ELSA = 1.66
      real, parameter :: REI = 15.
      real, parameter :: REC_REI = 1. / REI
      real, parameter :: KI = .0003 + 1.290 * REC_REI

      real, pointer, contiguous :: znt(:)
      real, pointer, dimension(:,:), contiguous :: zlwc

      !----------------------------------------------------------------

      MKPTR1D(znt, nt, fbus)
      MKPTR2D(zlwc, lwc, fbus)

      call init2nan(lwcth, vliqwcin)

      rec_grav = 1./GRAV
      hascond_L = (stcond /= 'NIL')

      icewcin = 0.0

      ! bug below: should be dp(i,nkm1) = 1. - 0.5*(sigma(i,nkm1) + sigma(i,nkm1-1))
      do i = 1,ni
         dp(i,1) = 0.5*(sigma(i,1)+sigma(i,2))
         dp(i,1) = max(dp(i,1)*ps(i),0.)
         dp(i,nkm1) = 0.5*(1.-sigma(i,nkm1))
         dp(i,nkm1) = max(dp(i,nkm1)*ps(i),0.)
      end do

      do k = 2,nkm1-1
         do i = 1,ni
            dp(i,k) = 0.5*(sigma(i,k+1)-sigma(i,k-1))
            dp(i,k) = max(dp(i,k)*ps(i),0.)
         end do
      end do

      call liqwc(lwcth, sigma, tm, ps, ni, nkm1, ni, satuco)    
      DO_K: do k = 1,nkm1
         do i = 1,ni
            liqwcin(i,k) = max(zlwc(i,k),0.)

            if (hascond_L) then
               if ((liqwcin(i,k)) > 1.e-6) then
                  cloud(i,k) = max(cloud(i,k) ,0.01)
               else
                  cloud(i,k) = 0.0
               endif
            endif

            if (cloud(i,k) < 0.01) then
               liqwcin(i,k) = 0.
            endif

            cloud(i,k) = min(cloud(i,k),1.)
            cloud(i,k) = max(cloud(i,k),0.)

            if (hascond_L) then

               !     Normalize water contents to get in-cloud values

               zz = max(cloud(i,k),0.05)
               lwcm1 = liqwcin(i,k)/zz

               !     Impose diabatic lifting limit when Sundquist scheme only

               liqwcin(i,k) = min(lwcm1,lwcth(i,k))

            endif

            !       calculation of argument for call vsexp
            
            !     liquid/solid water partition 
            !     Rockel et al, Beitr. Atmos. Phys, 1991,
            !     p.10 (depends on T only [frac = .0059+.9941*Exp(-.003102 * tcel*tcel)]            
            
            tcel = tm(i,k)-TCDK
            frac = exp(-.003102*tcel*tcel)
            if (tcel >= 0.) then
               frac = 1.0
            else
               frac = .0059+.9941*frac
            endif
            if (frac < 0.01) frac = 0.

            icewcin(i,k) = (1.-frac)*liqwcin(i,k)
            liqwcin(i,k) = frac*liqwcin(i,k)
         enddo
      enddo DO_K

      !     calculate in-cloud liquid and ice water paths in each layer

      do k = 1,nkm1
         do i = 1,ni
            icewpin(i,k) = icewcin(i,k)*dp(i,k)*rec_grav*1000.
            liqwpin(i,k) = liqwcin(i,k)*dp(i,k)*rec_grav*1000.
         end do
      end do

      do k=1,nkm1
         do I=1,ni
            EMW(i,k) = 1. - exp(-0.087 * ELSA * liqwpin(i,k))
            EMI(i,k) = 1. - exp(-ELSA * KI * icewpin(i,k))
            EMISS   = 1. - (1.-EMI(i,k)) * (1.-EMW(i,k))
            ENEB(i,k) = cloud(i,k)*emiss
         end do
      end do

      !...  maximum random overlap
      do i=1,ni
         ff(i,1) = 1.
         tmem(i) = 1.
         trmin(i) = 1.
      enddo
      do k=2,nk
         kind = k-2
         kind = max0(kind,1)
         do i=1,ni
            xnu = 1.-eneb(i,k-1)
            if (cloud(i,kind) < 0.01) then
               tmem(i) = ff(i,k-1)
               trmin(i) = xnu
            else
               trmin(i) = min(trmin(i), xnu)
            endif
            ff(i,k) = tmem(i) * trmin(i)
         enddo
      enddo

      do  i=1,ni
         znt(i) = 1.-ff(i,nk)
      enddo


      call series_xst(znt, 'nt', trnch)

      return
   end subroutine vintage_nt1

end module vintage_nt
