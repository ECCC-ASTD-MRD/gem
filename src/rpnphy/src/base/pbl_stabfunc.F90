!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module pbl_stabfunc
   use phy_status, only: PHY_OK, PHY_ERROR
   ! Calculation and manipulation of stability functions in the PBL
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   include "surface.cdk"
   private

   ! Internal parameters
   integer, parameter :: CONVERGE_ERR=2               !Return value for a convergence error
   real, parameter :: ASX = 5.                        !PBL stability function coefficient (Delage 1997)
   real, parameter :: F_MIN = 1e-6                    !Minimum absolute value for stability functions

   ! API subprograms
   public :: psf_stabfunc

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function psf_stabfunc(rig,z,fm,fh,blend_bottom,blend_top) result(stat)
      ! PBL stability function matching the surface layer
      use phy_options, only: pbl_func_stab,pbl_func_unstab

      ! Argument declaration
      real, dimension(:,:), intent(in) :: rig            !gradient Richardson number
      real, dimension(:,:), intent(in) :: z              !height above ground (m)
      real, dimension(:,:), intent(out) :: fm            !stability function for momentum
      real, dimension(:,:), intent(out) :: fh            !stability function for heat
      real, intent(in), optional :: blend_bottom         !height to begin blending of PBL functions [40m]
      real, intent(in), optional :: blend_top            !height to end blending of SL functions [200m]
      integer :: stat                                    !return status

      ! Common block
!!!#include <arch_specific.hf>
      include "clefcon.cdk"

      ! Local variables
      integer :: istat
      real :: myBlend_bottom,myBlend_top
      real, dimension(size(rig,dim=1),size(rig,dim=2)) :: fm_upper,fh_upper,fm_sfc,fh_sfc,w
      logical, dimension(size(rig,dim=1),size(rig,dim=2)) :: upper_layer,near_sfc,stable

      ! Set return status
      stat = PHY_ERROR

      ! Set default values
      myBlend_bottom = 40.
      if (present(blend_bottom)) myBlend_bottom = blend_bottom
      myBlend_top = 200.
      if (present(blend_top)) myBlend_top = blend_top

      ! Initializations
      stable = .true.
      where (rig < 0.) stable = .false.
      upper_layer = (z >= myBlend_bottom)
      near_sfc = (z <= myBlend_top)
      fm_upper = 1.; fh_upper = 1.
      fm_sfc = 1.; fh_sfc = 1.

      ! Stable case for the upper layers
      istat = PHY_OK
      UPPER_LAYER_STABLE: select case (pbl_func_stab)
      case ('LOUIS79')
         where (upper_layer .and. stable)
            fm_upper = 1. + 5.*rig
            fh_upper = beta*fm_upper
         endwhere
      case ('DERBY97')
         where (upper_layer .and. stable)
            w = 0.5*( 1. - sign(1., rig - 0.1))
            w = max(0., min(1.,w))
            fm_upper = w*(1./(1.-5.*rig)) + (1.-w)*(20.*rig)
            fh_upper = beta*fm_upper
         endwhere
      case ('BELJAARS99')
         where (upper_layer .and. stable)
            fm_upper = sqrt(1.+ 10.*rig/sqrt(1.+rig))
            fh_upper = beta*(1.+ 10.*rig*sqrt(1.+rig))/fm_upper
         endwhere
      case ('BELJAARS91')
         istat = sf_calc(rig,upper_layer.and.stable,bh91_fg, &
              bh91_fm,bh91_dfm,bh91_fh,bh91_dfh,fm_upper,fh_upper)
      case ('LOCK07')
         istat = sf_calc(rig,upper_layer.and.stable,l07_fg, &
              l07_fm,l07_dfm,l07_fh,l07_dfh,fm_upper,fh_upper)
      case ('DELAGE97')
         where (upper_layer .and. stable)
            fm_upper = min(1.+d97_as*rig,1./max(epsilon(rig),1.-ASX*rig))
            fh_upper = beta*fm_upper
         endwhere
      case DEFAULT
         istat = PHY_ERROR
      end select UPPER_LAYER_STABLE
      if (istat /= PHY_OK) then
         if (istat == CONVERGE_ERR) then
            call msg_toall(MSG_WARNING, '(psf_stabfunc) PBL convergence failure for '// &
                 trim(pbl_func_stab))
         else
            call msg_toall(MSG_ERROR, '(psf_stabfunc) Error computing stable PBL functions for '// &
                 trim(pbl_func_stab))
            stat = istat
            return
         endif
      endif

      ! Unstable case for the upper layers
      istat = PHY_OK
      UPPER_LAYER_UNSTABLE: select case (pbl_func_unstab)
      case ('DYER74')
         where (upper_layer .and. .not.stable)
            fm_upper = (1.-16.*rig)**(-0.25)
            fh_upper = beta*fm_upper**2
         endwhere
      case ('DELAGE92')
         where (upper_layer .and. .not.stable)
            fm_upper = (1.-dg92_ci*rig)**(-1./6.)
            fh_upper = beta*fm_upper**2
         endwhere
      case DEFAULT
         istat = PHY_ERROR
      end select UPPER_LAYER_UNSTABLE
      if (istat /= PHY_OK) then
         call msg_toall(MSG_ERROR, '(psf_stabfunc) Error computing unstable PBL functions for '// &
              trim(pbl_func_unstab))
         stat = istat
         return
      endif

      ! Stable case for the near-surface
      istat = PHY_OK
      NEAR_SFC_STABLE: select case (sl_func_stab)
      case ('BELJAARS91')
         istat = sf_calc(rig,near_sfc.and.stable,bh91_fg, &
              bh91_fm,bh91_dfm,bh91_fh,bh91_dfh,fm_sfc,fh_sfc)
      case ('LOCK07')
         istat = sf_calc(rig,near_sfc.and.stable,l07_fg, &
              l07_fm,l07_dfm,l07_fh,l07_dfh,fm_sfc,fh_sfc)
      case ('DELAGE97')
         where (near_sfc .and. stable)
            fm_sfc = min(1.+d97_as*rig,1./max(epsilon(rig),1.-ASX*rig))
            fh_sfc = beta*fm_sfc
         endwhere
      case DEFAULT
         istat = PHY_ERROR
      end select NEAR_SFC_STABLE
      if (istat /= PHY_OK) then
         if (istat == CONVERGE_ERR) then
            call msg_toall(MSG_WARNING, '(psf_stabfunc) SL convergence failure for '// &
                 trim(sl_func_stab))
         else
            call msg_toall(MSG_ERROR, '(psf_stabfunc) Error computing stable SL functions for '// &
                 trim(sl_func_stab))
            stat = istat
            return
         endif
      endif

      ! Unstable case for the near-surface
      istat = PHY_OK
      NEAR_SFC_UNSTABLE: select case (sl_func_unstab)
      case ('DYER74')
         where (near_sfc .and. .not.stable)
            fm_sfc = (1.-16.*rig)**(-0.25)
            fh_sfc = beta*fm_sfc**2
         endwhere
      case ('DELAGE92')
         where (near_sfc .and. .not.stable)
            fm_sfc = (1.-dg92_ci*rig)**(-1./6.)
            fh_sfc = beta*fm_sfc**2
         endwhere
      case DEFAULT
         istat = PHY_ERROR
      end select NEAR_SFC_UNSTABLE
      if (istat /= PHY_OK) then
         call msg_toall(MSG_ERROR, '(psf_stabfunc) Error computing unstable SL functions for '// &
              trim(sl_func_unstab))
         stat = istat
         return
      endif

      ! Blend linearly between SL (s) and PBL (u) functions
      w = max(0., min(1., (z - myBlend_bottom)/max(myBlend_top - myBlend_bottom,1.)))
      fm = w*fm_upper + (1.-w)*fm_sfc
      fh = w*fh_upper + (1.-w)*fh_sfc

      ! Make sure that the absolute value of the stability functions are reasonable
      fm = sign(max(abs(fm),F_MIN),fm)
      fh = sign(max(abs(fh),F_MIN),fh)

      ! Successful end of subprogram
      stat = PHY_OK
      return
   end function psf_stabfunc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function sf_calc(rig,mask,fg,sfm,dsfm,sfh,dsfh,fm,fh) result(stat)
      ! Compute stability function values by Newton-Raphson iterations

      ! Argument declaration
      real, dimension(:,:), intent(in) :: rig     !gradient Richardson number
      logical, dimension(:,:), intent(in) :: mask !calculation mask
      real, external :: fg                        !first-guess stability function (callback)
      real, external :: sfm                       !stability function for momentum (callback)
      real, external :: dsfm                      !d(sfm)/d(RI) for momentum (callback)
      real, external :: sfh                       !stability function for heat (callback)
      real, external :: dsfh                      !d(sfh)/d(RI) for heat (callback)
      real, dimension(:,:), intent(out) :: fm     !stability function values for momentum
      real, dimension(:,:), intent(out) :: fh     !stability function values for heat
      integer :: stat                             !return status

      ! Local parameters
      real, parameter :: EPSLN=1e-5               !minimum slope for derivative
      integer, parameter :: ITERMAX=50            !maximum number of interations for solution
      real, parameter :: TOLERANCE=0.001          !tolerance for iterations


      ! Local variables and common blocks
      integer :: i,n,k,nk,iter
      real :: x,g,dg

      ! Set return status
      stat = PHY_ERROR

      ! Initialization
      n  = size(rig,dim=1)
      nk = size(rig,dim=2)
      fm  = 1.
      fh  = 1.

      do k=1,nk
         do i=1,n
            if (.not.mask(i,k)) cycle

            ! First guess
            x = fg(rig(i,k))
            fm(i,k)  = sfm(x)
            fh(i,k)  = sfh(x)
            g = rig(i,k) - x*fh(i,k)/(fm(i,k)**2)

            ! Newton-Raphson iterations
            iter = 1
            do while (abs(g) > TOLERANCE .and. iter <= ITERMAX)
               dg  = - fh(i,k)/(fm(i,k)**2)*(1. + x*dsfh(x)/fh(i,k) &
                    - 2.*x*dsfm(x)/fm(i,k))
               dg  = sign(max(abs(dg),EPSLN),dg)
               x   = x - g/dg
               fm(i,k)  = sfm(x)
               fh(i,k)  = sfh(x)
               g = rig(i,k) - x*fh(i,k)/(fm(i,k)**2)
               iter = iter+1
            enddo
            if (iter > ITERMAX) stat = CONVERGE_ERR

         enddo
      enddo

      ! End of subprogram
      if (.not.stat==CONVERGE_ERR) stat = PHY_OK
      return

   end function sf_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Stable-case stability functions proposed by Beljaars and Holtslag (1991)

   real function bh91_fg(rig)
      ! First guess calculation (B-J 1991)
      real, intent(in) :: rig
      if ( rig <= 0.1 ) then
         bh91_fg = rig
      else
         bh91_fg = (3.*bh91_a/2.)*rig**2
      endif
   end function bh91_fg

   real function bh91_fm(v)
      ! Stability function for momentum (B-J 1991)
      real, intent(in) :: v
      bh91_fm = 1. + v*( bh91_a + bh91_b*(1. + bh91_c - bh91_d*v)*exp(-bh91_d*v) )
   end function bh91_fm

   real function bh91_dfm(v)
      ! Derivative of stability function for momentum wrt Ri (B-J 1991)
      real, intent(in) :: v
      bh91_dfm = bh91_a + &
           bh91_b*(1.+bh91_c -(3+bh91_c)*bh91_d*v + bh91_d**2*v**2)*exp(-bh91_d*v)
   end function bh91_dfm

   real function bh91_fh(v)
      ! Stability function for heat (B-J 1991)
      real, intent(in) :: v
      bh91_fh  = beta * ( 1. + &
           v*( bh91_a*sqrt(1. + 2./3.*bh91_a*v) + bh91_b*(1. + bh91_c - bh91_d*v)*exp(-bh91_d*v) ) )
   end function bh91_fh

   real function bh91_dfh(v)
      ! Derivative of stability function for heat wrt Ri (B-J 1991)
      real, intent(in) :: v
      bh91_dfh = beta * ( bh91_a*(1+bh91_a*v)/sqrt(1. + 2./3.*bh91_a*v) + &
           bh91_b*(1.+bh91_c -(3+bh91_c)*bh91_d*v + bh91_d**2*v**2)*exp(-bh91_d*v) )
   end function bh91_dfh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Stable-case Stability function by Lock (2007)

   real function l07_fg(rig)
      ! First guess calculation (Lock 2007)
      real, intent(in) :: rig
      if ( rig <= 0.1 ) then
         l07_fg = rig
      else
         l07_fg = 2.*rig
      endif
   end function l07_fg

   real function l07_fm(v)
      ! Stability function for momentum (Lock 2007)
      real, intent(in) :: v
      l07_fm  = 1. + l07_am*v
   end function l07_fm

   real function l07_dfm(v)
      ! Derivative of stability function for momentum wrt Ri (Lock 2007)
      real, intent(in) :: v  !Not used
      if (zua == -99999.9999) print *,'l07_dfm',v  !Avoid compiler complain of not used
      l07_dfm = l07_am
   end function l07_dfm

   real function l07_fh(v)
      ! Stability function for heat (Lock 2007)
      real, intent(in) :: v
      l07_fh = beta * ( 1. + l07_ah*v + 0.5*l07_ah**2*v**2 )
   end function l07_fh

   real function l07_dfh(v)
      ! Derivative of stability function for heat wrt Ri (Lock 2007)
      real, intent(in) :: v
      l07_dfh = beta * ( l07_ah + l07_ah**2*v )
   end function l07_dfh

end module pbl_stabfunc


