!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------

subroutine tkealg2(estar,en,zn,ze,dvdz2,buoy_flux,diss,shr,buoy,tau,n,nk)

   implicit none
!!!#include <arch_specific.hf>

   ! Argument declaration
   integer, intent(in) :: n                          !horizontal dimension
   integer, intent(in) :: nk                         !vertical dimension
   real, intent(in) :: tau                           !time step (s)
   real, dimension(n,nk), intent(in) :: en           !turbulent kinetic energy (m2/s2)
   real, dimension(n,nk), intent(in) :: zn           !mixing length scale (m)
   real, dimension(n,nk), intent(in) :: ze           !dissipation length scale (m)
   real, dimension(n,nk), intent(in) :: dvdz2        !square of vertical shear (/s2)
   real, dimension(n,nk), intent(in) :: buoy_flux    !buoyancy flux (/s2)
   real, dimension(n,nk), intent(out) :: estar       !updated turbulent kinetic energy (time *; m2/s2)
   real, dimension(n,nk), intent(out) :: diss        !viscous dissipation term (m2/s3)
   real, dimension(n,nk), intent(out) :: shr         !shear generation term (m2/s3)
   real, dimension(n,nk), intent(out) :: buoy        !buoyancy generation/suppression term (m2/s3)

   !Author
   !          J. Mailhot (August 1999)
   !Revision
   ! 001      A.-M. Leduc (Oct 2001) Automatic arrays
   ! 002      B. Bilodeau (Mar 2003) Include machcon.cdk
   ! 003      A. Plante   (May 2003) IBM conversion
   !             - calls to vsqrt routine (from massvp4 library)
   !             - calls to vstan routine (from massvp4 library)
   !             - calls to vslog routine (from massvp4 library)
   ! 004      a. Plante   (juillet 2003) correction of bugs in IBM conversion
   !             - optimization of tan removed
   !             - optimization of log removed
   !
   !Object
   !          Solve the algebraic part of the TKE equation
   !
   !Arguments
   !
   !
   !IMPLICITS
#include "machcon.cdk"
#include "clefcon.cdk"

   ! Local variable declaration
   integer :: j,k
   real :: clamda,sqrt_tke
   real, dimension(n,nk) :: b,c,b_over_c
   real(kind=8) :: r8tmp

   ! Compute coefficients (dble precision) and solve the algebraic TKE equation
   do k=1,nk-1
      do j=1,n
         b(j,k) = BLCONST_CK*zn(j,k)*(dvdz2(j,k) - buoy_flux(j,k))
         c(j,k) = BLCONST_CE/ze(j,k)
         r8tmp = (abs(b(j,k)/c(j,k)))
         b_over_c(j,k) = sqrt(r8tmp)
         r8tmp = en(j,k)
         sqrt_tke = sqrt(r8tmp)
         if( abs(b(j,k)) < epsilon(b) ) then
            estar(j,k) = sqrt_tke / (1.+0.5*sqrt_tke*c(j,k)*tau)
         elseif( b(j,k) > 0. ) then
            estar(j,k) = b_over_c(j,k)*( -1.0+2.0/(1.0+exp(-b_over_c(j,k)*c(j,k)*tau) * &
                 (-1.0+2.0/(1.0+sqrt_tke/b_over_c(j,k)))) )
         else
            estar(j,k)=b_over_c(j,k)*tan(min(tanlim,max(atan(sqrt_tke/b_over_c(j,k)) - &
                 0.5*b_over_c(j,k)*c(j,k)*tau , 0.0 ) ) )
         endif
         estar(j,k) = max(ETRMIN,estar(j,k)**2)
      enddo
   enddo

   ! Compute TKE equation terms for output purposes (time series)
   do k=1,nk-1
      do j=1,n
         b_over_c(j,k) = b(j,k)/c(j,k)
         b(j,k)=0.0
         if( en(j,k) /= b_over_c(j,k) .and. estar(j,k) /= b_over_c(j,k) ) then
            b(j,k)=alog( abs( (estar(j,k)-b_over_c(j,k)) / (en(j,k)-b_over_c(j,k)) ) )
         endif
         clamda = blconst_ck*zn(j,k)*b(j,k)/(c(j,k)*tau)
         shr(j,k) = -clamda * dvdz2(j,k)
         buoy(j,k) = clamda * buoy_flux(j,k)
         diss(j,k) = (estar(j,k)-en(j,k))/tau-(shr(j,k)+buoy(j,k))
      enddo
   enddo

   ! Set surface values
   estar(:,nk) = 0.
   shr(:,nk) = 0.
   buoy(:,nk) = 0.
   diss(:,nk) = 0.

   ! End of subprogram
   return
end subroutine tkealg2
