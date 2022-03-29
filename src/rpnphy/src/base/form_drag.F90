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

module form_drag
   implicit none
   private
   public :: form_drag1

contains

   !/@*
   function form_drag1(utend,vtend,uwind,vwind,gzmom,sigma_s,z0,tau,n,nk) result(status)
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      !@Arguments
      integer, intent(in) :: n                          !horizontal dimension
      integer, intent(in) :: nk                         !vertical dimension
      real, intent(in) :: tau                           !time step length (s)
      real, dimension(n), intent(in) :: sigma_s         !standard deviation of subgrid orography at small scales (m)
      real, dimension(n), intent(in) :: z0              !momentum roughness length (m)
      real, dimension(n,nk), intent(in) :: uwind       !x-component wind at time-minus (m/s)
      real, dimension(n,nk), intent(in) :: vwind       !y-component wind at time-minus (m/s)
      real, dimension(n,nk), intent(in) :: gzmom        !momentum level heights (m)
      real, dimension(n,nk), intent(out) :: utend       !x-component TOFD tedency (m/s^2)
      real, dimension(n,nk), intent(out) :: vtend       !x-component TOFD tedency (m/s^2)
      integer :: status                                 !completion status (RMN_OK or RMN_ERR)

      !@Author R. McTaggart-Cowan and A. Zadra (Fall 2016)

      !@Object
      !          Calculate the turbulent orographic form drag following the
      !          approach proposed by Beljaars, Brown and Wood (2004; QJRMS).
      !*@/


      ! Local parameter definitions
      real, parameter :: N1=-1.9,N2=-2.8,K1=3e-3,KFLT=3.5e-4,IH=1.02e-3, &
           CM=0.1,BETA=1.,CMD=5e-3,CCORR=0.6

      ! Local variable declarations
      integer :: i,k
      real :: fac,wspd
      real, dimension(n) :: a2

      ! Initialize status
      status = RMN_ERR

      ! Set constants
      a2 = (sigma_s)**2 / (IH*KFLT**N1) * K1**(N1-N2)

      ! Implicit calculation of TOFD wind tendencies (Eq. 16 of Beljaars et al. (2004))
      do k=1,nk-1
         do i=1,n
            fac = tofd_alpha*BETA*CMD*CCORR*2.019*exp(-((gzmom(i,k)+z0(i))/1500.)**1.5)*a2(i)*(gzmom(i,k)+z0(i))**(-1.2)
            wspd = sqrt(uwind(i,k)*uwind(i,k)+vwind(i,k)*vwind(i,k))
            utend(i,k) = -fac*wspd*uwind(i,k) / (1.+fac*wspd*tau)
            vtend(i,k) = -fac*wspd*vwind(i,k) / (1.+fac*wspd*tau)
         enddo
      enddo
      utend(:,nk) = 0.
      vtend(:,nk) = 0.

      ! Successful completion
      status = RMN_OK
      return
   end function form_drag1

end module form_drag
