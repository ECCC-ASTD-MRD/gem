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
 
function pbl_simple(km,kt,uwind,vwind,gzmom,gztherm,z0,n,nk) result(status)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   !Arguments
   integer, intent(in) :: n                          !horizontal dimension
   integer, intent(in) :: nk                         !vertical dimension
   real, dimension(n), intent(in) :: z0              !momentum roughness length (m)
   real, dimension(n,nk), intent(in) :: uwind        !x-component wind at time-minus (m/s)
   real, dimension(n,nk), intent(in) :: vwind        !y-component wind at time-minus (m/s)
   real, dimension(n,nk), intent(in) :: gzmom        !momentum level heights (m)
   real, dimension(n,nk), intent(in) :: gztherm      !thermodynamic level heights (m)
   real, dimension(n,nk), intent(out) :: km          !momentum-level diffusion coefficients
   real, dimension(n,nk), intent(out) :: kt          !thermodynamic-level diffusion coefficients
   integer :: status                                 !completion status (RMN_OK or RMN_ERR)

   !Author
   !          R. McTaggart-Cowan (Fall 2016)
   
   !Revisions

   !Object
   !          Establish simple diffusion coefficients for a neutral profile as
   !          proposed by Beljaars, Brown and Wood (2004; QJRMS).
   
   !Notes

   ! Local parameter definitions
   real, parameter :: LAMBDA=150.

   ! Local variable declarations
   integer :: i,k
   real :: mlenm,mlent
   real, dimension(n,nk) :: shearmag

   ! Initialize status
   status = RMN_ERR

   ! Compute the profile of vector shear magnitude
   shearmag(:,1) = sqrt( &
        (uwind(:,2)-uwind(:,1) / (gzmom(:,2)-gzmom(:,1)))**2 + &
        (vwind(:,2)-vwind(:,1) / (gzmom(:,2)-gzmom(:,1)))**2)
   do k=2,nk-1
      shearmag(:,k) = sqrt( &
           (uwind(:,k+1)-uwind(:,k-1) / (gzmom(:,k+1)-gzmom(:,k-1)))**2 + &
           (vwind(:,k+1)-vwind(:,k-1) / (gzmom(:,k+1)-gzmom(:,k-1)))**2)
   enddo

   ! Compute coefficients using length scales from Blackadar formulation
   do k=1,nk-1
      do i=1,n
         mlenm = 1./((1./(KARMAN*(gzmom(i,k)+z0(i))))+(1./LAMBDA))
         km(i,k) = mlenm**2 * shearmag(i,k)
         mlent = 1./((1./(KARMAN*(gztherm(i,k)+z0(i))))+(1./LAMBDA))
         kt(i,k) = mlent**2 * shearmag(i,k)
      enddo
   enddo
   km(:,nk) = km(:,nk-1)
   kt(:,nk) = kt(:,nk-1)

   ! Successful completion
   status = RMN_OK
   return
end function pbl_simple
