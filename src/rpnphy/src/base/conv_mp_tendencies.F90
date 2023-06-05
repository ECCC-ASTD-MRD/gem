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

module conv_mp_tendencies
   implicit none
   private
   public :: conv_mp_tendencies1

contains

   subroutine conv_mp_tendencies1(liq_tend, ice_tend, ttp, qcp, ncp, qip, nip, qti1p, nti1p, tdmask, ni, nk, nkm1)
      use phy_options, only: stcond,my_ccntype,delt
      use tdpack_const, only: TRPL
      use tendency, only: apply_tendencies
      implicit none
      ! Argument declaration
      integer, intent(in) :: ni                                            !Row length
      integer, intent(in) :: nk                                            !Total number of vertical levels
      integer, intent(in) :: nkm1                                          !Number of prognostic vertical levels
      real, dimension(:,:), pointer, contiguous :: liq_tend                !Tendency for liquid condensate mixing ratio (kg/kg/s)
      real, dimension(:,:), pointer, contiguous :: ice_tend                !Tendency for solid condensate mixing ratio (kg/kg/s)
      real, dimension(:,:), pointer, contiguous :: ttp                     !Time-plus dry air temperature (K)
      real, dimension(:,:), pointer, contiguous :: qcp                     !Liquid condensate mixing ratio (kg/kg)
      real, dimension(:,:), pointer, contiguous :: ncp                     !Liquid droplet number concentration (#/m3)
      real, dimension(:,:), pointer, contiguous :: qip                     !Ice condensate mixing ratio (kg/kg)
      real, dimension(:,:), pointer, contiguous :: nip                     !Ice number concentration (#/m3)
      real, dimension(:,:), pointer, contiguous :: qti1p                   !Ice condensate mixing ratio (kg/kg)
      real, dimension(:,:), pointer, contiguous :: nti1p                   !Ice number concentration (#/m3)
      real, dimension(:), pointer, contiguous :: tdmask                    !Tendency mask

#include <rmn/msg.h>
#include "phymkptr.hf"

      ! Local variables
      integer :: i,k
      real :: ncp_prescribed
      real, dimension(ni,nk) :: total_tend

      ! Set aerosol number concentration
      if (my_ccntype == 1) then
         ncp_prescribed = 0.8e+8  !maritime aerosols
      else
         ncp_prescribed = 2.0e+8  !continential aerosols
      endif

      !Update ice and liquid independently if possible (deep == kfc .or. bkf and stcond == MY)
      IF_KFBE_MP: if (associated(liq_tend) .and. associated(ice_tend) .and. stcond(1:3) == 'MP_') then

         IF_MP_MY2: if (stcond(1:6) == 'MP_MY2') then !note: this includes 'MP_MY2'

            !cloud droplets (number):
            !#TODO: use apply tendencies for clarity? Conditional outside loop (where?)
            do k = 1,nkm1
               do i = 1,ni
                  if (liq_tend(i,k) > 0.) then
                     if (qcp(i,k) > 1.e-9 .and. ncp(i,k) > 1.e-3) then
                        !assume mean size of the detrained droplets is the same size as those in the pre-existing clouds
                        ncp(i,k) = ncp(i,k) + (ncp(i,k)/qcp(i,k))*tdmask(i)*liq_tend(i,k)*delt
                        ncp(i,k) = min(ncp(i,k), 1.e+9)
                     else
                        !initialize cloud number mixing ratios based on specified aerosol type
                        ncp(i,k) = ncp_prescribed  !maritime aerosols
                     endif
                  endif
               enddo
            enddo

            !cloud droplets (mass)
            call apply_tendencies(qcp, liq_tend, tdmask, ni, nk, nkm1)

            !ice crystals (number):
            !#TODO: use apply tendencies for clarity? Conditional outside loop (where?)
            do k = 1,nkm1
               do i = 1,ni
                  if (ice_tend(i,k) > 0.) then
                     if (qip(i,k) > 1.e-9 .and. nip(i,k) > 1.e-3) then
                        !assume mean size of the detrained ice is the same size as those in the pre-existing "anvil"
                        nip(i,k) = nip(i,k) + (nip(i,k)/qip(i,k))*tdmask(i)*ice_tend(i,k)*delt
                        nip(i,k) = min(nip(i,k), 1.e+7)
                     else
                        !initialize ice number mixing ratio based on Cooper (1986)
                        nip(i,k) = 5.*exp(0.304*(TRPL-max(233.,ttp(i,k))))
                        nip(i,k) = min(nip(i,k), 1.e+7)
                     endif
                  endif
               enddo
            enddo

            !ice crystals (mass):
            call apply_tendencies(qip, ice_tend, tdmask, ni, nk, nkm1)

         elseif (stcond(1:5) == 'MP_P3') then

            !cloud droplets (number):
            ! note:  In the current version of P3, cloud droplet number is prescribed (there is no 'ncp').
            !        For future 2-moment cloud configuration, 'ncp' should be updated here (as above for MP_MY2)

            !cloud droplets (mass):
            call apply_tendencies(qcp, liq_tend, tdmask, ni, nk, nkm1)

            !ice crystals (number):
            ! note:  The initialization of ice number is the same for P3 as for MY2, but the two schemes have
            !         different variables.
            !#TODO: use apply tendencies for clarity? Conditional outside loop (where?)
            do k = 1,nkm1
               do i = 1,ni
                  if (ice_tend(i,k) > 0.) then
                     if (qti1p(i,k) > 1.e-9 .and. nti1p(i,k) > 1.e-3) then
                        !assume mean size of the detrained ice is the same size as in the pre-existing "anvil"
                        nti1p(i,k) = nti1p(i,k) + (nti1p(i,k)/qti1p(i,k))*tdmask(i)*ice_tend(i,k)*delt
                        nti1p(i,k) = min(nti1p(i,k), 1.e+7)
                     else
                        !initialize ice number mixing ratio based on Cooper (1986)
                        nti1p(i,k) = 5.*exp(0.304*(TRPL-max(233.,ttp(i,k))))
                        nti1p(i,k) = min(nti1p(i,k), 1.e+7)
                     endif
                  endif
               enddo
            enddo

            !ice crystals (mass):
            call apply_tendencies(qti1p, ice_tend, tdmask, ni, nk, nkm1)

         endif IF_MP_MY2

      elseif (associated(qcp)) then! not using combination kfc-bkf and MY

         total_tend = 0.
         if (associated(liq_tend)) total_tend(:,1:nkm1) = total_tend(:,1:nkm1) + liq_tend(:,1:nkm1)
         if (associated(ice_tend)) total_tend(:,1:nkm1) = total_tend(:,1:nkm1) + ice_tend(:,1:nkm1)
         call apply_tendencies(qcp, total_tend, tdmask, ni, nk, nkm1)

      endif IF_KFBE_MP

   end subroutine conv_mp_tendencies1

end module conv_mp_tendencies
