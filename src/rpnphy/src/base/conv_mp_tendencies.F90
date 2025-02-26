
module conv_mp_tendencies
   implicit none
   private
   public :: conv_mp_tendencies1

contains

   subroutine conv_mp_tendencies1(liq_tend, ice_tend, ttp, qcp, ncp, qip, nip, qti1p, nti1p, zti1p, tdmaskxdt, ni, nk, nkm1)
      use phy_options, only: stcond,my_ccntype
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
      real, dimension(:,:), pointer, contiguous :: nti1p                   !Ice number mixing ratio (#/kg)
      real, dimension(:,:), pointer, contiguous :: zti1p                   !Ice 6th-moment mixing ratio (m6/kg)
      real, dimension(:), pointer, contiguous :: tdmaskxdt                 !Tendency mask * delt

#include <rmn/msg.h>
      
      ! Local variables
      logical :: is_p3v5
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

            !#TODO: use apply tendencies for clarity? Conditional outside loop (where?)
            do k = 1,nkm1
               do i = 1,ni
                  
                  !cloud droplets (number):
                  if (liq_tend(i,k) > 0.) then
                     if (qcp(i,k) > 1.e-9 .and. ncp(i,k) > 1.e-3) then
                        !assume mean size of the detrained droplets is the same size as those in the pre-existing clouds
                        ncp(i,k) = ncp(i,k) + (ncp(i,k)/qcp(i,k))*tdmaskxdt(i)*liq_tend(i,k)
                        ncp(i,k) = min(ncp(i,k), 1.e+9)
                     else
                        !initialize cloud number mixing ratios based on specified aerosol type
                        ncp(i,k) = ncp_prescribed  !maritime aerosols
                     endif
                  endif

                  !ice crystals (number):
                  if (ice_tend(i,k) > 0.) then
                     if (qip(i,k) > 1.e-9 .and. nip(i,k) > 1.e-3) then
                        !assume mean size of the detrained ice is the same size as those in the pre-existing "anvil"
                        nip(i,k) = nip(i,k) + (nip(i,k)/qip(i,k))*tdmaskxdt(i)*ice_tend(i,k)
                        nip(i,k) = min(nip(i,k), 1.e+7)
                     else
                        !initialize ice number mixing ratio based on Cooper (1986)
                        nip(i,k) = 5.*exp(0.304*(TRPL-max(233.,ttp(i,k))))
                        nip(i,k) = min(nip(i,k), 1.e+7)
                     endif
                  endif
                  
               enddo
            enddo

            !cloud droplets (mass) and ice crystals (mass)
            call apply_tendencies(qcp,      qip, &
                 &                liq_tend, ice_tend, tdmaskxdt, ni, nk, nkm1)

         elseif (stcond(1:5) == 'MP_P3') then

            !cloud droplets (number):
            ! note:  In the P3v3, cloud droplet number is prescribed (there is no 'ncp').
            is_p3v5 = (stcond == 'MP_P3')
            if (is_p3v5) then  !#Note: Only for P3v5 (not P3v3)
               do k = 1,nkm1
                  do i = 1,ni
                     if (liq_tend(i,k) > 0.) then
                        if (qcp(i,k) > 1.e-9 .and. ncp(i,k) > 1.e-3) then
                           !assume mean size of the detrained droplets is the same size as those in the pre-existing clouds
                           ncp(i,k) = ncp(i,k) + (ncp(i,k)/qcp(i,k))*tdmaskxdt(i)*liq_tend(i,k)
                           ncp(i,k) = min(ncp(i,k), 1.e+9)
                        else
                           !initialize cloud number mixing ratios based on specified aerosol type
                           ncp(i,k) = ncp_prescribed  !maritime aerosols
                        endif
                     endif
                  enddo
               enddo
            endif

            !ice crystals (number):
            ! note:  The initialization of ice number is the same for P3 as for MY2, but the two schemes have
            !         different variables.
            !#TODO: use apply tendencies for clarity? Conditional outside loop (where?)
            do k = 1,nkm1
               do i = 1,ni
                  if (ice_tend(i,k) > 0.) then
                     if (qti1p(i,k) > 1.e-9 .and. nti1p(i,k) > 1.e-3) then
                        !assume mean size of the detrained ice is the same size as in the pre-existing "anvil"
                        nti1p(i,k) = nti1p(i,k) + (nti1p(i,k)/qti1p(i,k))*tdmaskxdt(i)*ice_tend(i,k)
                        nti1p(i,k) = min(nti1p(i,k), 1.e+7)
                        if (associated(zti1p) .and. is_p3v5) zti1p(i,k) = min(zti1p(i,k), 1.e-18)
                     else
                        !initialize ice number mixing ratio based on Cooper (1986)
                        nti1p(i,k) = 5.*exp(0.304*(TRPL-max(233.,ttp(i,k))))
                        nti1p(i,k) = min(nti1p(i,k), 1.e+7)
                        if (associated(zti1p) .and. is_p3v5) zti1p(i,k) = min(zti1p(i,k), 1.e-18)
                     endif
                  endif
               enddo
            enddo

            !cloud droplets (mass) and ice crystals (mass):
            !
            call apply_tendencies(qcp,      qti1p, &
                 &                liq_tend, ice_tend, tdmaskxdt, ni, nk, nkm1)

         endif IF_MP_MY2

      elseif (associated(qcp)) then! not using combination kfc-bkf and MY

         total_tend = 0.
         if (associated(liq_tend)) total_tend(:,1:nkm1) = total_tend(:,1:nkm1) + liq_tend(:,1:nkm1)
         if (associated(ice_tend)) total_tend(:,1:nkm1) = total_tend(:,1:nkm1) + ice_tend(:,1:nkm1)
         call apply_tendencies(qcp, total_tend, tdmaskxdt, ni, nk, nkm1)

      endif IF_KFBE_MP

   end subroutine conv_mp_tendencies1

end module conv_mp_tendencies
