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

module shallconv
   implicit none
   private
   public :: shallconv5

contains

   !/@*
   subroutine shallconv5(ztplus, zhumoins, &
        zsigw, zgztherm, zpmoins, &
        ztshal, zhushal, ztlcs, ztscs, zqlsc, zqssc, zkshal, &
        zfsc, zqcz, zqdifv, cdt1, ni, nk, nkm1)
      use debug_mod, only: init2nan
      use tdpack_const, only: GRAV, RAUW
      use phy_options
      implicit none
!!!#include <arch_specific.hf>

      integer, intent(in) :: ni, nk, nkm1
      real, intent(in) :: cdt1
      real, dimension(:), pointer :: zpmoins, ztlcs, ztscs, zkshal
      real, dimension(:,:), pointer :: &
           ztplus, zhumoins, zsigw, &
           zgztherm, ztshal, zhushal, zqlsc, zqssc, zfsc, zqcz, zqdifv

      !@Author A-M Leduc. (March 2002)
      !@Object This is the interface routine for shallow convection schemes
      !@Arguments
      ! cdt1     timestep
      ! ni       horizontal running length
      ! nk       vertical dimension
      !*@/
#include <rmn/msg.h>

      integer :: i
      integer, dimension(ni,nk) :: zilab
      real, dimension(ni) :: zdbdt
      real, dimension(ni,nk) :: geop
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'shallconv [BEGIN]')

      call init2nan(zdbdt)

      geop = zgztherm*GRAV
      zilab = 0.

      if (conv_shal == 'KTRSNT')then

         !# Kuo transient shallow convection - it does generate precipitation
         call ktrsnt4(ztshal, zhushal, zilab, zfsc, zqlsc, &
              zqssc, zdbdt,ztlcs,ztscs,zqcz, &
              ztplus, zhumoins, &
              geop, zqdifv, zpmoins, &
              zsigw, cdt1, zkshal, ni, nkm1)

         !# Transformation en hauteur d'eau (m) - diviser par densite eau
         do i=1,ni
            ztscs(i) = ztscs(i) / RAUW
            ztlcs(i) = ztlcs(i) / RAUW
         end do

      endif
      call msg_toall(MSG_DEBUG, 'shallconv [END]')
      !----------------------------------------------------------------
      return
   end subroutine shallconv5

end module shallconv

