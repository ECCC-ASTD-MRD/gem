!-------------------------------------- LICENCE BEGIN -----------------------
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
!-------------------------------------- LICENCE END -------------------------

module water_integrated
   implicit none
   private
   public :: water_integrated1

contains

   !/@*
   subroutine water_integrated1(tt,qq,qc,qi,qnp,sigma,ps, &
        zicw,ziwv,ziwv700,ziwp,zlwp2,zslwp,zslwp2,zslwp3,zslwp4, &
        ni,nk)
!#TODO: never used qr,qgp
      use debug_mod, only: init2nan
      use tdpack_const, only: TCDK
      use phy_options
      implicit none
!!!#include <arch_specific.hf>
      !@Author L.Spacek, November 2011
      !@Object Calculate integrated quantities of some variables
      !@Arguments

      integer               :: ni,nk
      real,dimension(ni,nk) :: tt,qq,qc,qi,qnp,sigma
      real,dimension(ni) :: ps,zicw,ziwv,ziwv700,ziwp,zlwp2, &
           zslwp,zslwp2,zslwp3,zslwp4

      !          - Input -
      ! fbus     historic variables for the physics
      ! fsiz     dimension of fbus
      ! vsiz     dimension of vbus
      ! tt       temperature
      ! qq       humidity
      ! qc       total condensate mixing ratio at t+dT
      ! qi       ice mixing ratio (M-Y, K-Y) at t+dT
      ! qnp      snow    mixing ratio (M-Y) at t+dT
      ! sigma    vertical coordinate
      ! ni       horizontal running length
      ! nk       vertical dimension
      !          - Input/Output -
      ! vbus     physics tendencies and other output fields from the physics
      !          icw    - integrated cloud water/ice
      !          iwv    - integrated water vapor
      !          iwv700 - integrated water vapor (0-700 mb)
      !          iwp    - integrated ice water
      !          lwp2   - liquid water path (Sundqvist)
      !          slwp   - integrated SLW (supercooled liquid water)
      !          slwp2  - integrated SLW (bottom to s2)
      !          slwp3  - integrated SLW (s2 to s3)
      !          slwp4  - integrated SLW (s3 to s4)
      !
      !*@/
#include <msg.h>

      logical               :: integrate=.false.
      integer               :: i,k
      real                  :: tcel,frac
      real,dimension(ni)    :: temp1,temp2
      real,dimension(ni,nk) :: liquid, solid
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'water_integrated [BEGIN]')

      call init2nan(temp1, temp2)
      call init2nan(liquid, solid)

      ! Arrays 'liquid' and 'solid' are passed to intwat3  and used
      ! for diagnostic calculations only.
      ! Computaion of lwc and iwc used by radiation code is done in prep_cw.

      if (stcond.eq.'NEWSUND'.or. stcond.eq.'CONSUN') then
         integrate=.true.
         do k=1,nk
            do i=1,ni
               tcel = min(0.,tt(i,k) - TCDK)
               temp1(i) = -.003102 * tcel*tcel
            end do
            call vsexp(temp2,temp1,ni )
            do i=1,ni
               if (tt(i,k) .ge. TCDK) then
                  liquid(i,k) =  qc(i,k)
                  solid(i,k)  =  0.

               else
                  frac = .0059 + .9941 * temp2(i)
                  liquid(i,k) = frac*qc(i,k)
                  solid(i,k)  = (1.-frac)*qc(i,k)
               end if
            end do
         end do
      end if

      if (stcond(1:2)=='MP') then
         integrate=.true.
         do k=1,nk
            do i=1,ni
               liquid(i,k) = qc(i,k)
               !note: for stcond=mp_p3, qnp is passed in zero and qi is the sum of all qitot for all ice categories
               solid(i,k)  = qi(i,k)+qnp(i,k)
            end do
         end do
      endif

      !     calcul de quantites integrees
      if (integrate) &
           call intwat3(zicw,ziwv,ziwv700,ziwp,zlwp2, &
           zslwp,zslwp2,zslwp3,zslwp4, &
           tt,qq,liquid,solid,sigma,ps,ni,nk)

      call msg_toall(MSG_DEBUG, 'water_integrated [END]')
      !----------------------------------------------------------------
      return
   end subroutine water_integrated1

end module water_integrated
