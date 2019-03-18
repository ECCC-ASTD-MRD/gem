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

subroutine lightning2(foudre_rt, zp0_plus, zsigm, ztplus, zwplus, q_grpl, iiwc, ni, nk)

   use phy_options
   use phybus

   implicit none

#include <arch_specific.hf>

   !@Arguments
   integer, intent(in)                  :: ni,nk      !horizontal/vertical dimensions
   real, dimension(ni), intent(out)     :: foudre_rt  !lightning flash rate density,            # m-2 s-1
   real, dimension(ni), intent(in)      :: zp0_plus  !surface pressure,                        Pa
   real, dimension(ni,nk), intent(in)   :: zsigm      !sigma
   real, dimension(ni,nk), intent(in)   :: ztplus     !air temperature                          K
   real, dimension(ni,nk), intent(in)   :: zwplus     !vertical motion
   real, dimension(ni,nk-1), intent(in) :: q_grpl     !mixing ratio of graupel+hail,            kg kg-1
   real, dimension(ni,nk-1), intent(in) :: iiwc       !mixing ratio of sum of all ice species,  kg kg-1

   !@Author Anna Glazer   (Sept 2014)
   !@Revisions
   ! Anna Glazer (Avril 2015) - Pointers for version 4.8-a1
   ! Anna Glazer (Mai 2015)   - Adapted to 4.8.a5
   !@Object
   !  Compute the lightning threat expressed as
   !  flash-rate density (number of flashes per sec and m2)
   !  f3 = r1*f1 + r2*f2   ! f3 = f(foudre)
   !
   !  ref: McCaul et al., Wea. Forecasting 2009, vol. 40, pp. 709-729

#include "tdpack_const.hf"

  !Local variables and parameters:
   integer, dimension(ni,nk-1) :: array_t15
   integer, dimension(ni)      :: index_t15,k_t15
   integer :: i, k, ik !, nkm1

   real, dimension(ni)     :: iiwp, f1
   real :: t15, r1, r2, k1, k2, dsg, dpsg, i_conv, i_grav, tempo

   real, parameter :: conv  = 3.e+8    !to convert f3 unit from 5minlcn_p3i1)
   real, parameter :: fmin  = 1.e-12   !s-1 m-2, flash-rate density threshold for output
   real, parameter :: seuil = 1.e-6    !kg kg-1, mixing ratio threshold for graupel flux calculation

   t15 = tcdk - 15.
   r1  = 0.95
   r2  = 0.05
   k1  = 0.042
   k2  = 0.20
   i_conv = 1./conv
   i_grav = 1./grav

   !  Graupel flux (at -15 deg C level) contribution
   !  f1 = k1 * (wplus*QG); vertical velocity and graupel mixing ratio , both at -15 deg C level
   !
   !  Find level (k=k_t15(i))whose temperature most closely matches -15 deg C
   !  if there is none then index_t15(i)=0 and k_t15(i)=0

   i_loop: do i = 1,ni

      index_t15(i) = 0
      k_t15(i)     = 0
      f1(i)        = 0.
      iiwp(i)      = 0.

      do k = 1,nk-1
         if ((ztplus(i,k) .ge. t15 .and. ztplus(i,k+1) .lt. t15) .or. &
             (ztplus(i,k) .le. t15 .and. ztplus(i,k+1) .gt. t15)) then
            index_t15(i)=index_t15(i)+1
            if(abs(ztplus(i,k)-t15) .lt. abs(ztplus(i,k+1)-t15)) then
               array_t15(i,index_t15(i)) = k
            else
               array_t15(i,index_t15(i)) = k+1
            endif
         endif
      enddo

      if (index_t15(i) .ge. 1) k_t15(i) = array_t15(i,index_t15(i))

      if (index_t15(i) .ge. 1 .and. zwplus(i,k_t15(i)) .gt. 0. .and.  &
          q_grpl(i,k_t15(i)) .gt. seuil) then
         f1(i) = 1.e+3*k1*zwplus(i,k_t15(i))*q_grpl(i,k_t15(i))
      endif

      dsg  = 0.5 * ( zsigm(i,2) - zsigm(i,1) )
      dpsg = zp0_plus(i)*dsg*i_grav
      iiwp(i) = max(iiwc(i,1) , 0. ) * dpsg

      do k = 2,nk-1
         dsg  = 0.5 * ( zsigm(i,k+1) - zsigm(i,k-1) )
         dpsg = zp0_plus(i)*dsg*i_grav
         iiwp(i) = iiwp(i) + max(iiwc(i,k) , 0. ) * dpsg
      enddo

!       dsg  = 0.5 * ( zsigm(i,nk) - zsigm(i,nk-1) )
!       dpsg = zp0_plus(i)*dsg*i_grav
!       iiwp(i) = iiwp(i) + max(iiwc(i,nk) , 0. ) * dpsg

   !  Combine all together and find lighting threat f3 (variable foudre_rt)

      tempo = r1*f1(i) + r2*k2*iiwp(i)
      if (i_conv*tempo .gt. fmin) then
         foudre_rt(i) = i_conv * tempo
      else
         foudre_rt(i) = 0.
      endif

   enddo i_loop

   return
end subroutine lightning2

