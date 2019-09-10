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

subroutine radcons2(ni, trnch)
   use tdpack_const
   use phybus, only: entbus, perbus, sla, fsa, c1slop, c2slop, c3slop, c4slop, c5slop
   implicit none
!!!#include <arch_specific.hf> 

   integer ni, trnch

   !@Author  J.Mailhot, A. Erfani, J. Toviessi (July 2009)
   !@Object Calculate constants required to define the radiation along the 
   !        slopes. 
   !Arguments
   ! ni            horizontal dimension 

   integer i, n
   real sum
   real pi180 

   real cos_n    ! North 
   real cos_e    ! East
   real cos_s    ! South 
   real cos_w    ! West 
   real sin_n    ! North 
   real sin_e    ! East 
   real sin_s    ! South 
   real sin_w    ! West 

   real, dimension(ni,4) :: sla_cos
   real, dimension(ni,4) :: sla_sin
   real, dimension(ni,4) :: sla_2cos
   real, dimension(ni,4) :: sla_2sin
   real, dimension(ni,4) :: sla_rad

#define MKPTR1D(NAME1,NAME2,BUS) nullify(NAME1); if (NAME2 > 0) NAME1(1:ni) => BUS(NAME2:,trnch)
#define MKPTR2D(NAME1,NAME2,N3,BUS) nullify(NAME1); if (NAME2 > 0) NAME1(1:ni,1:N3) => BUS(NAME2:,trnch)

   real, pointer, dimension(:) :: zc1slop, zc2slop, zc3slop, zc4slop, zc5slop
   real, pointer, dimension(:,:) :: zsla, zfsa

   MKPTR2D(zsla, sla, 4, entbus)
   MKPTR2D(zfsa, fsa, 4, entbus)
   MKPTR1D(zc1slop, c1slop, perbus)
   MKPTR1D(zc2slop, c2slop, perbus)
   MKPTR1D(zc3slop, c3slop, perbus)
   MKPTR1D(zc4slop, c4slop, perbus)
   MKPTR1D(zc5slop, c5slop, perbus)

   cos_n  = cos(pi)         ! North 
   cos_e  = cos(1.5*pi)     ! East 
   cos_s  = cos(0.0)        ! South
   cos_w  = cos(0.5*pi)     ! West
   sin_n  = sin(pi)         ! North
   sin_e  = sin(1.5*pi)     ! East
   sin_s  = sin(0.0)        ! South
   sin_w  = sin(0.5*pi)     ! West

   ! facteur de conversion d'angle de degre a radian
   pi180 =   pi/180.0

   do n = 1, 4
      do i = 1, ni
         sla_rad(i,n)  = zsla(i,n)*pi180
         sla_2cos(i,n) = sla_rad(i,n)*0.50
      enddo
   enddo

   call  vscos(sla_cos, sla_rad,   4*ni)
   call  vssin(sla_sin, sla_rad,   4*ni)
   call  vssin(sla_2sin,sla_2cos,  4*ni)
   call  vscos(sla_2cos,sla_2cos,  4*ni)

   do i = 1, ni

      sum =  zfsa(i,1) + zfsa(i,2) + zfsa(i,3) + zfsa(i,4) 

      zc1slop(i) = zfsa(i,1)*sla_cos(i,1) +   & 
           zfsa(i,2)*sla_cos(i,2) +          &
           zfsa(i,3)*sla_cos(i,3) +         &
           zfsa(i,4)*sla_cos(i,4)
      zc1slop(i) = zc1slop(i) + (1-sum)

      zc2slop(i) = zfsa(i,1)*sla_sin(i,1)*cos_n +   & 
           zfsa(i,2)*sla_sin(i,2)*cos_e +           &
           zfsa(i,3)*sla_sin(i,3)*cos_s +          & 
           zfsa(i,4)*sla_sin(i,4)*cos_w

      zc3slop(i) = zfsa(i,1)*sla_sin(i,1)*sin_n +   & 
           zfsa(i,2)*sla_sin(i,2)*sin_e +           &
           zfsa(i,3)*sla_sin(i,3)*sin_s +          & 
           zfsa(i,4)*sla_sin(i,4)*sin_w

      zc4slop(i) = zfsa(i,1)*sla_2cos(i,1)*sla_2cos(i,1) + & 
           zfsa(i,2)*sla_2cos(i,2)*sla_2cos(i,2) +       &
           zfsa(i,3)*sla_2cos(i,3)*sla_2cos(i,3) +      & 
           zfsa(i,4)*sla_2cos(i,4)*sla_2cos(i,4)
      zc4slop(i) = zc4slop(i)+ (1-sum)

      zc5slop(i) = zfsa(i,1)*sla_2sin(i,1)*sla_2sin(i,1) +  & 
           zfsa(i,2)*sla_2sin(i,2)*sla_2sin(i,2) +        &
           zfsa(i,3)*sla_2sin(i,3)*sla_2sin(i,3) +       & 
           zfsa(i,4)*sla_2sin(i,4)*sla_2sin(i,4)

   end do

   return
end subroutine radcons2
