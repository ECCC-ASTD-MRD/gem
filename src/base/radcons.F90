
module radcons
   implicit none
   private
   
   public :: radcons2
   
contains

subroutine radcons2(pvars, ni)
   use tdpack_const
   use phybusidx, only: sla, fsa, c1slop, c2slop, c3slop, c4slop, c5slop
   use phymem, only: phyvar
   implicit none
!!!#include <arch_specific.hf> 

   type(phyvar), pointer, contiguous :: pvars(:)
   integer, intent(in) :: ni

   !@Author  J.Mailhot, A. Erfani, J. Toviessi (July 2009)
   !@Object Calculate constants required to define the radiation along the 
   !        slopes. 
   !Arguments
   ! pvars         list of all phy vars (meta + slab data)
   ! ni            horizontal dimension 

#include "phymkptr.hf"
  
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

   real, dimension(4) :: sla_cos, sla_sin, sla_2cos, sla_2sin
   real :: sla_rad

   real, pointer, contiguous, dimension(:) :: zc1slop, zc2slop, zc3slop, zc4slop, zc5slop
   real, pointer, contiguous, dimension(:,:) :: zsla, zfsa  !, entbus, perbus

   MKPTR2DNL(zsla, sla, pvars)
   MKPTR2DNL(zfsa, fsa, pvars)
   MKPTR1D(zc1slop, c1slop, pvars)
   MKPTR1D(zc2slop, c2slop, pvars)
   MKPTR1D(zc3slop, c3slop, pvars)
   MKPTR1D(zc4slop, c4slop, pvars)
   MKPTR1D(zc5slop, c5slop, pvars)

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

   do i = 1, ni

      do n = 1, 4
         sla_rad  = zsla(i,n)*pi180
         sla_cos(n) = cos(sla_rad)
         sla_sin(n) = sin(sla_rad)
         sla_2sin(n) = sin(sla_rad*0.50)
         sla_2cos(n) = cos(sla_rad*0.50)
      enddo

      sum =  zfsa(i,1) + zfsa(i,2) + zfsa(i,3) + zfsa(i,4) 

      zc1slop(i) = zfsa(i,1)*sla_cos(1) +   & 
           zfsa(i,2)*sla_cos(2) +          &
           zfsa(i,3)*sla_cos(3) +         &
           zfsa(i,4)*sla_cos(4)
      zc1slop(i) = zc1slop(i) + (1-sum)

      zc2slop(i) = zfsa(i,1)*sla_sin(1)*cos_n +   & 
           zfsa(i,2)*sla_sin(2)*cos_e +           &
           zfsa(i,3)*sla_sin(3)*cos_s +          & 
           zfsa(i,4)*sla_sin(4)*cos_w

      zc3slop(i) = zfsa(i,1)*sla_sin(1)*sin_n +   & 
           zfsa(i,2)*sla_sin(2)*sin_e +           &
           zfsa(i,3)*sla_sin(3)*sin_s +          & 
           zfsa(i,4)*sla_sin(4)*sin_w

      zc4slop(i) = zfsa(i,1)*sla_2cos(1)*sla_2cos(1) + & 
           zfsa(i,2)*sla_2cos(2)*sla_2cos(2) +       &
           zfsa(i,3)*sla_2cos(3)*sla_2cos(3) +      & 
           zfsa(i,4)*sla_2cos(4)*sla_2cos(4)
      zc4slop(i) = zc4slop(i)+ (1-sum)

      zc5slop(i) = zfsa(i,1)*sla_2sin(1)*sla_2sin(1) +  & 
           zfsa(i,2)*sla_2sin(2)*sla_2sin(2) +        &
           zfsa(i,3)*sla_2sin(3)*sla_2sin(3) +       & 
           zfsa(i,4)*sla_2sin(4)*sla_2sin(4)

   end do

   return
end subroutine radcons2

end module radcons
