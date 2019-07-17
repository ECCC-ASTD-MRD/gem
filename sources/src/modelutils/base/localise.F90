!---------------------------------- LICENCE BEGIN -----------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END -------------------------------

!/@*
subroutine localise(imx,imy,x,y,departx,departy, &
     hdepartx,hdeparty,Ndx,Ndy)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   include "rmnlib_basics.inc"
   !@object find the I,J that corresponds to the other Yin/Yang grid
   !@ author Abdessamad Qaddouri - October 2009
   !@arguments
   integer, intent(in) :: Ndx, Ndy
   integer, intent(out):: imx(Ndx,Ndy),imy(Ndx,Ndy)
   real(REAL64), intent(in) :: departx(*),departy(*),hdepartx,hdeparty
   real(REAL64), intent(in) :: x(Ndx,Ndy),y(Ndx,Ndy)
   !*@/
   integer i,j

   !     The chosen point must be less than the working point
   !

   do j=1,Ndy
      do i=1,Ndx
         imy(i,j)=floor((y(i,j)-departy(1))/hdeparty)+1
         imx(i,j)=floor((x(i,j)-departx(1))/hdepartx)+1
         !if (imx(i,j).le.0.or.imy(i,j).le.0) print *,i,j,'imx=',imx(i,j),'imy=',imy(i,j)
         if (y(i,j) < departy(imy(i,j))) imy(i,j)=imy(i,j)-1
         if (x(i,j) < departx(imx(i,j))) imx(i,j)=imx(i,j)-1
         if (y(i,j) > departy(imy(i,j)+1)) imy(i,j)=imy(i,j)+1
         if (x(i,j) > departx(imx(i,j)+1)) imx(i,j)=imx(i,j)+1
      end do
   end do
   return
end subroutine localise

subroutine localise2(imx,imy,x,y,departx,departy, &
     hdepartx,hdeparty,Ndx,Ndy,NI,NJ)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   include "rmnlib_basics.inc"
   !@object find the I,J that corresponds to the other Yin/Yang grid
   !@ author Abdessamad Qaddouri - October 2009
   !@arguments
   integer, intent(in) :: Ndx, Ndy, NI,NJ
   integer, intent(out):: imx(Ndx,Ndy),imy(Ndx,Ndy)
   real(REAL64), intent(in) :: departx(1-NI:2*NI),departy(1-NJ:2*NJ),hdepartx,hdeparty
   real(REAL64), intent(in) :: x(Ndx,Ndy),y(Ndx,Ndy)
   !*@/
   integer i,j

   !     The chosen point must be less than the working point
   !

   do j=1,Ndy
      do i=1,Ndx
         imy(i,j)=floor((y(i,j)-departy(1))/hdeparty)+1
         imx(i,j)=floor((x(i,j)-departx(1))/hdepartx)+1
         !if (imx(i,j).le.0.or.imy(i,j).le.0) print *,i,j,'imx=',imx(i,j),'imy=',imy(i,j)
         if (y(i,j) < departy(imy(i,j))) imy(i,j)=imy(i,j)-1
         if (x(i,j) < departx(imx(i,j))) imx(i,j)=imx(i,j)-1
         if (y(i,j) > departy(imy(i,j)+1)) imy(i,j)=imy(i,j)+1
         if (x(i,j) > departx(imx(i,j)+1)) imx(i,j)=imx(i,j)+1
      end do
   end do
   return
end subroutine localise2


!/@*
subroutine localise_blend(imx,imy,x,y,departx,departy, &
     xmin,xmax,ymin,ymax,hdepartx,hdeparty)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   include "rmnlib_basics.inc"
   !@object find the I,J that corresponds to the other Yin/Yang grid
   !        in the blending zone.
   !@author Abdessamad Qaddouri - October 2009
   !@arguments
   integer, intent(out) :: imx,imy
   integer, intent(in) :: xmin,xmax,ymin,ymax
   real(REAL64), intent(in) :: departx(xmin:xmax),departy(ymin:ymax),hdepartx,hdeparty
   real(REAL64), intent(in) :: x,y
   !*@/

   !     The chosen point must be less than the working point
   !     Division can be a bit inaccurate. Point found could be negative
   !     and should exist in depart* vector:
   !     (X(1-Gni) <imx < X(2*Gni-1))
   !     (Y(1-Gnj) <imy < y(2*Gnj-1))
   !     Control of min and max of point is in the communication initialization
   !

   imy=floor((y-departy(1))/hdeparty)+1
   imx=floor((x-departx(1))/hdepartx)+1
   !if (imx.le.0.or.imy.le.0) print *,'imx=',imx,'imy=',imy
   if (y < departy(imy)) imy=imy-1
   if (x < departx(imx)) imx=imx-1
   if (y > departy(imy+1)) imy=imy+1
   if (x > departx(imx+1)) imx=imx+1
   return
end subroutine localise_blend
