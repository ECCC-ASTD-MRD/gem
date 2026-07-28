!---------------------------------- LICENCE BEGIN -------------------------------
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
!---------------------------------- LICENCE END ---------------------------------
!**s/r hzd_flt9pt  - 9 points explicit horizontal conservatif diffusion
!                                     Del2 operator

      subroutine hzd_CvDel2_flt9pt (rfd,F_sfd,rho,Minx,Maxx,Miny,Maxy,nk, &
                                    nu,m,n,i0,in,j0,jn)
      use hzd_mod
!
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,nk,m,n,i0,j0,in,jn 
      real, intent(INOUT) :: rfd (Minx:Maxx,Miny:Maxy,nk), rho(Minx:Maxx,Miny:Maxy,nk), &
                             F_sfd (Minx:Maxx,Miny:Maxy,nk)

      integer i,j,k
      real wk(Minx:Maxx,Miny:Maxy)!l_minx:l_maxx,l_miny:l_maxy)
      real(kind=REAL64) :: c1,c2,c3
      real(kind=REAL64), parameter :: one=1.d0, two=2.d0, four=4.d0
      real(kind=REAL64), parameter :: half=0.5d0, pt25=0.25d0
      real(kind=REAL64), parameter :: epsilon = 1.0d-12,zero=0.d0

      real(kind=REAL64) :: a,b,c,d,e,f,g,h,nu
!----------------------------------------------------------------------
     ! wk=zero
      c1 = nu*(one-two*nu)
      c2 = nu**2
      c3 = nu*four*(nu-one)

      if(m==n) then
              
!$omp do collapse(2)
         do k=1,nk
            do j=j0,jn
               do i=i0,in
                  a=half*(rho(i+1,j,k)+rho(i,j,k))
                  b=half*(rho(i,j,k)+rho(i-1,j,k))
                  c=half*(rho(i,j+1,k)+rho(i,j,k))
                  d=half*(rho(i,j,k)+rho(i,j-1,k))

                  e= pt25*(rho(i+1,j,k)+rho(i+1,j+1,k)+rho(i,j+1,k)+rho(i,j,k))
                  f= pt25*(rho(i+1,j,k)+rho(i+1,j-1,k)+rho(i,j-1,k)+rho(i,j,k))
                  g= pt25*(rho(i,j,k)+rho(i,j+1,k)+rho(i-1,j+1,k)+rho(i-1,j,k))
                  h= pt25*(rho(i,j,k)+rho(i-1,j,k)+rho(i-1,j-1,k)+rho(i,j-1,k))

                  rfd(i,j,k)=  rfd(i,j,k) + &
                     c1/rho(i,j,k)*(&
                  a*(F_sfd(i+1,j,k)-F_sfd(i,j,k))-b*(F_sfd(i,j,k)-F_sfd(i-1,j,k))+&
                  c*(F_sfd(i,j+1,k)-F_sfd(i,j,k))-d*(F_sfd(i,j,k)-F_sfd(i,j-1,k)) ) +&
                  c2/rho(i,j,k)*( &
                  e*(F_sfd(i+1,j+1,k)-F_sfd(i,j,k))-h*(F_sfd(i,j,k)-F_sfd(i-1,j-1,k))+&
                  g*(F_sfd(i-1,j+1,k)-F_sfd(i,j,k))-f*(F_sfd(i,j,k)-F_sfd(i+1,j-1,k)) )
               end do
            end do
         end do
!$omp end do
      else
!$omp do 
         do k=1,nk
            do j=j0,jn
               do i=i0,in
                  a=half*(rho(i+1,j,k)+rho(i,j,k))
                  b=half*(rho(i,j,k)+rho(i-1,j,k))
                  c=half*(rho(i,j+1,k)+rho(i,j,k))
                  d=half*(rho(i,j,k)+rho(i,j-1,k))

                  e= pt25*(rho(i+1,j,k)+rho(i+1,j+1,k)+rho(i,j+1,k)+rho(i,j,k))
                  f= pt25*(rho(i+1,j,k)+rho(i+1,j-1,k)+rho(i,j-1,k)+rho(i,j,k))
                  g= pt25*(rho(i,j,k)+rho(i,j+1,k)+rho(i-1,j+1,k)+rho(i-1,j,k))
                  h= pt25*(rho(i,j,k)+rho(i-1,j,k)+rho(i-1,j-1,k)+rho(i,j-1,k))

                  wk(i,j)=   &
                         c1/rho(i,j,k)*(&
                  a*(F_sfd(i+1,j,k)-F_sfd(i,j,k))-b*(F_sfd(i,j,k)-F_sfd(i-1,j,k))+&
                  c*(F_sfd(i,j+1,k)-F_sfd(i,j,k))-d*(F_sfd(i,j,k)-F_sfd(i,j-1,k)) ) +&
                  c2/rho(i,j,k)*( &
                  e*(F_sfd(i+1,j+1,k)-F_sfd(i,j,k))-h*(F_sfd(i,j,k)-F_sfd(i-1,j-1,k))+&
                  g*(F_sfd(i-1,j+1,k)-F_sfd(i,j,k))-f*(F_sfd(i,j,k)-F_sfd(i+1,j-1,k)) )
               end do
            end do
            do j=j0,jn
               do i=i0,in
                  F_sfd(i,j,k)= rfd(i,j,k) + wk(i,j)
               enddo
            enddo  
         end do
!$omp end do
      endif
!
!----------------------------------------------------------------------
!
      return
      end subroutine hzd_CvDel2_flt9pt



