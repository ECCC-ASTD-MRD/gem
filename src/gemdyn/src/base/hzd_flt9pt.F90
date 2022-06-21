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
!**s/r hzd_flt9pt  - 9 points explicit horizontal diffusion operator

      subroutine hzd_flt9pt (rfd,sfd,Minx,Maxx,Miny,Maxy,lnk,&
                             nu,m,n,i0,in,j0,jn)
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,lnk,m,n,i0,in,j0,jn
      real, intent(INOUT) :: rfd (Minx:Maxx,Miny:Maxy,lnk), &
                             sfd (Minx:Maxx,Miny:Maxy,lnk)
      real(kind=REAL64), intent(IN) :: nu

!OBJECT
!
!    *******************************************************************
!    *                                                                 *
!    *                      OPERATOR nu_DEL2_n                         *
!    *                                                                 *
!    *          n  successive calls to this subroutine produces        *
!    *                                                                 *
!    *                                           ----------------      *
!    *                                           |            n |      *
!    *          the equivalent of a diffusion:   | -(-nu_DEL2)  |      *
!    *                                           |              |      *
!    *                                           ----------------      *
!    *                                                                 *
!    *                                                                 *
!    *                                                                 *
!    *                 Each call applies a 9-point filter              *
!    *                                                                 *
!    *              [ Shuman, M.W.R. #57, p.357-361, eq #5. ]          *
!    *                                                                 *
!    *                 to the difference:   rfd - sfd.                 *
!    *                                                                 *
!    *                                                                 *
!    *******************************************************************

      integer i,j,k
      real wk(Minx:Maxx,Miny:Maxy)!l_minx:l_maxx,l_miny:l_maxy)
      real(kind=REAL64) :: c1,c2,c3
      real(kind=REAL64), parameter :: one=1.d0, two=2.d0, four=4.d0
!
!----------------------------------------------------------------------
!
      c1 = nu*(one-two*nu)
      c2 = nu**2
      c3 = nu*four*(nu-one)

      if (m == n) then
!$omp do collapse(2)
         do k=1,lnk
            do j=j0,jn
            do i=i0,in
               rfd(i,j,k)= rfd(i,j,k) +  &
                       c1*(sfd(i  ,j+1,k)+sfd(i+1,j  ,k) + &
                           sfd(i  ,j-1,k)+sfd(i-1,j  ,k))+ &
                       c2*(sfd(i-1,j+1,k)+sfd(i+1,j+1,k) + &
                           sfd(i-1,j-1,k)+sfd(i+1,j-1,k))+ &
                       c3* sfd(i  ,j  ,k)
            end do
            end do
         end do
!$omp end do
      else

!$omp do
         do k=1,lnk
            do j=j0,jn
            do i=i0,in
               wk(i,j)= c1*(sfd(i  ,j+1,k)+sfd(i+1,j  ,k) + &
                                   sfd(i  ,j-1,k)+sfd(i-1,j  ,k))+ &
                               c2*(sfd(i-1,j+1,k)+sfd(i+1,j+1,k) + &
                                   sfd(i-1,j-1,k)+sfd(i+1,j-1,k))+ &
                               c3* sfd(i  ,j  ,k)
            end do
            end do
            do j=j0,jn
            do i=i0,in
               sfd(i,j,k) = rfd(i,j,k) + wk(i,j)
            end do
            end do
         end do
!$omp end do
      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine hzd_flt9pt



