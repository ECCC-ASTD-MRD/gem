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
!**s/r hzd_nudeln

      subroutine hzd_nudeln2(rfd,sfd,Minx,Maxx,Miny,Maxy,lnk,nu,m,n)
      use HORgrid_options
      use glb_ld
      implicit none
#include <arch_specific.hf>

      integer Minx,Maxx,Miny,Maxy,lnk,m,n
      real rfd (Minx:Maxx,Miny:Maxy,lnk), &
           sfd (Minx:Maxx,Miny:Maxy,lnk)
      real*8 nu

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

      integer i,j,k,id,jd,iff,jf,i0,in,j0,jn
      real wk(l_minx:l_maxx,l_miny:l_maxy)
      real*8 c1,c2,c3,one,two,four
      parameter(one=1.d0,two=2.d0,four=4.d0)
!
!----------------------------------------------------------------------
!
      if (Grd_yinyang_L) then
         id = 1    + 2*west
         jd = 1    + 2*south
         iff= l_ni - 2*east
         jf = l_nj - 2*north
      else
         id = 1    + pil_w
         jd = 1    + pil_s
         iff= l_ni - pil_e
         jf = l_nj - pil_n
      end if
      c1 = nu*(one-two*nu)
      c2 = nu**2
      c3 = nu*four*(nu-one)

!$omp parallel private (i,j,k,wk,i0,j0,in,jn) &
!$omp shared(c1,c2,c3,id,jd,iff,jf)
!$omp do
      do k=1,lnk
         i0= id-1 ; in= iff+1
         j0= jd-1 ; jn= jf +1
         if (m == 1) then
            sfd(i0:in,j0:jn,k) = rfd(i0:in,j0:jn,k)
         else if (m == 2) then
            sfd(i0:in,j0:jn,k) = rfd(i0:in,j0:jn,k) - sfd(i0:in,j0:jn,k)
         else
            i0= i0+west  ; in= in-east
            j0= j0+south ; jn= jn-north
            sfd(i0:in,j0:jn,k) = rfd(i0:in,j0:jn,k) - sfd(i0:in,j0:jn,k)
         end if
         if (m == n) then
            do j=jd,jf
            do i=id,iff
               rfd(i,j,k)= rfd(i,j,k) +  &
                       c1*(sfd(i  ,j+1,k)+sfd(i+1,j  ,k) + &
                           sfd(i  ,j-1,k)+sfd(i-1,j  ,k))+ &
                       c2*(sfd(i-1,j+1,k)+sfd(i+1,j+1,k) + &
                           sfd(i-1,j-1,k)+sfd(i+1,j-1,k))+ &
                       c3* sfd(i  ,j  ,k)
            end do
            end do
         else
            do j=jd,jf
            do i=id,iff
               wk(i,j)  = c1*(sfd(i  ,j+1,k)+sfd(i+1,j  ,k) + &
                              sfd(i  ,j-1,k)+sfd(i-1,j  ,k))+ &
                          c2*(sfd(i-1,j+1,k)+sfd(i+1,j+1,k) + &
                              sfd(i-1,j-1,k)+sfd(i+1,j-1,k))+ &
                          c3* sfd(i  ,j  ,k)
            end do
            end do
            sfd(id:iff,jd:jf,k) = rfd(id:iff,jd:jf,k) + wk(id:iff,jd:jf)
         end if
      end do
!$omp enddo
!$omp end parallel
!
!----------------------------------------------------------------------
!
      return
      end



