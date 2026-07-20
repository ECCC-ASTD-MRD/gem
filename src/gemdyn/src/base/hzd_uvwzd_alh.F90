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
!**s/r hzd_alh

      subroutine hzd_uvwzd_alh(rfd,HzdLnr,HzdPwr,Minx,Maxx,Miny,Maxy,nk,F_v)
      use HORgrid_options
      use gem_options
      use glb_ld
      use glb_pil
      use hvdif_options 
      use hzd_mod
      use mem_tstp
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer :: Minx,Maxx,Miny,Maxy,nk,F_v
      integer :: i,j,k,dim,n,m,i0,in,j0,jn
      integer :: HzdPwr
      real    :: HzdLnr  
      real, dimension(Minx:Maxx,Miny:Maxy,nk), intent(inout) ::rfd
!
!----------------------------------------------------------------------
!
      if (Grd_yinyang_L) then
         i0 = 1    + 2*west
         j0 = 1    + 2*south
         in = l_ni - 2*east
         jn = l_nj - 2*north
      else
         i0 = 1    + pil_w
         j0 = 1    + pil_s
         in = l_ni - pil_e
         jn = l_nj - pil_n
      end if

      n= Hzdpwr/2.d0

!$omp do collapse(2)
      do k=1,nk
         do j=1-G_haloy, l_nj+G_haloy
            do i=1-G_halox, l_ni+G_halox
               sfd(i,j,k)= rfd(i,j,k)
            end do
         end do
      end do
!$omp end do

      do m=1,n 
         if (m == 2) then
!$omp do collapse(2)
            do k=1, nk
               do j=j0-1, jn+1
                  do i=i0-1, in+1
                     sfd1(i,j,k) = sfd(i,j,k)
                     sfd(i,j,k)  = rfd(i,j,k) - sfd(i,j,k)
                  end do
               end do
            end do
!$omp end do
         end if
         if (F_v==1) call hzd_u_alh (sfd,HzdLnr,minx,maxx,miny,maxy,nk)
         if (F_v==2) call hzd_v_alh (sfd,HzdLnr,minx,maxx,miny,maxy,nk)
         if (F_v==3) call hzd_scal_alh (sfd,HzdLnr,minx,maxx,miny,maxy,nk)
!$omp single         
         if (F_v==4) call hzd_theta_cons_alh (sfd,HzdLnr,minx,maxx,miny,maxy,nk,1)
!$omp end single                                         
         if(m==2) then
!$omp do collapse(2)
            do k=1,nk
               do j=Miny, Maxy
                  do i=Minx, Maxx
                      sfd(i,j,k)= sfd1(i,j,k)+sfd(i,j,k)
                  enddo
               enddo
            enddo
!$omp end do
         endif
      end do

!$omp do collapse(2)
      do k=1,nk
         do j=Miny, Maxy
            do i=Minx, Maxx
               rfd(i,j,k)= sfd(i,j,k)
            enddo
         enddo
      enddo
!$omp end do

      return
      end



