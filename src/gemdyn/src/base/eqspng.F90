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

!**s/r eqspng - apply vertical diffusion at and near
!               the model top on u and v

      subroutine eqspng (F_u,F_v,Minx,Maxx,Miny,Maxy,Nk)
      use hvdif_options
      use glb_ld
      use lun
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(INOUT) :: F_u,F_v

   !______________________________________________________________________
   !        |                                             |           |   |
   ! NAME   |             DESCRIPTION                     | DIMENSION |I/O|
   !--------|---------------------------------------------|-----------|---|
   ! F_u    | x component of velocity                     | 3D (Nk)   |i/o|
   ! F_v    | y component of velocity                     | 3D (Nk)   |i/o|
   !________|_____________________________________________|___________|___|
   !
   !    There are nlev levels that will be diffused.
   !    There are nlev-1 coefficent passed in namelist by user.
   !
   !    The boundary conditions are flux=0. This is achieved by
   !    putting coef=0 at boundaries (see drawing below).
   !    Also the index km and kp make the wind derivatives zero
   !    at boundary.
   !
   !    ---- this u equals u(1) (see km in loop)   \
   !                                               |
   !    ==== coef(1)=0, Top boundary.              | this derivative = 0
   !                                               |
   !    ---- u(1) first levels diffused            /
   !
   !    ==== coef(2)=first coef passed by user
   !
   !    ---- u(2)
   !
   !        ...
   !
   !    ---- u(nlev-1)
   !
   !    ==== coef(nlev)=last coef passed by user
   !
   !    ---- u(nlev) last level diffused           \
   !                                               |
   !    ==== coef(nlev+1)=0, Bottom boundary.      | this derivative = 0
   !                                               |
   !    ---- this u equal u(nlev) (see kp in loop) /
   !
   !_____________________________________________________________________
   !

      real, dimension(l_ni,eq_nlev) :: u,v
      integer :: i,j,k,km,kp
!
!-------------------------------------------------------------------
!
!$omp do
   do j=1,l_nj
      do k=1,eq_nlev
         kp=min(eq_nlev,k+1)
         km=max(1,k-1)
         do i=1,l_ni
            u(i,k)=F_u(i,j,k)+eponmod(i,j)*(cp(k)*(F_u(i,j,kp)-F_u(i,j,k )) &
                                        -cm(k)*(F_u(i,j,k )-F_u(i,j,km)))
            v(i,k)=F_v(i,j,k)+eponmod(i,j)*(cp(k)*(F_v(i,j,kp)-F_v(i,j,k )) &
                                        -cm(k)*(F_v(i,j,k )-F_v(i,j,km)))
         end do
      end do
      F_u(1:l_ni,j,1:eq_nlev)=u(1:l_ni,1:eq_nlev)
      F_v(1:l_ni,j,1:eq_nlev)=v(1:l_ni,1:eq_nlev)
   end do
!$omp end do
!
!-------------------------------------------------------------------
!
      return
      end subroutine eqspng
