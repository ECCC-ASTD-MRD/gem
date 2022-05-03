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

!*s/r nsyeigl - prepares projection matrix for the vertical modes
!               NK-th mode = computational mode
!               symmetric an nonsymmetric cases
!               nonsymmetric case is for staggered (LN(Z)) model
!
      subroutine nsyeigl ( F_eval_8, F_Levec_8, F_Revec_8, F_b_8,&
                           nn, nmax, nwork )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: nn, nmax, nwork
      real(kind=REAL64), dimension(nmax), intent(out) :: F_eval_8
      real(kind=REAL64), dimension(nmax,nmax), intent(inout) :: F_Levec_8
      real(kind=REAL64), dimension(nmax,nmax), intent(in) :: F_b_8
      real(kind=REAL64), dimension(nmax,nmax), intent(out) :: F_Revec_8

!author  Abdessamad Qaddouri - 2007
!
!arguments
!  Name           Description
!------------------------------------------------------------
! F_eval_8      - eigenvalues
! F_Revec_8     - Right eigenvector matrix
! F_levec_8     - left eigenvector matrix
!
!------------------------------------------------------------------

      real(kind=REAL64), external :: ddot
      integer  i,j,k,info
      real(kind=REAL64) alfr(nmax),alfi(nmax),beta(nmax),wk1(nwork)
      real(kind=REAL64) BVL(nmax,nmax),dotvl(nn),vl(nmax,nmax),b(nmax,nmax)
!
! --------------------------------------------------------------------
!
      INFO= -1 ; b= F_b_8

      call dggev('V','V',nn, F_Levec_8, nmax, b, nmax, &
                 alfr, alfi,beta,VL, nmax, F_Revec_8, nmax, &
                 wk1, nwork, info )
      F_Levec_8=VL

! (ALPHAR(j) + ALPHAI(j)*i)/BETA(j) is the jth eigenvalue
      do j=1,nmax
         F_eval_8(j)= ALfR(j)/beta(j)
      end do

! normalize F_Revec_8
      do j=1,nn
         do i=1,nmax
            BVL(i,j)=0.0
            do k=1,nn
               BVL(i,j)= BVL(i,j)+F_b_8(i,k)*F_Revec_8(k,j)
            end do
         end do
      end do
      do j=1,nn
         dotvl(j)= ddot(nmax,F_Revec_8(1,j),1,BVL(1,j),1)
      end do

      do j=1,nn
         do i=1,nmax
            F_Revec_8(i,j)= F_Revec_8(i,j)/sqrt(dotvl(j))
         end do
      end do

! normalize F_Levec_8
      do j=1,nn
         do i=1,nmax
            BVL(i,j)=0.0
            do k=1,nn
               BVL(i,j)= BVL(i,j)+F_b_8(i,k)*F_Levec_8(k,j)
            end do
         end do
      end do
      do j=1,nn
         dotvl(j)= ddot(NMAX,F_Levec_8(1,j),1,BVL(1,j),1)
      end do

      do j=1,nn
         do i=1,nmax
            F_Levec_8(i,j)= F_Levec_8(i,j)/sqrt(dotvl(j))
         end do
      end do

! normalize for non symmetric case
      do j=1,nn
         do i=1,nmax
            BVL(i,j)=0.0
            do k=1,nn
               BVL(i,j)= BVL(i,j)+F_b_8(i,k)*F_Revec_8(k,j)
            end do
         end do
      end do

      do j=1,nn
         dotvl(j)= ddot(NMAX,F_Levec_8(1,j),1,BVL(1,j),1)
      end do

      do j=1,nn
         do i=1,nmax
            F_Levec_8(i,j)= F_Levec_8(i,j)/dotvl(j)
         end do
      end do

      return
      end

