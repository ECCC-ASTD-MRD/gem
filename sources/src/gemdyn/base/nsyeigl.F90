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
                                                 NN, NMAX, NWORK )
      implicit none
#include <arch_specific.hf>

      integer NN, NMAX, NWORK
      real*8 F_eval_8(NMAX), F_Levec_8(NMAX,NMAX), F_b_8(NMAX,NMAX), &
             F_Revec_8(NMAX,NMAX)

!author  Abdessamad Qaddouri - 2007
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_eval_8     O    - eigenvalues
! F_Revec_8    O    - Right eigenvector matrix
! F_levec_8    O    - left eigenvector matrix
!
!------------------------------------------------------------------

      real*8, external :: ddot
      integer  i,j,k,INFO
      real*8 alfr(NMAX),alfi(NMAX),beta(NMAX),wk1(NWORK)
      real*8 BVL(NMAX,NMAX),dotvl(NN),VL(NMAX,NMAX),b(NMAX,NMAX)
!
! --------------------------------------------------------------------
!
      INFO= -1 ; b= F_b_8

      call dggev('V','V',NN, F_Levec_8, NMAX, b, NMAX, &
                 alfr, alfi,beta,VL, NMAX, F_Revec_8, NMAX, &
                 wk1, NWORK, INFO )
      F_Levec_8=VL

! (ALPHAR(j) + ALPHAI(j)*i)/BETA(j) is the jeme eigenvalue
      do j=1,NMAX
         F_eval_8(j)= ALfR(j)/beta(j)
      end do

! normalize F_Revec_8
      do j=1,NN
         do i=1,NMAX
            BVL(i,j)=0.0
            do k=1,NN
               BVL(i,j)= BVL(i,j)+F_b_8(i,k)*F_Revec_8(k,j)
            end do
         end do
      end do
      do j=1,NN
         dotvl(j)= ddot(NMAX,F_Revec_8(1,j),1,BVL(1,j),1)
      end do

      do j=1,NN
         do i=1,NMAX
            F_Revec_8(i,j)= F_Revec_8(i,j)/sqrt(dotvl(j))
         end do
      end do

! normalize F_Levec_8
      do j=1,NN
         do i=1,NMAX
            BVL(i,j)=0.0
            do k=1,NN
               BVL(i,j)= BVL(i,j)+F_b_8(i,k)*F_Levec_8(k,j)
            end do
         end do
      end do
      do j=1,NN
         dotvl(j)= ddot(NMAX,F_Levec_8(1,j),1,BVL(1,j),1)
      end do

      do j=1,NN
         do i=1,NMAX
            F_Levec_8(i,j)= F_Levec_8(i,j)/sqrt(dotvl(j))
         end do
      end do

! normalize for non symmetric case
      do j=1,NN
         do i=1,NMAX
            BVL(i,j)=0.0
            do k=1,NN
               BVL(i,j)= BVL(i,j)+F_b_8(i,k)*F_Revec_8(k,j)
            end do
         end do
      end do

      do j=1,NN
         dotvl(j)= ddot(NMAX,F_Levec_8(1,j),1,BVL(1,j),1)
      end do

      do j=1,NN
         do i=1,NMAX
            F_Levec_8(i,j)= F_Levec_8(i,j)/dotvl(j)
         end do
      end do
!
! --------------------------------------------------------------------
!
      return
      end

