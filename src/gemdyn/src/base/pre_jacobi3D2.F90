!i---------------------------------- LICENCE BEGIN -------------------------------
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

!**s/r  pre_jacobio3D -  BlocJacobi_3D additive-Schwarz preconditioner with
!                        local solver = "separable" approximate elliptic problem
!
      subroutine pre_jacobi3D2 ( Sols,Rhs_b,Ni,Nj,Nk)
      use dyn_fisl_options
      use glb_ld
      use opr
      use sol_mem
      use prec
      use mem_tstp
      use, intrinsic :: iso_fortran_env
      implicit none

      integer Ni,Nj,Nk
      real(kind=REAL64), dimension(Ni,Nj,Nk), intent(in) :: Rhs_b
      real(kind=REAL64), dimension(Ni,Nj,Nk), intent(out) :: Sols

!author
!       Abdessamad Qaddouri -  2013

      integer i,j,k,ii,jj
!
!     ---------------------------------------------------------------
!
!$omp do
      do j= 1, nj
            jj=Sol_jj0 + j+l_j0 - 2
         call dgemm ( 'N','N', Ni, nk, nk, 1.0D0,     &
                     rhs_b(1,j,1), Ni*Nj,Opr_lzevec_8,&
                     G_nk,0.0d0, w2_8 (1,j,1), Ni*Nj )
            do k=1,Nk
               do i = 1, ni
                  ii= Sol_ii0+i+l_i0 - 2
                  w2_8(i,j,k)= Opr_opsxp0_8(G_ni+ii) * &
                      Opr_opsyp0_8(G_nj+jj) * w2_8(i,j,k)
               end do
            end do
         end do

!$omp end do

!$omp do
      do k=1,Nk
         call dgemm ( 'T','N',Ni,Nj,Ni,1.0d0,Prec_xevec_8,Ni,&
                      w2_8(1,1,k),Ni,0.0d0,fdg(1,1,k),Ni)

         do j =2, Nj
            do i=1,Ni
               fdg(i,j,k) = fdg(i,j,k) - Prec_ai_8(i,j,k)*fdg(i,j-1,k)
            end do
         end do

         j = Nj
         do i=1,Ni
            fdg(i,j,k) = fdg(i,j,k)/Prec_bi_8(i,j,k)
         end do

         do j = Nj-1, 1, -1
            do i=1 , Ni
               fdg(i,j,k)=(fdg(i,j,k)-Prec_ci_8(i,j,k)*fdg(i,j+1,k))/Prec_bi_8(i,j,k)
            end do
         end do

         call dgemm ( 'N','N',Ni,Nj,Ni,1.0d0,Prec_xevec_8,Ni,&
                      fdg(1,1,k),Ni,0.d0,w2_8(1,1,k),Ni )
      end do
!$omp end do

!$omp do
      do j=1,Nj
         call dgemm ( 'N','T', Ni, nk, nk, 1.0D0,     &
                      w2_8(1,j,1), Ni*Nj, Opr_zevec_8,&
                      G_nk,0.0d0, w3_8 (1,j,1), Ni*Nj  )
      end do
!$omp end do

!$omp single
        do k=1,nk
           do j=1,nj
              do i=1,ni
                 Sols(i,j,k)= w3_8 (i,j,k)
              enddo
           enddo
        enddo
!$omp end single
!
!     ---------------------------------------------------------------
!
      return
      end

