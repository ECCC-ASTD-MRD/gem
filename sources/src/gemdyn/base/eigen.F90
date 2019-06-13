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

!**s/r geneigl - solves a generalised symmetric eigenproblem

      subroutine eigen ( F_eval_8, F_evec_8, F_b_8, NN, NMAX, NWORK,&
                            F_eigen_filename_S )
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: NN, NMAX, NWORK
      character(len=*), intent(in) :: F_eigen_filename_S
      real(kind=REAL64), dimension(NN), intent(out) :: F_eval_8
      real(kind=REAL64), dimension(NMAX,NN), intent(inout) :: F_evec_8, F_b_8

!object
!    To solve a generalised symmetric eigenproblem
!
!            a * x = lambda * b * x
!
!arguments:
!  Name        I/O                 Description
!----------------------------------------------------------------
!  F_eval_8    O     - eigenvalues (lambda)
!  F_evec_8    I/O   - input: matrix A
!                     output: eigenvectors (x)
!  F_b_8       I/O   - input: matrix B
!                     output:
!  NN          I     - order of problem
!  NMAX        I     - leading dimension of matrices in calling programm
!  NWORK       I     - size of F_work_8 >= max(1,3*n-1)
!
!      Only upper triangular parts of A and B need to be specified
!      A and  B are overwritten


      integer :: i, j, k, info, errop, unf, err
      real(kind=REAL64) :: faz_8, sav_8
      real(kind=REAL64), dimension(NWORK) :: wk1
      real(kind=REAL64), parameter :: one_8 = 1.0d0
!
!--------------------------------------------------------------------
!
      unf= 0

      F_eval_8 = 0.d0

      if (Ptopo_myproc == 0) then
         open ( unf,file=trim(F_eigen_filename_S),status='OLD', &
                form='unformatted',iostat=errop )

         if ( errop == 0 ) then
            write(6,1001) 'READING', trim(F_eigen_filename_S)
            read (unf) F_evec_8,F_eval_8
            close(unf)
         else
            info = -1
            call DSYGV( 1, 'V', 'U', NN, F_evec_8, NMAX, F_b_8, NMAX, &
                        F_eval_8, wk1, NWORK, info )
            do j=1,NN
               faz_8 = sign( one_8, F_evec_8(1,j) )
               do i= 1, NN
                  F_evec_8(i,j) = faz_8 * F_evec_8(i,j)
               end do
            end do
            do j= 1, NN/2
               k = NN - j + 1
               sav_8 = F_eval_8(j)
               F_eval_8(j) = F_eval_8(k)
               F_eval_8(k) = sav_8
               do i= 1, NN
                  sav_8 = F_evec_8(i,j)
                  F_evec_8(i,j) = F_evec_8(i,k)
                  F_evec_8(i,k) = sav_8
               end do
            end do
            if (Ptopo_couleur == 0) then
               open ( unf, file=trim(F_eigen_filename_S), &
                      form='unformatted',iostat=errop )
               if ( errop == 0 ) then
                  write(6,1001) 'WRITING', trim(F_eigen_filename_S)
                  write(unf) F_evec_8, F_eval_8
                  close(unf)
               end if
            end if
         end if
      end if

      call RPN_COMM_bcast (F_evec_8, NMAX*NN, "MPI_DOUBLE_PRECISION",0, "grid", err)
      call RPN_COMM_bcast (F_eval_8, NN, "MPI_DOUBLE_PRECISION", 0, "grid", err)

 1001 format (/' GENEIGL: ',a,' FILE ',a)
!
!--------------------------------------------------------------------
!
      return
      end
