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

      subroutine adz_cfl
      use adz_mem
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none

      integer :: cfl_i(3,3)
      real(kind=REAL64) :: cfl_8(3)
!
!     ---------------------------------------------------------------
!
      cfl_8= 0.d0 ; cfl_i= 0

      call adz_courant (Adz_pxyzm, Adz_i0,Adz_in,Adz_j0,Adz_jn, &
                        l_ni,l_nj,Adz_k0,l_nk,cfl_i,cfl_8)

      if (lun_out > 0) write (output_unit,99) 'x,y',cfl_i(1,1),cfl_i(2,1), &
                                          cfl_i(3,1),sngl(cfl_8(1))
      if (lun_out > 0) write (output_unit,99) 'z'  ,cfl_i(1,2),cfl_i(2,2), &
                                          cfl_i(3,2),sngl(cfl_8(2))
      if (lun_out > 0) write (output_unit,99) '3D' ,cfl_i(1,3),cfl_i(2,3), &
                                          cfl_i(3,3),sngl(cfl_8(3))

 99   format(' MAX COURANT NUMBER:  ', a3,&
             ': [(',i4,',',i4,',',i4,') ',f12.5,']')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_cfl

      subroutine adz_courant ( F_xyz, i0, in, j0, jn, &
                        F_ni, F_nj, k0, F_nk, F_cfl_i,F_cfl_8 )
      use glb_ld
      use ver
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_ni,F_nj,F_nk
      integer, intent(in) :: i0, in, j0, jn, k0 !Scope
      real, dimension(3,F_ni,F_nj,F_nk), intent(in) :: F_xyz
      integer, intent(out) :: F_cfl_i(3,3)
      real(kind=REAL64),  intent(out) :: F_cfl_8(3)

      integer, dimension(3,3,Ptopo_numproc) :: iwk
      real(kind=REAL64) , dimension(3  ,Ptopo_numproc) :: wk_8

      integer :: i, j, k, cfl_i(3,3), err, iproc, imax,jmax,kmax
      integer :: imaxH,jmaxH,kmaxH,imaxV,jmaxV,kmaxV
      real(kind=REAL64)  :: x_cfl, y_cfl, z_cfl, xy_cfl, xyz_cfl, Hmax_cfl_8, &
                 Vmax_cfl_8, max_cfl_8, cfl_8(3)
!
!     ---------------------------------------------------------------
!
      imaxH = 0
      jmaxH = 0
      kmaxH = 0
      Hmax_cfl_8 = 0.D0
      imaxV = 0
      jmaxV = 0
      kmaxV = 0
      Vmax_cfl_8 = 0.D0
      imax = 0
      jmax = 0
      kmax = 0
      max_cfl_8 = 0.D0

      do k=k0,F_nk
         do j=j0,jn
            do i=i0,in
                 x_cfl= abs(F_xyz(1,i,j,k)-l_i0+1-dble(i))
                 y_cfl= abs(F_xyz(2,i,j,k)-l_j0+1-dble(j))
                 z_cfl= abs(F_xyz(3,i,j,k)       -dble(k))
                xy_cfl= sqrt(x_cfl*x_cfl + y_cfl*y_cfl)
               xyz_cfl= sqrt(x_cfl*x_cfl + y_cfl*y_cfl + z_cfl*z_cfl)

               if (xy_cfl > Hmax_cfl_8) then
                  imaxH = i
                  jmaxH = j
                  kmaxH = k
                  Hmax_cfl_8 = xy_cfl
               end if
               if (z_cfl > Vmax_cfl_8) then
                  imaxV = i
                  jmaxV = j
                  kmaxV = k
                  Vmax_cfl_8 = z_cfl
               end if
               if (xyz_cfl>max_cfl_8) then
                  imax = i
                  jmax = j
                  kmax = k
                  max_cfl_8=xyz_cfl
               end if
            end do
         end do
      end do

      cfl_8(1)   = Hmax_cfl_8
      cfl_i(1,1) = imaxH + l_i0 - 1
      cfl_i(2,1) = jmaxH + l_j0 - 1
      cfl_i(3,1) = kmaxH
      cfl_8(2)   = Vmax_cfl_8
      cfl_i(1,2) = imaxV + l_i0 - 1
      cfl_i(2,2) = jmaxV + l_j0 - 1
      cfl_i(3,2) = kmaxV
      cfl_8(3)   = max_cfl_8
      cfl_i(1,3) = imax + l_i0 - 1
      cfl_i(2,3) = jmax + l_j0 - 1
      cfl_i(3,3) = kmax

      call RPN_COMM_gather(cfl_8,3,"MPI_DOUBLE_PRECISION",wk_8,3, &
                           "MPI_DOUBLE_PRECISION",0,"GRID", err)
      call RPN_COMM_gather(cfl_i,9,"MPI_INTEGER",iwk,9, &
                           "MPI_INTEGER",0,"GRID", err)

      if (Ptopo_myproc == 0) then
         imax = iwk(1,1,1)
         jmax = iwk(2,1,1)
         kmax = iwk(3,1,1)
         max_cfl_8 = wk_8(1,1)
         do iproc = 2, Ptopo_numproc
            if (wk_8(1,iproc)>max_cfl_8) then
               imax = iwk(1,1,iproc)
               jmax = iwk(2,1,iproc)
               kmax = iwk(3,1,iproc)
               max_cfl_8 = wk_8(1,iproc)
            end if
         end do
         F_cfl_8(1)   = max_cfl_8
         F_cfl_i(1,1) = imax
         F_cfl_i(2,1) = jmax
         F_cfl_i(3,1) = kmax

         imax = iwk(1,2,1)
         jmax = iwk(2,2,1)
         kmax = iwk(3,2,1)
         max_cfl_8 = wk_8(2,1)
         do iproc = 2, Ptopo_numproc
            if (wk_8(2,iproc) > max_cfl_8) then
               imax = iwk(1,2,iproc)
               jmax = iwk(2,2,iproc)
               kmax = iwk(3,2,iproc)
               max_cfl_8 = wk_8(2,iproc)
            end if
         end do
         F_cfl_8(2)   = max_cfl_8
         F_cfl_i(1,2) = imax
         F_cfl_i(2,2) = jmax
         F_cfl_i(3,2) = kmax

         imax = iwk(1,3,1)
         jmax = iwk(2,3,1)
         kmax = iwk(3,3,1)
         max_cfl_8 = wk_8(3,1)
         do iproc = 2, Ptopo_numproc
            if (wk_8(3,iproc)>max_cfl_8) then
               imax = iwk(1,3,iproc)
               jmax = iwk(2,3,iproc)
               kmax = iwk(3,3,iproc)
               max_cfl_8 = wk_8(3,iproc)
            end if
         end do
         F_cfl_8(3)   = max_cfl_8
         F_cfl_i(1,3) = imax
         F_cfl_i(2,3) = jmax
         F_cfl_i(3,3) = kmax
      end if
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_courant
