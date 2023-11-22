!----------------------------------LICENCE BEGIN -------------------------------
!     GEM - Library of kernel routines for the GEM numerical atmospheric model
!     Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!     Environnement Canada
!     This library is free software; you can redistribute it and/or modify it
!     under the terms of the GNU Lesser General Public License as published by
!     the Free Software Foundation, version 2.1 of the License. This library is
!     distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!     PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
!     You should have received a copy of the GNU Lesser General Public License
!     along with this library; if not, write to the Free Software Foundation, Inc.,
!     59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!----------------------------------LICENCE END ---------------------------------

!**s/r hzd_flt5pt  - 5 points explicit horizontal diffusion operator
!                    which includes geometric terms (hzd_geom_*)

      subroutine hzd_flt5pt ( F_champ, F_geom, Minx,Maxx,Miny,Maxy, NK,&
                              F_coef_8, i0,in,j0,jn )
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(IN) :: Minx, Maxx, Miny, Maxy, i0,in,j0,jn, NK
      real, dimension(Minx:Maxx,Miny:Maxy,NK), intent (INOUT) :: F_champ
      real(kind=REAL64), dimension(NK), intent(IN) :: F_coef_8
      real(kind=REAL64), dimension(Miny:Maxy,*), intent(IN) :: F_geom

      integer :: i,j,k
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy) :: w1
!
!---------------------------------------------------------------------
!
!$omp do
      do k=1,NK
         do j= j0, jn
            do i= i0, in
               w1(i,j)= F_geom(j,1) * F_champ(i  ,j  ,k) + &
                        F_geom(j,2) * F_champ(i-1,j  ,k) + &
                        F_geom(j,3) * F_champ(i+1,j  ,k) + &
                        F_geom(j,4) * F_champ(i  ,j-1,k) + &
                        F_geom(j,5) * F_champ(i  ,j+1,k)
            end do
          end do
          do j= j0, jn
             do i= i0, in
                F_champ(i,j,k) = F_champ(i,j,k) + F_coef_8(k) * w1(i,j)
             end do
          end do
       end do
!$omp end do
!
!----------------------------------------------------------------------
!
      return
      end subroutine hzd_flt5pt

