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
!
      subroutine levels (F_levels,F_nk,F_gnk,F_lid,F_nkequal,F_first)
      use, intrinsic :: iso_fortran_env
      use lun
      implicit none

      integer, intent(INOUT) :: F_gnk
      integer, intent(IN   ) :: F_nk,F_nkequal
      real, intent(IN   ) :: F_lid,F_first
      real, intent(INOUT) :: F_levels(F_nk)

      integer i,k,nk_s,err,nkequal
      real tmp
      real(kind=REAL64) r,r_saved,inc,prdz(F_gnk),ratio
!
!---------------------------------------------------------------------
!
      if (F_levels(2).lt.0.) then ! AUTOMATIC VERTICAL LAYERING

         nkequal = max(min(F_nkequal,F_gnk-3),0)
         nk_s = F_gnk - nkequal

         if (F_levels(1).lt.0.) then ! STRETCHED DEPTHS
    
! Iterative process to determine the stretching factor (r).
            inc  = 1.e-3
            r    = 1.d0 ; err= 0
            r_saved = r
            do k=1,50
               call cal_depth(prdz,F_lid,F_first,r,inc,nk_s,F_gnk)
               if (r<0.) then
                  do i=1,10
                     inc=2.d0*inc
                     call cal_depth(prdz,F_lid,F_first,r,inc,nk_s,F_gnk)
                     if (r>0) exit
                  end do
                  inc=inc/2.d0
               else
                  ratio= (r-r_saved)/r_saved
                  if (lun_out>0) write (Lun_out, '(a,1pe18.12,2x,a,1pe18.12)')&
               'Vertical stretching factor convergence: r=',r,'ratio= ',ratio
                  if (ratio<1.d-12) exit
                  inc=inc/10.
                  r_saved = r
               endif
            end do
            if (r<0) err= -1
            call gem_error(err,'levels','COULD NOT CONVERGE TO VERTICAL LAYERING SPECS')

         else                   ! EQUAL

            do k=1,nk_s
               prdz(k)=dble(F_lid)/dble(nk_s)
            end do

         endif

         F_levels(1)= prdz(1)
!         print*, 1,F_levels(1),prdz(1)
         do k=2,nk_s
            F_levels(k)= F_levels(k-1) + prdz(k)
!            print*, k,F_levels(k),prdz(k)
         end do

         do k=nk_s+1,F_gnk
            F_levels(k) = F_levels(k-1) + prdz(nk_s)
!            print*, k,F_levels(k),prdz(nk_s),'ff'
         end do

         do k=1,F_gnk/2
            tmp= F_levels(k)
            F_levels(k)= F_levels(F_gnk-k+1)
            F_levels(F_gnk-k+1)= tmp
         end do

      else                      ! MANUAL VERTICAL LAYERING
         
         F_gnk = 0
         do k = 1, size(F_levels)
            if ( F_levels(k) < 0.) exit
            F_gnk = k
         end do

      endif
!
!---------------------------------------------------------------------
!
      return
      end subroutine levels

      subroutine cal_depth (dz,lid,fh,r,inc,nk_s,gnk)
      use, intrinsic :: iso_fortran_env
      implicit none
      
      integer, intent(IN)  :: nk_s,gnk
      real, intent(IN   )  :: lid,fh
      real(kind=REAL64), intent(IN   ) :: inc
      real(kind=REAL64), intent(INOUT) :: r
      real(kind=REAL64), intent(OUT  ) :: dz(gnk)

      integer, parameter :: maxiter_gl = 500
      integer k,cnt
      real(kind=REAL64) :: zm
!
!---------------------------------------------------------------------
!
      cnt= 0
 987  r   = r + inc ; cnt= cnt+1
      dz(1) = fh ; zm= fh
      do k=2,nk_s
         dz(k) = dz(k-1) * r
         zm = zm + dz(k)
      end do
      do k=nk_s+1, gnk
         zm = zm + dz(nk_s)
      end do

      if (zm<lid) then
         if (cnt==maxiter_gl) then
            r= -1.
            return
         else
            goto 987
         endif
      endif
      r= r - inc
!     
!---------------------------------------------------------------------
!
      return
      end subroutine cal_depth
