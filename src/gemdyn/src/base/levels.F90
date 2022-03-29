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
      use VERgrid_options
      implicit none

      integer, intent(INOUT) :: F_gnk
      integer, intent(IN   ) :: F_nk,F_nkequal
      real, intent(IN   ) :: F_lid,F_first
      real, intent(INOUT) :: F_levels(F_nk)

      integer k,kk,nkequal,nk_s,nk0,nk1,nk2,j0,jn,k0,kn
      real tmp
      real(kind=REAL64) prdz(F_gnk), dz, dz0, ddz, top
!
!---------------------------------------------------------------------
!
      if (F_levels(2).lt.0.) then ! AUTOMATIC VERTICAL LAYERING

         nkequal = max(min(F_nkequal,F_gnk-3),0)
         nk_s = F_gnk - nkequal

         if (F_levels(1).lt.0.) then ! STRETCHED DEPTHS
    
! Iterative process to determine the stretching factor (r).
            prdz(1) = F_first
            prdz(2) = F_first*2.
            F_levels(1)= prdz(1)
            F_levels(2)= F_levels(1) + prdz(2)

            nk0= 3 ; k0= 1 ; kn= 2 ; top= F_levels(2)
            do k= 1, F_gnk
               nk1= min(Hyb_lin_depth(1,k),real(nk_s))
               if (nk1<1) exit
               if (nk1 > kn) then
                  if (Hyb_lin_depth(2,k) > top) then
                     top= min(Hyb_lin_depth(2,k),F_lid)
                     dz = dble(top-F_levels(nk0-1))/dble(nk1-nk0+1)
                     prdz(nk1) = dz
                     do kk= nk0,nk1
                        F_levels(kk)= F_levels(kk-1) + dz
                     end do
                     if (((nk0-3)>k0).and.((nk0 +3)<nk1)) then
                        j0 = nk0 -3
                        jn = nk0 +3
                        nk2 = jn-j0+1
                        dz0= F_levels(j0)-F_levels(j0-1)
                        ddz = dz-F_levels(j0)+F_levels(j0-1)
                        do kk=j0,jn
                           dz = dz0+ddz/nk2*(kk-j0+1)
                           F_levels(kk)= F_levels(kk-1) + dz
                        end do
                     endif
                     k0= nk0 ; kn= nk1 ; top= F_levels(nk1)
                     nk0= nk1+1
                  endif
               endif
            end do
            
            nk0 = nk0 -1
            nk1 = nk_s -nk0+1
            nk2 = F_gnk-nk0+1
            call cal_depth (prdz(nk0),F_lid,F_levels(nk0),nk1,nk2)
            
            do k=nk0,nk_s
               F_levels(k)= F_levels(k-1) + prdz(k)
            end do

            do k=nk_s+1,F_gnk
               F_levels(k) = F_levels(k-1) + prdz(nk_s)
            end do

         else                   ! EQUAL

            do k=1,nk_s
               prdz(k)=dble(F_lid)/dble(nk_s)
            end do

         endif
         
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

      subroutine cal_depth (dz,lid,fh,F_nks,F_gnk)
      use, intrinsic :: iso_fortran_env
      use lun
      implicit none
      
      integer, intent(IN) :: F_nks,F_gnk
      real   , intent(IN) :: lid,fh
      real(kind=REAL64), intent(OUT) :: dz(0:F_gnk)

      integer i,k,err
      real(kind=REAL64) :: r,r_saved ,inc, ratio
!
!---------------------------------------------------------------------
!
      inc  = 1.e-3
      r    = 1.d0 ; err= 0
      r_saved = r
      do k=1,50
         call find_r (dz,lid,fh,r,inc,F_nks,F_gnk)
         if (r<0.) then
            do i=1,10
               inc=2.d0*inc
               call find_r (dz,lid,fh,r,inc,F_nks,F_gnk)
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
!     
!---------------------------------------------------------------------
!
      return
      end subroutine cal_depth

      subroutine find_r (dz,lid,fh,r,inc,nk_s,gnk)
      use, intrinsic :: iso_fortran_env
      implicit none
      
      integer, intent(IN)  :: nk_s,gnk
      real, intent(IN   )  :: lid,fh
      real(kind=REAL64), intent(INOUT) :: inc
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
      zm= fh
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
      if (cnt==1) then
         r   = r - inc
         inc = inc/10.
         goto 987
      endif
      r= r - inc
!     
!---------------------------------------------------------------------
!
      return
      end subroutine find_r
