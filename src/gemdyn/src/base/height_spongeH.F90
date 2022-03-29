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

!**s/r height_spongeH

      subroutine height_spongeH()
      use gmm_vt0
      use theo_options
      use metric
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none

      integer i,j,k
      real :: zblen_bot, zblen_top, deltaZ, beta
      real(kind=REAL64) :: work1, fact
!
!----------------------------------------------------------------------
!
      zblen_top = Ver_z_8%m(0)
      zblen_bot=zblen_top-mtn_zblen_thk
      deltaZ= zblen_top-zblen_bot
      fact=1.d0
      if(Theo_case_S == 'MTN_SCHAR' .or. Theo_case_S == 'MTN_SCHAR2' ) then
         fact=sqrt(2.0*mtn_flo*Cstv_dt_8/mtn_dx)
      end if

      if(Theo_case_S /= 'MTN_SCHAR' ) then
!$omp do collapse(2)
      do k=1,l_nk
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               work1=GVM%zmom_8(i,j,k)-zblen_bot
               work1=min(1.d0,max(0.d0,work1/deltaZ))
               beta= work1*work1*min(1.d0,fact)
               ut0(i,j,k)=(1.-beta)*ut0(i,j,k)+beta*mtn_flo
            end do
         end do
      end do
!$omp enddo
      end if
!$omp do collapse(2)
      do k=1,l_nk
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               work1=GVM%ztht_8(i,j,k)-Zblen_bot
               work1=min(1.d0,max(0.d0,work1/deltaZ))
               beta=work1*work1*min(1.d0,fact)
               wt0(i,j,k)=(1.-beta)*wt0(i,j,k)
            end do
         end do
      end do
!$omp enddo
!
!----------------------------------------------------------------------
!
      return
      end subroutine height_spongeH
