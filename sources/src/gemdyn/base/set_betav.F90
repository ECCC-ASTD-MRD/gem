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



      subroutine set_betav(betav_m, betav_t, F_s, F_sl, Minx, Maxx, Miny, Maxy, Nk)

      use cstv
      use metric
      use mtn_options
      use tdpack
      use theo_options
      use ver

      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk

      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: betav_m, betav_t
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in) :: F_s, F_sl
!
      real :: htop, zblen_bot, zblen_top

      real*8 :: work1, work2, fact

      integer :: i, j, k

      zblen_top = Ver_z_8%m(0)

      fact=1.d0
      if(Theo_case_S == 'MTN_SCHAR' .or. Theo_case_S == 'MTN_SCHAR2' ) then
         fact=sqrt(2.0*mtn_flo*Cstv_dt_8/mtn_dx)
      end if

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then
         zblen_bot=zblen_top-mtn_zblen_thk
         do k=1,l_nk
            do j=Miny,Maxy
            do i=Minx,Maxx
               work1=zmom(i,j,k)-zblen_bot
               work2=zblen_top-zblen_bot
               work1=min(1.d0,max(0.d0,work1/work2))
               betav_m(i,j,k)=work1*work1*min(1.d0,fact)
               work1=ztht(i,j,k)-Zblen_bot
               work1=min(1.d0,max(0.d0,work1/work2))
               betav_t(i,j,k)=work1*work1*min(1.d0,fact)
            end do
            end do
         end do
      else
         htop = (mtn_nk+1)*mtn_dz
         zblen_bot=Cstv_Zsrf_8-grav_8*(htop-mtn_zblen_thk)/(287.*Cstv_tstr_8)
         do k=1,l_nk
            do j=Miny,Maxy
            do i=Minx,Maxx
               work1=zblen_bot-(Ver_a_8%m(k)+Ver_b_8%m(k)*F_s(i,j)+Ver_c_8%m(k)*F_sl(i,j))
               work2=zblen_bot-zblen_top
               work1=min(1.d0,max(0.d0,work1/work2))
               betav_m(i,j,k)=work1*work1*min(1.d0,fact)
               work1=zblen_bot-(Ver_a_8%t(k)+Ver_b_8%t(k)*F_s(i,j)+Ver_c_8%t(k)*F_sl(i,j))
               work1=min(1.d0,max(0.d0,work1/work2))
               betav_t(i,j,k)=work1*work1*min(1.d0,fact)
            end do
            end do
         end do
      end if

      return

      end
