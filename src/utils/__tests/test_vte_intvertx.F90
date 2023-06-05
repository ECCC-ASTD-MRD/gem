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

!/@
subroutine test_vte_intvertx()
   use, intrinsic :: iso_fortran_env, only: INT64
   use testutils
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-02
   !@/
#include <clib_interface_mu.hf>
#include <rmnlib_basics.hf>
   integer,parameter :: NI0=3,NJ0=2,NK0=3,NK2=NK0+4,LINBOT=1
   real,parameter :: dst_data_c(NK0) = (/1.343915,2.517857,3.458333/)
   real,parameter :: dst_data_l(NK0) = (/1.333333,2.5,3.5/)
   real,parameter :: dst_data_cl1(NK0) = (/1.333333,2.517857,3.458333/)
   integer :: i,j,k
   real :: dst_lev1d(NK0),src_lev1d(NK2)
   real,dimension(NI0,NJ0,NK2) :: src_lev,src_data
   real,dimension(NI0,NJ0,NK0) :: dst_lev,dst_data
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_vte_intvertx')
   call testutils_set_tolerence(1.d-5)

   dst_lev1d(1:NK0) = (/300.,500.,700./)*100.
   src_lev1d(1:NK2) = (/100.,250.,400.,600.,800.,900.,950./)*100.

   do j=1,NJ0
      do i=1,NI0
         dst_lev(i,j,1:NK0) = dst_lev1d(1:NK0)
         src_lev(i,j,1:NK2) = src_lev1d(1:NK2) + 10.*100.*float(i-2+(j-2)*NI0)/(NI0*NJ0)
      enddo
   enddo

   do k=1,NK2
      src_data(:,:,k) = float(k-1)
   enddo

!!$   print *,'',src_data(2,2,1:NK2)
!!$   print *,'',src_lev(2,2,1:NK2)
!!$   print *,'',dst_lev(2,2,1:NK0)

   call vte_intvertx3(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'','cubic')
   call testutils_assert_eq(dst_data(2,2,:),dst_data_c,'vte_intvertx3 c')

!!$   call vte_intvertx3b(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'','cubic')
!!$   call testutils_assert_eq(dst_data(2,2,:),dst_data_c,'vte_intvertx3b c')

   call vte_intvertx4(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'',0)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_c,'vte_intvertx4 c')

   call vte_intvertx_isodst(dst_data,src_data,src_lev,dst_lev1d,NI0*NJ0,NK2,NK0,'',0)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_c,'vte_intvertx_i c')

   call vte_intvertx3(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'','linear')
   call testutils_assert_eq(dst_data(2,2,:),dst_data_l,'vte_intvertx3 l')

!!$   call vte_intvertx3b(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'','linear')
!!$   call testutils_assert_eq(dst_data(2,2,:),dst_data_l,'vte_intvertx3b l')

   call vte_intvertx4(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'',NK0)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_l,'vte_intvertx4 l')

   call vte_intvertx_isodst(dst_data,src_data,src_lev,dst_lev1d,NI0*NJ0,NK2,NK0,'',NK0)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_l,'vte_intvertx_i l')
   
   call vte_intvertx4(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'',LINBOT)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_cl1,'vte_intvertx4 cl')

   call vte_intvertx_isodst(dst_data,src_data,src_lev,dst_lev1d,NI0*NJ0,NK2,NK0,'',LINBOT)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_cl1,'vte_intvertx_i cl')

   !- reverse order not working
   src_lev1d(1:NK2) = (/950.,900.,800.,600.,400.,250.,100./)*100.

   do j=1,NJ0
      do i=1,NI0
         src_lev(i,j,1:NK2) = src_lev1d(1:NK2) + 10.*100.*float(i-2+(j-2)*NI0)/(NI0*NJ0)
      enddo
   enddo
   do k=1,NK2
      src_data(:,:,NK2-k+1) = float(k-1)
   enddo

   call vte_intvertx4(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'',0)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_c,'vte_intvertx4 c sort down')

   call vte_intvertx_isodst(dst_data,src_data,src_lev,dst_lev1d,NI0*NJ0,NK2,NK0,'',0)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_c,'vte_intvertx_i c sort down')

   call vte_intvertx4(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'',NK0)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_l,'vte_intvertx4 l sort down')

   call vte_intvertx_isodst(dst_data,src_data,src_lev,dst_lev1d,NI0*NJ0,NK2,NK0,'',NK0)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_l,'vte_intvertx_i l sort down')

   call vte_intvertx4(dst_data,src_data,src_lev,dst_lev,NI0*NJ0,NK2,NK0,'',LINBOT)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_cl1,'vte_intvertx4 cl sort down')

   call vte_intvertx_isodst(dst_data,src_data,src_lev,dst_lev1d,NI0*NJ0,NK2,NK0,'',LINBOT)
   call testutils_assert_eq(dst_data(2,2,:),dst_data_cl1,'vte_intvertx_i cl sort down')

!!$   print *,'4i.lr',dst_data(2,2,:)
   ! ---------------------------------------------------------------------
   return
end subroutine test_vte_intvertx
