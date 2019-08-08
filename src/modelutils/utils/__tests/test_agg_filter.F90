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
subroutine test_agg_filter()
   use, intrinsic :: iso_fortran_env, only: INT64
   use testutils
   use agg_filter_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-02
!@/
#include <clib_interface_mu.hf>
#include <rmnlib_basics.hf>
   real, pointer :: data2d(:,:),data3d(:,:,:)
   integer :: i0,j0,in,jn,halo,ni,nj,nk,istat,i
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_agg_filter')

   halo=2 ; ni=4 ; nj=5 ; nk=2
   allocate(data2d(1-halo:ni+halo,1-halo:nj+halo),stat=istat)
   allocate(data3d(1-halo:ni+halo,1-halo:nj+halo,nk),stat=istat)

   data2d = -99.
   do i=1,ni
      data2d(i,:) = float(i)
   enddo
   i0=1;j0=1
   in=ubound(data2d,1)-halo
   jn=ubound(data2d,2)-halo
   call aggreg(data2d,i0,j0,in,jn,2,1)
   call testutils_assert_eq((/i0,j0,in,jn/),(/1,1,(ni+1)/2,nj/),'aggreg di=2 dims')
   call testutils_assert_eq(data2d(i0:in,j0),(/1.5,3.5/),'aggreg di=2 data')

   data2d = -99.
   do i=1,nj
      data2d(2:,i) = float(i)
   enddo
   i0=2;j0=1
   in=ubound(data2d,1)-halo
   jn=ubound(data2d,2)-halo
   call aggreg(data2d,i0,j0,in,jn,1,2)
   call testutils_assert_eq((/i0,j0,in,jn/),(/2,1,ni,(nj+1)/2/),'aggreg dj=2 dims')
   call testutils_assert_eq(data2d(i0,j0:jn),(/1.5,3.5,4.5/),'aggreg dj=2 data')

   data3d = -99.
   do i=2,nj
      data3d(:,i,nk) = float(i)
   enddo
   i0=1;j0=2
   in=ubound(data3d,1)-halo
   jn=ubound(data3d,2)-halo
   call aggreg(data3d,i0,j0,in,jn,1,2)
   call testutils_assert_eq((/i0,j0,in,jn/),(/1,2,ni,1+nj/2/),'aggreg3d dj=2 dims')
   call testutils_assert_eq(data3d(i0,j0:jn,nk),(/2.5,4.5/),'aggreg3d dj=2 data')
   ! ---------------------------------------------------------------------
   return
end subroutine test_agg_filter
