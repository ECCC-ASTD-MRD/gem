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
subroutine test_hgrid_wb()
   use, intrinsic :: iso_fortran_env, only: INT64
   use testutils
   use ezgrid_mod
   use hgrid_wb
   use rmn_gmm
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
   !@/
#include <rmnlib_basics.hf>
#include <clib_interface_mu.hf>
   integer,parameter :: NI0=11,NJ0=13
   integer :: istat,ig1,ig2,ig3,ig4, i,j,gid_z1,gid_z2,gid_z1b,i0,j0,lni,lnj,hx,hy
   real :: ax(NI0,1),ay(1,NJ0)
   character(len=2) :: grtyp_S,grref_S
   logical :: ok_L
   type(gmm_metadata) :: mymeta
   ! ---------------------------------------------------------------------
   call testutils_verbosity()

   ig1=900 ; ig2=0 ; ig3=43200 ; ig4=43100
   do i=1,NI0
      ax(i,1) = 10.+float(i)*0.25
   enddo
   do j=1,NJ0
      ay(1,j) = float(j)*0.25
   enddo
   gid_z1 = ezgdef_fmem(NI0,NJ0, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)
   ax = ax + 20.
   ay = ay + 5.
   gid_z1b = ezgdef_fmem(NI0,NJ0, 'Z', 'E', ig1,ig2,ig3,ig4, ax, ay)
   ig1=0 ; ig2=10 ; ig3=43200 ; ig4=43150
   gid_z2 = ezgdef_fmem(NI0,NJ0, 'Z',  'E', ig1,ig2,ig3,ig4, ax, ay)


   istat = hgrid_wb_put(' ',gid_z1)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_hgrid_wb','put no name')

   istat = hgrid_wb_put(' ',-1)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_hgrid_wb','invalid id')


   istat = hgrid_wb_put('z1',gid_z1)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_hgrid_wb','put')

   istat = hgrid_wb_put('z1',gid_z1)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_hgrid_wb','put duplicate')

   istat = hgrid_wb_put('z2',gid_z2,1,2,NI0-2,NJ0-3,0,0)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_hgrid_wb','put + params')

   istat = hgrid_wb_put('z2h',gid_z2,2,3,NI0-3,NJ0-4,2,3)
   call testutils_assert_ok(RMN_IS_OK(istat),'test_hgrid_wb','put + params h')


   istat = hgrid_wb_get(' ',gid_z1b)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_hgrid_wb','get no name')

   istat = hgrid_wb_get('z3',gid_z1b)
   call testutils_assert_ok(.not.RMN_IS_OK(istat),'test_hgrid_wb','get no such name')

   istat = hgrid_wb_get('z1',gid_z1b)
   ok_L = (RMN_IS_OK(istat) .and. gid_z1==gid_z1b)
   call testutils_assert_ok(ok_L,'test_hgrid_wb','get')
   if (.not.ok_L) then
      print *,gid_z1,':',istat,gid_z1b
   endif

   istat = hgrid_wb_get('z1',gid_z1b,i0,j0,lni,lnj)
   ok_L = (RMN_IS_OK(istat) .and. gid_z1==gid_z1b .and. &
        all((/i0,j0,lni,lnj/) == (/1,1,NI0,NJ0/)))
   call testutils_assert_ok(ok_L,'test_hgrid_wb','get + params')
   if (.not.ok_L) then
      print *,gid_z1,':',istat,gid_z1b
      print *,1,1,NI0,NJ0,':',i0,j0,lni,lnj
   endif

   istat = hgrid_wb_get('z2',gid_z1b,i0,j0,lni,lnj,hx,hy)
   ok_L = (RMN_IS_OK(istat) .and. gid_z2==gid_z1b .and. &
        all((/i0,j0,lni,lnj,hx,hy/) == (/1,2,NI0-2,NJ0-3,0,0/)))
   call testutils_assert_ok(ok_L,'test_hgrid_wb','get + params2')
   if (.not.ok_L) then
      print *,gid_z2,':',istat,gid_z1b
      print *,1,2,NI0-2,NJ0-3,0,0,':',i0,j0,lni,lnj,hx,hy
   endif

   istat = hgrid_wb_get('z2h',gid_z1b,i0,j0,lni,lnj,hx,hy)
   ok_L = (RMN_IS_OK(istat) .and. gid_z2==gid_z1b .and. &
        all((/i0,j0,lni,lnj,hx,hy/) == (/2,3,NI0-3,NJ0-4,2,3/)))
   call testutils_assert_ok(ok_L,'test_hgrid_wb','get + params2h')
   if (.not.ok_L) then
      print *,gid_z2,':',istat,gid_z1b
      print *,2,3,NI0-3,NJ0-4,2,3,':',i0,j0,lni,lnj,hx,hy
   endif

   istat = hgrid_wb_gmmmeta('z2h',mymeta)
   ok_L = (RMN_IS_OK(istat) .and. &
        all((/mymeta%l(1)%low,mymeta%l(1)%high,mymeta%l(1)%halo,mymeta%l(1)%n/) == (/1-hx,NI0-1,hx,NI0-hx-1/)) .and. &
        all((/mymeta%l(2)%low,mymeta%l(2)%high,mymeta%l(2)%halo,mymeta%l(2)%n/) == (/1-hy,NJ0-1,hy,NJ0-hy-1/)) .and. &
        all((/mymeta%l(3)%low,mymeta%l(3)%high,mymeta%l(3)%halo,mymeta%l(3)%n/) == (/0,0,0,0/)))
   call testutils_assert_ok(ok_L,'test_hgrid_wb','gmmmeta 2d')
   if (.not.ok_L) then
      print *,'l1',mymeta%l(1)
      print *,'  ',1-hx,NI0-1,hx,hx,NI0-hx-1
      print *,'l2',mymeta%l(2)
      print *,'  ',1-hy,NJ0-1,hy,hy,NJ0-hy-1
      print *,'l3',mymeta%l(3)
   endif

   istat = hgrid_wb_gmmmeta('z2h',mymeta,11)
   ok_L = (RMN_IS_OK(istat) .and. &
        all((/mymeta%l(1)%low,mymeta%l(1)%high,mymeta%l(1)%halo,mymeta%l(1)%n/) == (/1-hx,NI0-1,hx,NI0-hx-1/)) .and. &
        all((/mymeta%l(2)%low,mymeta%l(2)%high,mymeta%l(2)%halo,mymeta%l(2)%n/) == (/1-hy,NJ0-1,hy,NJ0-hy-1/)) .and. &
        all((/mymeta%l(3)%low,mymeta%l(3)%high,mymeta%l(3)%halo,mymeta%l(3)%n/) == (/1,11,0,11/)))
   call testutils_assert_ok(ok_L,'test_hgrid_wb','gmmmeta 3d')
   if (.not.ok_L) then
      print *,'l1',mymeta%l(1)
      print *,'  ',1-hx,NI0-1,hx,hx,NI0-hx-1
      print *,'l2',mymeta%l(2)
      print *,'  ',1-hy,NJ0-1,hy,hy,NJ0-hy-1
      print *,'l3',mymeta%l(3)
      print *,'  ',1,11,0,0,11
   endif

   ! ---------------------------------------------------------------------
   return
end subroutine test_hgrid_wb
