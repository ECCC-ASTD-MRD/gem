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
subroutine test_bus_fold()
   use testutils
   use bus_fold_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2012-02
   !@/
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   integer,parameter :: NI=4, NJ=3, NK=5, HX=1, HY=2
   integer,parameter :: BNI1=NI, BNI2=2*NI
   integer,parameter :: BHXY3 = -1,BHXY4 = 1
   integer,parameter :: NVAR2D=5,NVAR3D=6
   real,pointer,dimension(:,:,:) :: v3d1,v3d1o
   real,pointer,dimension(:,:) :: v3d1k,v2d1,v2d1o
   real,pointer,dimension(:,:) :: b1,b2,b3,b4,b5
   integer :: bni3,bni4,bnj1,bnj2,bnj3,bnj4,istat,i,j,k,i1,j1,k0,kn,size1
   integer :: idx_v3d1k,idx_v3d1,idx_v2d1,ni3,nj3,ni4,nj4,ij0n3(4),ij0n4(4),i0,j0,in,jn,n,row
   character(len=256) :: msg_S
   logical :: ok_L
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_bus_fold')
   call msg(MSG_DEBUG,'(test_bus_fold) BEGIN')

   bnj1 = ceiling(float(NI*NJ)/float(BNI1))
   bnj2 = ceiling(float(NI*NJ)/float(BNI2))

   ij0n3(:) = (/1-BHXY3,1-BHXY3,NI+BHXY3,NJ+BHXY3/)
   ni3 = ij0n3(3) - ij0n3(1) + 1
   nj3 = ij0n3(4) - ij0n3(2) + 1
   bni3 = BNI2
   bnj3 = ceiling(float(ni3*nj3)/float(bni3))

   ij0n4(:) = (/1-BHXY4,1-BHXY4,NI+BHXY4,NJ+BHXY4/)
   ni4 = ij0n4(3) - ij0n4(1) + 1
   nj4 = ij0n4(4) - ij0n4(2) + 1
   bni4 = BNI2
   bnj4 = ceiling(float(ni4*nj4)/float(bni4))

   write(msg_S,'(a,3i4,a,2i4,a,2i4)') '(test_bus_fold) allocate:',NI,NJ,NK,';',BNI1,bnj1,';',BNI2,bnj2
   call msg(MSG_DEBUG,msg_S)
   allocate( &
        v3d1(1-HX:NI+HX,1-HY:NJ+HY,NK), &
        v3d1o(1-HX:NI+HX,1-HY:NJ+HY,NK), &
        v2d1(1-HX:NI+HX,1-HY:NJ+HY), &
        v2d1o(1-HX:NI+HX,1-HY:NJ+HY), &
        v3d1k(1-HX:NI+HX,1-HY:NJ+HY), &
        b1(BNI1*(NVAR2D+NVAR3D*NK),bnj1), & !- same shape & no halo as 2d
        b2(BNI2*(NVAR2D+NVAR3D*NK),bnj2), & !- diff shape & no halo as 2d
        b3(bni3*(NVAR2D+NVAR3D*NK),bnj3), & !- diff shape & subset  as 2d
        b4(bni4*(NVAR2D+NVAR3D*NK),bnj4), & !- diff shape & w/ halo as 2d
        b5(bni4*(NVAR2D+NVAR3D*(NK-1)),bnj4), & !- diff shape & w/ halo as 2d
        stat=istat)
   write(msg_S,'(a,i4)') '(test_bus_fold) allocate status:',istat
   call msg(MSG_DEBUG,msg_S)

   v2d1o = -1.
   v3d1o = -1.
   do j=1,NJ
      j1 = j-1
      do i=1,NI
         i1 = i-1
         v2d1o = 6. + float(i1+j1*NI)
         do k=1,NK
            v3d1o = 8. + float(i1+j1*NI+k*NI*NJ)
         enddo
      enddo
   enddo

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_fold_set')
   b1 = -1.
   b2 = -1.
   b3 = -1.
   b4 = -1.
   istat = bus_fold_set('b1',b1,1,NI,1,NJ,1,NK,BNI1)
   istat = min(bus_fold_set('b2',b2,1,NI,1,NJ,1,NK,BNI2),istat)
   istat = min(bus_fold_set('b3',b3,ij0n3(1),ij0n3(3),ij0n3(2),ij0n3(4),1,NK,bni3),istat)
   istat = min(bus_fold_set('b4',b4,ij0n4(1),ij0n4(3),ij0n4(2),ij0n4(4),1,NK,bni4),istat)
   istat = min(bus_fold_set('b5',b5,ij0n4(1),ij0n4(3),ij0n4(2),ij0n4(4),2,NK,bni4),istat)
   call testutils_assert_eq(RMN_IS_OK(istat),.true.,'fold_set')

   k0 = 1
   kn = NK

   !---- bus with same shape & no halo as 2d/3d

   idx_v3d1k = 1
   size1 = BNI1
   idx_v3d1  = idx_v3d1k + size1
   size1 = BNI1*NK
   idx_v2d1  = idx_v3d1 + size1

   v2d1 = -1.
   v3d1 = -1.
   v3d1k = -1.

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_fold b1')
   istat = bus_fold(v2d1o,'b1',idx_v2d1)
   istat = bus_fold(v3d1o,'b1',idx_v3d1,k0,kn)
   istat = bus_fold(v3d1o,'b1',idx_v3d1k,kn,kn)

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_unfold b1')
   istat = bus_unfold(v2d1,'b1',idx_v2d1)
   istat = bus_unfold(v3d1,'b1',idx_v3d1,k0,kn)
   istat = bus_unfold(v3d1k,'b1',idx_v3d1k)

   ok_L = (all(v2d1(1:NI,1:NJ) == v2d1o(1:NI,1:NJ)))
   call testutils_assert_eq(ok_L,.true.,'b1 v2d1')
   if (.not.ok_L) then
      print *,'got:',v2d1(1:NI,1:NJ)
      print *,'exp:',v2d1o(1:NI,1:NJ)
   endif
   ok_L = (all(v3d1(1:NI,1:NJ,1:NK) == v3d1o(1:NI,1:NJ,1:NK)))
   call testutils_assert_eq(ok_L,.true.,'b1 v3d1')
   if (.not.ok_L) then
      print *,'got:',v3d1(1:NI,1:NJ,1)
      print *,'exp:',v3d1o(1:NI,1:NJ,1)
   endif
   ok_L = (all(v3d1k(1:NI,1:NJ) == v3d1o(1:NI,1:NJ,NK)))
   call testutils_assert_eq(ok_L,.true.,'b1 v3d1k')
   if (.not.ok_L) then
      print *,'got:',v3d1k(1:NI,1:NJ)
      print *,'exp:',v3d1o(1:NI,1:NJ,NK)
   endif

   !---- bus with diff shape & no halo as 2d/3d

   idx_v3d1k = 1
   size1 = BNI2
   idx_v3d1  = idx_v3d1k + size1
   size1 = BNI2*NK
   idx_v2d1  = idx_v3d1 + size1

   v2d1 = -1.
   v3d1 = -1.
   v3d1k = -1.

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_fold b2')
   istat = bus_fold(v2d1o,'b2',idx_v2d1)
   istat = bus_fold(v3d1o,'b2',idx_v3d1,k0,kn)
   istat = bus_fold(v3d1o,'b2',idx_v3d1k,kn,kn)

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_unfold b2')
   istat = bus_unfold(v2d1,'b2',idx_v2d1)
   istat = bus_unfold(v3d1,'b2',idx_v3d1,k0,kn)
   istat = bus_unfold(v3d1k,'b2',idx_v3d1k)

   ok_L = (all(v2d1(1:NI,1:NJ) == v2d1o(1:NI,1:NJ)))
   call testutils_assert_eq(ok_L,.true.,'b2 v2d1')
   if (.not.ok_L) then
      print *,'got:',v2d1(1:NI,1:NJ)
      print *,'exp:',v2d1o(1:NI,1:NJ)
   endif
   ok_L = (all(v3d1(1:NI,1:NJ,1:NK) == v3d1o(1:NI,1:NJ,1:NK)))
   call testutils_assert_eq(ok_L,.true.,'b2 v3d1')
   if (.not.ok_L) then
      print *,'got:',v3d1(1:NI,1:NJ,1)
      print *,'exp:',v3d1o(1:NI,1:NJ,1)
   endif
   ok_L = (all(v3d1k(1:NI,1:NJ) == v3d1o(1:NI,1:NJ,NK)))
   call testutils_assert_eq(ok_L,.true.,'b2 v3d1k')
   if (.not.ok_L) then
      print *,'got:',v3d1k(1:NI,1:NJ)
      print *,'exp:',v3d1o(1:NI,1:NJ,NK)
   endif

   !---- bus with diff shape & subset as 2d/3d
   i0 = ij0n3(1)
   in = ij0n3(3)
   j0 = ij0n3(2)
   jn = ij0n3(4)

   v2d1o = -1.
   v3d1o = -1.
   do j=j0,jn
      j1 = j-j0
      do i=i0,in
         i1 = i-i0
         v2d1o = 6. + float(i1+j1*ni3)
         do k=1,NK
            v3d1o = 8. + float(i1+j1*ni3+k*ni3*nj3)
         enddo
      enddo
   enddo

   idx_v3d1k = 1
   size1 = bni3
   idx_v3d1  = idx_v3d1k + size1
   size1 = bni3*NK
   idx_v2d1  = idx_v3d1 + size1

   v2d1 = -1.
   v3d1 = -1.
   v3d1k = -1.

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_fold b3')
   istat = bus_fold(v2d1o,'b3',idx_v2d1)
   istat = bus_fold(v3d1o,'b3',idx_v3d1,k0,kn)
   istat = bus_fold(v3d1o,'b3',idx_v3d1k,kn,kn)

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_unfold b3')
   istat = bus_unfold(v2d1,'b3',idx_v2d1)
   istat = bus_unfold(v3d1,'b3',idx_v3d1,k0,kn)
   istat = bus_unfold(v3d1k,'b3',idx_v3d1k)

   ok_L = (all(v2d1(i0:in,j0:jn) == v2d1o(i0:in,j0:jn)))
   call testutils_assert_eq(ok_L,.true.,'b3 v2d1')
   if (.not.ok_L) then
      print *,'got:',v2d1(i0:in,j0:jn)
      print *,'exp:',v2d1o(i0:in,j0:jn)
   endif
   ok_L = (all(v3d1(i0:in,j0:jn,1:NK) == v3d1o(i0:in,j0:jn,1:NK)))
   call testutils_assert_eq(ok_L,.true.,'b3 v3d1')
   if (.not.ok_L) then
      print *,'got:',v3d1(i0:in,j0:jn,1)
      print *,'exp:',v3d1o(i0:in,j0:jn,1)
   endif
   ok_L = (all(v3d1k(i0:in,j0:jn) == v3d1o(i0:in,j0:jn,NK)))
   call testutils_assert_eq(ok_L,.true.,'b3 v3d1k')
   if (.not.ok_L) then
      print *,'got:',v3d1k(i0:in,j0:jn)
      print *,'exp:',v3d1o(i0:in,j0:jn,NK)
   endif

   !---- bus with diff shape & halo as 2d/3d
   i0 = ij0n4(1)
   in = ij0n4(3)
   j0 = ij0n4(2)
   jn = ij0n4(4)

   v2d1o = -1.
   v3d1o = -1.
   do j=j0,jn
      j1 = j-j0
      do i=i0,in
         i1 = i-i0
         v2d1o = 6. + float(i1+j1*ni4)
         do k=1,NK
            v3d1o = 8. + float(i1+j1*ni4+k*ni4*nj4)
         enddo
      enddo
   enddo

   idx_v3d1k = 1
   size1 = bni4
   idx_v3d1  = idx_v3d1k + size1
   size1 = bni4*NK
   idx_v2d1  = idx_v3d1 + size1

   v2d1 = -1.
   v3d1 = -1.
   v3d1k = -1.

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_fold b4')
   istat = bus_fold(v2d1o,'b4',idx_v2d1)
   istat = bus_fold(v3d1o,'b4',idx_v3d1,k0,kn)
   istat = bus_fold(v3d1o,'b4',idx_v3d1k,kn,kn)

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_unfold b4')
   istat = bus_unfold(v2d1,'b4',idx_v2d1)
   istat = bus_unfold(v3d1,'b4',idx_v3d1,k0,kn)
   istat = bus_unfold(v3d1k,'b4',idx_v3d1k)

   ok_L = (all(v2d1(i0:in,j0:jn) == v2d1o(i0:in,j0:jn)))
   call testutils_assert_eq(ok_L,.true.,'b4 v2d1')
   if (.not.ok_L) then
      print *,'got:',v2d1(i0:in,j0:jn)
      print *,'exp:',v2d1o(i0:in,j0:jn)
   endif
   ok_L = (all(v3d1(i0:in,j0:jn,1:NK) == v3d1o(i0:in,j0:jn,1:NK)))
   call testutils_assert_eq(ok_L,.true.,'b4 v3d1')
   if (.not.ok_L) then
      print *,'got:',v3d1(i0:in,j0:jn,1)
      print *,'exp:',v3d1o(i0:in,j0:jn,1)
   endif
   ok_L = (all(v3d1k(i0:in,j0:jn) == v3d1o(i0:in,j0:jn,NK)))
   call testutils_assert_eq(ok_L,.true.,'b4 v3d1k')
   if (.not.ok_L) then
      print *,'got:',v3d1k(i0:in,j0:jn)
      print *,'exp:',v3d1o(i0:in,j0:jn,NK)
   endif

   ok_L = .true.
   n = 0
   row = 1
   do k=1,NK
      do j=j0,jn
         do i=i0,in
            n = n + 1
            if (n > bni4*NK) then
               row = row+1
               n = 0
            endif
            if (v3d1o(i,j,k) /= b4(idx_v3d1+n-1,row)) then
               print *,i,j,k,v3d1o(i,j,k),' != ',idx_v3d1+n-1,row,b4(idx_v3d1+n-1,row)
               ok_L = .false.
            endif
         enddo
      enddo
   enddo
   call testutils_assert_eq(ok_L,.true.,'b4 v3d1 bus values')

   !---- bus with diff shape & halo as 2d/3d & notop (k0=2)
   i0 = ij0n4(1)
   in = ij0n4(3)
   j0 = ij0n4(2)
   jn = ij0n4(4)

   k0 = 2

   v2d1o = -1.
   v3d1o = -1.
   do j=j0,jn
      j1 = j-j0
      do i=i0,in
         i1 = i-i0
         v2d1o = 6. + float(i1+j1*ni4)
         do k=2,NK
            v3d1o = 8. + float(i1+j1*ni4+k*ni4*nj4)
         enddo
      enddo
   enddo

   idx_v3d1k = 1
   size1 = bni4
   idx_v3d1  = idx_v3d1k + size1
   size1 = bni4*(kn-k0+1)
   idx_v2d1  = idx_v3d1 + size1

   v2d1 = -1.
   v3d1 = -1.
   v3d1k = -1.

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_fold b5')
   istat = bus_fold(v2d1o,'b5',idx_v2d1)
   istat = bus_fold(v3d1o,'b5',idx_v3d1,k0,kn)
   istat = bus_fold(v3d1o,'b5',idx_v3d1k,kn,kn)

   call msg(MSG_DEBUG,'(test_bus_fold) to bus_unfold b5')
   istat = bus_unfold(v2d1,'b5',idx_v2d1)
   istat = bus_unfold(v3d1,'b5',idx_v3d1,k0,kn)
   istat = bus_unfold(v3d1k,'b5',idx_v3d1k)

   ok_L = (all(v2d1(i0:in,j0:jn) == v2d1o(i0:in,j0:jn)))
   call testutils_assert_eq(ok_L,.true.,'b5 v2d1')
   if (.not.ok_L) then
      print *,'got:',v2d1(i0:in,j0:jn)
      print *,'exp:',v2d1o(i0:in,j0:jn)
   endif
   ok_L = (all(v3d1(i0:in,j0:jn,k0:kn) == v3d1o(i0:in,j0:jn,k0:kn)))
   call testutils_assert_eq(ok_L,.true.,'b5 v3d1')
   if (.not.ok_L) then
      print *,'got:',v3d1(i0:in,j0:jn,k0)
      print *,'exp:',v3d1o(i0:in,j0:jn,k0)
   endif
   ok_L = (all(v3d1k(i0:in,j0:jn) == v3d1o(i0:in,j0:jn,NK)))
   call testutils_assert_eq(ok_L,.true.,'b5 v3d1k')
   if (.not.ok_L) then
      print *,'got:',v3d1k(i0:in,j0:jn)
      print *,'exp:',v3d1o(i0:in,j0:jn,NK)
   endif

   ok_L = .true.
   n = 0
   row = 1
   do k=k0,NK
      do j=j0,jn
         do i=i0,in
            n = n + 1
            if (n > bni4*(kn-k0+1)) then
               row = row+1
               n = 0
            endif
            if (v3d1o(i,j,k) /= b5(idx_v3d1+n-1,row)) then
               print *,i,j,k,v3d1o(i,j,k),' != ',idx_v3d1+n-1,row,b5(idx_v3d1+n-1,row)
               ok_L = .false.
            endif
         enddo
      enddo
   enddo
   call testutils_assert_eq(ok_L,.true.,'b5 v3d1 bus values')


   call msg(MSG_DEBUG,'(test_bus_fold) END')
   ! ---------------------------------------------------------------------
   return
end subroutine test_bus_fold
