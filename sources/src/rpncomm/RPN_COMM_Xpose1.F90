!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */

	SUBROUTINE RPN_COMM_Xpose1(n3partiel,npe,pecomm, &
      n1partiel,za,min1,max1,n1g,n2,min3,max3,n3g,zb,size)
!
	use rpn_comm
!	forward transpose, za to zb
!	gather first dimension into processor,
!	distribute last dimension
!	result array has gathered index as last dimension
!	(last dimension of arrays is always in-processor)
!
	implicit none

	integer n3partiel,n1partiel,npe,pecomm,size
	integer min1,max1,n1g,n2,min3,max3,n3g
	integer za(size,min1:max1,n2,n3g)
	integer zb(size,n2,min3:max3,n1g)
!	integer, allocatable :: za(:,:,:,:)
!	integer, allocatable :: zb(:,:,:,:)
!
	integer *8 za8(size/2,min1:max1,n2,n3g)
	integer *8 zb8(size/2,n2,min3:max3,n1g)
!	integer*8, allocatable ::  za8(:,:,:,:)
!	integer*8, allocatable ::  zb8(:,:,:,:)
	pointer (za8_,za8)
	pointer (zb8_,zb8)

!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
!	integer, dimension(size,n2,min3:max3,n1partiel,npe) :: ta
!	integer*8, dimension(size/2,n2,min3:max3,n1partiel,npe) :: ta8
	integer, allocatable :: ta(:,:,:,:,:)
	integer*8 :: ta8(size/2,n2,min3:max3,n1partiel,npe)
!
	pointer (ta8_,ta8)
	
	integer i,j,k,iter,i0,ierr,isz,isz2,n3w
	logical odd
!
	za8_ = loc(za)
	zb8_ = loc(zb)
!	allocate(za(size,min1:max1,n2,n3g))
!	allocate(zb(size,n2,min3:max3,n1g))
!	allocate(za8(size/2,min1:max1,n2,n3g))
!	allocate(zb8(size/2,n2,min3:max3,n1g))
!	za8 => za
!	zb8 => zb
	isz2 = size/2
	odd = isz2*2 .ne. size
	if(odd) isz2=size
	n3w = max3-min3+1
!
	if(npe.eq.1)then
	  do isz=1,isz2
	  do k=1,n3g
	  do j=1,n2
	    if(odd) then
	      do i=1,n1g
	        zb(isz,j,k,i)=za(isz,i,j,k)
	      enddo
	    else
!VDIR NODEP
	      do i=1,n1g
	        zb8(isz,j,k,i)=za8(isz,i,j,k)
	      enddo
	    endif
	  enddo
	  enddo
	  enddo
	  return
	endif

	allocate(ta(size,n2,n3w,n1partiel,npe))
!	allocate(ta8(size/2,n2,n3w,n1partiel,npe))
	ta8_=loc(ta)

	do isz=1,isz2
	i0 = 0
	do k=1,n3g,n3partiel
	 i0 = 1+i0
	 do i=1,n1partiel
	  if(odd) then
!VDIR NODEP
	    do j=1,n2*min(n3partiel,n3g+1-k)
	      ta(isz,j,min3,i,i0)=za(isz,i,j,k)
	    enddo
	  else
!VDIR NODEP
	    do j=1,n2*min(n3partiel,n3g+1-k)
	      ta8(isz,j,min3,i,i0)=za8(isz,i,j,k)
	    enddo
	  endif
	 enddo
	enddo
	enddo
!
!	call tmg_start(97,'COMM XPOSE1')
	call MPI_ALLTOALL( &
      ta,size*n1partiel*n2*n3w,MPI_INTEGER, &
      zb,size*n1partiel*n2*n3w,MPI_INTEGER, &
      pecomm,ierr)
!	call tmg_stop(97)
!
	deallocate(ta)
!
	return
	end
