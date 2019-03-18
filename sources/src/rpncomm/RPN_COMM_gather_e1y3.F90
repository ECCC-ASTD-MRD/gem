!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/

	SUBROUTINE RPN_COMM_gather_e1y3(srca,srcc,srcr,&
     &        min1,max1,min2,max2, n1,n2,n3,dest,size)
	use rpn_comm
	implicit none
	integer min1,max1,min2,max2,n1,n2,n3,size
	real dest(size,6*(2*(max(max(n2*n3,n1*n3),n1*n2)/2)+1)*size)
	real srca(size,min1:max1,min2:max2,n3)
	real srcc(size,min1:max1,min2:max2,n3)
	real srcr(size,min1:max1,min2:max2,n3)
!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
	real, allocatable :: zxij(:,:,:)
!	real, dimension(size,2*((n2*n3)/2)+1,6) :: zxij
!	pointer (zij_,zxij)
	real *8, dimension(size/2,2*((n2*n3)/2)+1,6) :: zxij8
	pointer (zij8_,zxij8)
!
	real, dimension(size,2*((n1*n3)/2)+1,6) :: zyij
	pointer (zyij_,zyij)
	real *8, dimension(size/2,2*((n1*n3)/2)+1,6) :: zyij8
	pointer (zyij8_,zyij8)
!
	real *8 srca8(size/2,min1:max1,min2:max2,n3)
	real *8 srcc8(size/2,min1:max1,min2:max2,n3)
	real *8 srcr8(size/2,min1:max1,min2:max2,n3)
	pointer (srca8_,srca8)
	pointer (srcc8_,srcc8)
	pointer (srcr8_,srcr8)
!
	integer i,j,k,ie,isz,incr,ierr,nwords,mode,direction
	integer size2
!
	nwords=6*(2*((n2*n3)/2)+1)*size
	allocate(zxij(size,2*((n2*n3)/2)+1,6))
	mode=1
	direction=pe_mycol
!
 1	continue

	zij8_=loc(zxij)
	zyij_=zij8_
	zyij8_=zij8_
	srca8_=loc(srca)
	srcc8_=loc(srcc)
	srcr8_=loc(srcr)
!
	size2 = size
	if(mod(size,2).eq.0) size2 = size/2

	if(mode.eq.1)then    ! edge 1 (along x or y )
	  do isz=1,size2
	   incr=0
	   do j=1,n2
	    if(size.eq.size2)then
	     do k=1,n3
	      zxij(isz,incr+k,1) = srca(isz, 1,j,k)
	      zxij(isz,incr+k,4) = srca(isz,n1,j,k)
	      zxij(isz,incr+k,2) = srcc(isz, 1,j,k)
	      zxij(isz,incr+k,5) = srcc(isz,n1,j,k)
	      zxij(isz,incr+k,3) = srcr(isz, 1,j,k)
	      zxij(isz,incr+k,6) = srcr(isz,n1,j,k)
	     enddo
	    else
	     do k=1,n3
	      zxij8(isz,incr+k,1) = srca8(isz, 1,j,k)
	      zxij8(isz,incr+k,4) = srca8(isz,n1,j,k)
	      zxij8(isz,incr+k,2) = srcc8(isz, 1,j,k)
	      zxij8(isz,incr+k,5) = srcc8(isz,n1,j,k)
	      zxij8(isz,incr+k,3) = srcr8(isz, 1,j,k)
	      zxij8(isz,incr+k,6) = srcr8(isz,n1,j,k)
	     enddo
	    endif
	    incr=incr+n3
	   enddo
	  enddo
	else                  !  edge 2 (along x or y)
	  do isz=1,size2
	   incr=0
	   do i=1,n1
	    if(size.eq.size2)then
	     do k=1,n3
	      zyij(isz,incr+k,1) = srca(isz,i, 1,k)
	      zyij(isz,incr+k,4) = srca(isz,i,n2,k)
	      zyij(isz,incr+k,2) = srcc(isz,i, 1,k)
	      zyij(isz,incr+k,5) = srcc(isz,i,n2,k)
	      zyij(isz,incr+k,3) = srcr(isz,i, 1,k)
	      zyij(isz,incr+k,6) = srcr(isz,i,n2,k)
	     enddo
	    else
	     do k=1,n3
	      zyij8(isz,incr+k,1) = srca8(isz,i, 1,k)
	      zyij8(isz,incr+k,4) = srca8(isz,i,n2,k)
	      zyij8(isz,incr+k,2) = srcc8(isz,i, 1,k)
	      zyij8(isz,incr+k,5) = srcc8(isz,i,n2,k)
	      zyij8(isz,incr+k,3) = srcr8(isz,i, 1,k)
	      zyij8(isz,incr+k,6) = srcr8(isz,i,n2,k)
	     enddo
	    endif
	    incr=incr+n3
	   enddo
	  enddo
	endif

	call MPI_ALLGATHER(zxij,nwords,MPI_REAL,&
     &                     dest,nwords,MPI_REAL,&
     &                     direction,ierr)
	deallocate(zxij)
	return

	entry RPN_COMM_gather_e1x3(srca,srcc,srcr,&
     &        min1,max1,min2,max2, n1,n2,n3,dest,size)

	nwords=6*(2*((n2*n3)/2)+1)*size
	allocate(zxij(size,2*((n2*n3)/2)+1,6))
	mode=1
	direction=pe_myrow
	goto 1

	entry RPN_COMM_gather_e2y3(srca,srcc,srcr,&
     &        min1,max1,min2,max2, n1,n2,n3,dest,size)

	nwords=6*(2*((n1*n3)/2)+1)*size
	allocate(zxij(size,2*((n1*n3)/2)+1,6))
	mode=2
	direction=pe_mycol
	goto 1

	entry RPN_COMM_gather_e2x3(srca,srcc,srcr,&
     &        min1,max1,min2,max2, n1,n2,n3,dest,size)

	nwords=6*(2*((n1*n3)/2)+1)*size
	allocate(zxij(size,2*((n1*n3)/2)+1,6))
	mode=2
	direction=pe_myrow
	goto 1

	end
