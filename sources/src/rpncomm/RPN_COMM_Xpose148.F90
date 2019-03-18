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

	SUBROUTINE RPN_COMM_Xpose148(npe,pecomm, &
        za,min1,max1,n1g,min2,max2,n2l,min3,max3, &
        n3g,zb_8,a_8,b_8)
!
	use rpn_comm
!	forward transpose, za to zb
!	gather first dimension into processor,
!	distribute last dimension
!	result array has gathered index as last dimension
!	(last dimension of arrays is always in-processor)
!
	implicit none

	integer npe,pecomm
	integer min1,max1,n1g,n2l,min2,max2,min3,max3,n3g
	real za(min1:max1,min2:max2,n3g)
	real*8 zb_8(min2:max2,min3:max3,n1g)
	real*8 a_8,b_8
!	integer, allocatable :: za(:,:,:,:)
!	integer, allocatable :: zb(:,:,:,:)
!
!	real*8, allocatable ::  za8(:,:,:,:)
!	real*8, allocatable ::  zb8(:,:,:,:)


!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
!	integer, dimension(size,n2,min3:max3,n1partiel,npe) :: ta
!	real*8, dimension(size/2,n2,min3:max3,n1partiel,npe) :: ta8
	real, allocatable, dimension(:) :: ta,tb
	integer, dimension(npe) :: icount,idispl,scount,sdispl
	integer, dimension(npe) :: kcount,kdispl,rcount,rdispl
	integer istart,iend,kstart,kend,level,proc
	integer index0, index1, index2, index3
!
	
	integer i,j,k,ierr, myproc,i0
	integer nrecv,nsend, nkl, nil
	integer RPN_COMM_limit
!

!	allocate(za(size,min1:max1,n2,n3g))
!	allocate(zb(size,n2,min3:max3,n1g))
!	allocate(za8(size/2,min1:max1,n2,n3g))
!	allocate(zb8(size/2,n2,min3:max3,n1g))
!	za8 => za
!	zb8 => zb

!
	if(npe.eq.1)then
	  do k=1,n3g
	  do j=1,n2l
	     do i=1,n1g
	        zb_8(j,k,i)=a_8*dble(za(i,j,k))+b_8
	     enddo
	  enddo
	  enddo
	  return
	endif

	if(pecomm.eq.pe_mycol) then
	   myproc = pe_mey
	else
	   myproc = pe_mex
	endif


	ierr = RPN_COMM_limit(myproc,npe,1,n1g,istart,iend,icount,idispl)
	ierr = RPN_COMM_limit(myproc,npe,1,n3g,kstart,kend,kcount,kdispl)
	nil = icount(myproc+1)
	nkl = kcount(myproc+1)

	sdispl(1) = 0
	scount(1) = nil*n2l*kcount(1)
	rdispl(1) = 0
	rcount(1) = icount(1)*n2l*nkl

	do i = 2,npe
	   scount(i) = nil*n2l*kcount(i)
	   sdispl(i) = sdispl(i-1) + scount(i-1)
	   rcount(i) = icount(i)*n2l*nkl
	   rdispl(i) = rdispl(i-1) + rcount(i-1)
	enddo
	
	nsend = sdispl(npe)+scount(npe)
	nrecv = rdispl(npe)+rcount(npe)

	allocate(ta(nsend), tb(nrecv))


	i0=0
	do k=1,n3g
	do j=1,n2l
	do i=1,nil
	   i0=i0+1
!	   if(pe_me.eq.0) write(*,*) 'tata1',za(i,j,k),i,j,k
	   ta(i0) = za(i,j,k)
	enddo
	enddo
	enddo

	call mpi_alltoallv(ta,scount,sdispl, MPI_REAL, &
                          tb,rcount,rdispl, MPI_REAL, &
                          pecomm, ierr)

	level=0
	
	do proc = 1,npe
	   i0 = idispl(proc)
	   nil = icount(proc)
	   do k =1,nkl
	      do j=1,n2l-3,4
              do i=1,nil
		 index0=i+level+(j-1)*nil
		 index1=i+level+(j)*nil
		 index2=i+level+(j+1)*nil
		 index3=i+level+(j+2)*nil
		 zb_8(j,k,i+i0)   = a_8*dble(tb(index0))+b_8
		 zb_8(j+1,k,i+i0) = a_8*dble(tb(index1))+b_8
		 zb_8(j+2,k,i+i0) = a_8*dble(tb(index2))+b_8
		 zb_8(j+3,k,i+i0) = a_8*dble(tb(index3))+b_8
!	   if(pe_me.eq.0) write(*,*) 'tata2',zb_8(j,k,i+i0),j,k,i+i0
	      enddo
	      enddo
	      do j=j,n2l
	      do i=1,nil
		 index0= i + level + (j-1)*nil
		 zb_8(j,k,i+i0) = a_8*dble(tb(index0))+b_8
	      enddo
	      enddo
	      level = level + nil* n2l
	   enddo
	enddo
	deallocate(ta,tb)

	return
	end
