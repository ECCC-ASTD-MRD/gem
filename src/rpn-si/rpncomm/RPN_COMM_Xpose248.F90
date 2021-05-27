!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! *
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! *
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! */

	SUBROUTINE RPN_COMM_Xpose248(npe,pecomm, &
        za,min1,max1,n1g,min2,max2,n2l,min3,max3,n3g, &
        zb_8,a_8,b_8)
	use rpn_comm
!
!	backward transpose, zb to za
!	(last dimension of arrays is always in-processor)
!
	implicit none

	integer npe,pecomm
	integer min1,max1,n1g,n2l,min3,max3,n3g,min2,max2
	real za(min1:max1,min2:max2,n3g)
	real*8 zb_8(min2:max2,min3:max3,n1g)
	real*8 a_8,b_8

	real, allocatable, dimension(:) :: ta, tb
	integer nsend, nrecv, myproc
	integer index0, index1, index2, index3, njkl, level, k0, proc
	integer, dimension(npe) :: icount,idispl,scount,sdispl
	integer, dimension(npe) :: kcount,kdispl,rcount,rdispl
	integer RPN_COMM_limit
	integer nil,nkl,istart,iend,kstart,kend
!

!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!

	integer i,j,k,ierr,n3w

	n3w = max3-min3+1
!
!
	if(npe.eq.1)then
	  do k=1,n3g
	  do j=1,n2l
	      do i=1,n1g
	        za(i,j,k)=sngl(a_8*zb_8(j,k,i)+b_8)
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
	scount(1) = icount(1)*n2l*nkl
	rdispl(1) = 0
	rcount(1) = nil*n2l*kcount(1)

	do i = 2,npe
	   scount(i) = icount(i)*n2l*nkl
	   sdispl(i) = sdispl(i-1) + scount(i-1)
	   rcount(i) = nil*n2l*kcount(i)
	   rdispl(i) = rdispl(i-1) + rcount(i-1)
	enddo
!
!
	nsend = sdispl(npe)+scount(npe)
	nrecv = rdispl(npe)+rcount(npe)

	allocate(ta(nsend), tb(nrecv))

	level = 0

	do i = 1, n1g
	do k = 1, nkl
	do j = 1, n2l
	   level = level+1
	   ta(level) = sngl(zb_8(j,k,i))
	enddo
	enddo
	enddo

!	call tmg_start(99,'COMM XPOSE2')
	call MPI_ALLTOALLV( &
            ta,scount,sdispl,MPI_REAL, &
            tb,rcount,rdispl,MPI_REAL, &
            pecomm,ierr)
!	call tmg_stop(99)
!
	level = 0
	
	do proc = 1,npe
	   k0 = kdispl(proc)
	   nkl = kcount(proc)
	   njkl = n2l*nkl
	   do i = 1,nil-3,4
	      do k = 1,nkl
	      do j = 1,n2l
		 index0=j+level+(k-1)*n2l+(i-1)*njkl
		 index1=j+level+(k-1)*n2l+(i)*njkl
		 index2=j+level+(k-1)*n2l+(i+1)*njkl
		 index3=j+level+(k-1)*n2l+(i+2)*njkl
		 za(i   ,j,k+k0) = tb(index0)
		 za(i+1 ,j,k+k0) = tb(index1)
		 za(i+2 ,j,k+k0) = tb(index2)
		 za(i+3 ,j,k+k0) = tb(index3)
	      enddo
	      enddo
	   enddo
	   do i = i,nil
	      do k = 1,nkl
	      do j = 1,n2l
		 index0 = j + level + (k-1)*n2l +(i-1)*njkl
		 za(i,j,k+k0) = tb(index0)
	      enddo
	      enddo
	   enddo
	   level = level + nil*n2l*nkl
	enddo
	      
!
	deallocate(ta,tb)
!
	return
	end
