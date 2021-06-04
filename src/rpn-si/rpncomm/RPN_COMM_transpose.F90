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
        SUBROUTINE RPN_COMM_transpose_44(za,min1,max1,n1g,n2,min3,max3,n3g,zb,direction)
        use rpn_comm
        implicit none
        integer, intent(IN) :: min1,max1,n1g,n2,min3,max3,n3g
        integer, intent(IN) :: direction
        integer, intent(INOUT) :: za(min1:max1,n2,n3g)
        integer, intent(INOUT) :: zb(n2,min3:max3,n1g)
!
        integer n3p,n1p,npe,pecomm,ierr,mype,i,n23
        integer, dimension(:,:,:), allocatable :: ta
        integer, dimension(:,:), allocatable :: ndim
        integer, dimension(:), allocatable :: sdispls,rdispls
        integer, dimension(:), allocatable :: scnts,rcnts
!
        if (direction == RPN_COMM_FORWARD_X .or. direction == RPN_COMM_BACKWARD_X ) then
          npe=pe_nx
          pecomm=pe_myrow
          mype=pe_mex
        else if (direction == RPN_COMM_FORWARD_Y .or. direction == RPN_COMM_BACKWARD_Y ) then
          npe=pe_ny
          pecomm=pe_mycol
          mype=pe_mey
        else
          return
        endif
        if (npe == 1) then  ! only one PE in this direction
          call transpose_44(za,min1,max1,zb,n1g,n2*(max3-min3+1),128,direction)
        else
          allocate(ndim(2,npe))
          allocate(scnts(npe),sdispls(npe),rcnts(npe),rdispls(npe))
          call MPI_allgather( (/ n2,max3-min3+1 /), 2, MPI_INTEGER, ndim,2,MPI_INTEGER,pecomm,ierr)
          n1p = (n1g+npe-1)/npe
          n3p = (n3g+npe-1)/npe
          scnts = n3p
          scnts(npe) = n3g - (npe-1)*n3p
          if(scnts(npe) <= 0 )then
            n3p = n3g / npe
            scnts = n3p
            scnts(1+mod(n3g,npe):) = 1 + n3g / npe
          endif
!          if(scnts /= ndim(2,:)) error
          rcnts = scnts(mype)
          sdispls(1) = 0
          rdispls(1) = 0
          do i=2,npe
            sdispls(i)=sdispls(i-1)+scnts(i-1)
            rdispls(i)=rdispls(i-1)+rcnts(i-1)
          enddo
          n23 = n2*n3p
          allocate(ta(n23,n1p,npe))
          if(direction > 0 ) then  ! forward transpose
            do i=1,npe
              call transpose_44(za(min1,1,sdispls(i)),min1,max1,  &
                                ta(1,1,i),n1g,n2*n3p,128,direction)
            enddo
            call MPI_alltoallv(ta,scnts*n23,sdispls*n23,MPI_INTEGER, &
                                zb,rcnts*n23,rdispls*n23,MPI_INTEGER, &
                                pecomm,ierr)
          else                     ! backward transpose
            call MPI_alltoallv(zb,scnts*n23,sdispls*n23,MPI_INTEGER, &
                                ta,rcnts*n23,rdispls*n23,MPI_INTEGER, &
                                pecomm,ierr)
            do i=1,npe
              call transpose_44(za(min1,1,sdispls(i)),min1,max1, &
                                ta(1,1,i),n1g,n2*n3p,128,direction)
            enddo
          endif
          deallocate(ndim,sdispls,rdispls,scnts,rcnts,ta)
        endif
        return

        contains
!       forward direction : transpose a into b
!       backward direction : transpose b into a
        subroutine transpose_44(a,mini,maxi,b,ni,nj,block,direction)
        implicit none
        integer, intent(INOUT), dimension(mini:maxi,nj)  :: a
        integer, intent(INOUT), dimension(nj,ni) :: b
        integer, intent(IN) :: ni,nj,mini,maxi
        integer, intent(IN) :: block, direction
        !
        integer i,j,i0,in,j0,jn
        !
        if (direction == RPN_COMM_FORWARD_X .or. direction == RPN_COMM_FORWARD_Y ) then
          do i0=1,ni,block
            in=min(ni,i0+block-1)
            do j0=1,nj,block
              jn=min(nj,j0+block-1)
              if(iand(in,1)==0 .and. iand(jn,1)==0) then ! in and n are even
                do i=i0,in,2  ! unroll i loop by 2
                do j=j0,jn,2  ! unroll j loop by 2
                  b(j,i)=a(i,j)
                  b(j+1,i)=a(i,j+1)
                  b(j,i+1)=a(i+1,j)
                  b(j+1,i+1)=a(i+1,j+1)
                enddo
                enddo
              else                   ! in or jn is odd (or both)
                do i=i0,in
                do j=j0,jn
                  b(j,i)=a(i,j)
                enddo
                enddo
              endif
            enddo
          enddo
        else
          do j0=1,nj,block
            jn=min(nj,j0+block-1)
            do i0=1,ni,block
              in=min(ni,i0+block-1)
              if(iand(in,1)==0 .and. iand(jn,1)==0) then ! in and n are even
                do j=j0,jn,2  ! unroll j loop by 2
                do i=i0,in,2  ! unroll i loop by 2
                  a(i,j)=b(j,i)
                  a(i,j+1)=b(j+1,i)
                  a(i+1,j)=b(j,i+1)
                  a(i+1,j+1)=b(j+1,i+1)
                enddo
                enddo
              else                   ! in or jn is odd (or both)
                do j=j0,jn
                do i=i0,in
                  b(j,i)=a(i,j)
                enddo
                enddo
              endif
            enddo
          enddo
        endif
        return
        end subroutine transpose_44

        end SUBROUTINE RPN_COMM_transpose_44
!
      SUBROUTINE RPN_COMM_transpose(za,min1,max1,n1g,n2,min3,max3,n3g,zb,type,size)
      use rpn_comm
      implicit none
      integer min1,max1,n1g,n2,min3,max3,n3g,type,size
      integer za(size,min1:max1,n2,n3g)
      integer zb(size,n2,min3:max3,n1g)
!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
      integer n3partiel,n1partiel,npe,pecomm
!
      if(abs(type).eq.1) then   ! transpose along X
        npe=pe_nx
        pecomm=pe_myrow
      else                      ! transpose along Y
        npe=pe_ny
        pecomm=pe_mycol
      endif
      n3partiel=(n3g+npe-1)/npe
      n1partiel=(n1g+npe-1)/npe
!
!	check that min1>0, max1>=n1partiel
!	check that min3>0, max3>=n3partiel
!	check that size = 1 or 2 (integer/real, real*8)
!
      if(type.gt.0) then  ! forward transpose
!
!	  call tmg_start(96,'RPN_COMM_Xpose1')
        call RPN_COMM_Xpose1(n3partiel,npe,pecomm,n1partiel, za,min1,max1,n1g,n2,min3,max3,n3g,zb,size)
!	  call tmg_stop(96)
!
      else ! backward transpose
!
!	  call tmg_start(98,'RPN_COMM_Xpose2')
        call RPN_COMM_Xpose2(n3partiel,npe,pecomm,n1partiel,za,min1,max1,n1g,n2,min3,max3,n3g,zb,size)
!	  call tmg_stop(98)
!
      endif
!
      return
      end SUBROUTINE RPN_COMM_transpose
