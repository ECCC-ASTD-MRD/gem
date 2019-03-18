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
module RPN_COMM_transpose_utils
  use rpn_comm
  implicit none

  private :: forward, backward

  contains
!
  function forward(direction) result(yesno)
  logical :: yesno
  integer, intent(IN) :: direction
  yesno = direction == RPN_COMM_FORWARD_X .or. direction == RPN_COMM_FORWARD_Y
  return
  end function forward
  function backward(direction) result(yesno)
  logical :: yesno
  integer, intent(IN) :: direction
  yesno = direction == RPN_COMM_BACKWARD_X .or. direction == RPN_COMM_BACKWARD_Y
  return
  end function backward
!
!       forward direction : transpose a into b
!       backward direction : transpose b into a
  subroutine transpose_with_halo_a4_b4(a,mini,maxi,b,ni,nj,block,forward)
  implicit none
  integer, intent(INOUT), dimension(mini:maxi,nj)  :: a  ! a may have a halo along first dimension
  integer, intent(INOUT), dimension(nj,ni) :: b
  integer, intent(IN) :: ni,nj,mini,maxi
  integer, intent(IN) :: block
  logical, intent(IN) :: forward
  !
  integer i,j,i0,in,j0,jn
  !
  if (forward) then
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
  end subroutine transpose_with_halo_a4_b4

  SUBROUTINE RPN_COMM_tr_44(za,min1,max1,n1p,n1g,n2,n3p,n3g,zb,npe,pecomm,mype,direction,min3o,max3o,ierr)
  use rpn_comm
  implicit none
  integer, intent(IN) :: min1,max1,n1p,n1g,n2,n3p,n3g
  integer, intent(IN) :: direction,npe,pecomm,mype
  integer, intent(OUT) :: min3o,max3o,ierr
  integer, intent(INOUT) :: za(min1:max1,n2,n3g)
  integer, intent(INOUT) :: zb(n2,n3p,n1g)
!
  integer, dimension(npe) :: sdispls,rdispls,scnts,rcnts
  integer,dimension(n2,0:n1p*n3g-1) :: ta
  integer :: i, n12
!
  scnts = n3p   ! start with "old" decomposition (last tile shorter)
  scnts(npe) = n3g - (npe-1)*n3p
  n12 = n1p * n2
  if(scnts(npe) <= 0 )then  ! "old" failed, try "new" decomposition (last tiles one shorter)
    scnts = n3g / npe
    scnts(1:mod(n3g,npe)) = scnts(1:mod(n3g,npe)) + 1
  endif
  rcnts = scnts(mype)   ! receive counts are constant for a PE
  sdispls(1) = 0        ! displacements = integral of counts
  rdispls(1) = 0
  do i=2,npe
    sdispls(i)=sdispls(i-1)+scnts(i-1)
    rdispls(i)=rdispls(i-1)+rcnts(i-1)
  enddo
!  if(sdispls(npe)+scnts(npe) /= n3g) OUCH ! can't happen !!
  if (forward(direction)) then  ! forward data transpose, ZA(ni,n2,n3g) to TA(n2,n3p(pe),n1p,npe) to ZB(n2,n3p(pe),n1G)
    do i=1,npe  ! transpose sub block before exchange
      call transpose_with_halo_a4_b4(za(min1,1,1+sdispls(i)),min1,max1,ta(1,n1p*sdispls(i)),n1p,n2*scnts(i),128,.true.)
    enddo
    call MPI_alltoallv(ta,scnts*n12,sdispls*n12,MPI_INTEGER,zb,rcnts*n12,rdispls*n12,MPI_INTEGER,pecomm,ierr)
    min3o = sdispls(mype) + 1   ! return this PE's range in n3g
    max3o = sdispls(mype) + scnts(mype)
  else          ! backward data transpose, ZB(n2,n3p(pe),n1G) to TA to ZA(ni,n2,n3g)
    call MPI_alltoallv(zb,scnts*n12,sdispls*n12,MPI_INTEGER,ta,rcnts*n12,rdispls*n12,MPI_INTEGER,pecomm,ierr)
    do i=1,npe  ! transpose sub block after exchange
      call transpose_with_halo_a4_b4(za(min1,1,1+sdispls(i)),min1,max1,ta(1,n1p*sdispls(i)),n1p,n2*scnts(i),128,.false.)
    enddo
  endif
  return
  end SUBROUTINE RPN_COMM_tr_44

  SUBROUTINE RPN_COMM_transpose_44_d(za,min1,max1,n2,min3,max3,decomp1,zb,decomp3,forward,status)
  use rpn_comm
  implicit none
  integer, intent(IN) :: min1,max1,n2,min3,max3
  integer, intent(IN) :: decomp1   ! tag for decomposition along min1,max1
  integer, intent(IN) :: decomp3   ! tag for decomposition along min3,max3
  logical, intent(IN) :: forward
  integer, intent(OUT) :: status
  integer, intent(INOUT) :: za(min1:max1,n2,*)   ! * is n3g
  integer, intent(INOUT) :: zb(n2,min3:max3,*)   ! * is n1g
  status = -1
  end SUBROUTINE RPN_COMM_transpose_44_d

  SUBROUTINE RPN_COMM_transpose_44(za,min1,max1,n1g,n2,min3,max3,n3g,zb,direction,min3o,max3o,status)
  use rpn_comm
  implicit none
  integer, intent(IN) :: min1,max1,n1g,n2,min3,max3,n3g
  integer, intent(IN) :: direction
  integer, intent(OUT) :: status,min3o,max3o
  integer, intent(INOUT) :: za(min1:max1,n2,n3g)
  integer, intent(INOUT) :: zb(n2,min3:max3,n1g)
!
  integer n3p,n1p,npe,pecomm,ierr,mype,lstatus
!
  lstatus = RPN_COMM_OK
  if (forward(direction)) then
    npe=pe_nx
    pecomm=pe_myrow
    mype=pe_mex
  else if (backward(direction)) then
    npe=pe_ny
    pecomm=pe_mycol
    mype=pe_mey
  else
    lstatus = RPN_COMM_ERROR
    return
  endif
  call MPI_allreduce(lstatus,status,1,MPI_INTEGER,MPI_LOR,pecomm,ierr)
  if(status /= RPN_COMM_OK) return
  if (npe == 1) then  ! only one PE in this direction
    call transpose_with_halo_a4_b4(za,min1,max1,zb,n1g,n2*(max3-min3+1),128,forward(direction))
  else
    n1p = (n1g+npe-1)/npe
    n3p = (n3g+npe-1)/npe
    call RPN_COMM_tr_44(za,min1,max1,n1p,n1g,n2,n3p,n3g,zb,npe,pecomm,mype,direction,min3o,max3o,ierr)
  endif
  status = RPN_COMM_OK
  if(ierr /= MPI_SUCCESS) status = RPN_COMM_ERROR
  return
  end SUBROUTINE RPN_COMM_transpose_44
end module RPN_COMM_transpose_utils
