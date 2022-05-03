#define ARGUMENT_TYPE real
subroutine low_pass_dwt2d_r4(z,ni,nj,nx,ny)
! 2D low pass filter using 53 dwt linear prediction wavelet
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z

  call fwd_linear53_dwt2d_r4(z,ni,nj,nx,ny)   ! forward transform
  call zero_high_dwt2d_r4(z,ni,nj,nx,ny)      ! get rid of high frequency terms
  call inv_linear53_dwt2d_r4(z,ni,nj,nx,ny)   ! inverse transform
end subroutine low_pass_dwt2d_r4

subroutine low_pass_dwt2d_r4_4x4(z,ni,nj,nx,ny)
! 2D low pass filter using 53 dwt linear prediction wavelet
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z
  integer :: ni2, nj2

  ni2 = (ni+1)/2
  nj2 = (nj+1)/2
  call fwd_linear53_dwt2d_r4(z,ni,nj,ni,nj)   ! forward transform
  call zero_high_dwt2d_r4(z,ni,nj,ni,nj)      ! get rid of high frequency terms
  call fwd_linear53_dwt2d_r4(z,ni/2,nj/2,ni*2,nj/2)
  call zero_high_dwt2d_r4(z,ni/2,nj/2,ni*2,nj/2)      ! get rid of high frequency terms
  call inv_linear53_dwt2d_r4(z,ni2,nj2,ni*2,nj2)
  call inv_linear53_dwt2d_r4(z,ni,nj,ni,nj)   ! inverse transform
end subroutine low_pass_dwt2d_r4_4x4

subroutine low_pass_quant_dwt2d_r4(z,ni,nj,nx,ny)
! 2D low pass filter using 53 dwt linear prediction wavelet
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z

  call fwd_linear53_dwt2d_r4(z,ni,nj,nx,ny)   ! forward transform
  call zero_high_dwt2d_r4(z,ni,nj,nx,ny)      ! get rid of high frequency terms
  call quant_low_r4(z,ni,nj,nx,ny)            ! quantize by 2**16 non zero terms
  call inv_linear53_dwt2d_r4(z,ni,nj,nx,ny)   ! inverse transform
end subroutine low_pass_quant_dwt2d_r4

subroutine quant_low_r4(z,ni,nj,nx,ny)
! quantize (2**16 intervals) non zero terms terms in dwt transformed data
! this routine assumes the transform layout used by fwd_linear53_dwt2d_r4 / inv_linear53_dwt2d_r4
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z

  integer :: j, neven
  ARGUMENT_TYPE :: the_min, the_max, the_range, delta

  the_min = z(0,0)
  the_max = z(0,0)
  neven = (ni+1)/2            ! number of even (low frequency) terms in a row
  do j = 1 , nj-1 , 2         ! min and max loop over odd/even row pairs
    ! ignore odd rows
    the_max = max(the_max,maxval(z(0:neven-1,j-1))) ! lower half of even rows
    the_min = min(the_max,minval(z(0:neven-1,j-1))) ! lower half of even rows
  end do
  if(iand(nj,1) == 1) then  ! last row is an even row (odd number of rows)
    the_max = max(the_max,maxval(z(0:neven-1,nj-1))) ! lower half of even rows
    the_min = min(the_max,minval(z(0:neven-1,nj-1))) ! lower half of even rows
  endif

  the_range = the_max - the_min
  if(the_range == 0) the_range = 1.0
  delta = the_range / 65535.0      ! quantification interval
  print *,'DEBUG: the_min, the_max, the_range, delta', the_min, the_max, the_range, delta

  do j = 1 , nj-1 , 2              ! quantize loop over odd/even row pairs
    z(0:neven-1,j-1) = the_min + (delta * nint( (z(0:neven-1,j-1)-the_min) / delta) )
  end do
  if(iand(nj,1) == 1) z(0:neven-1,nj-1) = the_min + (delta * nint( (z(0:neven-1,nj-1)-the_min) / delta) )
!   print *,'DEBUG: end of quant_low_r4'
  return
end subroutine quant_low_r4

subroutine zero_high_dwt2d_r4(z,ni,nj,nx,ny)
! zero out all odd (high frequency) terms in dwt transformed data
! this routine assumes the transform layout used by fwd_linear53_dwt2d_r4 / inv_linear53_dwt2d_r4
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z

  integer :: j, neven

  neven = (ni+1)/2            ! number of even (low frequency) terms in a row
  do j = 1 , nj-1 , 2         ! loop over odd rows
    z(0:ni-1    ,j  ) = 0.0   ! full row zeroed in odd rows
    z(neven:ni-1,j-1) = 0.0   ! odd terms zeroed in even rows (upper half)
  end do
  if(iand(nj,1) == 1) z(neven:ni-1,nj-1) = 0.0  ! last row is an even row (odd number of rows)
  return
end subroutine zero_high_dwt2d_r4

subroutine inv_linear53_dwt2d_r4z(z,ni,nj,nx,ny)   ! 2D inverse transform , along j first, then along i
! in place INVERSE lifting transform using linear prediction wavelet for ARGUMENT_TYPE numbers
! all odd terms are asumed to be zero
  use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z

  ARGUMENT_TYPE, dimension(-1:ni) :: even
  integer :: j, j00, jm1, jm2, jp1
  integer :: nodd, neven

  nodd  = ishft(ni,-1)     ! number of odd terms
  neven = ishft(ni+1,-1)   ! number of even terms (nodd + 1 if ni is odd)
  j00 = 0

  do j = 0 , nj-1 , 2
    jm2 = j - 2           ! not used for j == 0
    jm1 = abs(j - 1)      ! mirror condition at lower boundary
    jp1 = j + 1
    if(jp1 == nj) jp1 = nj - 2                     ! upper mirror boundary condition
    z(0:neven-1,j) = z(0:neven-1,j) - .25 * (z(0:neven-1,jm1) + z(0:neven-1,jp1))  ! unupdate even row j
    if(j > 0) then
      z(0:neven-1,jm1) = z(0:neven-1,jm1) + .5 * (z(0:neven-1,jm2) + z(0:neven-1,j)) ! unpredict odd row jm1
      call inv_linear53_dwt1d_r4z(z(:,j-2))
      call inv_linear53_dwt1d_r4z(z(:,j-1))
      j00 = j
    endif
  end do
  if(iand(nj,1)==0) z(0:neven-1,nj-1) = z(0:neven-1,nj-1) + z(0:neven-1,nj-2) ! last row is odd, unpredict it
  do j = j00 , nj-1
    call inv_linear53_dwt1d_r4z(z(:,j))
  end do
contains
  subroutine inv_linear53_dwt1d_r4z(f)  ! 1D along i transform
    use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
    implicit none
    ARGUMENT_TYPE, intent(INOUT), dimension(0:ni-1) :: f

    integer :: i

    do i=0,nodd-1           ! split into even / odd (assumed zero)
      even(i) = f(i)
    end do
    if(iand(ni,1) /= 0) even(nodd) = f(nodd)      ! one more even values than odd values if ni is odd
    if(iand(ni,1) == 0) even(nodd) = even(nodd-1) ! mirror condition at upper boundary if ni is even
    do i = 0,nodd-1
      f(2*i)   = even(i)
      f(2*i+1) = .5 * (even(i) + even(i+1))
    end do
    if(iand(ni,1) /= 0) f(ni-1) = even(neven-1)   ! odd number of points, one extra even element
  end subroutine inv_linear53_dwt1d_r4z
end subroutine inv_linear53_dwt2d_r4z

subroutine inv_linear53_dwt2d_r4(z,ni,nj,nx,ny)   ! 2D inverse transform , along j first, then along i
! in place INVERSE lifting transform using linear prediction wavelet for ARGUMENT_TYPE numbers
  use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z

  ARGUMENT_TYPE, dimension(-1:ni) :: even, odd
  integer :: j, j00, jm1, jm2, jp1
  integer :: nodd, neven

  nodd  = ishft(ni,-1)     ! number of odd terms
  neven = ishft(ni+1,-1)   ! number of even terms (nodd + 1 if ni is odd)
  j00 = 0

  do j = 0 , nj-1 , 2
    jm2 = j - 2           ! not used for j == 0
    jm1 = abs(j - 1)      ! mirror condition at lower boundary
    jp1 = j + 1
    if(jp1 == nj) jp1 = nj - 2                     ! upper mirror boundary condition
    z(0:ni-1,j) = z(0:ni-1,j) - .25 * (z(0:ni-1,jm1) + z(0:ni-1,jp1))  ! unupdate even row j
    if(j > 0) then
      z(0:ni-1,jm1) = z(0:ni-1,jm1) + .5 * (z(0:ni-1,jm2) + z(0:ni-1,j)) ! unpredict odd row jm1
      call inv_linear53_dwt1d_r4(z(:,j-2))
      call inv_linear53_dwt1d_r4(z(:,j-1))
      j00 = j
    endif
  end do
  if(iand(nj,1)==0) z(0:ni-1,nj-1) = z(0:ni-1,nj-1) + z(0:ni-1,nj-2) ! last row is odd, unpredict it
  do j = j00 , nj-1
    call inv_linear53_dwt1d_r4(z(:,j))
  end do
contains
  subroutine inv_linear53_dwt1d_r4(f)  ! 1D along i transform
    use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
    implicit none
    ARGUMENT_TYPE, intent(INOUT), dimension(0:ni-1) :: f

    integer :: i

    do i=0,nodd-1           ! split into even / odd
      even(i) = f(i)
      odd(i)  = f(neven+i)
    end do
    if(iand(ni,1) /= 0) even(nodd) = f(nodd)   ! one more even values than odd values if ni is odd
    odd(-1)   = odd(0)         ! mirror condition at lower boundary
    odd(nodd) = odd(nodd-1)    ! mirror condition at upper boundary, used only of ni is odd

    even(0) = even(0) - odd(0)
    do i = 1,neven-1  ! unupdate even points
      even(i) = even(i) - .25 * (odd(i-1) + odd(i))
    end do

    if(iand(ni,1) == 0) even(nodd) = even(nodd-1) ! mirror condition at upper boundary if ni is even
    do i = 0,nodd-1  ! unpredict odd points
      odd(i) = odd(i) + .5 * (even(i) + even(i+1))
    end do
    do i = 0,nodd-1
      f(2*i)   = even(i)
      f(2*i+1) = odd(i)
    end do
    if(iand(ni,1) /= 0) f(ni-1) = even(neven-1)  ! odd number of points, one extra even element
  end subroutine inv_linear53_dwt1d_r4
end subroutine inv_linear53_dwt2d_r4

subroutine fwd_linear53_dwt2d_r4z(z,ni,nj,nx,ny)   ! 2D forward transform, along i first, then along j
! in place FORWARD lifting transform using linear prediction wavelet for ARGUMENT_TYPE numbers
! all odd terms in transform will be assumed 0 by inverse transform and are not stored
  use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z

  ARGUMENT_TYPE, dimension(-1:ni) :: even, odd
  integer :: j, j00, jm1, jm2, jp1, jj
  integer :: nodd, neven

  nodd  = ishft(ni,-1)     ! number of odd terms
  neven = ishft(ni+1,-1)   ! number of even terms (nodd + 1 if ni is odd)
  j00 = 0
  do j = 1,nj-1,2
    jm2 = abs(j - 2)
    jm1 = j - 1
    jp1 = j + 1
    if(jp1 == nj) jp1 = nj - 2                              ! upper mirror boundary condition
    do jj = j00,max(j,jp1)
      call fwd_linear53_dwt1d_r4z(z(:,jj))
    end do
    j00 = jp1 + 1   ! odd terms are ignored in the following 2 lines
    z(0:neven-1,  j) = z(0:neven-1,  j) -  .5 * (z(0:neven-1,jm1) + z(0:neven-1,jp1))       ! predict odd rows
    z(0:neven-1,jm1) = z(0:neven-1,jm1) + .25 * (z(0:neven-1,jm2) + z(0:neven-1,  j))       ! update even rows (below odd row)
  end do
  if(mod(nj,2)==1) z(0:neven-1,nj-1) = z(0:neven-1,nj-1) + .5 * z(0:neven-1,nj-2)   ! odd number of rows, last row is en even row
  return
contains
  subroutine fwd_linear53_dwt1d_r4z(f)   ! 1D along i transform
      use, intrinsic :: iso_fortran_env
    implicit none
    ARGUMENT_TYPE, intent(INOUT), dimension(0:ni-1) :: f

    integer :: i

    do i=0,nodd-1           ! split into even / odd
      even(i) = f(2*i)
      odd(i)  = f(2*i+1)
    end do
    if(iand(ni,1) /= 0) then  ! is ni odd ?
      even(nodd) = f(ni-1)      ! one more even values than odd values if ni is odd
    else
      even(nodd) = even(nodd-1) ! mirror condition at upper boundary if ni is even
    endif
    do i = 0, nodd-1           ! predict odd values (and copy updated values into f)
      odd(i) = odd(i) - .5 * (even(i) + even(i+1))
    end do
    odd(-1)   = odd(0)         ! mirror condition at lower boundary
    odd(nodd) = odd(nodd-1)    ! mirror condition at upper boundary, used only of ni is odd
    do  i=0,neven-1            ! update even values (and copy updated values into f)
      f(i) = even(i) + .25 * (odd(i) + odd(i-1))
    end do
  end subroutine fwd_linear53_dwt1d_r4z
end subroutine fwd_linear53_dwt2d_r4z

subroutine fwd_linear53_dwt2d_r4(z,ni,nj,nx,ny)   ! 2D forward transform, along i first, then along j
! in place FORWARD lifting transform using linear prediction wavelet for ARGUMENT_TYPE numbers
  use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
  implicit none
  integer, intent(IN) :: ni, nj, nx, ny
  ARGUMENT_TYPE, intent(INOUT), dimension(0:nx-1,0:ny-1) :: z

  ARGUMENT_TYPE, dimension(-1:ni) :: even, odd
  integer :: j, j00, jm1, jm2, jp1, jj
  integer :: nodd, neven

  nodd  = ishft(ni,-1)     ! number of odd terms
  neven = ishft(ni+1,-1)   ! number of even terms (nodd + 1 if ni is odd)
  j00 = 0
  do j = 1,nj-1,2
    jm2 = abs(j - 2)
    jm1 = j - 1
    jp1 = j + 1
    if(jp1 == nj) jp1 = nj - 2                              ! upper mirror boundary condition
    do jj = j00,max(j,jp1)
      call fwd_linear53_dwt1d_r4(z(:,jj))
    end do
    j00 = jp1 + 1
    z(0:ni-1,  j) = z(0:ni-1,  j) -  .5 * (z(0:ni-1,jm1) + z(0:ni-1,jp1))       ! predict odd rows
    z(0:ni-1,jm1) = z(0:ni-1,jm1) + .25 * (z(0:ni-1,jm2) + z(0:ni-1,  j))       ! update even rows (below odd row)
  end do
  if(mod(nj,2)==1) z(0:ni-1,nj-1) = z(0:ni-1,nj-1) + .5 * z(0:ni-1,nj-2)   ! odd number of rows, last row is en even row
  return
contains
  subroutine fwd_linear53_dwt1d_r4(f)   ! 1D along i transform
      use, intrinsic :: iso_fortran_env
    implicit none
    ARGUMENT_TYPE, intent(INOUT), dimension(0:ni-1) :: f

    integer :: i

    do i=0,nodd-1           ! split into even / odd
      even(i) = f(2*i)
      odd(i)  = f(2*i+1)
    end do
    if(iand(ni,1) /= 0) then  ! is ni odd ?
      even(nodd) = f(ni-1)      ! one more even values than odd values if ni is odd
    else
      even(nodd) = even(nodd-1) ! mirror condition at upper boundary if ni is even
    endif
    do i = 0, nodd-1           ! predict odd values (and copy updated values into f)
      odd(i) = odd(i) - .5 * (even(i) + even(i+1))
      f(neven+i) = odd(i)
    end do
    odd(-1)   = odd(0)         ! mirror condition at lower boundary
    odd(nodd) = odd(nodd-1)    ! mirror condition at upper boundary, used only of ni is odd
    do  i=0,neven-1            ! update even values (and copy updated values into f)
      f(i) = even(i) + .25 * (odd(i) + odd(i-1))
    end do
  end subroutine fwd_linear53_dwt1d_r4
end subroutine fwd_linear53_dwt2d_r4
#if defined(SELF_TEST)
#define NI 12000
#define NJ 8000
#define NREP 1
program test
  ARGUMENT_TYPE, dimension(0:NI-1,0:NJ-1) :: z, z0, z1
  integer :: i, j, ni, nj, ni2, nj2, ni4, nj4, ni8, nj8
  real(kind=REAL64), dimension(NREP) :: T0,T1,T2,T3
  real(kind=REAL64), external :: MPI_WTIME
  ni = NI
  nj = NJ
  do j = 0,nj-1
  do i = 0,ni-1
    z(i,j) = 1 + i*1.2 + j*1.7
    z0(i,j) = z(i,j)
  end do
  end do
  print *,'min,max z0',minval(z0),maxval(z0)
1 format(12F6.2)
  if(ni <=10 .and. nj <=10) then
    print *,'=========================================================='
    do j = nj-1,0,-1
      print 1,z(:,j)
    end do
  endif
  ni2 = (ni+1)/2
  nj2 = (nj+1)/2
  ni4 = (ni2+1)/2
  nj4 = (nj2+1)/2
  ni8 = (nj4+1)/2
  nj8 = (nj4+1)/2
  do irep = 1, NREP
  T0(irep) = MPI_WTIME()
  call fwd_linear53_dwt2d_r4(z,ni,nj,ni,nj)
  call fwd_linear53_dwt2d_r4(z,ni2,nj2,ni*2,nj2)
  call fwd_linear53_dwt2d_r4(z,ni4,nj4,ni*4,nj4)
  call fwd_linear53_dwt2d_r4(z,ni8,nj8,ni*8,nj8)
  T1(irep) = MPI_WTIME()
  if(ni <=10 .and. nj <=10) then
    print *,'=========================================================='
    do j = nj-1,0,-1
      print 1,z(:,j)
    end do
  endif
  T2(irep) = MPI_WTIME()
  call inv_linear53_dwt2d_r4(z,ni8,nj8,ni*8,nj8)
  call inv_linear53_dwt2d_r4(z,ni4,nj4,ni*4,nj4)
  call inv_linear53_dwt2d_r4(z,ni2,nj2,ni*2,nj2)
  call inv_linear53_dwt2d_r4(z,ni,nj,ni,nj)
  T3(irep) = MPI_WTIME()
  end do
  if(ni <=10 .and. nj <=10) then
    print *,'=========================================================='
    do j = nj-1,0,-1
      print 1,z(:,j)
    end do
  endif
  t1 = t1-t0
  t3 = t3-t2
  print *,'min transform time:',minval(t1),minval(t3),minval(t1)/ni/nj,minval(t3)/ni/nj
  print *,'max transform time:',maxval(t1),maxval(t3),maxval(t1)/ni/nj,maxval(t3)/ni/nj
  print *,'min,max z0',minval(z0),maxval(z0)
  print *,'min,max z',minval(z),maxval(z)
  z1 = abs(z0-z)/z0
  print *,'average,max,min per point relative error:',sum(z1)/ni/nj,maxval(z1),minval(z1)
  z1 = (z0-z)
  print *,'average,max,min per point absolute error:',sum(z1)/ni/nj,maxval(z1),minval(z1)
  print *,'=========================================================='
  print *,'after linear filtering'

  call low_pass_dwt2d_r4(z,ni,nj,ni,nj)
  z1 = abs(z0-z)/z0
  print *,'average,max,min per point relative error:',sum(z1)/ni/nj,maxval(z1),minval(z1)
  z1 = (z0-z)
  print *,'average,max,min per point absolute error:',sum(z1)/ni/nj,maxval(z1),minval(z1)

  print *,'=========================================================='
  print *,'after linear filtering and quantification'
  z = z0
  call low_pass_quant_dwt2d_r4(z,ni,nj,ni,nj)
  z1 = abs(z0-z)/z0
  print *,'z0 average:',sum(z0)/ni/nj
  print *,'average,max,min per point relative error:',sum(z1)/ni/nj,maxval(z1),minval(z1)
  z1 = (z0-z)
  print *,'average,max,min per point absolute error:',sum(z1)/ni/nj,maxval(z1),minval(z1)
  stop
end program test
#endif

