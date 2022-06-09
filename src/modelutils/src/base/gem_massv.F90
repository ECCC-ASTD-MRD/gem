
#ifndef __INTEL_COMPILER
#ifdef USE_GEM_MKLMASSV
#undef USE_GEM_MKLMASSV
#endif
#endif

#ifdef USE_GEM_MKLMASSV

!!include "mkl_vml.f90"
subroutine gem_vasin(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdasin( n, x, y )

   return
end subroutine gem_vasin

subroutine gem_vatan2(z,y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: y,x
   real(kind=8), dimension(n), intent(out) :: z

   call vdatan2( n, y, x, z )

   return
end subroutine gem_vatan2

subroutine gem_vcos(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdcos( n, x, y )

   return
end subroutine gem_vcos

subroutine gem_vcosisin(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   complex(kind=8), dimension(n), intent(out) :: y

   call vzcis( n, x, y )

   return
end subroutine gem_vcosisin

subroutine gem_vdint(y,x,n)
   !removes fractional part of the number.
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdtrunc( n, x, y )

   return
end subroutine gem_vdint

subroutine gem_vdiv(z,x,y,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x,y
   real(kind=8), dimension(n), intent(out) :: z      

   call vddiv( n, x, y, z )

   return
end subroutine gem_vdiv

subroutine gem_vdnint(y,x,n)
   !anint- rounds to nearest integer value
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdround( n, x, y )

   return
end subroutine gem_vdnint

subroutine gem_vexp(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdexp( n, x, y )

   return
end subroutine gem_vexp

subroutine gem_vlog(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdln( n, x, y )

   return
end subroutine gem_vlog

subroutine gem_vpow1n(r,x,y,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), intent(in) :: x
   real(kind=8), dimension(n) :: x1
   real(kind=8), dimension(n), intent(in) :: y
   real(kind=8), dimension(n), intent(out) :: r

   x1(:)=x
   call vdpowx( n, x1, y, r )

   return
end subroutine gem_vpow1n

subroutine gem_vpown1(r,x,y,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), intent(in) :: y
   real(kind=8), dimension(n) :: y1
   real(kind=8), dimension(n), intent(out) :: r

   y1(:)=y
   call vdpowx( n, x, y1, r )

   return
end subroutine gem_vpown1

subroutine gem_vpownn(r,x,y,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x,y
   real(kind=8), dimension(n), intent(out) :: r

   call vdpowx( n, x, y, r )

   return
end subroutine gem_vpownn

subroutine gem_vrec(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdinv( n, x, y )

   return
end subroutine gem_vrec

subroutine gem_vrsqrt(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdinvsqrt( n, x, y )

   return
end subroutine gem_vrsqrt

subroutine gem_vsatan2(z,y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: y,x
   real(kind=4), dimension(n), intent(out) :: z

   call vdatan2( n, y, x, z )

   return
end subroutine gem_vsatan2

subroutine gem_vscos(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y

   call vscos( n, x, y )

   return
end subroutine gem_vscos

subroutine gem_vscosisin(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   complex(kind=4), dimension(n), intent(out) :: y

   call vccis( n, x, y )

   return
end subroutine gem_vscosisin

subroutine gem_vsdiv(z,x,y,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x,y
   real(kind=4), dimension(n), intent(out) :: z

   call vsdiv( n, x, y, z )

   return
end subroutine gem_vsdiv

subroutine gem_vsexp(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y

   call vsexp( n, x, y )

   return
end subroutine gem_vsexp

subroutine gem_vsin(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdsin( n, x, y )

   return
end subroutine gem_vsin

subroutine gem_vsincos(x,y,z,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: z
   real(kind=8), dimension(n), intent(out) :: x,y

   call vdsincos( n, z, x, y )

   return
end subroutine gem_vsincos

subroutine gem_vslog(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y

   call vsln( n, x, y )

   return
end subroutine gem_vslog

subroutine gem_vspow1n(r,x,y,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), intent(in) :: x
   real(kind=4), dimension(n) :: x1
   real(kind=4), dimension(n), intent(in) :: y
   real(kind=4), dimension(n), intent(out) :: r

   x1(:)=x
   call vspowx( n, x1, y, r )

   return
end subroutine gem_vspow1n

subroutine gem_vspown1(r,x,y,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), intent(in) :: y
   real(kind=4), dimension(n) :: y1
   real(kind=4), dimension(n), intent(out) :: r

   y1(:)=y
   call vspowx( n, x, y1, r )

   return
end subroutine gem_vspown1

subroutine gem_vspownn(r,x,y,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x,y
   real(kind=4), dimension(n), intent(out) :: r

   call vdpowx( n, x, y, r )

   return
end subroutine gem_vspownn

subroutine gem_vsqrt(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdsqrt( n, x, y )

   return
end subroutine gem_vsqrt

subroutine gem_vsrec(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y

   call vsinv( n, x, y )

   return
end subroutine gem_vsrec

subroutine gem_vsrsqrt(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y

   call vsinvsqrt( n, x, y )

   return
end subroutine gem_vsrsqrt

subroutine gem_vssin(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y

   call vssin( n, x, y )

   return
end subroutine gem_vssin

subroutine gem_vssincos(x,y,z,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: z
   real(kind=4), dimension(n), intent(out) :: x,y

   call vssincos( n, z, x, y )

   return
end subroutine gem_vssincos

subroutine gem_vssqrt(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y

   call vssqrt( n, x, y )

   return
end subroutine gem_vssqrt

subroutine gem_vstan(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y

   call vstan( n, x, y )

   return
end subroutine gem_vstan

subroutine gem_vtan(y,x,n)
   implicit none
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y

   call vdtan( n, x, y )

   return
end subroutine gem_vtan


#else


subroutine gem_vasin(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=asin(x(j))
   enddo
   return
end subroutine gem_vasin

subroutine gem_vatan2(z,y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: y,x
   real(kind=8), dimension(n), intent(out) :: z
   do j = 1, n
      z(j)=atan2(y(j),x(j))
   enddo
   return
end subroutine gem_vatan2

subroutine gem_vcos(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=cos(x(j))
   enddo
   return
end subroutine gem_vcos

subroutine gem_vcosisin(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   complex(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=dcmplx(cos(x(j)),sin(x(j)))
   enddo
   return
end subroutine gem_vcosisin

subroutine gem_vdint(y,x,n)
   !removes fractional part of the number.
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=dint(x(j))
   enddo
   return
end subroutine gem_vdint

subroutine gem_vdiv(z,x,y,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x,y
   real(kind=8), dimension(n), intent(out) :: z      
   do j = 1, n
      z(j)=x(j)/y(j)
   enddo
   return
end subroutine gem_vdiv

subroutine gem_vdnint(y,x,n)
   !anint- rounds to nearest integer value
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=dnint(x(j))
   enddo
   return
end subroutine gem_vdnint

subroutine gem_vexp(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=exp(x(j))
   enddo
   return
end subroutine gem_vexp

subroutine gem_vlog(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=log(x(j))
   enddo
   return
end subroutine gem_vlog

subroutine gem_vpow1n(r,x,y,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), intent(in) :: x
   real(kind=8), dimension(n), intent(in) :: y
   real(kind=8), dimension(n), intent(out) :: r
   do j = 1, n
      r(j)=x**y(j)
   enddo
   return
end subroutine gem_vpow1n

subroutine gem_vpown1(r,x,y,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), intent(in) :: y
   real(kind=8), dimension(n), intent(out) :: r
   do j = 1, n
      r(j)=x(j)**y
   enddo
   return
end subroutine gem_vpown1

subroutine gem_vpownn(r,x,y,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x,y
   real(kind=8), dimension(n), intent(out) :: r
   do j = 1, n
      r(j)=x(j)**y(j)
   enddo
   return
end subroutine gem_vpownn

subroutine gem_vrec(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=1.d0/x(j)
   enddo
   return
end subroutine gem_vrec

subroutine gem_vrsqrt(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=1.d0/sqrt(x(j))
   enddo
   return
end subroutine gem_vrsqrt

subroutine gem_vsatan2(z,y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: y,x
   real(kind=4), dimension(n), intent(out) :: z
   do j = 1, n
      z(j)=atan2(y(j),x(j))
   enddo
   return
end subroutine gem_vsatan2

subroutine gem_vscos(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=cos(x(j))
   enddo
   return
end subroutine gem_vscos

subroutine gem_vscosisin(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   complex(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=cmplx(cos(x(j)),sin(x(j)))
   enddo
   return
end subroutine gem_vscosisin

subroutine gem_vsdiv(z,x,y,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x,y
   real(kind=4), dimension(n), intent(out) :: z
   do j = 1, n
      z(j)=x(j)/y(j)
   enddo
   return
end subroutine gem_vsdiv

subroutine gem_vsexp(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=exp(x(j))
   enddo
   return
end subroutine gem_vsexp

subroutine gem_vsin(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=sin(x(j))
   enddo
   return
end subroutine gem_vsin

subroutine gem_vsincos(x,y,z,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: z
   real(kind=8), dimension(n), intent(out) :: x,y
   do j = 1, n
      x(j)=sin(z(j))
      y(j)=cos(z(j))
   enddo
   return
end subroutine gem_vsincos

subroutine gem_vslog(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=log(x(j))
   enddo
   return
end subroutine gem_vslog

subroutine gem_vspow1n(r,x,y,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), intent(in) :: x
   real(kind=4), dimension(n), intent(in) :: y
   real(kind=4), dimension(n), intent(out) :: r
   do j = 1, n
      r(j)=x**y(j)
   enddo
   return
end subroutine gem_vspow1n

subroutine gem_vspown1(r,x,y,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), intent(in) :: y
   real(kind=4), dimension(n), intent(out) :: r
   do j = 1, n
      r(j)=x(j)**y
   enddo
   return
end subroutine gem_vspown1

subroutine gem_vspownn(r,x,y,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x,y
   real(kind=4), dimension(n), intent(out) :: r
   do j = 1, n
      r(j)=x(j)**y(j)
   enddo
   return
end subroutine gem_vspownn

subroutine gem_vsqrt(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=sqrt(x(j))
   enddo
   return
end subroutine gem_vsqrt

subroutine gem_vsrec(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=1.d0/x(j)
   enddo
   return
end subroutine gem_vsrec

subroutine gem_vsrsqrt(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=1.d0/sqrt(x(j))
   enddo
   return
end subroutine gem_vsrsqrt

subroutine gem_vssin(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=sin(x(j))
   enddo
   return
end subroutine gem_vssin

subroutine gem_vssincos(x,y,z,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: z
   real(kind=4), dimension(n), intent(out) :: x,y
   do j = 1, n
      x(j)=sin(z(j))
      y(j)=cos(z(j))
   enddo
   return
end subroutine gem_vssincos

subroutine gem_vssqrt(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=sqrt(x(j))
   enddo
   return
end subroutine gem_vssqrt

subroutine gem_vstan(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=4), dimension(n), intent(in) :: x
   real(kind=4), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=tan(x(j))
   enddo
   return
end subroutine gem_vstan

subroutine gem_vtan(y,x,n)
   implicit none
   integer :: j
   integer, intent(in) :: n
   real(kind=8), dimension(n), intent(in) :: x
   real(kind=8), dimension(n), intent(out) :: y
   do j = 1, n
      y(j)=tan(x(j))
   enddo
   return
end subroutine gem_vtan

#endif
