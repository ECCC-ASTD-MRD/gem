#if defined (AIX)
      subroutine vpownn(r,x,y,n)
      real*8 x(*),y(*),r(*)
      real*8 t(n)
      call vlog(t,x,n)
      do 10 j=1,n
      t(j)=t(j)*y(j)
   10 continue
      call vexp(r,t,n)
      return
      end
      subroutine vpow1n(r,x,y,n)
      real*8 x,y(*),r(*)
      real*8 t(n),tt
      tt = log(x)
      do 10 j=1,n
      t(j)=tt*y(j)
   10 continue
      call vexp(r,t,n)
      return
      end
      subroutine vpown1(r,x,y,n)
      real*8 x(*),y,r(*)
      real*8 t(n)
      call vlog(t,x,n)
      do 10 j=1,n
      t(j)=t(j)*y
   10 continue
      call vexp(r,t,n)
      return
      end
      subroutine vspownn(r,x,y,n)
      real*4 x(*),y(*),r(*)
      real (kind=8), dimension(n) :: t8,x8,y8,r8
      do i = 1,n
         x8(i) = x(i)
         y8(i) = y(i)
      enddo
      call vlog(t8,x8,n)
      do 10 j=1,n
      t8(j)=t8(j)*y8(j)
   10 continue
      call vexp(r8,t8,n)
      do i = 1,n
         r(i) = r8(i)
      enddo
      return
      end
      subroutine vspow1n(r,x,y,n)
      real*4 x,y(*),r(*)
      real*8 t8(n),tt,y8(n),r8(n)
      do i = 1,n
         y8(i) = y(i)
      enddo
      tt = log(x)
      do 10 j=1,n
      t8(j)=tt*y8(j)
   10 continue
      call vexp(r8,t8,n)
      do i = 1,n
         r(i) = r8(i)
      enddo
      return
      end
      subroutine vspown1(r,x,y,n)
      real*4 x(*),y,r(*)
      real*8 t8(n),x8(n),y8,r8(n)
      do i = 1,n
         x8(i) = x(i)
      enddo
      y8 = y
      call vlog(t8,x8,n)
      do 10 j=1,n
      t8(j)=t8(j)*y8
   10 continue
      call vexp(r8,t8,n)
      do i = 1,n
         r(i) = r8(i)
      enddo
      return
      end
#else
      subroutine vpow1n(r,x,y,n)
      real*8 x,y(*),r(*)
      do 10 j=1,n
      r(j)=x**y(j)
   10 continue
      return
      end
      subroutine vpown1(r,x,y,n)
      real*8 x(*),y,r(*)
      do 10 j=1,n
      r(j)=x(j)**y
   10 continue
      return
      end
      subroutine vpownn(r,x,y,n)
      real*8 x(*),y(*),r(*)
      do 10 j=1,n
      r(j)=x(j)**y(j)
   10 continue
      return
      end
      subroutine vspow1n(r,x,y,n)
      real*4 x,y(*),r(*)
      do 10 j=1,n
      r(j)=x**y(j)
   10 continue
      return
      end
      subroutine vspown1(r,x,y,n)
      real*4 x(*),y,r(*)
      do 10 j=1,n
      r(j)=x(j)**y
   10 continue
      return
      end
      subroutine vspownn(r,x,y,n)
      real*4 x(*),y(*),r(*)
      do 10 j=1,n
      r(j)=x(j)**y(j)
   10 continue
      return
      end
#endif