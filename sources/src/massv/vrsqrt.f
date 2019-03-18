      subroutine vrsqrt(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=1.d0/sqrt(x(j))
   10 continue
      return
      end
