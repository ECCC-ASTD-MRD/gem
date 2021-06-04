      subroutine vstan(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=tan(x(j))
   10 continue
      return
      end
