      subroutine vdnint(y,x,n)
      real*8 x(*),y(*)
      do 10 j=1,n
      y(j)=dnint(x(j))
   10 continue
      return
      end
