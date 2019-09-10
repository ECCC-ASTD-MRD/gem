      subroutine vpow1n(r,x,y,n)
      real*8 x,y(*),r(*)
      do 10 j=1,n
      r(j)=x**y(j)
   10 continue
      return
      end
