      subroutine vspown1(r,x,y,n)
      real*4 x(*),y,r(*)
      do 10 j=1,n
      r(j)=x(j)**y
   10 continue
      return
      end
