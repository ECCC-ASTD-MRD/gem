      subroutine vdiv(z,x,y,n)
      real*8 x(*),y(*),z(*)
      do 10 j=1,n
      z(j)=x(j)/y(j)
   10 continue
      return
      end
