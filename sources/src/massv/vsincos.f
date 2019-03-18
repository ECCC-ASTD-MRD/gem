      subroutine vsincos(x,y,z,n)
      real*8 x(*),y(*),z(*)
      do 10 j=1,n
      x(j)=sin(z(j))
      y(j)=cos(z(j))
   10 continue
      return
      end
