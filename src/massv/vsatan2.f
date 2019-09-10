      subroutine vsatan2(z,y,x,n)
      real*4 x(*),y(*),z(*)
      do 10 j=1,n
      z(j)=atan2(y(j),x(j))
   10 continue
      return
      end
