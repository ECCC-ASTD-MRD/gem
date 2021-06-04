      subroutine vscosisin(y,x,n)
      complex*8 y(*)
      real*4 x(*)
      do 10 j=1,n
      y(j)= cmplx(cos(x(j)),sin(x(j)))
   10 continue
      return
      end
