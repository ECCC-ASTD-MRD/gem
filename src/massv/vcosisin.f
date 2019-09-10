      subroutine vcosisin(y,x,n)
      complex*16 y(*)
      real*8 x(*)
      do 10 j=1,n
      y(j)=dcmplx(cos(x(j)),sin(x(j)))
   10 continue
      return
      end
