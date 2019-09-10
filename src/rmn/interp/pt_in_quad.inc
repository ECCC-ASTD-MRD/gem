   logical function pt_in_quad(x,y,x1,y1,x2,y2,x3,y3,x4,y4)
   implicit none

   real x,y,x1,y1,x2,y2,x3,y3,x4,y4
   logical pt_in_triangle
   external pt_in_triangle

   pt_in_quad = .false.
   if (pt_in_triangle(x,y,x1,y1,x2,y2,x3,y3)) then
      pt_in_quad = .true.
   else
      pt_in_quad = pt_in_triangle(x,y,x1,y1,x3,y3,x4,y4)
   endif

   return
   end function pt_in_quad

   logical function pt_in_triangle(x,y,x1,y1,x2,y2,x3,y3)
   implicit none

   real x,y,x1,y1,x2,y2,x3,y3
   real a,b,c,d, det, lambda1, lambda2, lambda3
   pt_in_triangle = .false.

   a = x1 - x3
   b = x2 - x3
   c = y1 - y3
   d = y2 - y3

   det = 1.0 / (a * d - b * c)

   lambda1 = (d * (x - x3) - b * (y - y3)) * det
   lambda2 = (a * (y - y3) - c * (x - x3)) * det
   lambda3 = 1.0 - lambda1 - lambda2

   if (lambda1 < 0.0 .or. lambda1 > 1.0) then
      return
   else if (lambda2 < 0.0 .or. lambda2 > 1.0) then
      return
   else if (lambda3 < 0.0 .or. lambda3 > 1.0) then
      return
   else
      pt_in_triangle =  .true.
   endif

   return
   end function pt_in_triangle

